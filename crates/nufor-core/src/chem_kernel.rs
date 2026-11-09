//! runtime-width chemistry kernel adapter: advance one cell's
//! composition over a timestep at fixed total energy, whichever
//! mechanism the case selected.
//!
//! two arms sit behind one call:
//!
//! - the lithium/hydrogen arm reproduces the legacy nine-species
//!   kernel's arithmetic operation for operation: the same closure
//!   solve, the same `omega` call, the same rate-scaled positivity
//!   handling, the same accumulation order, and no sanitising
//!   afterwards. the numbers the legacy march produced do not move,
//!   and the parity test holds the two to each other bit for bit.
//! - the GRI-Mech arm steps the fifty-three-species composition
//!   implicitly: backward euler with a newton solve on the mass
//!   fractions, the temperature riding the fixed-energy constraint. a
//!   solve that does not converge leaves the cell alone and defers it.
//!
//! the widths are decided at runtime from [`ChemLane`], so a caller
//! never has to special-case the mechanism.

use crate::chem_lane::ChemLane;
use crate::combustion as comb;
use crate::combustion_closure as closure9;
use crate::combustion_thermo as ct;
use crate::error::Error;

/// what one cell's chemistry step did.
#[derive(Clone, Copy, Debug)]
pub struct CellOutcome {
    /// the chemistry timescale violation: the requested step measured
    /// in units of the adaptive sub-step the arm took. zero on the
    /// methane arm, which is unconditionally stable.
    pub oversub: f64,
    /// whether the cell was deferred, meaning its chemistry did not
    /// converge and its composition came through untouched.
    pub deferred: bool,
}

/// advance one cell's composition over `dt` at fixed total energy.
///
/// `e_int` is the cell's specific internal energy in J/kg, held fixed
/// across the step; the temperature follows from the composition and
/// that energy. `rho` is the density in kg/m3 and `y` the mass
/// fractions, one per species of the lane, updated in place. the return
/// value is the chemistry timescale violation, the march's stiffness
/// measure.
pub fn advance_cell(
    lane: ChemLane,
    e_int: f64,
    rho: f64,
    y: &mut [f64],
    dt: f64,
) -> Result<f64, Error> {
    Ok(advance_cell_outcome(lane, e_int, rho, y, dt)?.oversub)
}

/// the same cell step, with the deferral flag alongside the violation.
pub fn advance_cell_outcome(
    lane: ChemLane,
    e_int: f64,
    rho: f64,
    y: &mut [f64],
    dt: f64,
) -> Result<CellOutcome, Error> {
    if y.len() != lane.n_species() || dt <= 0.0 || !dt.is_finite() {
        return Err(Error::InvalidArgs);
    }
    if !e_int.is_finite() || !rho.is_finite() {
        return Err(Error::InvalidArgs);
    }
    match lane {
        ChemLane::Li2004 => li_arm(e_int, rho, y, dt),
        ChemLane::Gri30 => gri_arm(e_int, rho, y, dt),
    }
}

/// what one field pass did.
#[derive(Clone, Copy, Debug, Default)]
pub struct FieldReport {
    /// the largest timescale violation over the field.
    pub worst_oversub: f64,
    /// the cells whose chemistry stepped cleanly.
    pub cells_stepped: usize,
    /// the cells whose implicit solve did not converge and were left
    /// unchanged.
    pub cells_deferred: usize,
}

/// advance a whole field of cells, each independently, over `dt` at
/// fixed total energy.
///
/// `y` is one row per species (`y[s][cell]`), `density` and
/// `internal_energy` are per cell. each cell reads only its own
/// composition, density and energy, so the pass holds no shared mutable
/// state and a slab-split run lands on the same bits as a serial one.
/// the return value is the largest timescale violation over the field.
pub fn advance_field(
    lane: ChemLane,
    y: &mut [Vec<f64>],
    density: &[f64],
    internal_energy: &[f64],
    dt: f64,
) -> Result<f64, Error> {
    Ok(advance_field_report(lane, y, density, internal_energy, dt)?.worst_oversub)
}

/// the same field pass, with the tallies alongside the violation.
pub fn advance_field_report(
    lane: ChemLane,
    y: &mut [Vec<f64>],
    density: &[f64],
    internal_energy: &[f64],
    dt: f64,
) -> Result<FieldReport, Error> {
    let ncells = y.first().map_or(0, |row| row.len());
    if y.len() != lane.n_species() {
        return Err(Error::InvalidArgs);
    }
    if y.iter().any(|row| row.len() != ncells) {
        return Err(Error::InvalidArgs);
    }
    if density.len() != ncells || internal_energy.len() != ncells {
        return Err(Error::InvalidArgs);
    }
    if dt <= 0.0 || !dt.is_finite() {
        return Err(Error::InvalidArgs);
    }
    // one buffer for the cell being worked on, refilled per cell, so the
    // pass allocates once and shares nothing between cells.
    let mut parcel = vec![0.0f64; lane.n_species()];
    let mut report = FieldReport::default();
    for cell in 0..ncells {
        for (species, row) in y.iter().enumerate() {
            parcel[species] = row[cell];
        }
        let rho = density[cell];
        let e_int = internal_energy[cell];
        if !rho.is_finite() || !e_int.is_finite() || rho <= 0.0 {
            report.cells_deferred += 1;
            continue;
        }
        match advance_cell_outcome(lane, e_int, rho, &mut parcel, dt) {
            Ok(outcome) => {
                if outcome.deferred {
                    report.cells_deferred += 1;
                } else {
                    for (species, row) in y.iter_mut().enumerate() {
                        row[cell] = parcel[species];
                    }
                    report.worst_oversub = report.worst_oversub.max(outcome.oversub);
                    report.cells_stepped += 1;
                }
            }
            Err(_) => {
                report.cells_deferred += 1;
            }
        }
    }
    Ok(report)
}

/// the nine-species hydrogen kernel.
///
/// this mirrors the legacy march's cell step operation for operation:
/// the same closure solve, the same `omega` call, the same rate scaling
/// for positivity, the same accumulation order, and no sanitising
/// anywhere afterwards. any change here moves the legacy numbers, so
/// the parity test in `tests/chem_kernel_parity.rs` is the gate.
fn li_arm(e_int: f64, rho: f64, y: &mut [f64], dt: f64) -> Result<CellOutcome, Error> {
    let mut a = [0.0f64; ct::NS];
    a.copy_from_slice(y);

    let mut worst_oversub = 0.0f64;
    {
        let mut t = closure9::temperature_at(e_int, &a);
        let mut remaining = dt;
        let mut steps = 0;
        while remaining > 1e-14 && steps < 512 {
            let (dydt, _) = comb::omega(rho, &a, t);
            let rate_scale = dydt.iter().fold(0.0f64, |acc, d| acc.max(d.abs()));
            let tau = if rate_scale > 1e-30 {
                1.0 / rate_scale
            } else {
                remaining
            };
            let h = remaining.min(tau * 0.2);
            // positivity by RATE SCALING, not by shrinking the step:
            // the reaction vector is scaled down by the largest factor
            // that keeps every species non-negative over this step. the
            // stoichiometry is untouched, so the species sum and the
            // element pools stay exactly conserved.
            let mut scale = 1.0f64;
            for (s, &d) in dydt.iter().enumerate() {
                if d < 0.0 && a[s] > 0.0 {
                    let allowed = a[s] / (-d * h);
                    if allowed < scale {
                        scale = allowed;
                    }
                }
            }
            for (ys, d) in a.iter_mut().zip(dydt.iter()) {
                let d = d * h * scale;
                *ys += d;
            }
            t = closure9::temperature_at(e_int, &a);
            remaining -= h;
            steps += 1;
            if steps == 1 {
                worst_oversub = worst_oversub.max(dt / (tau * 0.2));
            }
        }
    }
    y.copy_from_slice(&a);
    Ok(CellOutcome {
        oversub: worst_oversub,
        deferred: false,
    })
}

/// the fifty-three-species methane kernel.
///
/// the temperature rides the cell's fixed specific internal energy and
/// the composition is stepped implicitly by backward euler, the same
/// step the zero-dimensional parcel reactor takes. the energy handed to
/// the solver is the caller's per-mass value on the mechanism's own
/// thermodynamic basis, which is the same basis the lane's closure
/// reports. a solve that does not converge leaves the cell's
/// composition exactly as it came in.
fn gri_arm(e_int: f64, rho: f64, y: &mut [f64], dt: f64) -> Result<CellOutcome, Error> {
    use crate::gri30_reactor;

    let mut a = [0.0f64; gri30_reactor::NS];
    a.copy_from_slice(y);

    let mut iterations = 0u32;
    match gri30_reactor::backward_euler_step(&a, rho, e_int, dt, &mut iterations) {
        Ok(Some(stepped)) => {
            // the reactor's own output, copied straight across: the
            // kernel adds no cleaning of its own, so its composition is
            // the reactor's composition bit for bit. whoever consumes it
            // decides the hygiene policy.
            for (dst, src) in y.iter_mut().zip(stepped.iter()) {
                *dst = *src;
            }
            Ok(CellOutcome {
                oversub: 0.0,
                deferred: false,
            })
        }
        Ok(None) => {
            // no convergence: the cell keeps its incoming composition.
            Ok(CellOutcome {
                oversub: 0.0,
                deferred: true,
            })
        }
        Err(_message) => Err(Error::KernelFailure),
    }
}

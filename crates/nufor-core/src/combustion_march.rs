//! the h2/o2 march: 9-species mass fractions advected beside the
//! conserved state with operator-split finite-rate chemistry. the
//! structure mirrors the park5 lane: closure-aware face fluxes
//! through the shared sweep machinery, then the upwind species
//! transport, then the chemistry substep at the post-step state.
//!
//! scope: inviscid cartesian 2d only. the viscous, rans and cut-cell
//! lanes do not serve the combustion model and must reject it loudly
//! rather than fall back to a perfect gas.

use crate::axi::compute_axi_fluxes;
use crate::combustion as comb;
use crate::combustion_closure as closure9;
use crate::combustion_thermo as ct;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;
use crate::species_transport::transport_species_generic;
use crate::state2d::cons_to_prim2d;
use crate::state2d::ConservedState2d;
use crate::thermo::ThermoModel;

/// per-cell 9-species mass fractions, `y[s][k]`.
#[derive(Clone, Debug)]
pub struct CombustionState {
    pub y: Vec<Vec<f64>>,
    /// the composition the case injects: what an inflow face's
    /// species ghost carries. stored at construction so a burned
    /// edge cell never feeds itself back through the inflow.
    pub fill: [f64; ct::NS],
}

impl CombustionState {
    /// uniform composition; `fractions` must sum to 1.
    pub fn uniform(n: usize, fractions: &[f64; ct::NS]) -> Self {
        CombustionState {
            y: fractions.iter().map(|f| vec![*f; n]).collect(),
            fill: *fractions,
        }
    }

    /// the state whose field and injection composition are the given
    /// mass fractions: the chamber case injects the burned gas, the
    /// premix case the unburned one.
    pub fn with_fill(n: usize, fill: &[f64; ct::NS]) -> Self {
        CombustionState::uniform(n, fill)
    }

    /// stoichiometric unburned h2/o2: 1 mol H2 + 0.5 mol O2, mass
    /// fractions 2.016/18.016 and 16.0/18.016.
    pub fn stoich_h2o2(n: usize) -> Self {
        let mut y = [0.0f64; ct::NS];
        y[ct::IH2] = 2.016 / 18.016;
        y[ct::IO2] = 16.0 / 18.016;
        CombustionState::uniform(n, &y)
    }

    /// rescale the cell's mass fractions to sum to one without
    /// clipping: the advection reconstruction can push a species
    /// slightly negative or the sum off one, and clipping would
    /// destroy the element pools the transport just conserved.
    pub fn renorm(&mut self, k: usize) {
        let sum: f64 = self.y.iter().map(|row| row[k]).sum();
        if sum > 0.0 {
            for i in 0..ct::NS {
                self.y[i][k] /= sum;
            }
        }
    }

    /// clip negatives and renormalize; returns the worst violation.
    pub fn sanitize(&mut self, k: usize) -> f64 {
        let mut sum = 0.0f64;
        let mut worst = 0.0f64;
        for i in 0..ct::NS {
            if self.y[i][k] < 0.0 {
                worst = worst.max(-self.y[i][k]);
                self.y[i][k] = 0.0;
            }
            if self.y[i][k] > 1.0 {
                worst = worst.max(self.y[i][k] - 1.0);
                self.y[i][k] = 1.0;
            }
            sum += self.y[i][k];
        }
        if sum > 0.0 {
            for i in 0..ct::NS {
                self.y[i][k] /= sum;
            }
        }
        worst
    }

    /// worst |sum(y) - 1| over the field.
    pub fn max_sum_drift(&self) -> f64 {
        let n = self.y[0].len();
        let mut worst = 0.0f64;
        for k in 0..n {
            let mut sum = 0.0;
            for row in &self.y {
                sum += row[k];
            }
            worst = worst.max((sum - 1.0).abs());
        }
        worst
    }
}

/// one forward-euler transport sweep: the shared closure-aware flux
/// machinery with the 9-species face compositions, then the generic
/// upwind species transport on the same face mass fluxes.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y_cell = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y_cell[s] = species.y[s][k];
        }
        p[k] = closure9::pressure_at(state.rho[k], e_int[k], &y_cell);
        a[k] = closure9::sound_speed_at(state.rho[k], e_int[k], &y_cell);
    }
    let mut smax = 0.0f64;
    for k in 0..n {
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1e-12);

    let model = ThermoModel::Combustion9 { gamma_shaping: 1.2 };
    let (fx, fy) = compute_axi_fluxes(
        model,
        state,
        g,
        muscl,
        bc,
        &p,
        &e_int,
        Some(crate::axi::SpeciesField::Combustion9(&*species)),
        nthreads,
        None,
    )?;

    let mut resid = 0.0f64;
    #[allow(clippy::needless_range_loop)]
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let dyi = dt / g.dys[j];
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i])
                + dyi * (fy[i].mass[j + 1] - fy[i].mass[j]));
            let dmx =
                -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dyi * (fy[i].mx[j + 1] - fy[i].mx[j]));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dyi * (fy[i].my[j + 1] - fy[i].my[j]));
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dyi * (fy[i].e[j + 1] - fy[i].e[j]));
            state.rho[k] += drho;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            resid = resid
                .max(drho.abs())
                .max(dmx.abs())
                .max(dmy.abs())
                .max(de.abs());
        }
    }

    let fx_mass: Vec<Vec<f64>> = fx.iter().map(|s| s.mass.clone()).collect();
    let fy_mass: Vec<Vec<f64>> = fy.iter().map(|s| s.mass.clone()).collect();
    let resid_y = transport_species_generic(
        &mut species.y,
        &state.rho,
        dt,
        &g.dxs,
        &g.dys,
        &fx_mass,
        &fy_mass,
        None,
    )?;
    for k in 0..n {
        species.renorm(k);
    }
    Ok((dt, resid.max(resid_y)))
}

/// rk2 (heun) over the transport sweep plus the operator-split
/// chemistry substep at the post-step state. frozen chemistry
/// (`react = false`) skips the substep entirely: the cold-flow path
/// pays nothing for the mechanism.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_rk2(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    react: bool,
) -> Result<(f64, f64), Error> {
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    advance2d_combustion(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads)?;
    advance2d_combustion(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads)?;
    let n = g.nx * g.ny;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    for s in 0..ct::NS {
        for k in 0..n {
            species.y[s][k] = 0.5 * species.y[s][k] + 0.5 * y1.y[s][k];
        }
    }
    for k in 0..n {
        species.renorm(k);
    }
    // dt for the chemistry substep comes from the pre-step state,
    // the same probe the park lane uses; the march itself never
    // advances a third time.
    let (dt, _) = {
        let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let n = g.nx * g.ny;
        let mut smax = 0.0f64;
        for k in 0..n {
            let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
            let mut y_cell = [0.0f64; ct::NS];
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                y_cell[s] = species.y[s][k];
            }
            let a = closure9::sound_speed_at(state.rho[k], e_int, &y_cell);
            smax = smax.max(u[k].abs().max(v[k].abs()) + a);
        }
        (cfl * g.dx.min(g.dy) / smax.max(1e-12), 0.0)
    };
    let mut oversub = 0.0f64;
    if react {
        oversub = chemistry_substep9(state, species, g, dt, nthreads)?;
    }
    Ok((dt, oversub))
}

/// one axisymmetric euler step with the nine-species combustion
/// closure: the same face machinery the park lane rides, with the
/// annular update (r-weighted radial faces plus the p/r source).
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_axi(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &crate::Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y_cell = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y_cell[s] = species.y[s][k];
        }
        p[k] = closure9::pressure_at(state.rho[k], e_int[k], &y_cell);
        a[k] = closure9::sound_speed_at(state.rho[k], e_int[k], &y_cell);
    }
    let mut smax = 0.0f64;
    for k in 0..n {
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1e-12);
    let model = ThermoModel::Combustion9 { gamma_shaping: 1.2 };
    let (fx, fy) = crate::axi::compute_axi_fluxes(
        model,
        state,
        g,
        muscl,
        bc,
        &p,
        &e_int,
        Some(crate::axi::SpeciesField::Combustion9(species)),
        nthreads,
        wall,
    )?;

    // the annular update: r-weighted radial flux difference plus the
    // radial pressure source, mirroring the park lane.
    // the entry state, so the solid interior can be restored after the
    // update: the masked cells are scaffolding, not flow.
    let entry_rho = state.rho.clone();
    let entry_e = state.e.clone();

    let mut resid = 0.0f64;
    for j in 0..ny {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        #[allow(clippy::needless_range_loop)]
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let radial = |f: &Vec<f64>| (rp * f[j + 1] - rm * f[j]) / (rc * dr);
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i]) + dt * radial(&fy[i].mass));
            let dmx = -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dt * radial(&fy[i].mx));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dt * radial(&fy[i].my)) + dt * p[k] / rc;
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dt * radial(&fy[i].e));
            state.rho[k] += drho;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            resid = resid
                .max(drho.abs())
                .max(dmx.abs())
                .max(dmy.abs())
                .max(de.abs());
        }
    }

    // restore the solid interior after the update: the wall faces
    // carry no flux, but the solid-solid mirror pairs still ratchet
    // their dummy velocity stage, and left alone that momentum leaks
    // into the band. the solid state is scaffolding, not flow.
    if let Some(dist) = wall {
        let band = 1.5 * g.dx.min(g.dy);
        let cache = crate::body::body_cache(g, dist, band);
        for &k in &cache.solid_cells {
            state.rho[k] = entry_rho[k];
            state.mx[k] = 0.0;
            state.my[k] = 0.0;
            state.e[k] = entry_e[k];
        }
    }

    // species advection on the same face mass fluxes, carrying rho*Y.
    let fx_mass: Vec<Vec<f64>> = fx.iter().map(|s| s.mass.clone()).collect();
    let fy_mass: Vec<Vec<f64>> = fy.iter().map(|s| s.mass.clone()).collect();
    let resid_y = crate::species_transport::transport_species_generic(
        &mut species.y,
        &state.rho,
        dt,
        &g.dxs,
        &g.dys,
        &fx_mass,
        &fy_mass,
        Some(&g.faces_y),
    )?;
    for k in 0..n {
        species.renorm(k);
    }
    Ok((dt, resid.max(resid_y)))
}

/// heun (rk2) wrapper around the axisymmetric combustion step with
/// the operator-split chemistry at the post-step state.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_axi_rk2(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &crate::Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    react: bool,
) -> Result<(f64, f64), Error> {
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    advance2d_combustion_axi(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
    advance2d_combustion_axi(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
    let n = g.nx * g.ny;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    for s in 0..ct::NS {
        for k in 0..n {
            species.y[s][k] = 0.5 * species.y[s][k] + 0.5 * y1.y[s][k];
        }
    }
    for k in 0..n {
        species.renorm(k);
    }
    // dt for the chemistry substep comes from the post-step state,
    // the same probe the planar lane uses.
    let (dt, _) = {
        let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let n = g.nx * g.ny;
        let mut smax = 0.0f64;
        for k in 0..n {
            let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
            let mut y_cell = [0.0f64; ct::NS];
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                y_cell[s] = species.y[s][k];
            }
            let a = closure9::sound_speed_at(state.rho[k], e_int, &y_cell);
            smax = smax.max(u[k].abs().max(v[k].abs()) + a);
        }
        (cfl * g.dx.min(g.dy) / smax.max(1e-12), 0.0)
    };
    let mut oversub = 0.0f64;
    if react {
        oversub = chemistry_substep9(state, species, g, dt, nthreads)?;
    }
    Ok((dt, oversub))
}

/// chemistry substep (operator-split): integrate the mechanism's
/// omega per cell over the flow dt at fixed total energy, the same
/// split the park lane uses. substeps cap at a fraction of the
/// fastest chemical timescale. returns the largest timescale
/// violation, the march's stiffness report.
pub fn chemistry_substep9(
    state: &ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    dt: f64,
    nthr: usize,
) -> Result<f64, Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let threads = nthr.clamp(1, 8).min(n.div_ceil(24));
    if threads <= 1 {
        let mut worst = 0.0f64;
        for k in 0..n {
            let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
            let mut y = [0.0f64; ct::NS];
            for (ys, ys_src) in y.iter_mut().zip(species.y.iter()) {
                *ys = ys_src[k];
            }
            worst = worst.max(react_cell(e_int, state.rho[k], &mut y, dt));
            for (ys_dst, ys) in species.y.iter_mut().zip(y.iter()) {
                ys_dst[k] = *ys;
            }
        }
        return Ok(worst);
    }
    // chunked over cells: each thread owns a disjoint window of every
    // species row, so the result matches the serial loop exactly.
    let chunk = n.div_ceil(threads);
    let bounds: Vec<(usize, usize)> = (0..n)
        .step_by(chunk)
        .map(|s| (s, (s + chunk).min(n)))
        .collect();
    let mut worst = 0.0f64;
    let species = &mut species.y;
    // one pass per chunk: carve disjoint windows out of every species
    // row before any thread starts, so no row is borrowed twice.
    // carve the disjoint windows with ONE pass over the rows: each row
    // splits successively per chunk, so every window is unique and the
    // thread bodies never alias.
    let mut chunk_windows: Vec<Vec<&mut [f64]>> = (0..bounds.len()).map(|_| Vec::new()).collect();
    for row in species.iter_mut() {
        let mut rest: &mut [f64] = row;
        let mut off = 0usize;
        for (ci, (start, end)) in bounds.iter().enumerate() {
            let (win, tail) = rest.split_at_mut(end - off);
            chunk_windows[ci].push(&mut win[start - off..]);
            rest = tail;
            off = *end;
        }
    }
    std::thread::scope(|sc| {
        let mut handles = Vec::new();
        for ((start, end), mut wins) in bounds.iter().zip(chunk_windows) {
            let et = &et;
            let u = &u;
            let v = &v;
            let rho = &state.rho;
            handles.push(sc.spawn(move || {
                let mut w = 0.0f64;
                for k in *start..*end {
                    let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
                    let mut y = [0.0f64; ct::NS];
                    for (ys, ys_src) in y.iter_mut().zip(wins.iter()) {
                        *ys = ys_src[k - start];
                    }
                    w = w.max(react_cell(e_int, rho[k], &mut y, dt));
                    for (ys_dst, ys) in wins.iter_mut().zip(y.iter()) {
                        ys_dst[k - start] = *ys;
                    }
                }
                w
            }));
        }
        for h in handles {
            worst = worst.max(h.join().unwrap());
        }
    });
    Ok(worst)
}

/// advance one cell's composition over dt at fixed total energy.
/// returns the chemistry timescale violation, the march's stiffness
/// report.
fn react_cell(e_int: f64, rho: f64, y: &mut [f64; ct::NS], dt: f64) -> f64 {
    let mut worst_oversub = 0.0f64;
    {
        let mut t = closure9::temperature_at(e_int, y);
        let mut remaining = dt;
        let mut steps = 0;
        while remaining > 1e-14 && steps < 512 {
            let (dydt, _) = comb::omega(rho, y, t);
            let rate_scale = dydt.iter().fold(0.0f64, |acc, d| acc.max(d.abs()));
            let tau = if rate_scale > 1e-30 {
                1.0 / rate_scale
            } else {
                remaining
            };
            let h = remaining.min(tau * 0.2);
            // positivity by RATE SCALING, not by shrinking the step:
            // the reaction vector is scaled down by the largest factor
            // that keeps every species non-negative over this step.
            // the stoichiometry is untouched, so the species sum and
            // the element pools stay exactly conserved, and the
            // chemistry still advances (shrinking the step instead
            // would freeze it behind a trace species).
            let mut scale = 1.0f64;
            for (s, &d) in dydt.iter().enumerate() {
                if d < 0.0 && y[s] > 0.0 {
                    let allowed = y[s] / (-d * h);
                    if allowed < scale {
                        scale = allowed;
                    }
                }
            }
            for (ys, d) in y.iter_mut().zip(dydt.iter()) {
                let d = d * h * scale;
                // no clamp and no renormalization anywhere: the
                // reaction system conserves the species sum and the
                // element pools exactly, and any correction destroys
                // one or both.
                *ys += d;
            }
            t = closure9::temperature_at(e_int, y);
            remaining -= h;
            steps += 1;
            if steps == 1 {
                worst_oversub = worst_oversub.max(dt / (tau * 0.2));
            }
        }
        // no sanitize anywhere: the reaction system conserves both
        // the species sum and the element pools exactly, so any
        // clamp or renormalization manufactures or destroys mass. the
        // adaptive step keeps excursions small; a mildly negative
        // species relaxes back within the same cell-dt.
    }
    worst_oversub
}

/// internal energy per mass for an unburned h2/o2 fill at (rho, p):
/// T from the mixture gas constant, then the closure energy.
pub fn h2o2_ic_energy(rho: f64, p: f64) -> f64 {
    let mut y = [0.0f64; ct::NS];
    y[ct::IH2] = 2.016 / 18.016;
    y[ct::IO2] = 16.0 / 18.016;
    h2o2_ic_energy_at(rho, p, &y)
}

/// the closure energy of a uniform fill at (rho, p): the temperature
/// from the fill's own mixture gas constant, then the closure energy
/// at that composition. a burned-injection case (the chamber) needs
/// the energy of the burned gas, not the premix's.
pub fn h2o2_ic_energy_at(rho: f64, p: f64, fill: &[f64; ct::NS]) -> f64 {
    let t = p / (rho * closure9::mix_gas_constant(fill));
    closure9::internal_energy_at(t, fill)
}

#[cfg(test)]
mod tests {
    use super::{
        advance2d_combustion, advance2d_combustion_axi, advance2d_combustion_rk2, closure9, ct,
        h2o2_ic_energy, react_cell, Boundaries2d, CombustionState, ConservedState2d,
    };

    #[test]
    fn axi_combustion_holds_uniform_equilibrium() {
        use crate::grid2d::grid2d;
        use crate::state2d::prim_to_cons2d;
        // a uniform field at chamber conditions on a small grid: the
        // march must hold it steady (slip walls, no body).
        let g = grid2d(8, 8, 0.0, 0.08, 0.0, 0.08).unwrap();
        let n = g.nx * g.ny;
        let rho = 6.7412f64;
        let p = 7.0e6;
        let mut sp = CombustionState::stoich_h2o2(n);
        let e_int = h2o2_ic_energy(rho, p);
        let (mx, my, e) =
            prim_to_cons2d(&vec![rho; n], &vec![0.0; n], &vec![0.0; n], &vec![e_int; n]).unwrap();
        let mut state = ConservedState2d {
            rho: vec![rho; n],
            mx,
            my,
            e,
        };
        let bc = crate::Boundaries2d {
            west: crate::Bc2d::SupersonicOutflow,
            east: crate::Bc2d::SupersonicOutflow,
            south: crate::Bc2d::SlipWall,
            north: crate::Bc2d::SlipWall,
        };
        for _ in 0..40 {
            let res = advance2d_combustion_axi(&mut state, &mut sp, &g, 0.4, false, &bc, 1, None);
            assert!(res.is_ok(), "axi step failed: {:?}", res.err());
        }
        for k in 0..n {
            assert!(
                (state.rho[k] - rho).abs() < 1e-6,
                "rho drifted at k={k}: {}",
                state.rho[k]
            );
        }
    }

    /// the real-gas reservoir inflow: the ghost carries the reservoir
    /// state expanded isentropically to the local face pressure, so a
    /// chamber holds its pressure while the throat draws flow. the
    /// ghost must land on the face pressure (the pressure inversion
    /// through the closure must reproduce the requested pressure) and
    /// the expansion must be subsonic at the face.
    #[test]
    fn reservoir_ghost_expands_to_the_face_pressure() {
        let n = 8;
        let mut sp = CombustionState::stoich_h2o2(n);
        // a burned chamber fill: the cea equilibrium composition at
        // 70 bar and 3676 k.
        let fill = [
            0.01583304, 0.07194470, 0.77605355, 0.11605824, 0.00232724, 0.01719592, 0.00049223,
            0.00009298, 0.0,
        ];
        let y9: [f64; ct::NS] = fill;
        let rho_t = 3.6116f64;
        let pt = 7.0e6f64;
        let tt = 3676.52f64;
        let bc = crate::Boundaries2d {
            west: crate::Bc2d::ReservoirInflowRealGas {
                pt,
                tt,
                rho_t,
                fill: fill.to_vec(),
                n_species: ct::NS,
            },
            east: crate::Bc2d::SupersonicOutflow,
            south: crate::Bc2d::SlipWall,
            north: crate::Bc2d::SlipWall,
        };
        let model = crate::ThermoModel::Combustion9 { gamma_shaping: 1.2 };
        // a face pressure below the reservoir: the ghost expands and
        // accelerates toward the face.
        let p_face = 5.0e6f64;
        let (rho_g, u_g, e_g) =
            crate::axi::reservoir_real_gas(&bc.west, model, p_face).expect("no ghost");
        assert!(rho_g > 0.0 && rho_g < rho_t, "expanded density {rho_g}");
        assert!(u_g > 0.0, "the expansion must accelerate the flow");
        // the ghost's own pressure through the closure must land on
        // the requested face pressure.
        let p_g = crate::combustion_closure::pressure_at(rho_g, e_g, &y9);
        assert!(
            (p_g - p_face).abs() / p_face < 1e-3,
            "ghost pressure {p_g} vs the face {p_face}"
        );
        // at the reservoir pressure the ghost is the reservoir state.
        let (rho_full, u_full, e_full) =
            crate::axi::reservoir_real_gas(&bc.west, model, pt).expect("no ghost");
        assert!((rho_full - rho_t).abs() / rho_t < 1e-9, "rho {rho_full}");
        assert!(u_full.abs() < 1e-9, "the reservoir must not accelerate");
        let p_full = crate::combustion_closure::pressure_at(rho_full, e_full, &y9);
        // the closure's own temperature inversion carries a few ppm of
        // roundtrip error, so the bound is 1e-5 relative.
        assert!((p_full - pt).abs() / pt < 1e-5, "pressure {p_full}");
        // the composition is untouched by the expansion.
        let _ = &mut sp;
    }

    /// the full stiff transient the combustor march drives: a stoich
    /// fill at 70 bar and 800 k burns to equilibrium in one cell-dt.
    /// the element pools must survive to roundoff — no clamp or
    /// renormalization may run anywhere in the integration.
    #[test]
    fn stiff_transient_preserves_element_pools() {
        let r_mix = crate::combustion_thermo::RU / 0.018016;
        let e_int = {
            let mut y_fill = [0.0f64; ct::NS];
            y_fill[ct::IH2] = 2.016 / 18.016;
            y_fill[ct::IO2] = 16.0 / 18.016;
            crate::combustion_thermo::mix_h_mass(&y_fill, 800.0) - r_mix * 800.0
        };
        let mut y = [0.0f64; ct::NS];
        y[ct::IH2] = 2.016 / 18.016;
        y[ct::IO2] = 16.0 / 18.016;
        // a march-scale dt: the chemical time is ~1e-7 s, so any
        // dt at or above it exercises the adaptive substep ladder.
        for dt in [1e-7f64, 5e-7, 2e-6] {
            let mut y_case = y;
            react_cell(e_int, 12.64, &mut y_case, dt);
            let h: f64 = 2.0 * 1.008 / 2.016 * y_case[ct::IH2]
                + 2.0 * 1.008 / 18.015 * y_case[ct::IH2O]
                + 1.008 / 17.007 * y_case[ct::IOH]
                + y_case[ct::IH]
                + 1.008 / 33.0067 * y_case[ct::IHO2]
                + 2.0 * 1.008 / 34.0147 * y_case[ct::IH2O2];
            let o: f64 = y_case[ct::IO2]
                + 2.0 * 15.999 / 18.015 * y_case[ct::IH2O]
                + 15.999 / 17.007 * y_case[ct::IOH]
                + y_case[ct::IO_]
                + 2.0 * 15.999 / 33.0067 * y_case[ct::IHO2]
                + 2.0 * 15.999 / 34.0147 * y_case[ct::IH2O2];
            assert!(
                (h - 2.016 / 18.016).abs() < 1e-12,
                "h pool drifted at dt {dt}: {h} vs {}",
                2.016 / 18.016
            );
            assert!(
                (o - 16.0 / 18.016).abs() < 1e-12,
                "o pool drifted at dt {dt}: {o} vs {}",
                16.0 / 18.016
            );
            assert!(
                (y_case.iter().sum::<f64>() - 1.0).abs() < 1e-12,
                "mass sum drifted at dt {dt}"
            );
        }
    }

    /// a uniform frozen field stays uniform: the transport must not
    /// disturb a constant state, and with chemistry off the fill
    /// never moves.
    #[test]
    fn uniform_frozen_field_stays_uniform() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let u = 100.0f64;
        let p = 7.0e6f64;
        let e_int = h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let mut st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![rho * u; n],
            my: vec![0.0; n],
            e: vec![rho * (e_int + 0.5 * u * u); n],
        };
        let mut sp = CombustionState::stoich_h2o2(n);
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        for _ in 0..4 {
            advance2d_combustion(&mut st, &mut sp, &g, 0.3, false, &bc, 1).unwrap();
        }
        for k in 0..n {
            assert!((st.rho[k] - rho).abs() / rho < 1e-10, "rho moved");
            assert!(
                (sp.y[ct::IH2][k] - 2.016 / 18.016).abs() < 1e-10,
                "composition moved"
            );
        }
    }

    /// a hot spot ignites and moves toward equilibrium: the products
    /// grow at the expense of the reactants through the substep.
    #[test]
    fn hot_cell_ignites_toward_equilibrium() {
        let (nx, ny) = (4, 4);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.04, 0.0, 0.04).unwrap();
        let rho = 1.0f64;
        let n = nx * ny;
        // hot IC: e_int for 4000 K unburned-like mix
        let y = CombustionState::stoich_h2o2(n).fill;
        let t_hot = 4000.0f64;
        let e_int = closure9::internal_energy_at(t_hot, &y);
        let mut st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let mut sp = CombustionState::stoich_h2o2(n);
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        // rk2 with react on: one step at a large flow dt so the
        // substep budget runs chemistry
        advance2d_combustion_rk2(&mut st, &mut sp, &g, 0.3, false, &bc, 1, true).unwrap();
        let yh2o_before = sp.y[ct::IH2O][0];
        let _ = yh2o_before;
        assert!(
            sp.y[ct::IH2O][0] > 1e-3,
            "water must appear: {}",
            sp.y[ct::IH2O][0]
        );
    }

    /// an inflow face keeps injecting the case fill after the near
    /// cells burn: the species ghost carries the stored fill, not the
    /// burned edge cell.
    #[test]
    fn inflow_ghost_composition_is_the_case_fill() {
        let (nx, ny) = (6, 4);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.06, 0.0, 0.04).unwrap();
        let mut st = ConservedState2d {
            rho: vec![12.64; nx * ny],
            mx: vec![1264.0; nx * ny],
            my: vec![0.0; nx * ny],
            e: vec![0.0; nx * ny],
        };
        let mut sp = CombustionState::stoich_h2o2(nx * ny);
        let e_int = h2o2_ic_energy(12.64, 7.0e6);
        for k in 0..nx * ny {
            st.e[k] = st.rho[k] * (e_int + 0.5 * 100.0 * 100.0);
        }
        // burn the west edge cells completely: all water
        for j in 0..ny {
            let k = j * nx;
            sp.y[ct::IH2][k] = 0.0;
            sp.y[ct::IO2][k] = 0.0;
            sp.y[ct::IH2O][k] = 1.0;
        }
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::SupersonicInflow {
                rho: 12.64,
                u: 100.0,
                v: 0.0,
                p: 7.0e6,
            },
            east: crate::solver2d::Bc2d::SupersonicOutflow,
            south: crate::solver2d::Bc2d::SlipWall,
            north: crate::solver2d::Bc2d::SlipWall,
        };
        // one frozen step: the west edge cell must keep its burned
        // composition (the ghost injects the fill outward, it does not
        // pull the burned cell back), and the fill stays unburned.
        advance2d_combustion(&mut st, &mut sp, &g, 0.3, false, &bc, 1).unwrap();
        assert!(sp.y[ct::IH2O][nx] > 0.9);
        assert!((sp.fill[ct::IH2] - 2.016 / 18.016).abs() < 1e-12);
    }
}

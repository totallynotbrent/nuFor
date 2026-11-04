//! finite-rate GRI-Mech source terms applied to a spatial composition.
//!
//! the transport seam (`chem_spatial`) moves a scheme-width composition
//! through a flow. this module closes the other half of spatial
//! chemistry: the reaction source. each cell is a closed parcel for the
//! duration of the source substep, so the chemistry each cell runs is
//! exactly the zero-dimensional reactor's — the same backward-euler
//! step on the fifty-three mass fractions with the temperature following
//! from the cell's fixed internal energy.
//!
//! the caller supplies each cell's density and internal energy per mass.
//! the composition is updated in place, clipping only rounding-negative
//! fractions and renormalizing, the same hygiene the parcel march uses.
//!
//! no reaction subset is selected here. every cell runs the full
//! mechanism. adaptive selection is a separate question and belongs to
//! its own seam.

use crate::chem_spatial::SpatialComposition;
use crate::error::Error;
use crate::gri30_reactor::{backward_euler_step, NS};

/// what one source pass did.
#[derive(Clone, Copy, Debug, Default)]
pub struct ReactionReport {
    /// the cells whose chemistry stepped successfully.
    pub cells_reacted: usize,
    /// the cells whose solve did not converge and were left unchanged.
    pub cells_deferred: usize,
    /// the largest number of newton iterations any cell needed.
    pub worst_iterations: u32,
    /// the worst per-cell composition-sum error after the pass.
    pub max_sum_error: f64,
}

/// apply one finite-rate GRI-Mech source step to every cell.
///
/// `density` and `internal_energy` are per cell, in kg/m3 and J/kg. the
/// composition must carry the mechanism's fifty-three species. a cell
/// whose implicit solve does not converge is left at its incoming
/// composition and counted in the report, so a stiff cell pauses its
/// chemistry rather than corrupting the field.
pub fn apply_gri30_source(
    composition: &mut SpatialComposition,
    density: &[f64],
    internal_energy: &[f64],
    dt: f64,
) -> Result<ReactionReport, Error> {
    let ncells = composition.ncells();
    if composition.n_species() != NS
        || density.len() != ncells
        || internal_energy.len() != ncells
        || dt.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
    {
        return Err(Error::InvalidArgs);
    }
    let mut report = ReactionReport::default();
    for cell in 0..ncells {
        let mut parcel = [0.0f64; NS];
        for (species, row) in composition.y.iter().enumerate() {
            parcel[species] = row[cell];
        }
        let rho = density[cell];
        if rho.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) || !rho.is_finite() {
            report.cells_deferred += 1;
            continue;
        }
        // the parcel's energy is reconstructed from its own composition
        // and the caller-supplied internal energy, so the chemistry sees
        // the same constraint the reactor marches under.
        let energy = internal_energy[cell];
        if !energy.is_finite() {
            report.cells_deferred += 1;
            continue;
        }
        let mut iterations = 0u32;
        match backward_euler_step(&parcel, rho, energy, dt, &mut iterations) {
            Ok(Some(stepped)) => {
                for (species, value) in stepped.iter().enumerate() {
                    composition.y[species][cell] = if *value < 0.0 { 0.0 } else { *value };
                }
                report.worst_iterations = report.worst_iterations.max(iterations);
                report.cells_reacted += 1;
            }
            _ => {
                report.cells_deferred += 1;
            }
        }
        // the step preserves the sum analytically; renormalizing clears
        // the solver's rounding dust the way the march does.
        let _ = composition.renormalize(cell);
    }
    report.max_sum_error = composition.max_sum_error();
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chem_scheme::Gri30;
    use crate::gri30;
    use crate::gri30_reactor;
    use crate::gri30_reactor::internal_energy_mass;

    /// a stoichiometric methane-air parcel at a temperature that is hot
    /// enough to react but gentle enough for one implicit step.
    fn premix(fraction_ch4: f64) -> ([f64; NS], f64, f64) {
        let air = gri30_reactor::stoich_air("CH4", 1.0).expect("stoichiometric air");
        let mut y = [0.0f64; NS];
        let ch4 = gri30::species_index("CH4").unwrap();
        // dilute the stoichiometric mixture toward air so the cell does
        // not explode in a single step.
        for s in 0..NS {
            y[s] = fraction_ch4 * {
                // stoich_air puts the fuel in the CH4 slot already.
                if s == ch4 {
                    1.0
                } else {
                    0.0
                }
            } + (1.0 - fraction_ch4) * air[s];
        }
        let temperature = 1_200.0f64;
        let rho = 1.0f64;
        (y, rho, internal_energy_mass(&y, temperature))
    }

    /// one cell advanced by the spatial pass must equal the reactor's
    /// own step on the same parcel, to rounding.
    #[test]
    fn one_spatial_cell_matches_the_reactor_step() {
        let (y, rho, energy) = premix(0.055);
        let fractions = y;
        let mut composition = SpatialComposition::uniform(&Gri30, 1, &fractions).unwrap();
        let dt = 1.0e-6;
        let report = apply_gri30_source(&mut composition, &[rho], &[energy], dt).unwrap();
        assert_eq!(report.cells_reacted, 1);
        assert_eq!(report.cells_deferred, 0);

        let mut iterations = 0u32;
        let expected = backward_euler_step(&y, rho, energy, dt, &mut iterations)
            .unwrap()
            .unwrap();
        for (s, &want_raw) in expected.iter().enumerate() {
            let got = composition.y[s][0];
            let want = want_raw.max(0.0);
            assert!(
                (got - want).abs() < 1.0e-9,
                "species {s}: spatial {got:e} against reactor {want:e}"
            );
        }
    }

    /// a reacting parcel burns down its fuel under the source step, and
    /// the sum stays a partition of unity.
    #[test]
    fn a_reacting_cell_burns_its_fuel_and_stays_normalized() {
        let air = gri30_reactor::stoich_air("CH4", 1.0).expect("stoichiometric air");
        let ch4 = gri30::species_index("CH4").unwrap();
        let fraction = 0.055f64;
        let mut y = [0.0f64; NS];
        for s in 0..NS {
            y[s] = if s == ch4 {
                fraction
            } else {
                (1.0 - fraction) * air[s]
            };
        }
        let rho = 1.0f64;
        // normalize so the composition constructor accepts the mixture.
        let total: f64 = y.iter().sum();
        for value in y.iter_mut() {
            *value /= total;
        }
        // hot enough to ignite well inside twenty steps of ten microseconds.
        let temperature = 1_800.0f64;
        let energy = internal_energy_mass(&y, temperature);
        let mut composition = SpatialComposition::uniform(&Gri30, 1, &y).unwrap();
        let mut dt = 1.0e-4f64;
        let mut reacted_steps = 0usize;
        // subdivide the way the parcel march subdivides: a step whose
        // solve does not converge is retried at half the size.
        for _ in 0..40 {
            let report = apply_gri30_source(&mut composition, &[rho], &[energy], dt).unwrap();
            if report.cells_reacted == 1 {
                reacted_steps += 1;
                dt *= 1.5;
            } else {
                dt *= 0.5;
            }
        }
        assert!(reacted_steps > 0, "the cell's solve never converged");
        let mut parcel = [0.0f64; NS];
        for (s, row) in composition.y.iter().enumerate() {
            parcel[s] = row[0];
        }
        let fuel_left = parcel[ch4];
        let water_formed =
            parcel[gri30::species_index("H2O").unwrap()] - y[gri30::species_index("H2O").unwrap()];
        assert!(
            fuel_left < y[ch4],
            "the fuel did not move at all: {fuel_left:e}"
        );
        assert!(
            water_formed > 0.0,
            "burning makes water: {water_formed:e} kg/kg of H2O formed"
        );
        assert!(composition.max_sum_error() < 1.0e-12);
        for row in &composition.y {
            assert!(row.iter().all(|v| *v >= 0.0 && v.is_finite()));
        }
    }

    /// a cell of pure nitrogen has no chemistry: the source step leaves
    /// it bit for bit unchanged.
    #[test]
    fn an_inert_cell_is_untouched() {
        let mut nitrogen = [0.0f64; NS];
        nitrogen[gri30::species_index("N2").unwrap()] = 1.0;
        let mut composition = SpatialComposition::uniform(&Gri30, 2, &nitrogen).unwrap();
        let rho = 1.2f64;
        let e = internal_energy_mass(&nitrogen, 300.0);
        let report = apply_gri30_source(&mut composition, &[rho, rho], &[e, e], 1.0e-5).unwrap();
        assert_eq!(report.cells_reacted, 2);
        for (s, &initial) in nitrogen.iter().enumerate() {
            assert!((composition.y[s][0] - initial).abs() < 1.0e-15);
        }
    }
}

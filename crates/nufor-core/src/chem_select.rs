//! per-cell reaction-set selection against the dense mechanism.
//!
//! every cell currently runs the full 325-reaction GRI-Mech. that is
//! correct but expensive, and most of a cell's species are trace
//! passengers whose evolution is numerically irrelevant over one
//! step. this module picks, per cell, which SPECIES a reduced step
//! evolves: rank by production activity, keep the loudest above a
//! threshold, freeze the rest for the step.
//!
//! freezing is a property of the STATE VECTOR, not of the kinetics.
//! the rate evaluation still sees all fifty-three species, so
//! third bodies and colliders keep their bath and the mechanism keeps
//! its reversibility. only the evolved set shrinks. that separation is
//! what makes a reduced step defensible against the dense reference
//! instead of a silently different chemistry.
//!
//! selection is recomputed per call from the cell's current state, so a
//! cell that quenches sheds its radical set on the next pass.

use crate::chem_spatial::SpatialComposition;
use crate::error::Error;
use crate::gri30_reactor::{backward_euler_step_selected, mass_fraction_rates, NS};

/// how a cell's reaction set is chosen.
#[derive(Clone, Copy, Debug)]
pub struct SelectionPolicy {
    /// a species joins the evolved set when its absolute mass-fraction
    /// rate exceeds this fraction of the largest rate in the cell.
    pub relative_activity: f64,
    /// a floor under which a species is considered dormant regardless
    /// of the ranking, in kg/kg/s.
    pub activity_floor: f64,
    /// always evolve at least this many of the most active species, so
    /// a quiescent cell still relaxes its fastest pair.
    pub minimum_active: usize,
    /// never evolve fewer than the major carriers: fuel, oxidiser, and
    /// the diluent keep their degrees of freedom.
    pub protected: &'static [&'static str],
}

impl Default for SelectionPolicy {
    fn default() -> Self {
        Self {
            relative_activity: 1.0e-3,
            activity_floor: 1.0e-9,
            minimum_active: 4,
            protected: &["CH4", "O2", "N2", "CO2", "H2O"],
        }
    }
}

/// what one selection pass chose.
#[derive(Clone, Copy, Debug, Default)]
pub struct SelectionReport {
    /// the average number of species each cell evolved.
    pub mean_active_species: f64,
    /// the smallest evolved set any cell got.
    pub min_active_species: usize,
    /// the largest evolved set any cell got.
    pub max_active_species: usize,
    /// the cells whose step did not converge and stayed unchanged.
    pub cells_deferred: usize,
}

/// choose the evolved species for one cell from its current state.
///
/// ranks by `|dy/dt|` and keeps everything above the policy's relative
/// threshold plus floor, always including the protected majors and at
/// least `minimum_active` species.
pub fn select_mask(
    parcel: &[f64; NS],
    rho: f64,
    energy: f64,
    policy: &SelectionPolicy,
) -> [bool; NS] {
    let (dy, _) = mass_fraction_rates(parcel, rho, energy);
    let loudest = dy.iter().fold(0.0f64, |acc, rate| acc.max(rate.abs()));
    let cutoff = policy
        .activity_floor
        .max(policy.relative_activity * loudest);
    let mut ranked: Vec<(usize, f64)> = dy
        .iter()
        .enumerate()
        .map(|(s, rate)| (s, rate.abs()))
        .collect();
    ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut mask = [false; NS];
    for &(s, activity) in &ranked {
        if activity >= cutoff {
            mask[s] = true;
        }
    }
    // the protected majors always evolve: they carry the mass and the
    // heat release, and freezing them would decouple the step from the
    // parcel's energy.
    for name in policy.protected {
        if let Some(slot) = crate::gri30::species_index(name) {
            mask[slot] = true;
        }
    }
    // honour the floor on the count, taking the loudest extras if the
    // threshold marked too few.
    let active_count = mask.iter().filter(|&&m| m).count();
    if active_count < policy.minimum_active {
        for &(s, _) in ranked.iter().take(policy.minimum_active) {
            mask[s] = true;
        }
    }
    mask
}

/// apply one implicitly integrated source step per cell, each with its
/// own selected reaction set.
///
/// like [`crate::chem_react::apply_gri30_source`] but the evolved
/// species per cell come from [`select_mask`] at the cell's incoming
/// state. a cell whose solve does not converge is left unchanged and
/// counted.
pub fn apply_selected_source(
    composition: &mut SpatialComposition,
    density: &[f64],
    internal_energy: &[f64],
    dt: f64,
    policy: &SelectionPolicy,
) -> Result<SelectionReport, Error> {
    let ncells = composition.ncells();
    if composition.n_species() != NS
        || density.len() != ncells
        || internal_energy.len() != ncells
        || dt.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
    {
        return Err(Error::InvalidArgs);
    }
    let mut report = SelectionReport::default();
    let mut total_active = 0usize;
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
        let energy = internal_energy[cell];
        if !energy.is_finite() {
            report.cells_deferred += 1;
            continue;
        }
        let mask = select_mask(&parcel, rho, energy, policy);
        let active = mask.iter().filter(|&&m| m).count();
        total_active += active;
        report.min_active_species = if report.min_active_species == 0 {
            active
        } else {
            report.min_active_species.min(active)
        };
        report.max_active_species = report.max_active_species.max(active);
        let mut iterations = 0u32;
        match backward_euler_step_selected(&parcel, rho, energy, dt, &Some(mask), &mut iterations) {
            Ok(Some(stepped)) => {
                for (species, value) in stepped.iter().enumerate() {
                    composition.y[species][cell] = if *value < 0.0 { 0.0 } else { *value };
                }
            }
            _ => {
                report.cells_deferred += 1;
            }
        }
        let _ = composition.renormalize(cell);
    }
    report.mean_active_species = if ncells == 0 {
        0.0
    } else {
        total_active as f64 / ncells as f64
    };
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chem_scheme::Gri30;
    use crate::gri30;
    use crate::gri30_reactor::{internal_energy_mass, temperature_at_energy};

    fn ch4_air(dilution: f64, temperature: f64) -> ([f64; NS], f64, f64) {
        let air = crate::gri30_reactor::stoich_air("CH4", 1.0).unwrap();
        let ch4 = gri30::species_index("CH4").unwrap();
        let mut y = [0.0f64; NS];
        for s in 0..NS {
            y[s] = if s == ch4 {
                dilution
            } else {
                (1.0 - dilution) * air[s]
            };
        }
        let total: f64 = y.iter().sum();
        for value in y.iter_mut() {
            *value /= total;
        }
        let rho = 1.0f64;
        let energy = internal_energy_mass(&y, temperature);
        (y, rho, energy)
    }

    /// a masked step must reproduce the unmasked step when every species
    /// is selected: the selection parameter must reduce to the plain
    /// step, proving the mask touches only WHICH DOFs move.
    #[test]
    fn a_full_mask_matches_the_plain_step() {
        let (y, rho, energy) = ch4_air(0.055, 1200.0);
        let dt = 1.0e-6;
        let mut it_plain = 0u32;
        let mut it_masked = 0u32;
        let plain = backward_euler_step_selected(&y, rho, energy, dt, &None, &mut it_plain)
            .unwrap()
            .unwrap();
        let full = [true; NS];
        let masked = backward_euler_step_selected(&y, rho, energy, dt, &Some(full), &mut it_masked)
            .unwrap()
            .unwrap();
        for (s, (&a, &b)) in plain.iter().zip(masked.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1.0e-12,
                "species {s}: plain {a:e} masked {b:e}"
            );
        }
    }

    /// a masked step must hold the frozen species EXACTLY at their
    /// incoming values, while the selected ones move.
    #[test]
    fn a_mask_freezes_the_unselected_species() {
        let (y, rho, energy) = ch4_air(0.055, 1200.0);
        let dt = 1.0e-6;
        // evolve only the majors: freeze everything else.
        let mut mask = [false; NS];
        for name in ["CH4", "O2", "N2", "CO2", "H2O"] {
            mask[gri30::species_index(name).unwrap()] = true;
        }
        let mut it = 0u32;
        let stepped = backward_euler_step_selected(&y, rho, energy, dt, &Some(mask), &mut it)
            .unwrap()
            .unwrap();
        for s in 0..NS {
            if !mask[s] {
                assert!(
                    (stepped[s] - y[s]).abs() < 1.0e-18,
                    "frozen species {s} moved: {} to {}",
                    y[s],
                    stepped[s]
                );
            }
        }
        // the fuel must have MOVED, proving the step did something.
        let ch4 = gri30::species_index("CH4").unwrap();
        assert!(stepped[ch4] != y[ch4]);
    }

    /// selection on a hot reacting parcel keeps far fewer than the full
    /// fifty-three species, and always protects the majors.
    #[test]
    fn selection_shrinks_the_set_and_protects_majors() {
        let (y, rho, energy) = ch4_air(0.055, 1350.0);
        let policy = SelectionPolicy::default();
        let mask = select_mask(&y, rho, energy, &policy);
        let active = mask.iter().filter(|&&m| m).count();
        assert!(active < NS, "nothing was screened out: {active}");
        assert!(active >= policy.minimum_active);
        for name in policy.protected {
            let slot = gri30::species_index(name).unwrap();
            assert!(mask[slot], "protected species {name} was frozen");
        }
    }

    /// the selected pass on a field defers nothing on mild parcels and
    /// shrinks the mean evolved set below the full mechanism.
    #[test]
    fn the_selected_pass_reports_a_smaller_set() {
        let (y, _, _) = ch4_air(0.055, 1200.0);
        let (_, rho, energy) = ch4_air(0.055, 1200.0);
        let mut composition = SpatialComposition::uniform(&Gri30, 3, &y).unwrap();
        let report = apply_selected_source(
            &mut composition,
            &[rho, rho, rho],
            &[energy, energy, energy],
            1.0e-6,
            &SelectionPolicy::default(),
        )
        .unwrap();
        assert_eq!(report.cells_deferred, 0);
        assert!(report.mean_active_species < NS as f64);
        assert!(report.min_active_species <= report.max_active_species);
        // sanity: temperature is still recoverable from the result.
        let mut parcel = [0.0f64; NS];
        for (s, row) in composition.y.iter().enumerate() {
            parcel[s] = row[0];
        }
        let t = temperature_at_energy(&parcel, energy);
        assert!(t.is_finite() && t > 0.0);
    }
}

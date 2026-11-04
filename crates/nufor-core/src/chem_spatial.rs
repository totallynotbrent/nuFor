//! scheme-width composition state for spatial chemistry.
//!
//! the fixed nine-species combustion lane predates the GRI-Mech
//! mechanism. this module supplies the small, scheme-independent piece
//! that a wider lane needs first: a conserved composition whose width
//! follows its chemistry scheme and a transport pass that treats every
//! species row equally.
//!
//! it deliberately does not choose an active reaction subset. dropping
//! a reaction participant requires preserving the other branch of a
//! reversible reaction, preserving a broadened falloff reaction, and
//! accounting for third-body colliders. the conservative composition
//! and its movement do not depend on any of those decisions.

use crate::chem_closure;
use crate::chem_scheme::ChemScheme;
use crate::error::Error;
use crate::species_transport::transport_species_generic;

/// a spatial composition with one mass-fraction row per species.
#[derive(Clone, Debug)]
pub struct SpatialComposition {
    /// the species mass fractions, `y[species][cell]`.
    pub y: Vec<Vec<f64>>,
}

impl SpatialComposition {
    /// build a uniform composition from a scheme-sized mass-fraction row.
    ///
    /// the row must have one entry per scheme species and its entries
    /// must sum to one within `1.0e-12`. a mismatch is rejected rather
    /// than normalized silently because normalization hides an indexing
    /// mistake between a mechanism and its case file.
    pub fn uniform<S: ChemScheme>(
        scheme: &S,
        ncells: usize,
        fractions: &[f64],
    ) -> Result<Self, Error> {
        if fractions.len() != scheme.n_species() {
            return Err(Error::InvalidArgs);
        }
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > 1.0e-12 {
            return Err(Error::InvalidArgs);
        }
        Ok(Self {
            y: fractions
                .iter()
                .map(|fraction| vec![*fraction; ncells])
                .collect(),
        })
    }

    /// the number of species carried by the composition.
    pub fn n_species(&self) -> usize {
        self.y.len()
    }

    /// the number of cells carried by the composition.
    pub fn ncells(&self) -> usize {
        self.y.first().map_or(0, |row| row.len())
    }

    /// rescale one cell's mass fractions to sum to one.
    ///
    /// the divisor has to be positive. a nonpositive cell has lost the
    /// meaning of a composition and renormalizing it would fabricate
    /// species mass.
    pub fn renormalize(&mut self, cell: usize) -> Result<(), Error> {
        let mut sum = 0.0;
        for row in &self.y {
            sum += row[cell];
        }
        if !sum.is_finite() || sum <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        for row in &mut self.y {
            row[cell] /= sum;
        }
        Ok(())
    }

    /// the worst composition-sum error over all cells.
    #[allow(clippy::needless_range_loop)]
    pub fn max_sum_error(&self) -> f64 {
        let mut worst = 0.0f64;
        for cell in 0..self.ncells() {
            let sum: f64 = self.y.iter().map(|row| row[cell]).sum();
            worst = worst.max((sum - 1.0).abs());
        }
        worst
    }
    /// the mixture state in one cell.
    pub fn mixture_state<S: ChemScheme>(
        &self,
        scheme: &S,
        cell: usize,
        density: f64,
        internal_energy: f64,
    ) -> Result<chem_closure::MixtureState, Error> {
        if self.y.len() != scheme.n_species() {
            return Err(Error::InvalidArgs);
        }
        let fractions: Vec<f64> = self.y.iter().map(|row| row[cell]).collect();
        Ok(chem_closure::mixture_state(
            scheme,
            density,
            internal_energy,
            &fractions,
        ))
    }
}

/// advance a composition with the shared conservative species transport.
///
/// the face mass fluxes use the same layout as the planar and annular
/// flow lanes: `fx[row][face]` and `fy[column][face]`. the caller's
/// radius weights are optional, so the planar and axisymmetric paths
/// share this function.
#[allow(clippy::too_many_arguments)]
pub fn advance_composition(
    composition: &mut SpatialComposition,
    density: &[f64],
    dt: f64,
    dxs: &[f64],
    dys: &[f64],
    fx: &[Vec<f64>],
    fy: &[Vec<f64>],
    radial_weights: Option<&[f64]>,
) -> Result<f64, Error> {
    let ncells = dxs.len() * dys.len();
    if composition.ncells() != ncells
        || density.len() != ncells
        || fx.len() != dys.len()
        || fy.len() != dxs.len()
    {
        return Err(Error::InvalidArgs);
    }
    let residual = transport_species_generic(
        &mut composition.y,
        density,
        dt,
        dxs,
        dys,
        fx,
        fy,
        radial_weights,
    )?;
    for cell in 0..ncells {
        composition.renormalize(cell)?;
    }
    Ok(residual)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chem_scheme::{Gri30, Li2004};

    /// the wider mechanism can carry all of its species into the spatial
    /// composition without borrowing a fixed-width buffer.
    #[test]
    fn gri30_fifty_three_species_form_a_valid_composition() {
        let mut fractions = vec![0.0; Gri30.n_species()];
        fractions[Gri30.n_species() - 1] = 1.0;
        let composition = SpatialComposition::uniform(&Gri30, 12, &fractions).unwrap();
        assert_eq!(composition.n_species(), 53);
        assert_eq!(composition.ncells(), 12);
        assert!(composition.max_sum_error() < 1.0e-15);
    }

    /// a mismatched composition row is rejected rather than resized.
    #[test]
    fn a_species_count_mismatch_is_rejected() {
        let fractions = vec![1.0; Gri30.n_species() - 1];
        assert!(SpatialComposition::uniform(&Gri30, 1, &fractions).is_err());
    }

    /// a generic closure reads the wide composition without a
    /// nine-species translation.
    #[test]
    fn the_generic_closure_accepts_the_wide_composition() {
        let mut fractions = vec![0.0; Gri30.n_species()];
        let methane = crate::gri30::species_index("CH4").unwrap();
        let oxygen = crate::gri30::species_index("O2").unwrap();
        let nitrogen = crate::gri30::species_index("N2").unwrap();
        fractions[methane] = 0.1;
        fractions[oxygen] = 0.2;
        fractions[nitrogen] = 0.7;
        let composition = SpatialComposition::uniform(&Gri30, 1, &fractions).unwrap();
        let energy = chem_closure::internal_energy_at(&Gri30, 1_500.0, &fractions);
        let state = composition.mixture_state(&Gri30, 0, 1.0, energy).unwrap();
        assert!((state.temperature - 1_500.0).abs() < 1.0);
        assert!(state.pressure > 0.0);
        assert!(state.sound_speed > 100.0);
    }

    /// a sixty-five species contact moves every row with the flow. this
    /// exercises the generic transport path used by the nine- and
    /// fifty-three-species lanes without selecting a reaction subset.
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn a_sixty_five_species_contact_moves_every_row() {
        let ns = 65;
        let nx = 6;
        let ny = 3;
        let ncells = nx * ny;
        let mut y = vec![vec![0.0; ncells]; ns];
        for (cell, row) in y.iter_mut().enumerate().take(ncells) {
            let left_column = cell % nx == 0;
            let witness = if left_column { 0.9 } else { 0.005 };
            let remainder = 1.0 - witness;
            for (species, species_value) in row.iter_mut().enumerate() {
                if species == 4 {
                    *species_value = witness;
                } else {
                    *species_value = remainder / (ns - 1) as f64;
                }
            }
        }
        let mut composition = SpatialComposition { y };
        let density = vec![1.0; ncells];
        let dxs = vec![0.1; nx];
        let dys = vec![0.1; ny];
        let mut fx = vec![vec![0.0; nx + 1]; ny];
        let fy = vec![vec![0.0; ny + 1]; nx];
        for row in &mut fx {
            for face in row.iter_mut() {
                *face = 0.2;
            }
        }
        let dt = 0.05;
        advance_composition(&mut composition, &density, dt, &dxs, &dys, &fx, &fy, None).unwrap();
        let receiver = 1;
        assert!(composition.y[4][receiver] > 0.005);
        for species in 0..ns {
            assert!(composition.y[species].iter().all(|value| value.is_finite()));
            assert!(composition.y[species].iter().all(|value| *value >= 0.0));
        }
        assert!(composition.max_sum_error() < 1.0e-12);
    }

    /// the original nine-species scheme uses the same interface.
    #[test]
    fn the_nine_species_scheme_uses_the_same_interface() {
        let fractions = vec![0.0; Li2004.n_species()];
        let fractions = [
            fractions[0],
            fractions[1],
            fractions[2],
            fractions[3],
            fractions[4],
            fractions[5],
            fractions[6],
            fractions[7],
            fractions[8],
        ];
        assert!(SpatialComposition::uniform(&Li2004, 1, &fractions).is_err());
    }
}

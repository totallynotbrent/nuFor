//! wide-composition chamber state: the same dual-lane structure as
//! [`super::CombustionState`] but carrying an arbitrary scheme width, so
//! the GRI-Mech 53-species mechanism can ride the chamber march through
//! the exact transport and boundary seams the nine-species lane uses.
//!
//! the nine-species lane keeps its fixed-width arrays untouched; this is
//! the additive twin, selected at the case level, not a refactor of the
//! shipped path.

use crate::chem_scheme::ChemScheme;
use crate::error::Error;

/// per-cell composition over a runtime-width scheme, plus the injected
/// fill. mirrors [`super::CombustionState`] field for field.
#[derive(Clone, Debug)]
pub struct WideCombustionState {
    /// `y[species][cell]` mass fractions.
    pub y: Vec<Vec<f64>>,
    /// the composition the case injects at inflow faces.
    pub fill: Vec<f64>,
    /// the scheme this state is laid out against.
    pub n_species: usize,
}

impl WideCombustionState {
    /// uniform field of `ncells` cells at the given fractions, which must
    /// sum to one; the fill is the same vector.
    pub fn uniform(ncells: usize, fractions: &[f64]) -> Result<Self, Error> {
        let n_species = fractions.len();
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > 1.0e-12 || fractions.iter().any(|f| !f.is_finite() || *f < 0.0) {
            return Err(Error::InvalidArgs);
        }
        Ok(WideCombustionState {
            y: fractions.iter().map(|f| vec![*f; ncells]).collect(),
            fill: fractions.to_vec(),
            n_species,
        })
    }

    /// uniform field with a DIFFERENT injection fill (burned gas in, unburned
    /// everywhere else, etc.). both vectors must sum to one.
    pub fn with_fill(ncells: usize, fractions: &[f64], fill: &[f64]) -> Result<Self, Error> {
        if fractions.len() != fill.len() {
            return Err(Error::InvalidArgs);
        }
        let mut state = Self::uniform(ncells, fractions)?;
        let sum: f64 = fill.iter().sum();
        if (sum - 1.0).abs() > 1.0e-12 {
            return Err(Error::InvalidArgs);
        }
        state.fill = fill.to_vec();
        Ok(state)
    }

    /// convenience constructors tied to schemes, so callers stay off raw
    /// indices where possible.
    pub fn gri30_uniform(
        ncells: usize,
        scheme: &dyn ChemScheme,
        fractions: &[f64],
    ) -> Result<Self, Error> {
        if fractions.len() != scheme.n_species() {
            return Err(Error::InvalidArgs);
        }
        Self::uniform(ncells, fractions)
    }

    pub fn ncells(&self) -> usize {
        if self.y.is_empty() {
            0
        } else {
            self.y[0].len()
        }
    }

    /// copy cell `k`'s composition into a scratch buffer.
    pub fn cell(&self, k: usize, out: &mut [f64]) {
        for (s, row) in self.y.iter().enumerate() {
            out[s] = row[k];
        }
    }

    /// write a scratch buffer back into cell `k`.
    pub fn set_cell(&mut self, k: usize, src: &[f64]) {
        for (s, row) in self.y.iter_mut().enumerate() {
            row[k] = src[s];
        }
    }

    /// rescale cell `k` to sum to one without clipping negatives.
    pub fn renorm(&mut self, k: usize) {
        let sum: f64 = self.y.iter().map(|row| row[k]).sum();
        if sum > 0.0 {
            for row in self.y.iter_mut() {
                row[k] /= sum;
            }
        }
    }

    /// worst per-cell deviation of the composition sum from one.
    pub fn max_sum_drift(&self) -> f64 {
        let ncells = self.ncells();
        let mut worst = 0.0f64;
        for k in 0..ncells {
            let sum: f64 = self.y.iter().map(|row| row[k]).sum();
            worst = worst.max((sum - 1.0).abs());
        }
        worst
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chem_scheme::Gri30;

    /// a uniform gri30 fill constructs with the scheme's width and sums to one.
    #[test]
    fn gri30_uniform_width_is_53() {
        let scheme = Gri30;
        let fractions = vec![1.0 / 53.0; 53];
        let state = WideCombustionState::gri30_uniform(8, &scheme, &fractions).unwrap();
        assert_eq!(state.n_species, 53);
        assert_eq!(state.ncells(), 8);
        assert!(state.max_sum_drift() < 1.0e-12);
        let sum: f64 = state.fill.iter().sum();
        assert!((sum - 1.0).abs() < 1.0e-12);
    }

    /// mismatched lengths and bad sums are rejected, not normalized silently.
    #[test]
    fn malformed_vectors_are_rejected() {
        assert!(WideCombustionState::uniform(4, &[0.5, 0.4]).is_err());
        assert!(WideCombustionState::with_fill(4, &[0.5, 0.5], &[0.9, 0.2]).is_err());
        assert!(WideCombustionState::uniform(4, &[0.5, 0.5, -0.1]).is_err());
    }

    /// renorm rescales a drifted cell back to one.
    #[test]
    fn renorm_restores_unit_sum() {
        let mut state = WideCombustionState::uniform(2, &[0.25, 0.75]).unwrap();
        state.y[0][0] = 0.5;
        state.y[1][0] = 0.75;
        assert!(state.max_sum_drift() > 0.1);
        state.renorm(0);
        assert!(state.max_sum_drift() < 1.0e-15);
    }
}

//! 5-species air state and mixture thermodynamics (rung 8).
//!
//! species ride beside the conserved state as mass fractions `y`,
//! one array per species. the perfect-gas and eqair lanes never
//! allocate any of this: the species machinery is only built when a
//! case asks for `eos = "park5"`.
//!
//! species order is fixed across the crate (`park::N_N2 ..= N_O`):
//! N2, O2, NO, N, O.

use crate::error::Error;
use crate::park::{M, NS};

/// per-cell species mass fractions, `y[i][k]` = fraction of species i in cell k.
#[derive(Clone, Debug)]
pub struct SpeciesState {
    pub y: Vec<Vec<f64>>,
}

impl SpeciesState {
    /// uniform composition: `fractions` must sum to 1.
    pub fn uniform(n: usize, fractions: &[f64]) -> Result<Self, Error> {
        if fractions.len() != NS {
            return Err(Error::InvalidArgs);
        }
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > 1e-9 {
            return Err(Error::InvalidArgs);
        }
        Ok(SpeciesState {
            y: fractions.iter().map(|f| vec![*f; n]).collect(),
        })
    }

    /// standard air: 0.767 N2 / 0.233 O2 by mass (the 79/21 molar split).
    pub fn air(n: usize) -> Self {
        SpeciesState {
            y: vec![
                vec![0.767; n],
                vec![0.233; n],
                vec![0.0; n],
                vec![0.0; n],
                vec![0.0; n],
            ],
        }
    }

    /// mixture molar mass, kg per mol of mixture, per cell.
    pub fn mix_molar_mass(&self, k: usize) -> f64 {
        let mut inv = 0.0;
        for (yi, mi) in self.y.iter().zip(M.iter()) {
            if yi[k] > 0.0 {
                inv += yi[k] / mi;
            }
        }
        if inv <= 0.0 {
            M[0]
        } else {
            1.0 / inv
        }
    }

    /// clip negatives, renormalize to sum 1. returns the largest
    /// pre-clip violation so the caller can log how hard the
    /// positivity floor was hit.
    pub fn sanitize(&mut self, k: usize) -> f64 {
        let mut sum = 0.0f64;
        let mut worst = 0.0f64;
        for i in 0..NS {
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
            for i in 0..NS {
                self.y[i][k] /= sum;
            }
        }
        worst
    }

    /// worst |sum(y) - 1| over the field, the drift diagnostic for the run log.
    pub fn max_sum_drift(&self) -> f64 {
        let n = self.y[0].len();
        let mut worst = 0.0f64;
        for k in 0..n {
            let mut sum = 0.0;
            for i in 0..NS {
                sum += self.y[i][k];
            }
            worst = worst.max((sum - 1.0).abs());
        }
        worst
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_requires_unit_sum() {
        assert!(SpeciesState::uniform(4, &[0.5, 0.5, 0.0, 0.0, 0.0]).is_ok());
        assert!(SpeciesState::uniform(4, &[0.5, 0.6, 0.0, 0.0, 0.0]).is_err());
        assert!(SpeciesState::uniform(4, &[0.767, 0.233]).is_err());
    }

    #[test]
    fn mix_molar_mass_of_air() {
        let s = SpeciesState::air(1);
        // 79/21 molar -> mm = 1 / (0.79/28.016 + 0.21/32) e-3 kg/mol ~ 28.85 g/mol.
        // the mass-fraction split 0.767/0.233 is the same composition.
        let mm = s.mix_molar_mass(0) * 1e3;
        assert!((mm - 28.85).abs() < 0.05, "mm = {mm}");
    }

    #[test]
    fn pure_species_molar_mass() {
        let mut s = SpeciesState::air(1);
        for i in 0..NS {
            s.y[i][0] = 0.0;
        }
        s.y[crate::park::N_O2][0] = 1.0;
        assert!((s.mix_molar_mass(0) - M[1]).abs() < 1e-15);
    }

    #[test]
    fn sanitize_renormalizes_and_reports() {
        let mut s = SpeciesState::air(1);
        s.y[0][0] = 1.2;
        s.y[1][0] = -0.1;
        let worst = s.sanitize(0);
        assert!((worst - 0.2).abs() < 1e-12);
        let mut sum = 0.0;
        for i in 0..NS {
            sum += s.y[i][0];
            assert!(s.y[i][0] >= 0.0 && s.y[i][0] <= 1.0);
        }
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(s.max_sum_drift() < 1e-12);
    }

    #[test]
    fn frozen_state_has_zero_drift() {
        let s = SpeciesState::air(16);
        assert_eq!(s.max_sum_drift(), 0.0);
    }
}

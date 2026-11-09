//! runtime-selected chemistry mechanism for the reacting chamber.
//!
//! a case names its mechanism in the case file, and every consumer that
//! touches a composition (the march, the faces, the boundaries, the
//! writers) reads the width and the mixture properties from this enum
//! instead of the fixed nine-species constants.
//!
//! the hydrogen arm evaluates the same polynomials and the same closure
//! arithmetic the legacy nine-species lane used, so its numbers do not
//! move; a numerics test holds the two to each other bit for bit.

use crate::chem_scheme::{ChemScheme, Gri30, Li2004};
use crate::combustion_thermo as li;

/// the chemistry mechanism a case runs: the nine-species hydrogen set or
/// the fifty-three-species gri-mech methane set.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ChemLane {
    /// the nine-species hydrogen lane.
    #[default]
    Li2004,
    /// the fifty-three-species methane lane.
    Gri30,
}

impl ChemLane {
    /// the mechanism named in a case file: `"li2004-h2"` or `"gri30"`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "li2004-h2" => Some(ChemLane::Li2004),
            "gri30" => Some(ChemLane::Gri30),
            _ => None,
        }
    }

    /// the mechanism identifier a case file writes.
    pub fn name(&self) -> &'static str {
        match self {
            ChemLane::Li2004 => "li2004-h2",
            ChemLane::Gri30 => "gri30",
        }
    }

    /// the lane's species count.
    pub fn n_species(&self) -> usize {
        match self {
            ChemLane::Li2004 => Li2004.n_species(),
            ChemLane::Gri30 => Gri30.n_species(),
        }
    }

    /// species names in index order, for the writers and the case parser.
    pub fn species_names(&self) -> &'static [&'static str] {
        match self {
            ChemLane::Li2004 => &li::SPECIES_NAMES,
            ChemLane::Gri30 => &crate::gri30_thermo::SPECIES_NAMES,
        }
    }

    /// a species' index by its symbol, or none when the mechanism has no
    /// such species.
    pub fn species_index(&self, name: &str) -> Option<usize> {
        self.species_names()
            .iter()
            .position(|s| s.eq_ignore_ascii_case(name))
    }

    /// a species' molar mass, kg/mol.
    pub fn molar_mass(&self, s: usize) -> f64 {
        match self {
            ChemLane::Li2004 => li::M[s],
            ChemLane::Gri30 => Gri30.molar_mass(s),
        }
    }

    /// the mixture gas constant, j/(kg k).
    ///
    /// the hydrogen arm sums species-wise, which is the form the legacy
    /// closure has always used, so its numbers do not move.
    pub fn mix_gas_constant(&self, y: &[f64]) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::mix_gas_constant(&Li2004, y),
            ChemLane::Gri30 => crate::chem_closure::mix_gas_constant(&Gri30, y),
        }
    }

    /// the mixture's mean molar mass, kg/mol.
    ///
    /// a cell carrying no positive mass falls back to nitrogen on the
    /// hydrogen lane, which is the legacy closure's choice, and to the
    /// scheme's first species on the methane lane.
    pub fn mix_molar_mass(&self, y: &[f64]) -> f64 {
        match self {
            ChemLane::Li2004 => {
                let mut inv = 0.0f64;
                for (s, ys) in y.iter().enumerate() {
                    if *ys > 0.0 {
                        inv += ys / li::M[s];
                    }
                }
                if inv <= 0.0 {
                    li::M[li::IN2]
                } else {
                    1.0 / inv
                }
            }
            ChemLane::Gri30 => crate::chem_closure::mix_molar_mass(&Gri30, y),
        }
    }

    /// the mixture specific heat at constant pressure, j/(kg k).
    pub fn mix_cp_mass(&self, y: &[f64], t: f64) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::mix_cp_mass(&Li2004, y, t),
            ChemLane::Gri30 => crate::chem_closure::mix_cp_mass(&Gri30, y, t),
        }
    }

    /// the mixture static enthalpy per mass, j/kg.
    pub fn mix_h_mass(&self, y: &[f64], t: f64) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::mix_h_mass(&Li2004, y, t),
            ChemLane::Gri30 => crate::chem_closure::mix_h_mass(&Gri30, y, t),
        }
    }

    /// internal energy per mass at temperature `t`, j/kg.
    pub fn internal_energy_at(&self, t: f64, y: &[f64]) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::internal_energy_at(&Li2004, t, y),
            ChemLane::Gri30 => crate::chem_closure::internal_energy_at(&Gri30, t, y),
        }
    }

    /// temperature from internal energy per mass, kelvin.
    pub fn temperature_at(&self, e_int: f64, y: &[f64]) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::temperature_at(&Li2004, e_int, y),
            ChemLane::Gri30 => crate::chem_closure::temperature_at(&Gri30, e_int, y),
        }
    }

    /// pressure from conserved data, pa.
    pub fn pressure_at(&self, rho: f64, e_int: f64, y: &[f64]) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::pressure_at(&Li2004, rho, e_int, y),
            ChemLane::Gri30 => crate::chem_closure::pressure_at(&Gri30, rho, e_int, y),
        }
    }

    /// frozen-composition sound speed, m/s.
    pub fn sound_speed_at(&self, rho: f64, e_int: f64, y: &[f64]) -> f64 {
        match self {
            ChemLane::Li2004 => crate::chem_closure::sound_speed_at(&Li2004, rho, e_int, y),
            ChemLane::Gri30 => crate::chem_closure::sound_speed_at(&Gri30, rho, e_int, y),
        }
    }

    /// the energy a cell may fall to before the closure reads it as
    /// below its temperature range: the composition's energy at the
    /// fits' lower temperature bound, and not zero.
    pub fn floor_energy_at(&self, y: &[f64]) -> f64 {
        self.internal_energy_at(200.0, y)
    }

    /// a zero composition in the lane's width.
    pub fn zero_fill(&self) -> Vec<f64> {
        vec![0.0f64; self.n_species()]
    }

    /// the hydrogen lane's stoichiometric premix, h2/o2 by mass.
    pub fn li2004_stoich_premix(&self) -> Vec<f64> {
        let mut y = self.zero_fill();
        y[li::IH2] = 2.016 / 18.016;
        y[li::IO2] = 16.0 / 18.016;
        y
    }
}

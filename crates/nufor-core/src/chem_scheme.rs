//! the chemistry schemes the combustion lane can run.
//!
//! two mechanisms sit beside each other today: the nine-species li-2004
//! hydrogen set and the fifty-three-species gri-mech 3.0 methane set.
//! both bring the same seven things to the solver — a species count,
//! molar masses, heat capacities, enthalpies, and transport properties —
//! so they meet the solver through one interface.
//!
//! the interface is deliberately small. a scheme answers questions
//! about individual species; the mixture rules live once in
//! `chem_closure`, which works on slices and therefore cares nothing
//! about how many species there are.
//!
//! the li-2004 numbers must not move, so `Li2004` forwards to the
//! tables that were already there, and `Gri30` forwards to the
//! generated mechanism. the equivalence test in the crate's test
//! directory checks the shared closure against the li-2004 closure it
//! grew out of.

use crate::combustion_thermo as li;
use crate::gri30;
use crate::gri30_transport_generated::{Shape, TRANSPORT_ROWS};

/// what the solver needs to know about one chemistry scheme.
pub trait ChemScheme {
    /// the scheme's name, as it appears in a case file.
    fn name(&self) -> &'static str;

    /// the number of species the scheme carries.
    fn n_species(&self) -> usize;

    /// the molar mass of species `s`, kg/mol.
    fn molar_mass(&self, s: usize) -> f64;

    /// the specific heat of species `s` at temperature `t`, j/(kg k).
    fn cp_mass(&self, s: usize, t: f64) -> f64;

    /// the enthalpy of species `s` at temperature `t`, j/kg, including
    /// its formation enthalpy.
    fn h_mass(&self, s: usize, t: f64) -> f64;

    /// the lennard-jones well depth of species `s`, kelvin.
    fn well_depth(&self, s: usize) -> f64;

    /// the lennard-jones collision diameter of species `s`, angstrom.
    fn diameter(&self, s: usize) -> f64;

    /// the molecular shape of species `s`.
    fn shape(&self, s: usize) -> Shape;

    /// the universal gas constant in j/(mol k), matching the molar
    /// mass basis the scheme reports.
    fn ru(&self) -> f64;
}

/// the nine-species li-2004 hydrogen scheme.
#[derive(Clone, Copy, Debug, Default)]
pub struct Li2004;

/// the fifty-three-species gri-mech 3.0 methane scheme.
#[derive(Clone, Copy, Debug, Default)]
pub struct Gri30;

impl ChemScheme for Li2004 {
    fn name(&self) -> &'static str {
        "li2004-h2"
    }

    fn n_species(&self) -> usize {
        li::NS
    }

    fn molar_mass(&self, s: usize) -> f64 {
        li::M[s]
    }

    fn cp_mass(&self, s: usize, t: f64) -> f64 {
        li::cp_mass(s, t)
    }

    fn h_mass(&self, s: usize, t: f64) -> f64 {
        li::h_mass(s, t)
    }

    fn well_depth(&self, s: usize) -> f64 {
        crate::combustion_transport::SPECIES_WELL_DEPTH[s]
    }

    fn diameter(&self, s: usize) -> f64 {
        crate::combustion_transport::SPECIES_DIAMETER[s]
    }

    fn shape(&self, s: usize) -> Shape {
        // the nine-species set is h2, o2, h2o, oh, h, o, ho2, h2o2, n2;
        // only water and the two radicals are bent, the rest are
        // diatomic or atomic.
        match s {
            li::IH2O | li::IHO2 => Shape::Nonlinear,
            li::IH2 | li::IO2 | li::IOH | li::IN2 => Shape::Linear,
            _ => Shape::Atom,
        }
    }

    fn ru(&self) -> f64 {
        li::RU
    }
}

impl ChemScheme for Gri30 {
    fn name(&self) -> &'static str {
        "gri30"
    }

    fn n_species(&self) -> usize {
        crate::gri30_thermo::NS_GRI30
    }

    fn molar_mass(&self, s: usize) -> f64 {
        // the mechanism's masses are kg/kmol; the solver's basis is
        // kg/mol.
        gri30::molar_mass(s) * 1.0e-3
    }

    fn cp_mass(&self, s: usize, t: f64) -> f64 {
        gri30::cp_mass(s, t)
    }

    fn h_mass(&self, s: usize, t: f64) -> f64 {
        gri30::h_mass(s, t)
    }

    fn well_depth(&self, s: usize) -> f64 {
        TRANSPORT_ROWS[s].well_depth
    }

    fn diameter(&self, s: usize) -> f64 {
        TRANSPORT_ROWS[s].diameter
    }

    fn shape(&self, s: usize) -> Shape {
        TRANSPORT_ROWS[s].shape
    }

    fn ru(&self) -> f64 {
        gri30::RU_KMOL * 1.0e-3
    }
}

/// the transport rows of a scheme, ready for the transport formulas.
///
/// the rows carry the mass in kg/mol, the depth in kelvin, and the
/// diameter in angstrom, which is what `transport_general` expects.
pub fn transport_rows<S: ChemScheme>(
    scheme: &S,
) -> Vec<crate::transport_general::SpeciesTransport> {
    (0..scheme.n_species())
        .map(|s| crate::transport_general::SpeciesTransport {
            molar_mass: scheme.molar_mass(s),
            well_depth: scheme.well_depth(s),
            diameter: scheme.diameter(s),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the two schemes report the species counts the mechanisms claim.
    #[test]
    fn scheme_counts_are_right() {
        assert_eq!(Li2004.n_species(), 9);
        assert_eq!(Gri30.n_species(), 53);
    }

    /// the hydrogen lane's masses come straight from its own table.
    #[test]
    fn li2004_forwards_to_its_own_tables() {
        for s in 0..Li2004.n_species() {
            assert_eq!(Li2004.molar_mass(s), li::M[s]);
            assert_eq!(Li2004.ru(), li::RU);
        }
    }

    /// the methane lane's masses convert from the mechanism's kmol
    /// basis to the solver's mol basis without losing precision.
    #[test]
    fn gri30_masses_convert_cleanly() {
        assert!((Gri30.molar_mass(0) - 2.016e-3).abs() < 1e-12);
        assert!((Gri30.molar_mass(3) - 31.998e-3).abs() < 1e-12);
        assert!((Gri30.ru() - 8.314462618).abs() < 1e-9);
    }

    /// every species of both schemes carries transport data.
    #[test]
    fn every_species_has_transport_data() {
        for s in 0..Li2004.n_species() {
            assert!(Li2004.well_depth(s) > 0.0, "hydrogen species {s}");
            assert!(Li2004.diameter(s) > 0.0, "hydrogen species {s}");
        }
        for s in 0..Gri30.n_species() {
            assert!(Gri30.well_depth(s) > 0.0, "methane species {s}");
            assert!(Gri30.diameter(s) > 0.0, "methane species {s}");
        }
    }
}

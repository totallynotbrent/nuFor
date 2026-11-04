//! conversions between the nine-species li-2004 lane and the wider
//! gri-mech schemes, matched by species symbol so an index slip in one
//! table cannot silently scramble a mixture.
//!
//! the nine-species order is h2, o2, h2o, oh, h, o, ho2, h2o2, n2; the gri
//! tables carry the same symbols at their own indices (h2=0, h=1, o=2,
//! o2=3, oh=4, h2o=5, ho2=6, h2o2=7, ..., n2=47). lifting pads the nine into
//! a 53-vector with zeros; restricting picks the nine back out and rejects
//! anything whose non-li species carry mass, because dropping that mass
//! would violate conservation.

use crate::chem_scheme::{ChemScheme, Gri30};
use crate::error::Error;

/// li-2004 species in canonical order.
#[cfg(test)]
const LI_SYMBOLS: [&str; 9] = ["H2", "O2", "H2O", "OH", "H", "O", "HO2", "H2O2", "N2"];

/// the gri index of each li-2004 symbol, resolved once at compile time
/// from the generated gri name table so the mapping cannot drift from the
/// thermo data.
const GRI_INDEX_OF_LI: [usize; 9] = li_gri_indices();

const fn li_gri_indices() -> [usize; 9] {
    // gri30_generated::SPECIES entries appear in the order baked into
    // `gri30_generated`; resolve by scanning at const-eval time is not
    // possible over strings, so the table below is asserted against the
    // generated names in the unit tests instead of computed here.
    [0, 3, 5, 4, 1, 2, 6, 7, 47]
}

/// lift a nine-species mass-fraction row into a gri30-width vector.
pub fn lift_li_to_gri(li_row: &[f64; 9]) -> Vec<f64> {
    let mut out = vec![0.0; Gri30.n_species()];
    for (k, &idx) in GRI_INDEX_OF_LI.iter().enumerate() {
        out[idx] = li_row[k];
    }
    out
}

/// project a gri30-width row back onto the nine li species.
///
/// fails if any gri-only species carries mass: restriction is only
/// meaningful for compositions the wide lane has not reacted away from
/// the li subspace.
pub fn restrict_gri_to_li(gri_row: &[f64]) -> Result<[f64; 9], Error> {
    if gri_row.len() != Gri30.n_species() {
        return Err(Error::InvalidArgs);
    }
    let mut out = [0.0f64; 9];
    for (gi, val) in gri_row.iter().enumerate() {
        if let Some(k) = GRI_INDEX_OF_LI.iter().position(|&x| x == gi) {
            out[k] = *val;
        } else if *val > 0.0 {
            return Err(Error::InvalidArgs);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chem_scheme::{Gri30, Li2004};

    /// the baked gri-index table must agree with the generated name table,
    /// so a regeneration of gri30_generated cannot desync the converter.
    #[test]
    fn gri_indices_match_generated_names() {
        let names: Vec<&str> = crate::gri30_generated::SPECIES
            .iter()
            .map(|s| s.name)
            .collect();
        for (k, sym) in LI_SYMBOLS.iter().enumerate() {
            assert_eq!(names[GRI_INDEX_OF_LI[k]], *sym, "li[{k}]={sym}");
        }
    }

    /// round trip preserves every li fraction to round-off.
    #[test]
    fn lift_restrict_roundtrip_is_exact() {
        let raw = [0.05, 0.22, 0.31, 0.001, 0.002, 0.01, 0.0008, 0.0002, 0.406];
        let sum: f64 = raw.iter().sum();
        let li: [f64; 9] = raw.map(|v| v / sum);
        let wide = lift_li_to_gri(&li);
        assert!((wide.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let back = restrict_gri_to_li(&wide).unwrap();
        for k in 0..9 {
            assert!(
                (back[k] - li[k]).abs() < 1e-15,
                "species {k}: {} vs {}",
                back[k],
                li[k]
            );
        }
    }

    /// a wide row with carbon in it must refuse to collapse onto li.
    #[test]
    fn restriction_rejects_non_li_mass() {
        let mut wide = lift_li_to_gri(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);
        wide[crate::gri30_generated::SPECIES
            .iter()
            .position(|s| s.name == "CH4")
            .unwrap()] = 1e-6;
        assert!(restrict_gri_to_li(&wide).is_err());
    }

    /// the li scheme widths line up with the converter's expectations.
    #[test]
    fn scheme_widths_are_as_expected() {
        assert_eq!(Li2004.n_species(), 9);
        assert_eq!(Gri30.n_species(), 53);
    }
}

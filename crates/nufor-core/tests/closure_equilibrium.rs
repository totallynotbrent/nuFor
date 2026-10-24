//! the closure's equilibrium constants against cea's chamber table.
//!
//! the reaction quotients at cea's own composition must match the
//! closure's equilibrium constants, which is the condition for the
//! mechanism to hold that composition. every reaction's quotient is
//! checked, so a wrong equilibrium constant in one branch fails here
//! rather than quietly moving the chemistry.

use nufor_core::combustion::NU_PAIRS;
use nufor_core::combustion_thermo as ct;

/// cea's 70 bar chamber composition in mass fractions.
const CEA_CHAMBER: [f64; ct::NS] = [
    0.01583304, 0.07194470, 0.77605355, 0.11605824, 0.00232724, 0.01719592, 0.00049223, 0.00009298,
    0.0,
];

/// the equilibrium constant's form: the closure's Kc is the
/// concentration form, so the mole-fraction quotient picks up the
/// total concentration.
fn quotient_residual(t: f64, p: f64, y: &[f64; ct::NS], r: usize) -> f64 {
    let mut x = [0.0f64; ct::NS];
    for s in 0..ct::NS {
        x[s] = y[s] / ct::M[s];
    }
    let xsum: f64 = x.iter().sum::<f64>();
    for v in x.iter_mut() {
        *v /= xsum;
    }
    let mut q = 0.0f64;
    let mut dnu = 0.0f64;
    for &(s, nu) in NU_PAIRS[r].iter() {
        if nu == 0.0 {
            continue;
        }
        q += nu * x[s].max(1e-300).ln();
        dnu += nu;
    }
    // the concentration quotient: n_s = x_s * p / (R T).
    let ln_q = q + dnu * (p / (ct::RU * t)).ln();
    ln_q - ct::ln_kc(&NU_PAIRS[r], t)
}

/// every reaction's quotient at cea's chamber composition must match
/// the closure's equilibrium constant.
#[test]
fn the_equilibrium_constants_match_ceas_chamber() {
    let t = 3676.53f64;
    let p = 70.0e5f64;
    let mut worst = 0.0f64;
    let mut worst_reaction = 0usize;
    for r in 0..nufor_core::combustion::NR {
        let res = quotient_residual(t, p, &CEA_CHAMBER, r).abs();
        if res > worst {
            worst = res;
            worst_reaction = r;
        }
    }
    // the residuals are the log of the quotient's error: a tenth means
    // the constant is within about ten percent of cea's.
    assert!(
        worst < 0.15,
        "the equilibrium constants drift from cea's: reaction {worst_reaction} at {worst:.4}"
    );
}

/// the dissociation branch carries the tightest constants: the
/// hydrogen and oxygen dissociation and the hydroxyl formation must
/// all agree with cea's composition.
#[test]
fn the_dissociation_branches_agree_with_cea() {
    let t = 3676.53f64;
    let p = 70.0e5f64;
    // the dissociation and formation reactions: H2 -> 2H, 2O -> O2,
    // O + H -> OH, H2O -> 2OH, and the peroxide branches.
    for r in [4usize, 5, 6, 14, 15, 16, 17, 18] {
        let res = quotient_residual(t, p, &CEA_CHAMBER, r).abs();
        assert!(
            res < 0.15,
            "reaction {r}'s equilibrium constant drifts from cea's: {res:.4}"
        );
    }
}

/// the chamber's own closure must reproduce cea's chamber pressure and
/// sound speed at cea's composition, temperature and density.
#[test]
fn the_chamber_state_matches_cea() {
    use nufor_core::combustion_closure as cl;
    let rho = 3.6116f64;
    let t = 3676.53f64;
    let e = cl::internal_energy_at(t, &CEA_CHAMBER);
    let p = cl::pressure_at(rho, e, &CEA_CHAMBER);
    assert!(
        (p - 70.0e5).abs() / 70.0e5 < 0.02,
        "the chamber pressure {p} against cea's 70 bar"
    );
    let a = cl::sound_speed_at(rho, e, &CEA_CHAMBER);
    assert!(
        (a - 1480.49).abs() / 1480.49 < 0.05,
        "the chamber sound speed {a} against cea's 1480.49"
    );
}

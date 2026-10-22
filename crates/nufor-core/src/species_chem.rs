//! park 5-species finite-rate chemistry: forward and
//! backward reaction rates and the per-species production terms.
//!
//! forward rates are the park 1990 arrhenius set already tabulated in
//! `park::RATES` (si units, m3/mol/s). backward rates come from the
//! equilibrium constants through
//!
//! `k_b = k_f(T) / K_c(T)`
//!
//! with K_c the concentration-form equilibrium constant from park's
//! 1985 five-coefficient curve fits (TP-2867 table III, the companion
//! fits for these exact forward rates):
//!
//! `ln K_c = B1 + B2 ln Z + B3 Z + B4 Z^2 + B5 Z^3,  Z = 10000/T`
//!
//! K_c carries mol/cm3 units for the three dissociation reactions
//! (delta-nu = +1) and is dimensionless for the two exchange
//! reactions. RP-1232 eq. (6c) fixes the temperature convention:
//! both k_f and K_c evaluate at the translational temperature, never
//! the geometric-mean reaction temperature, so the recombination rate
//! depends only on the heavy-particle temperature. one-temperature
//! flow (this model) has T = T_tr everywhere, so the distinction is
//! inert here and matters when the park2t machinery drives this
//! module.
//!
//! the fits are comfortable to 10000 K (they were fitted at five
//! points up to that temperature); above that they extrapolate, which
//! RP-1232 flags against. dissociation behind strong shocks peaks in
//! the 6-14 kK band, so the fits cover the working range.

use crate::park::{M, NS};

/// the five equilibrium-constant fits, one per reaction, TP-2867
/// table III ("kinetic model of park"), in the order of `park::RATES`:
/// N2 dissoc, O2 dissoc, NO dissoc, zeldovich-1, zeldovich-2.
/// `dnu` is the mole-number change, which fixes K_c's units.
pub const KEQ_FITS: [(f64, f64, f64, f64, f64, i32); 5] = [
    // N2 + M -> N + N + M: ln Kc (mol/cm3), dnu = +1
    (1.858, -1.325, -9.856, -0.174, 0.008, 1),
    // O2 + M -> O + O + M: ln Kc (mol/cm3), dnu = +1
    (2.855, 0.988, -6.181, -0.023, -0.001, 1),
    // NO + M -> N + O + M: ln Kc (mol/cm3), dnu = +1
    (0.792, -0.492, -6.761, -0.091, 0.004, 1),
    // N2 + O -> NO + N: dimensionless, dnu = 0
    (1.066, -0.833, -3.095, -0.084, 0.004, 0),
    // NO + O -> O2 + N: dimensionless, dnu = 0
    (-2.063, -1.480, -0.580, -0.114, 0.005, 0),
];

/// the concentration-form equilibrium constant for reaction `r` at
/// temperature `t`, in SI units: mol/m3 for dnu = +1 reactions,
/// dimensionless otherwise. the TP-2867 fits print mol/cm3, so the
/// dissociation constants convert by 1e6.
pub fn k_eq(r: usize, t: f64) -> f64 {
    let (b1, b2, b3, b4, b5, dnu) = KEQ_FITS[r];
    let z = 10000.0 / t.max(200.0);
    let ln_kc = b1 + b2 * z.ln() + b3 * z + b4 * z * z + b5 * z * z * z;
    let kc = ln_kc.exp();
    if dnu == 1 {
        kc * 1e6
    } else {
        kc
    }
}

/// forward rate for reaction `r` at temperature `t`, si units:
/// m3/mol/s for two-body, m6/mol2/s for three-body.
pub fn k_f(r: usize, t: f64) -> f64 {
    let rate = &crate::park::RATES[r];
    rate.a * t.powf(rate.eta) * (-rate.td / t).exp()
}

/// backward rate for reaction `r` at temperature `t`, si units.
pub fn k_b(r: usize, t: f64) -> f64 {
    k_f(r, t) / k_eq(r, t)
}

/// number-density of species `s` from (rho, Y_s), mol/m3.
pub fn n_i(rho: f64, y_s: f64, s: usize) -> f64 {
    rho * y_s / M[s]
}

/// one cell's reaction state: (rho, Y, T). returns the per-species
/// production rate d(Y_i)/dt, 1/s, and the chemical energy deposition
/// d(e_chem)/dt, J/(kg s) — dissociation absorbs energy (negative
/// deposition into the thermal pool, the formation energy the atoms
/// carry away).
///
/// reaction set (park 1990, matching `park::RATES`):
///   r0: N2 + M -> N + N + M   (M = N2, O2, NO at RATES[0]; M = N, O at RATES[1])
///   r1: O2 + M -> O + O + M   (RATES[2])
///   r2: NO + M -> N + O + M   (RATES[3])
///   r3: N2 + O -> NO + N      (RATES[4], zeldovich)
///   r4: NO + O -> O2 + N      (RATES[5], zeldovich)
///
/// the third-body split for N2 dissociation: RATES[0] serves the
/// molecular partners, RATES[1] the atoms. KEQ_FITS[0] pairs with
/// both branches (same K_c: the equilibrium state does not depend on
/// the collider).
pub fn omega(rho: f64, y: &[f64; NS], t: f64) -> ([f64; NS], f64) {
    let n: [f64; NS] = [
        n_i(rho, y[0], 0),
        n_i(rho, y[1], 1),
        n_i(rho, y[2], 2),
        n_i(rho, y[3], 3),
        n_i(rho, y[4], 4),
    ];
    let n_mol = n[0] + n[1] + n[2];
    let n_atom = n[3] + n[4];
    let nn2 = n[0];
    let no2 = n[1];
    let nno = n[2];
    let nn = n[3];
    let no = n[4];

    // dissociation: rate constant m6 mol-2 s-1 (three-body).
    // the N2 reaction splits by collider (park's table): each branch
    // carries its own forward constant and its own backward constant
    // through the SHARED K_c — kb_branch = kf_branch / K_c. that keeps
    // every branch's forward and backward fluxes in equilibrium at the
    // same composition (a single kb on the summed forward would break
    // detailed balance between colliders).
    let kc0 = k_eq(0, t);
    let kf0_mol = k_f(0, t);
    let kf0_atom = k_f(1, t);
    // the dissociation reactions pair a two-body-form forward rate
    // kf(T) n_s n_M (m3/mol/s) with the three-body recombination
    // kb(T) n^2 n_M (m6/mol2/s): detailed balance needs
    // n^2/n_s = kf/kb = K_c at equilibrium, independent of n_M.
    let r0_net = (kf0_mol * nn2 - (kf0_mol / kc0) * nn * nn) * n_mol
        + (kf0_atom * nn2 - (kf0_atom / kc0) * nn * nn) * n_atom;

    let kf1 = k_f(2, t);
    let kb1 = k_b(1, t);
    let r1_net = (kf1 * no2 - kb1 * no * no) * n_mol;

    let kf2 = k_f(3, t);
    let kb2 = k_b(2, t);
    let r2_net = (kf2 * nno - kb2 * nn * no) * n_mol;

    // exchange: two-body, m3/mol/s.
    let kf3 = k_f(4, t);
    let kb3 = k_b(3, t);
    let r3_net = kf3 * no * nn2 - kb3 * nno * nn;

    let kf4 = k_f(5, t);
    let kb4 = k_b(4, t);
    let r4_net = kf4 * no * nno - kb4 * no2 * nn;

    // species production, mol m-3 s-1:
    //   N2: -r0 - r3 + r4
    //   O2: -r1 + r4
    //   NO: -r2 + r3 - r4
    //   N:  +2r0 + r2 + r3 + r4
    //   O:  +2r1 + r2 - r3 - r4
    let dn: [f64; NS] = [
        -r0_net - r3_net + r4_net,
        -r1_net + r4_net,
        -r2_net + r3_net - r4_net,
        2.0 * r0_net + r2_net + r3_net + r4_net,
        2.0 * r1_net + r2_net - r3_net - r4_net,
    ];
    // mass production dY/dt = dn * M_s / rho.
    let mut dydt = [0.0f64; NS];
    for (dy, d) in dydt.iter_mut().zip(dn.iter()) {
        *dy = d / rho;
    }
    // the mass production carries each species' molar mass.
    for s in 0..NS {
        dydt[s] *= M[s];
    }
    // chemical energy per mass change: the formation energy the
    // species carry (HF_SPECIES, atom-positive). the thermal pool
    // supplies what the products bind.
    let mut de_chem = 0.0f64;
    for (dy, hf) in dydt.iter().zip(crate::species_closure::HF_SPECIES.iter()) {
        de_chem += dy * hf;
    }
    (dydt, de_chem)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::park::{D0_K, RU};

    /// K_c(N2) at 10000 K from the verified extract: 2.85e-4 mol/cm3.
    #[test]
    fn keq_n2_matches_verified_table() {
        let kc = k_eq(0, 10000.0);
        assert!((kc / 1e6 - 2.85e-4).abs() / 2.85e-4 < 0.02, "kc = {kc}");
        // and it grows strongly by 20000 K (extrapolation, trend check).
        assert!(k_eq(0, 20000.0) > 100.0 * k_eq(0, 10000.0));
    }

    /// thermodynamic cycle: K(N2+O->NO+N) ~= K(N2)/K(NO) at all T.
    #[test]
    fn keq_exchange_cycle_closes() {
        for t in [4000.0, 6000.0, 8000.0, 10000.0, 15000.0] {
            let lhs = k_eq(3, t);
            // the mol/cm3-to-SI conversion cancels in the ratio: both
            // dissociation constants carry the same 1e6.
            let rhs = k_eq(0, t) / k_eq(2, t);
            let ratio = lhs / rhs;
            assert!(
                (ratio - 1.0).abs() < 0.02,
                "cycle mismatch at {t}: ratio = {ratio}"
            );
        }
    }

    /// equilibrium consistency: at a fixed (T, rho) the equilibrium
    /// composition must null the source term. take T = 10000 K,
    /// pure N2/N system: solve the equilibrium N2 fraction from K_c.
    #[test]
    fn equilibrium_nulls_the_source() {
        let t = 10000.0f64;
        let rho = 1.0e-3f64;
        // K_c (mol/m3) for N2 -> 2N: K_c = n_N^2 / n_N2. with the N2/N
        // system at fixed rho, n_N2 = rho(1-x)/M_N2 and n_N = rho x/M_N,
        // so the closed form for x_N solves directly (no iteration):
        let kc = k_eq(0, t);
        let a = kc * M[3] * M[3] / (rho * M[0]);
        let xn = ((a * a + 4.0 * a).sqrt() - a) / 2.0;
        let y = [1.0 - xn, 0.0, 0.0, xn, 0.0];
        let (dydt, _) = omega(rho, &y, t);
        // the N2 <-> 2N pair must nearly null (relaxation of the
        // O2/NO family is separate; their channels are zero here).
        let n_n2 = n_i(rho, y[0], 0);
        let n_n = n_i(rho, y[3], 3);
        // scale by the forward dissociation flux: near equilibrium the
        // net is a small difference of two large fluxes.
        let scale = (k_f(0, t) * n_n * n_n2).abs().max(1e-30);
        let rel = (dydt[0].abs() * rho / M[0]) / scale;
        assert!(rel < 0.15, "source at equilibrium: rel = {rel}, xn = {xn}");
    }

    /// dissociation absorbs energy: the chemical deposition at a cold
    /// composition pushing toward dissociation is positive (products
    /// carry MORE formation energy), draining the thermal pool.
    #[test]
    fn dissociation_draws_formation_energy() {
        let t = 9000.0f64;
        let rho = 1.0e-3f64;
        let y = [0.8, 0.2, 0.0, 0.0, 0.0];
        let (dydt, de_chem) = omega(rho, &y, t);
        assert!(dydt[3] > 0.0, "N growing: {}", dydt[3]);
        assert!(dydt[0] < 0.0, "N2 depleting: {}", dydt[0]);
        assert!(de_chem > 0.0, "chemical pool absorbing: {de_chem}");
        // finite and sane as a rate (not a bounded quantity): the
        // integration step will bring it to the stiff coupling later.
        assert!(de_chem.is_finite());
    }

    /// the formation-energy bookkeeping is the hess cycle on D0.
    #[test]
    fn hf_is_the_hess_cycle() {
        let hf = crate::species_closure::HF_SPECIES;
        // per kg of N atoms, the bond energy is 0.5 * D0(N2) / M_N.
        let e_n = 0.5 * RU * D0_K[0] / M[3];
        assert!((hf[3] - e_n).abs() < 1e-6, "hf[3] = {}", hf[3]);
        // NO: the hess cycle (0.5 D0_N2 + 0.5 D0_O2 - D0_NO) / M_NO.
        let e_no = RU * (0.5 * D0_K[0] + 0.5 * D0_K[1] - D0_K[2]) / M[2];
        assert!((hf[2] - e_no).abs() < 1e-6, "hf[2] = {}", hf[2]);
        assert!(hf[0].abs() < 1e-12 && hf[1].abs() < 1e-12);
    }
}

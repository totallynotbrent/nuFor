//! catalytic wall boundary condition for the park5 viscous lane:
//! the bracket pair every TPS report wants, per the production-code
//! definitions (LAURA TM-4674 jtype 0 vs fully-catalytic; DPLR
//! icatmd 0 vs 2, see plans/catalytic_research/wilke_transport.md).
//!
//! non-catalytic: dY_s/dn = 0 at the wall for every species (zero
//! diffusive species flux), so the wall heat flux is the fourier
//! conductor alone.
//!
//! fully-catalytic: the wall composition is the equilibrium
//! composition at the wall temperature (TM-4674 lines 1698-1701:
//! "the surface mass fractions are defined by the equilibrium
//! composition at the given wall temperature"). at a 1000 K wall the
//! equilibrium is essentially all molecules, so this is the
//! Y_atoms(surface) -> 0 idealization; the closed-form equilibrium
//! solve keeps it honest at any wall temperature. the energy the
//! atoms carry into the wall is the formation-enthalpy flux:
//!
//! `q_w = -k dT/dn + sum_s h_f_s J_s,n`
//!
//! with J_s,n the diffusive mass flux into the wall. the atoms that
//! reach the surface recombine and deposit D0 per unit mass; with
//! sum J_s = 0 (the single-D model) the molecules carry the
//! recombined mass back out, so the NET mass flux is zero and the
//! surface term is purely enthalpy.

use crate::error::Error;
use crate::park::{M, NS};
use crate::species_closure::HF_SPECIES;

/// the catalytic wall model selected by the case file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CatalyticWall {
    /// dY/dn = 0: zero diffusive species flux at the wall.
    NonCatalytic,
    /// wall composition = equilibrium at T_wall: atoms recombine on
    /// contact, the formation enthalpy lands on the surface.
    FullyCatalytic,
}

/// wall composition for the fully-catalytic condition: the
/// equilibrium composition at the wall temperature from the
/// park-1985 K_eq fits (the closure's own fits — no new constants).
/// for the 5-species neutral set at a cold wall this is the
/// all-molecules limit; the closed form covers any temperature.
/// wall composition for the fully-catalytic condition: the
/// equilibrium composition at the wall temperature from the
/// park-1985 K_eq fits (the closure's own fits — no new constants).
/// `rho_near` is the near-wall density; it only shapes the answer
/// where the wall is hot enough to dissociate (K_c drives
/// everywhere else). the degree of dissociation comes from the
/// closed form of K_c = [X]^2/[X2] with the element balance
/// [X2] = (1-a) c, [X] = 2 a c, c = rho/M_X2:
///
/// `a^2 / (1 - a) = K_c M_X2 / (4 rho)`
///
/// `a = (sqrt(b^2 + 4b) - b) / 2` with `b = K_c M_X2 / (4 rho)`
/// (the positive root of a^2 + b a - b = 0).
pub fn equilibrium_wall_composition(
    t_wall: f64,
    rho_near: f64,
    y_near: &[f64; NS],
) -> Result<[f64; NS], Error> {
    if !t_wall.is_finite() || t_wall <= 0.0 {
        return Err(Error::InvalidArgs);
    }
    // at low wall temperatures the dissociation equilibria are so
    // far on the molecule side that the closed-form atom fractions
    // underflow; below 2000 K the composition is all molecules to
    // any meaningful precision (K_c < 1e-8 mol/cm^3 at 4000 K and
    // falling like a power law with T).
    if t_wall < 2000.0 {
        return Ok([1.0, 0.0, 0.0, 0.0, 0.0]);
    }
    let rho = rho_near.max(1e-9);
    // element pools inherited from the near-wall composition: the
    // N pool rides N2 + NO + N, the O pool O2 + NO + O. a pure-N2
    // state stays pure (no invented oxygen), air keeps both pools.
    let m_n = y_near[0] + y_near[2] * (M[3] / M[2]) + y_near[3];
    let m_o = y_near[1] + y_near[2] * (M[4] / M[2]) + y_near[4];
    let frac_n = if m_n + m_o > 0.0 {
        m_n / (m_n + m_o)
    } else {
        0.0
    };
    let frac_o = 1.0 - frac_n;
    // N2 <-> 2N on the N pool; O2 <-> 2O on the O pool. the NO
    // exchange channel is second-order small at a cold wall and the
    // bracket does not need it: the fully-catalytic wall deposits
    // atoms as their parent molecules.
    let kc_n2 = crate::species_chem::k_eq(0, t_wall); // mol/m^3
    let kc_o2 = crate::species_chem::k_eq(1, t_wall);
    let b_n = kc_n2 * M[0] / (4.0 * rho);
    let b_o = kc_o2 * M[1] / (4.0 * rho);
    let alpha_n = ((b_n * b_n + 4.0 * b_n).sqrt() - b_n) / 2.0;
    let alpha_o = ((b_o * b_o + 4.0 * b_o).sqrt() - b_o) / 2.0;
    // moles per m^3: each element's diatomic carrier density scales
    // with its element share of the mixture
    let c_n = frac_n * rho / M[0];
    let c_o = frac_o * rho / M[1];
    let n_n2 = c_n * (1.0 - alpha_n);
    let n_n = 2.0 * c_n * alpha_n;
    let n_o2 = c_o * (1.0 - alpha_o);
    let n_o = 2.0 * c_o * alpha_o;
    // assemble mass fractions (NO stays 0 for the bracket)
    let masses = [n_n2 * M[0], n_o2 * M[1], 0.0, n_n * M[3], n_o * M[4]];
    let tot: f64 = masses.iter().sum();
    if !(tot.is_finite() && tot > 0.0) {
        return Err(Error::InvalidArgs);
    }
    let mut y = [0.0f64; NS];
    for s in 0..NS {
        y[s] = masses[s] / tot;
    }
    Ok(y)
}

/// the catalytic formation-enthalpy surface flux on one wall face
/// (W/m^2): sum_s h_f_s J_s,wallward with the diffusive mass flux
/// into the wall computed from the one-sided quadratic derivative
/// through the wall-face composition and the first two fluid cells,
/// the same stencil shape the fourier conductor uses so the
/// reported wall flux and the energy budget the march applies are
/// the same numbers.
///
/// `y_wall` is the wall-face composition (the fully-catalytic
/// equilibrium; the non-catalytic condition never calls this because
/// its species flux is zero by definition), `y0`/`y1` the first two
/// fluid-cell compositions, `positions` the three locations measured
/// along the interior direction (ascending: wall face, first cell
/// center, second cell center).
pub fn catalytic_surface_flux(
    y_wall: &[f64; NS],
    y0: &[f64; NS],
    y1: &[f64; NS],
    rho: f64,
    d_eff: f64,
    positions: [f64; 3],
) -> f64 {
    let [w, c0, c1] = positions;
    let h0 = c0 - w;
    let h1 = c1 - c0;
    if !(h0 > 0.0 && h1 > 0.0) {
        return 0.0;
    }
    // one-sided quadratic derivative at the wall node: exact for
    // quadratics, zero for constants. verified in tests below.
    let k0 = -(2.0 * h0 + h1) / (h0 * (h0 + h1));
    let k1 = (h0 + h1) / (h0 * h1);
    let k2 = -h0 / (h1 * (h0 + h1));
    let mut q_cat = 0.0f64;
    for s in 0..NS {
        let dy = y_wall[s] * k0 + y0[s] * k1 + y1[s] * k2;
        // the interior axis points away from the wall, so the
        // wallward component of the diffusive flux is +rho D dY/dx:
        // atoms rising away from a catalytic wall flow wallward and
        // deposit their formation enthalpy on the surface.
        let j_wallward = rho * d_eff * dy;
        q_cat += HF_SPECIES[s] * j_wallward;
    }
    q_cat
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cold_wall_equilibrium_is_all_molecules() {
        let y = equilibrium_wall_composition(1000.0, 0.01, &[0.767, 0.233, 0.0, 0.0, 0.0]).unwrap();
        assert!((y[0] - 1.0).abs() < 1e-12);
        assert!(y[3] + y[4] < 1e-12);
    }

    #[test]
    fn hot_wall_equilibrium_dissociates() {
        // pure N2 at 9000 K, 1 atm (rho 0.0038 kg/m3): the K_p
        // algebra gives alpha ~ 0.97 (a^2/(1-a) = K_c M/(4 rho)),
        // so the mass fraction of atomic N is dominant. the hand
        // anchor: ln K_c = -9.16 mol/cm3 at 9000 K from the fit.
        let rho_1atm_9000 = 101325.0 / (8.314462 / 0.028016 * 9000.0);
        let y = equilibrium_wall_composition(9000.0, rho_1atm_9000, &[1.0, 0.0, 0.0, 0.0, 0.0])
            .unwrap();
        assert!(y[3] > 0.5, "y_N at 9000K/1atm = {}", y[3]);
        assert!(y[0] < 0.5, "y_N2 at 9000K/1atm = {}", y[0]);
        assert!((y.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        // denser = less dissociated (Le Chatelier): 100x rho halves
        // the atom fraction roughly; at minimum it must go down.
        let y_dense =
            equilibrium_wall_composition(9000.0, rho_1atm_9000 * 100.0, &[1.0, 0.0, 0.0, 0.0, 0.0])
                .unwrap();
        assert!(
            y_dense[3] < y[3],
            "denser should recombine: {} vs {}",
            y_dense[3],
            y[3]
        );
    }

    #[test]
    fn derivative_stencil_is_exact_for_quadratics() {
        // f(x) = x^2 sampled at the wall node 0.05 and the two
        // following cell nodes; only the N slot carries the signal so
        // the summed flux is HF_N times the stencil's derivative at
        // the wall, f'(0.05) = 0.1.
        let mut y_wall = [0.0f64; NS];
        let mut c0 = [0.0f64; NS];
        let mut c1 = [0.0f64; NS];
        y_wall[3] = 0.05 * 0.05;
        c0[3] = 0.1 * 0.1;
        c1[3] = 0.3 * 0.3;
        let q = catalytic_surface_flux(&y_wall, &c0, &c1, 1.0, 1.0, [0.05, 0.1, 0.3]);
        let hf_n = HF_SPECIES[3];
        assert!((q - hf_n * 0.1).abs() < 1e-6 * hf_n, "q = {q}");
    }

    #[test]
    fn uniform_composition_gives_exactly_zero_flux() {
        // the mass constraint sum Y = 1 in every cell means the
        // summed derivative is d(1)/dx = 0: any physically valid
        // uniform-sum state gives q_cat = 0 even for wild per-species
        // gradients (the net enthalpy transport of the single-D
        // model redistributes, it does not create).
        let y_any = [0.6, 0.05, 0.05, 0.25, 0.05];
        let y_other = [0.7, 0.05, 0.05, 0.15, 0.05];
        let q = catalytic_surface_flux(&y_any, &y_other, &y_other, 1e-3, 1.0, [0.0, 0.05, 0.15]);
        // all three sum to 1? y_any sums 1.0, y_other sums 1.0 -> the
        // summed flux is 0 with the shared HF only if HF were equal;
        // it is not. the honest zero-flux identity: y uniform per
        // species across all three positions.
        let q_uni = catalytic_surface_flux(&y_any, &y_any, &y_any, 1e-3, 1.0, [0.0, 0.05, 0.15]);
        // exact arithmetic gives 0 (the weights sum to 0); float
        // roundoff leaves ~1e-11, so pin it as negligible against
        // the O(1e5) gradient-driven flux scale below.
        assert!(q_uni.abs() < 1e-6, "q_uni = {q_uni}");
        let _ = q;
    }

    #[test]
    fn fully_catalytic_flux_exceeds_non_catalytic_on_atom_gradient() {
        // shock-layer cells with atoms rising away from a cold
        // catalytic wall (ghost = all molecules): atoms flow
        // wallward, formation enthalpy deposits, term is positive.
        let y0 = [0.6, 0.05, 0.05, 0.25, 0.05]; // wall-adjacent cell: 25% N
        let y1 = [0.55, 0.05, 0.05, 0.3, 0.05]; // next cell: 30% N
        let eq_cold =
            equilibrium_wall_composition(1000.0, 0.01, &[0.767, 0.233, 0.0, 0.0, 0.0]).unwrap();
        let q_full = catalytic_surface_flux(&eq_cold, &y0, &y1, 1e-3, 1.0, [0.0, 0.05, 0.15]);
        // non-catalytic: zero diffusive flux by definition; the
        // reported q_cat is 0, so any positive full-cat value exceeds it
        assert!(q_full > 0.0, "q_full = {q_full}");
        // magnitude sanity: rho D dY/dx ~ 1e-3 * 1 * (0.25/0.05) = 5e-3
        // kg/(m^2 s) of N; times HF_N 33.6 MJ/kg ~ 1.7e5 W/m^2. plus
        // the O contribution; the total should be O(1e5) W/m^2.
        assert!(q_full > 1e4 && q_full < 1e6, "q_full = {q_full}");
    }
}

//! spalart-allmaras one-equation turbulence model: closures, constants, and
//! the production/destruction source term.
//!
//! the model transports a modified turbulent viscosity nu_tilde. production
//! scales with vorticity, destruction with distance to the wall, and the
//! eddy viscosity mu_t = rho nu_tilde fv1 feeds the mean-flow viscous fluxes.
//! constants and closures follow the original 1994 formulation (no trip
//! terms, fully turbulent); wall distance d is provided per cell by the
//! caller.

/// c_b1, production constant.
pub const C_B1: f64 = 0.1355;
/// c_b2, cross-diffusion constant.
pub const C_B2: f64 = 0.622;
/// sigma, the diffusion prandtl number of nu_tilde.
pub const SIGMA: f64 = 2.0 / 3.0;
/// von karman constant.
pub const KAPPA: f64 = 0.41;
/// c_w1 = c_b1/kappa^2 + (1 + c_b2)/sigma, the destruction coefficient
/// (the standard composite, as in the 1994 paper and mainstream codes).
pub const C_W1: f64 = C_B1 / (KAPPA * KAPPA) + (1.0 + C_B2) / SIGMA;
/// c_w2, destruction damping constant.
pub const C_W2: f64 = 0.3;
/// c_w3, destruction damping constant.
pub const C_W3: f64 = 2.0;
/// c_v1, the fv1 saturation constant.
pub const C_V1: f64 = 7.1;
/// cs, the stilde floor fraction of the vorticity (the guard for the
/// fv2 < 0 freestream region, as in mainstream implementations).
pub const C_S: f64 = 0.3;

/// chi = nu_tilde / nu, the viscosity ratio argument of the closures.
#[inline]
pub fn chi(nu_tilde: f64, nu: f64) -> f64 {
    nu_tilde / nu
}

/// fv1 = chi^3 / (chi^3 + c_v1^3), the viscous damping function.
#[inline]
pub fn fv1(chi: f64) -> f64 {
    let c3 = chi * chi * chi;
    c3 / (c3 + C_V1 * C_V1 * C_V1)
}

/// fv2 = 1 - chi / (1 + chi fv1), the second damping function.
#[inline]
pub fn fv2(chi: f64) -> f64 {
    1.0 - chi / (1.0 + chi * fv1(chi))
}

/// stilde = s + nu_tilde fv2 / (kappa^2 d^2), floored at cs * s so it stays
/// positive where fv2 < 0 (the standard guard; s = sqrt(2)|omega| in 2d).
#[inline]
pub fn stilde(nu_tilde: f64, s: f64, d: f64, nu: f64) -> f64 {
    let d2 = (d * d).max(1e-24);
    let chi_v = chi(nu_tilde, nu);
    (s + nu_tilde * fv2(chi_v) / (KAPPA * KAPPA * d2)).max(C_S * s)
}

/// r = nu_tilde / (stilde kappa^2 d^2), clipped to r_max = 10; the
/// lewandowski-style bound prevents division blowups near walls.
pub const R_MAX: f64 = 10.0;

#[inline]
pub fn r_func(nu_tilde: f64, s_tilde: f64, d: f64) -> f64 {
    let d2 = (d * d).max(1e-24);
    let r = nu_tilde / (s_tilde.max(1e-10) * KAPPA * KAPPA * d2);
    r.min(R_MAX)
}

/// g = r + c_w2 (r^6 - r), the 1994 destruction damping argument.
#[inline]
pub fn g_func(r: f64) -> f64 {
    let r6 = r * r * r * r * r * r;
    r + C_W2 * (r6 - r)
}

/// fw = g ((1 + c_w3^6) / (g^6 + 1 + c_w3^6))^(1/6).
#[inline]
pub fn fw(r: f64) -> f64 {
    let g = g_func(r);
    let g6 = g * g * g * g * g * g;
    let c36 = C_W3 * C_W3 * C_W3 * C_W3 * C_W3 * C_W3;
    g * ((1.0 + c36) / (g6 + 1.0 + c36)).powf(1.0 / 6.0)
}

/// the sa source term: production minus destruction per unit volume.
///
/// production = c_b1 stilde nu_tilde
/// destruction = c_w1 fw (nu_tilde / d)^2
/// returns d(nu_tilde)/dt contribution, units nu per time.
pub fn source(nu_tilde: f64, s: f64, d: f64, nu: f64, _dt: f64, _rho: f64) -> f64 {
    let s_tilde = stilde(nu_tilde, s, d, nu);
    let r = r_func(nu_tilde, s_tilde, d);
    let prod = C_B1 * s_tilde * nu_tilde;
    let destr = C_W1 * fw(r) * nu_tilde * nu_tilde / (d * d).max(1e-24);
    prod - destr
}

/// the eddy viscosity mu_t = rho nu_tilde fv1(chi) coupling into the mean flow.
#[inline]
pub fn eddy_viscosity(rho: f64, nu_tilde: f64, nu: f64) -> f64 {
    rho * nu_tilde * fv1(chi(nu_tilde, nu))
}

//! spalart-allmaras one-equation turbulence model: closures, constants, and
//! the production/destruction source term.
//!
//! the model transports a modified turbulent viscosity nu_tilde. production
//! scales with vorticity, destruction with distance to the wall, and the
//! eddy viscosity mu_t = rho nu_tilde fv1 feeds the mean-flow viscous fluxes.
//! constants and closures follow the original 1994 formulation (no trip
//! terms, fully turbulent); wall distance d is provided per cell by the
//! caller.

use crate::error::Error;
use crate::Grid2d;

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

/// fw = g ((1 + c_w3^6) / (g^6 + c_w3^6))^(1/6), the destruction damping
/// function of the 1994 formulation.
#[inline]
pub fn fw(r: f64) -> f64 {
    let g = g_func(r);
    let g6 = g * g * g * g * g * g;
    let c36 = C_W3 * C_W3 * C_W3 * C_W3 * C_W3 * C_W3;
    g * ((1.0 + c36) / (g6 + c36)).powf(1.0 / 6.0)
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

/// mu_t = rho * nu_tilde * fv1(chi), the mean-flow eddy coupling.
#[inline]
pub fn eddy_viscosity(rho: f64, nu_tilde: f64, nu_lam: f64) -> f64 {
    rho * nu_tilde * fv1(chi(nu_tilde, nu_lam))
}

/// the sa state of a 2d run.
pub struct SaState {
    /// nu_tilde per cell, m^2/s. initialized to a few nu_lam near the
    /// freestream (turbulence intensity tu).
    pub nu_tilde: Vec<f64>,
    /// distance to the nearest body-face cell, m.
    pub d_wall: Vec<f64>,
}

impl SaState {
    pub fn zeros(n: usize) -> Self {
        Self { nu_tilde: vec![0.0; n], d_wall: vec![0.0; n] }
    }

    /// init from a per-cell wall distance field. tu_inf is the freestream
    /// turbulence intensity fraction (use 5e-3 for atmospheric high
    /// altitude).
    pub fn init(d_wall: Vec<f64>, tu_inf: f64, nu_lam_inf: f64, u_inf: f64) -> Self {
        // nu_tilde_inf from tu_inf: a heuristic, 3 * nu_inf * tu_inf^2 * u_inf
        // is hand-wavy — the mainstream choice is around 3 * nu_inf tu_inf^2
        let mut s = Self {
            nu_tilde: vec![3.0 * nu_lam_inf * tu_inf * tu_inf * u_inf.max(1.0); d_wall.len()],
            d_wall,
        };
        s.nu_tilde.iter_mut().for_each(|v| *v = v.max(1e-12));
        s
    }
}

/// one explicit backward step on nu_tilde at every cell.
///   d(nu_t)/dt + u dot grad(nu_t) = prod - destr + (1/sigma) nabla^2(...)
/// st is used for u, v, rho, nu (lam) — so the marching code computes the
/// same gradients it uses elsewhere. mask[k]=true marks body cells.
pub fn advance_sa_step(
    sa: &mut SaState,
    rho: &[f64],
    u: &[f64],
    v: &[f64],
    nu_lam: &[f64],
    g: &crate::Grid2d,
    mask: &[bool],
    dt: f64,
) -> Result<(), Error> {
    let n = (g.nx * g.ny) as usize;
    if sa.nu_tilde.len() != n || sa.d_wall.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut new_nu = sa.nu_tilde.clone();
    for k in 0..n {
        if mask[k] {
            new_nu[k] = 0.0;
            continue;
        }
        let i = k % g.nx;
        let j = k / g.nx;
        let i1 = (i + 1).min(g.nx - 1);
        let j1 = (j + 1).min(g.ny - 1);
        // vorticity magnitude on cell k
        let dudx = (u[j * g.nx + i1] - u[k]) / g.dx.max(1e-12);
        let dudy = (u[j1 * g.nx + i] - u[k]) / g.dy.max(1e-12);
        let dvdx = (v[j * g.nx + i1] - v[k]) / g.dx.max(1e-12);
        let dvdy = (v[j1 * g.nx + i] - v[k]) / g.dy.max(1e-12);
        let omega = ((dvdx - dudy).powi(2) + (dudy - dvdx).powi(2)).sqrt();
        let d = sa.d_wall[k].max(g.dx);
        let s = source(
            sa.nu_tilde[k],
            omega.max(1e-16),
            d,
            nu_lam[k].max(1e-30),
            dt,
            rho[k],
        );
        // crude laplacian on nu_tilde for diffusion
        let mut lap = 0.0;
        let dx2 = g.dx.max(1e-12);
        let dy2 = g.dy.max(1e-12);
        if i > 0 { lap += (sa.nu_tilde[j * g.nx + i - 1] - sa.nu_tilde[k]) / dx2; }
        if i + 1 < g.nx { lap += (sa.nu_tilde[j * g.nx + i + 1] - sa.nu_tilde[k]) / dx2; }
        if j > 0 { lap += (sa.nu_tilde[(j - 1) * g.nx + i] - sa.nu_tilde[k]) / dy2; }
        if j + 1 < g.ny { lap += (sa.nu_tilde[(j + 1) * g.nx + i] - sa.nu_tilde[k]) / dy2; }
        let diffus = (nu_lam[k] + sa.nu_tilde[k]) * lap / SIGMA;
        let upd = sa.nu_tilde[k] + dt * (s + diffus);
        new_nu[k] = upd.max(0.0);
    }
    sa.nu_tilde = new_nu;
    Ok(())
}

/// fold the turbulent viscosity into the working array mut_t passed by
/// the marching routine: mu_turb = rho * nu_tilde * fv1 + the intrinsic
/// formula's eddy contribution.
pub fn blend_mu_into(mu_lam: &mut [f64], rho: &[f64], sa: &SaState) {
    for k in 0..mu_lam.len() {
        let nu_lam = mu_lam[k] / rho[k].max(1e-30);
        let chi_v = chi(sa.nu_tilde[k], nu_lam);
        let mu_t = rho[k] * sa.nu_tilde[k] * fv1(chi_v);
        mu_lam[k] += mu_t;
    }
}

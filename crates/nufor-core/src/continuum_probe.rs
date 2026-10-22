//! continuum-validity probe: bird's gradient-length local knudsen number.
//!
//! the navier-stokes/fourier closure assumes the flow is locally near
//! equilibrium. bird's criterion quantifies that with a knudsen number
//! built on the local gradient length instead of a global body scale:
//! the mean free path over the macroscopic length over which the
//! quantity actually changes,
//!
//!   kn_gll(q) = (lambda / q) * |grad q|
//!
//! with q the density (bird's original form) or the temperature. the
//! velocity-gradient variant needs the shear scale and stays with the
//! dsmc. continuum holds where kn_gll < 0.05; boyd et al. (phys.
//! fluids 7, 1995) and lofthouse & boyd (aiaa 2007-3905) use that
//! threshold to switch hybrid cfd/dsmc domains, and the same number is
//! reported here so a run states where its continuum assumption holds
//! and where it does not.
//!
//! the mean free path comes from the dynamic viscosity through the
//! hard-sphere relation: mu = rho * lambda * cbar / 2 with
//! cbar = sqrt(8 r t / pi), so
//!
//!   lambda = 2 mu / (rho * sqrt(8 r t / pi))
//!
//! which stays well-defined for any gas whose viscosity is known
//! (sutherland here) and needs no collision-diameter data.

use crate::grid2d::Grid2d;

/// hard-sphere mean free path from viscosity, density, temperature.
///
/// mu in pa s, rho in kg/m3, t in k, gas constant r in j/(kg k).
/// returns lambda in meters.
pub fn mean_free_path(mu: f64, rho: f64, t: f64, gas_const: f64) -> f64 {
    let cbar = (8.0 * gas_const * t.max(1e-12) / std::f64::consts::PI).sqrt();
    (2.0 * mu) / (rho.max(1e-30) * cbar)
}

/// bird's breakdown threshold: continuum holds below this.
pub const BREAKDOWN_KN: f64 = 0.05;

/// the probe result over the flow field.
#[derive(Debug, Clone)]
pub struct ContinuumProbe {
    /// gradient-length knudsen per cell, density form (bird's original).
    pub kn_rho: Vec<f64>,
    /// gradient-length knudsen per cell, temperature form.
    pub kn_t: Vec<f64>,
    /// mean free path per cell, m.
    pub lambda: Vec<f64>,
    /// global knudsen number on the reference length, freestream based.
    pub kn_inf: f64,
    /// max kn_gll(rho) over the probed cells.
    pub kn_rho_max: f64,
    /// max kn_gll(t) over the probed cells.
    pub kn_t_max: f64,
    /// fraction of probed cells above the breakdown threshold in
    /// either form.
    pub frac_breakdown: f64,
    /// cell index (row-major) of the worst rho-form cell.
    pub worst_cell: usize,
}

/// the center-to-center span (m) between two cells along one axis: the
/// sum of the cell widths between the two centers, so skipped solid
/// cells contribute their true width to the quotient.
fn skip_span(g: &Grid2d, lo: usize, hi: usize, axis: usize) -> f64 {
    let nx = g.nx;
    if axis == 0 {
        let (i0, i1) = (lo % nx, hi % nx);
        g.dxs[i0..=i1].iter().sum::<f64>() - 0.5 * (g.dxs[i0] + g.dxs[i1])
    } else {
        let (j0, j1) = (lo / nx, hi / nx);
        g.dys[j0..=j1].iter().sum::<f64>() - 0.5 * (g.dys[j0] + g.dys[j1])
    }
}

/// evaluate the gradient-length knudsen field over the flow cells.
///
/// inputs are the primitive fields on the grid (row-major), the cell
/// temperature, the dynamic viscosity (pa s, one value per cell), the
/// gas constant, the solid-cell mask, and the freestream reference
/// (rho, temperature, mu) for the global knudsen on `l_ref` (nose
/// radius, throat diameter, body length: the run's stated scale).
#[allow(clippy::too_many_arguments)]
pub fn continuum_probe(
    g: &Grid2d,
    rho: &[f64],
    t: &[f64],
    mu_cell: &[f64],
    gas_const: f64,
    solid: &[bool],
    rho_inf: f64,
    t_inf: f64,
    mu_inf: f64,
    l_ref: f64,
) -> ContinuumProbe {
    let (nx, ny) = (g.nx, g.ny);
    let n = nx * ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let mut kn_rho = vec![0.0f64; n];
    let mut kn_t = vec![0.0f64; n];
    let mut lambda = vec![0.0f64; n];
    let mut kn_rho_max = 0.0f64;
    let mut kn_t_max = 0.0f64;
    let mut worst_cell = 0usize;
    let mut n_probe = 0usize;
    let mut n_break = 0usize;
    // the neighbor pair along an axis, skipping solid cells; a
    // degenerate (lo == hi) pair means no usable gradient that way.
    let neighbors = |i: usize, j: usize, axis: usize| -> (usize, usize) {
        let k = idx(i, j);
        let mut lo = k;
        let mut hi = k;
        if axis == 0 {
            if i > 0 && !solid[idx(i - 1, j)] {
                lo = idx(i - 1, j);
            }
            if i + 1 < nx && !solid[idx(i + 1, j)] {
                hi = idx(i + 1, j);
            }
        } else {
            if j > 0 && !solid[idx(i, j - 1)] {
                lo = idx(i, j - 1);
            }
            if j + 1 < ny && !solid[idx(i, j + 1)] {
                hi = idx(i, j + 1);
            }
        }
        (lo, hi)
    };
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            if solid[k] {
                continue;
            }
            let grad_mag = |f: &[f64], axis: usize| -> f64 {
                let (lo, hi) = neighbors(i, j, axis);
                if lo == hi {
                    return 0.0;
                }
                let span = skip_span(g, lo, hi, axis).max(1e-30);
                ((f[hi] - f[lo]) / span).abs()
            };
            let lam = mean_free_path(mu_cell[k], rho[k], t[k], gas_const);
            lambda[k] = lam;
            let g_rho = (grad_mag(rho, 0).powi(2) + grad_mag(rho, 1).powi(2)).sqrt();
            let g_t = (grad_mag(t, 0).powi(2) + grad_mag(t, 1).powi(2)).sqrt();
            kn_rho[k] = lam * g_rho / rho[k].max(1e-30);
            kn_t[k] = lam * g_t / t[k].max(1.0);
            n_probe += 1;
            if kn_rho[k] > kn_rho_max {
                kn_rho_max = kn_rho[k];
                worst_cell = k;
            }
            if kn_t[k] > kn_t_max {
                kn_t_max = kn_t[k];
            }
            if kn_rho[k] > BREAKDOWN_KN || kn_t[k] > BREAKDOWN_KN {
                n_break += 1;
            }
        }
    }
    let lam_inf = mean_free_path(mu_inf, rho_inf, t_inf, gas_const);
    ContinuumProbe {
        kn_rho,
        kn_t,
        lambda,
        kn_inf: lam_inf / l_ref.max(1e-30),
        kn_rho_max,
        kn_t_max,
        frac_breakdown: if n_probe > 0 {
            n_break as f64 / n_probe as f64
        } else {
            0.0
        },
        worst_cell,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid2d::grid2d;

    fn test_grid(nx: usize, ny: usize, x1: f64, y1: f64) -> Grid2d {
        grid2d(nx, ny, 0.0, x1, 0.0, y1).expect("test grid")
    }

    /// uniform field: zero gradients, zero local kn, finite lambda.
    #[test]
    fn uniform_field_has_zero_kn() {
        let g = test_grid(8, 8, 1.0, 1.0);
        let n = 64;
        let rho = vec![0.01; n];
        let t = vec![500.0; n];
        let mu = vec![1.7e-5; n];
        let solid = vec![false; n];
        let p = continuum_probe(&g, &rho, &t, &mu, 287.0, &solid, 0.01, 500.0, 1.7e-5, 1.0);
        assert!(p.kn_rho.iter().all(|v| v.abs() < 1e-14));
        assert!(p.kn_t.iter().all(|v| v.abs() < 1e-14));
        assert!(p.lambda[0] > 0.0 && p.lambda[0].is_finite());
        assert!(p.kn_inf > 0.0 && p.kn_inf < 1.0);
        assert_eq!(p.frac_breakdown, 0.0);
    }

    /// a linear density ramp: kn_gll matches the hand value on an
    /// interior cell (gradient exactly the ramp slope).
    #[test]
    fn linear_ramp_matches_hand_value() {
        let g = test_grid(10, 4, 1.0, 1.0);
        let (nx, ny) = (g.nx, g.ny);
        let n = nx * ny;
        let slope = 1000.0;
        let mut rho = vec![0.0f64; n];
        for j in 0..ny {
            for i in 0..nx {
                rho[j * nx + i] = 100.0 + slope * (i as f64 + 0.5) * 0.1;
            }
        }
        let t = vec![300.0; n];
        let mu = vec![1.8e-5; n];
        let solid = vec![false; n];
        let p = continuum_probe(&g, &rho, &t, &mu, 287.0, &solid, 100.0, 300.0, 1.8e-5, 1.0);
        let k = 2 * nx + 5;
        let expected = p.lambda[k] * slope / rho[k];
        assert!((p.kn_rho[k] - expected).abs() / expected < 1e-12);
    }

    /// solid neighbors: the gradient skips them instead of reading the
    /// pinned solid value.
    #[test]
    fn solid_neighbor_is_skipped() {
        let g = test_grid(6, 4, 1.0, 1.0);
        let (nx, ny) = (g.nx, g.ny);
        let n = nx * ny;
        let mut rho = vec![1.0f64; n];
        for j in 0..ny {
            rho[j * nx] = 500.0; // left column solid-pinned junk values
        }
        let t = vec![300.0; n];
        let mu = vec![1.8e-5; n];
        let mut solid = vec![false; n];
        for j in 0..ny {
            solid[j * nx] = true;
        }
        let p = continuum_probe(&g, &rho, &t, &mu, 287.0, &solid, 1.0, 300.0, 1.8e-5, 1.0);
        // cell (1, j): west neighbor solid, east is 1.0 same value -> zero.
        let k = 2 * nx + 1;
        assert!(p.kn_rho[k].abs() < 1e-14);
        // and the solid cells themselves report zero.
        let ks = 2 * nx;
        assert_eq!(p.kn_rho[ks], 0.0);
    }

    /// the mean free path scales inversely with density at fixed t
    /// and as sqrt(t) at fixed density (through cbar).
    #[test]
    fn lambda_scales_with_density_and_sqrt_t() {
        let a = mean_free_path(1.7e-5, 1.0, 300.0, 287.0);
        let b = mean_free_path(1.7e-5, 2.0, 300.0, 287.0);
        assert!((a / b - 2.0).abs() < 1e-12);
        let c = mean_free_path(1.7e-5, 1.0, 1200.0, 287.0);
        // at fixed mu the molecular speed grows as sqrt(t), so the
        // mean free path shrinks as 1/sqrt(t).
        assert!((c / a - 0.5).abs() < 1e-12);
    }
}

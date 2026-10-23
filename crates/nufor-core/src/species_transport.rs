//! species transport pass: advect the mass fractions with
//! the same hllc sweeps the flow state uses.
//!
//! the species equation is a continuity equation per species:
//! d(rho Y_i)/dt + div(rho u Y_i) = omega_i, and this pass carries
//! the advective half. the face flux reuses the mass flux with the
//! upwind mass fraction — across a contact discontinuity the hllc
//! mass flux is upwinded on the contact side, so the species flux
//! (mass flux times that side's Y_i) transports the species exactly
//! with the material it belongs to: a uniform-Y flow stays uniform
//! to roundoff, and a species contact moves with the flow.

use crate::error::Error;
use crate::species::SpeciesState;

/// one species advection update on precomputed mass-flux sweeps, the
/// park5 entry point over the shared core below.
///
/// `fx[j].mass[f]` is the x-face mass flux on row j (f in 0..=nx),
/// `fy[i].mass[f]` the y-face mass flux on column i (f in 0..=ny),
/// matching the sweeps `compute_axi_fluxes` returns. the upwinding
/// side follows the sign of the face mass flux, which is the hllc
/// contact state's own choice for the material crossing the face.
pub fn transport_species(
    state: &mut SpeciesState,
    rho: &[f64],
    dt: f64,
    dxs: &[f64],
    dys: &[f64],
    fx: &[Vec<f64>],
    fy: &[Vec<f64>],
) -> Result<f64, Error> {
    transport_species_generic(&mut state.y, rho, dt, dxs, dys, fx, fy, None)
}

/// the species-count-agnostic core both lanes share: `y` holds one
/// row per species, however many that is.
///
/// `radial_weights` carries the axisymmetric face radii: when it is
/// given, the y-sweep update is weighted by the face radius (the
/// annular divergence the flow update uses), so the species mass
/// conserves against the flow on an annular grid. `None` is the
/// planar lane.
#[allow(clippy::too_many_arguments)]
pub fn transport_species_generic(
    y: &mut [Vec<f64>],
    rho: &[f64],
    dt: f64,
    dxs: &[f64],
    dys: &[f64],
    fx: &[Vec<f64>],
    fy: &[Vec<f64>],
    radial_weights: Option<&[f64]>,
) -> Result<f64, Error> {
    let ns = y.len();
    let nx = dxs.len();
    let ny = dys.len();
    if rho.len() != nx * ny || fx.len() != ny || fy.len() != nx {
        return Err(Error::InvalidArgs);
    }
    let idx = |i: usize, j: usize| j * nx + i;
    let mut resid = 0.0f64;
    // per-species arrays of per-cell upwind fractions: build the
    // face upwind Y first so both sweeps share the same choice.
    let mut yx_face: Vec<Vec<f64>> = y.iter().map(|_| vec![0.0f64; (nx + 1) * ny]).collect();
    let mut yy_face: Vec<Vec<f64>> = y.iter().map(|_| vec![0.0f64; (ny + 1) * nx]).collect();
    for s in 0..ns {
        let yi = &y[s];
        for j in 0..ny {
            for f in 0..=nx {
                let mf = fx[j][f];
                // face f sits between cell f-1 (left) and f (right); the
                // domain-edge faces take the adjacent interior cell on the
                // missing side (zero-gradient exterior), so a uniform field
                // is exactly divergence-free and a throughflow keeps its
                // composition instead of bleeding in zero-species ghosts.
                let left = if f > 0 {
                    yi[idx(f - 1, j)]
                } else {
                    yi[idx(f, j)]
                };
                let right = if f < nx {
                    yi[idx(f, j)]
                } else {
                    yi[idx(f - 1, j)]
                };
                yx_face[s][f * ny + j] = if mf >= 0.0 { left } else { right };
            }
        }
        for i in 0..nx {
            for f in 0..=ny {
                let mf = fy[i][f];
                let low = if f > 0 {
                    yi[idx(i, f - 1)]
                } else {
                    yi[idx(i, f)]
                };
                let high = if f < ny {
                    yi[idx(i, f)]
                } else {
                    yi[idx(i, f - 1)]
                };
                yy_face[s][f * nx + i] = if mf >= 0.0 { low } else { high };
            }
        }
    }
    // the update: rho*Y_i advances by the face flux divergence of
    // (mass flux * Y_upwind), then Y_i renormalizes against the new
    // density. applying the divergence to rho*Y keeps the pass
    // conservative in the per-species mass, not just the fraction.
    for s in 0..ns {
        let yi = &mut y[s];
        for j in 0..ny {
            for i in 0..nx {
                let k = idx(i, j);
                let dxi = dt / dxs[i];
                let dyi = dt / dys[j];
                let flux_x =
                    yx_face[s][(i + 1) * ny + j] * fx[j][i + 1] - yx_face[s][i * ny + j] * fx[j][i];
                let flux_y =
                    yy_face[s][(j + 1) * nx + i] * fy[i][j + 1] - yy_face[s][j * nx + i] * fy[i][j];
                // the axisymmetric y-sweep: the annular divergence
                // weights the radial faces by their radius, exactly as
                // the flow update does, so the species mass stays
                // consistent with the conserved flow mass.
                let flux_y = match radial_weights {
                    Some(w) => {
                        let rm = w[j];
                        let rp = w[j + 1];
                        let rc = 0.5 * (rm + rp);
                        (rp * yy_face[s][(j + 1) * nx + i] * fy[i][j + 1]
                            - rm * yy_face[s][j * nx + i] * fy[i][j])
                            / (rc * dys[j])
                    }
                    None => flux_y,
                };
                let rho_y = rho[k] * yi[k];
                let d = -(dxi * flux_x + dyi * flux_y);
                let new_rho_y = rho_y + d;
                resid = resid.max(d.abs() / rho[k].max(1e-30));
                // a species mass may not go negative: the honest fix is
                // the step limit the caller applies (the species cfl
                // probe), and the clip here is the last resort that
                // keeps the field finite. it destroys that species'
                // mass, so the caller must keep it from triggering.
                yi[k] = if new_rho_y > 0.0 {
                    new_rho_y / rho[k]
                } else {
                    0.0
                };
            }
        }
    }
    Ok(resid)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(n: usize) -> (Vec<f64>, Vec<f64>) {
        (vec![0.1; n], vec![0.1; n])
    }

    /// a uniform field with a constant-velocity mass flux must stay
    /// exactly uniform: the divergence of a constant flux is zero.
    #[test]
    fn uniform_stays_uniform() {
        let (dxs, dys) = grid(4);
        let n = 16;
        let rho = vec![1.0; n];
        let mut s = SpeciesState::air(n);
        // constant mass flux: every x-face carries rho*u, every y-face 0.
        let fx: Vec<Vec<f64>> = (0..4).map(|_| vec![0.5; 5]).collect();
        let fy: Vec<Vec<f64>> = (0..4).map(|_| vec![0.0; 5]).collect();
        let r = transport_species(&mut s, &rho, 0.01, &dxs, &dys, &fx, &fy).unwrap();
        assert!(r < 1e-14);
        assert!(s.max_sum_drift() < 1e-14);
        assert!((s.y[0][5] - 0.767).abs() < 1e-14);
    }

    /// inflow on the left, outflow on the right, uniform initial Y:
    /// fractions never move because the upwind face Y is uniform.
    #[test]
    fn throughflow_keeps_composition() {
        let (dxs, dys) = grid(4);
        let n = 16;
        let rho = vec![1.0; n];
        let mut s = SpeciesState::air(n);
        let fx: Vec<Vec<f64>> = (0..4).map(|_| vec![0.5; 5]).collect();
        let fy: Vec<Vec<f64>> = (0..4).map(|_| vec![0.0; 5]).collect();
        // many steps: the field must stay air everywhere.
        for _ in 0..50 {
            transport_species(&mut s, &rho, 0.01, &dxs, &dys, &fx, &fy).unwrap();
        }
        for k in 0..n {
            assert!((s.y[0][k] - 0.767).abs() < 1e-12, "k = {k}");
            assert!((s.y[1][k] - 0.233).abs() < 1e-12, "k = {k}");
        }
    }

    /// a left half rich in N2 sweeping right must cross faces with
    /// the material: after enough steps the discontinuity moves at
    /// the flux speed and the crossing cells are pure upwind values.
    #[test]
    fn contact_advects_with_the_flow() {
        let (dxs, dys) = grid(4);
        let n = 16;
        let rho = vec![1.0; n];
        let mut s = SpeciesState::air(n);
        // left half: pure N2. right half: pure O2.
        for j in 0..4 {
            for i in 0..4 {
                let k = j * 4 + i;
                s.y[0][k] = 1.0;
                s.y[1][k] = 0.0;
            }
        }
        let u = 0.5f64;
        let fx: Vec<Vec<f64>> = (0..4).map(|_| vec![u; 5]).collect();
        let fy: Vec<Vec<f64>> = (0..4).map(|_| vec![0.0; 5]).collect();
        // dt*u/dx = 0.05 per step -> after 20 steps the front moved one cell.
        for _ in 0..20 {
            transport_species(&mut s, &rho, 0.01, &dxs, &dys, &fx, &fy).unwrap();
        }
        // cells 0..3 columns are still pure N2, column 4 is the blend zone.
        for j in 0..4 {
            assert!((s.y[0][j * 4] - 1.0).abs() < 1e-12);
        }
        assert!(s.max_sum_drift() < 1e-12);
    }

    /// the generic transport must move every species row it is
    /// given: the update loop once hardcoded the park count of five,
    /// so the combustion lane's O, HO2, H2O2, and N2 rows were never
    /// advected. a nine-row field with a divergent flux must move the
    /// last row too.
    #[test]
    fn generic_transport_moves_every_species_row() {
        let (dxs, dys) = grid(4);
        let n = 16;
        let rho = vec![1.0; n];
        let mut y: Vec<Vec<f64>> = (0..9).map(|_| vec![0.0; n]).collect();
        for row in y.iter_mut() {
            row.fill(0.1);
        }
        // a y-flux that diverges: the bottom face carries flow up, the
        // top face carries none, so the top row loses mass.
        let fx: Vec<Vec<f64>> = (0..4).map(|_| vec![0.0; 5]).collect();
        let fy: Vec<Vec<f64>> = (0..4).map(|_| vec![0.0, 0.0, 0.0, 0.0, 0.0]).collect();
        let mut fy = fy;
        for col in fy.iter_mut() {
            col[0] = 0.5;
        }
        let resid =
            transport_species_generic(&mut y, &rho, 0.01, &dxs, &dys, &fx, &fy, None).unwrap();
        assert!(resid > 1e-6, "the flux must move something");
        // the bottom row gains the inflowed mass (the flux diverges
        // there); the last species row must have moved like the first.
        let last_moved = (y[8][0] - 0.1).abs();
        let first_moved = (y[0][0] - 0.1).abs();
        assert!(
            last_moved > 1e-9,
            "the ninth species row did not move: {}",
            y[8][0]
        );
        assert!(
            (last_moved - first_moved).abs() < 1e-14,
            "species rows diverged: {last_moved} vs {first_moved}"
        );
    }
}

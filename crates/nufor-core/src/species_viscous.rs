//! the park5 laminar viscous pass: species-aware stresses, heat
//! conduction and species diffusion, plus the catalytic wall energy
//! flux on body faces. separate from the shared conductor
//! (viscous2d) because every transport coefficient there is a
//! scalar while park5 needs the per-cell composition through
//! wilke's rule; the stencil shapes (one-sided quadratic at walls,
//! centered face interpolation in the interior) are the same so the
//! two lanes report comparable wall fluxes.
//!
//! structure per rk2 stage:
//! 1. per-cell (p, T, mu_wilke, k_wilke, D_eff) from the closure
//!    temperature and the composition.
//! 2. x/y face sweeps: stress + fourier + species-diffusion fluxes;
//!    body wall faces get the no-slip one-sided stencil, the
//!    isothermal conductor and the catalytic formation-enthalpy term.
//! 3. the divergence update on the flow state and the composition.
//!
//! the species diffusion uses the single effective diffusivity
//! (wilke_transport.md): J_s = -rho D dY_s/dx_i with sum_s J_s = 0
//! identically. the energy flux carries sum_s h_f_s J_s
//! (formation enthalpy of the diffusive transport), and the wall
//! face adds the catalytic term: the atoms arriving at a
//! fully-catalytic surface recombine and deposit their formation
//! enthalpy, which is exactly the same enthalpy flux evaluated at
//! the wall ghost composition.
//!
//! sources: plans/catalytic_research/wilke_transport.md (RP-1232
//! tables, LAURA TM-4674 / DPLR wall definitions).

use crate::catalytic::{catalytic_surface_flux, equilibrium_wall_composition, CatalyticWall};
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::park::NS;
use crate::species::SpeciesState;
use crate::species_closure as closure5;
use crate::species_props::{effective_diffusivity, mixture_conductivity, mixture_viscosity};
use crate::state2d::ConservedState2d;

/// per-cell transport snapshot: (T, mu, k, D) at the cell's own state.
pub(crate) struct CellTransport {
    pub t: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub d: Vec<f64>,
}

/// the per-cell closure temperature and wilke transport coefficients.
/// temperature comes from the closure's own inversion (e depends on
/// Y); mu and k are wilke mixtures, D the effective binary value.
pub(crate) fn cell_transport(
    state: &ConservedState2d,
    species: &SpeciesState,
    g: &Grid2d,
) -> Result<CellTransport, Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut t = vec![0.0f64; n];
    let mut mu = vec![0.0f64; n];
    let mut k = vec![0.0f64; n];
    let mut d = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    for k_idx in 0..n {
        let e_int = et[k_idx] - 0.5 * (u[k_idx] * u[k_idx] + v[k_idx] * v[k_idx]);
        let mut y: [f64; NS] = [0.0; NS];
        for (ys, ysrc) in y.iter_mut().zip(species.y.iter()) {
            *ys = ysrc[k_idx];
        }
        t[k_idx] = closure5::temperature_at(e_int, &y);
        p[k_idx] = closure5::pressure_at(state.rho[k_idx], e_int, &y);
        mu[k_idx] = mixture_viscosity(&y, t[k_idx]);
        k[k_idx] = mixture_conductivity(&y, t[k_idx]);
        d[k_idx] = effective_diffusivity(&y, t[k_idx], p[k_idx].max(1.0));
    }
    Ok(CellTransport { t, mu, k, d })
}

/// the wall treatment selected for body faces.
#[derive(Clone, Copy)]
pub(crate) struct WallOpts {
    pub temperature: f64,
    pub catalytic: CatalyticWall,
}

/// one viscous application on (state, species) with the given dt.
/// `solid` is the body SDF (empty = body-free); wall faces are fluid
/// cells with a solid neighbor, the same classification the shared
/// conductor uses.
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_species_viscous(
    state: &mut ConservedState2d,
    species: &mut SpeciesState,
    g: &Grid2d,
    ct: &CellTransport,
    dt: f64,
    bc: &crate::solver2d::Boundaries2d,
    wall_opts: WallOpts,
    solid: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    // solid classification, same cache as the euler zeroing uses.
    let solid_mask: Vec<bool> = match solid {
        Some(dist) => crate::body::solid_mask(g, dist).to_vec(),
        None => Vec::new(),
    };
    let has_body = !solid_mask.is_empty();
    let is_solid =
        |i: usize, j: usize| -> bool { i < nx && j < ny && has_body && solid_mask[j * nx + i] };
    let _ = bc;
    // interior face fluxes, then the divergence update. face arrays
    // follow the conductor's layout: x faces (nx+1)*ny face-major,
    // y faces nx*(ny+1) column-major.
    let mut fxm = vec![0.0f64; (nx + 1) * ny];
    let mut fym = vec![0.0f64; (nx + 1) * ny];
    let mut fe = vec![0.0f64; (nx + 1) * ny];
    // species face fluxes: rho*D*dY/dn per species per x face.
    let mut fdx: Vec<Vec<f64>> = species
        .y
        .iter()
        .map(|_| vec![0.0f64; (nx + 1) * ny])
        .collect();
    let mut gym_x = vec![0.0f64; nx * (ny + 1)];
    let mut gym_y = vec![0.0f64; nx * (ny + 1)];
    let mut gye = vec![0.0f64; nx * (ny + 1)];
    let mut fdy: Vec<Vec<f64>> = species
        .y
        .iter()
        .map(|_| vec![0.0f64; nx * (ny + 1)])
        .collect();
    // the wall composition arrays, filled during the face sweeps.

    // x faces
    for j in 0..ny {
        for f in 1..nx {
            let (a, b) = (idx(f - 1, j), idx(f, j));
            let k = f * ny + j;
            if is_solid(f - 1, j) || is_solid(f, j) {
                // body wall face along x
                let (left_solid, right_solid) = (is_solid(f - 1, j), is_solid(f, j));
                if left_solid && right_solid {
                    continue;
                }
                let vals = wall_face_x(
                    left_solid, f, j, &u, &v, ct, species, &state.rho, nx, g, wall_opts,
                );
                let (txx, txy, q_cond, q_cat, _j_wall) = vals;
                fxm[k] = txx;
                fym[k] = txy;
                fe[k] = q_cond + q_cat;
                continue;
            }
            // centered interior face: average coefficients, then the
            // standard stress + conduction + species diffusion.
            let mu_f = 0.5 * (ct.mu[a] + ct.mu[b]);
            let k_f = 0.5 * (ct.k[a] + ct.k[b]);
            let d_f = 0.5 * (ct.d[a] + ct.d[b]);
            let (ca, cb) = (g.centers_x[a], g.centers_x[b]);
            let dx = cb - ca;
            let dudx = (u[b] - u[a]) / dx;
            let dvdx = (v[b] - v[a]) / dx;
            let dtdx = (ct.t[b] - ct.t[a]) / dx;
            let div = dudx + 0.0;
            let txx = mu_f * (2.0 * dudx - 2.0 / 3.0 * div);
            let txy = mu_f * dvdx;
            // the rhs flux stores +k grad(T) (the conductor's
            // convention: the divergence update subtracts it, so a
            // hot peak cools).
            let qx = k_f * dtdx;
            // species diffusion flux at the face (W/m^2 enthalpy part
            // folded into fe below)
            let rho_f = 0.5 * (state.rho[a] + state.rho[b]);
            let mut q_enth = 0.0f64;
            for (s, ys) in species.y.iter().enumerate() {
                let dydx = (ys[b] - ys[a]) / dx;
                let js = -rho_f * d_f * dydx;
                fdx[s][k] = js;
                q_enth += closure5::HF_SPECIES[s] * js;
            }
            // energy rhs flux: stress work + the folded heat/enthalpy
            // terms (each the negative of the physical flux, matching
            // the update's subtracted difference).
            fe[k] = 0.5 * (u[a] + u[b]) * txx + 0.5 * (v[a] + v[b]) * txy + qx - q_enth;
            fxm[k] = txx;
            fym[k] = txy;
        }
    }
    // y faces; the f=0 south and f=ny north domain faces carry no-slip
    // wall fluxes when the boundary asks for them (the shared
    // conductor's treatment, with the catalytic term added).
    for i in 0..nx {
        for f in 0..=ny {
            if f == 0 || f == ny {
                let (a, b, flip) = if f == 0 {
                    (idx(i, 0), idx(i, 1), 1.0)
                } else {
                    (idx(i, ny - 1), idx(i, ny - 2), -1.0)
                };
                let no_slip = if f == 0 {
                    matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall)
                } else {
                    matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall)
                };
                if !no_slip {
                    continue; // free face: zero viscous flux
                }
                let k = i * (ny + 1) + f;
                let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
                let face_y = if f == 0 { g.faces_y[0] } else { g.faces_y[ny] };
                let mu_f = ct.mu[a];
                let dudy = flip * quad_deriv_at_wall(0.0, u[a], u[b], face_y, c0, c1);
                let dvdy = flip * quad_deriv_at_wall(0.0, v[a], v[b], face_y, c0, c1);
                let txy = mu_f * dudy;
                let tyy = mu_f * dvdy;
                let mut q = 0.0;
                if wall_opts.temperature > 0.0 {
                    // rhs-flux convention: +k grad(T) on true
                    // positions (the derivative is direction-free and
                    // the divergence update supplies the orientation).
                    q = ct.k[a]
                        * quad_deriv_at_wall(
                            wall_opts.temperature,
                            ct.t[a],
                            ct.t[b],
                            face_y,
                            c0,
                            c1,
                        );
                }
                if wall_opts.catalytic == CatalyticWall::FullyCatalytic {
                    let mut y0 = [0.0f64; NS];
                    let mut y1v = [0.0f64; NS];
                    for s in 0..NS {
                        y0[s] = species.y[s][a];
                        y1v[s] = species.y[s][b];
                    }
                    let y_wall = crate::catalytic::equilibrium_wall_composition(
                        wall_opts.temperature.max(200.0),
                        state.rho[a],
                        &y0,
                    )
                    .unwrap_or(y0);
                    // positions ascend away from the wall into the
                    // interior on both faces: the south wall's cells
                    // sit above the face, the north's below it, so the
                    // offsets take the absolute interior distance and
                    // the cell order follows.
                    let d0 = (c0 - face_y).abs();
                    let d1 = (c1 - face_y).abs();
                    let (pa, pb, ya, yb) = if d0 <= d1 {
                        (d0, d1, y0, y1v)
                    } else {
                        (d1, d0, y1v, y0)
                    };
                    // the surface flux is the wallward heating (positive
                    // when atoms drain into the wall); the stored rhs
                    // flux is the +y-direction entry, which is the
                    // wallward direction on the north face and its
                    // opposite on the south.
                    let q_wall = crate::catalytic::catalytic_surface_flux(
                        &y_wall,
                        &ya,
                        &yb,
                        state.rho[a],
                        ct.d[a],
                        [0.0, pa, pb],
                    );
                    q += -flip * q_wall;
                }
                gym_x[k] = txy;
                gym_y[k] = tyy;
                gye[k] = u[a] * txy + v[a] * tyy + q;
                // the species wall-face flux uses the wall composition
                // (fully-catalytic: the equilibrium ghost; non-cat:
                // zero-gradient, so the plain central difference
                // between the cells gives zero net flux through the
                // wall face into the ghost-free boundary).
                let y_face_wall: Option<[f64; NS]> =
                    if wall_opts.catalytic == CatalyticWall::FullyCatalytic {
                        let mut y0 = [0.0f64; NS];
                        for (ys, ysrc) in y0.iter_mut().zip(species.y.iter()) {
                            *ys = ysrc[a];
                        }
                        crate::catalytic::equilibrium_wall_composition(
                            wall_opts.temperature.max(200.0),
                            state.rho[a],
                            &y0,
                        )
                        .ok()
                    } else {
                        None
                    };
                let mut q_enth = 0.0;
                for s in 0..NS {
                    let js = match &y_face_wall {
                        Some(yw) => {
                            // one-sided wallward flux with the ghost
                            // composition, ascending interior offsets
                            let d0 = (c0 - face_y).abs();
                            let d1 = (c1 - face_y).abs();
                            let (pa, pb, ya, yb) = if d0 <= d1 {
                                (d0, d1, species.y[s][a], species.y[s][b])
                            } else {
                                (d1, d0, species.y[s][b], species.y[s][a])
                            };
                            // derivative at the wall node (0.0)
                            let h0 = pa;
                            let h1 = pb - pa;
                            let k0 = -(2.0 * h0 + h1) / (h0 * (h0 + h1));
                            let k1 = (h0 + h1) / (h0 * h1);
                            let k2 = -h0 / (h1 * (h0 + h1));
                            let dy = yw[s] * k0 + ya * k1 + yb * k2;
                            // the stored face flux is the +y-direction
                            // component of the diffusive transport.
                            // the ascending interior axis points away
                            // from the wall (up for the south face,
                            // down for the north), so the +y sign
                            // relative to that axis flips with the
                            // face: south (interior above) J = -rho D
                            // dy; north (interior below) J = +rho D dy.
                            // atoms draining wallward leave the domain
                            // in both cases.
                            let sign = if f == 0 { -1.0 } else { 1.0 };
                            sign * state.rho[a] * ct.d[a] * dy
                        }
                        None => 0.0,
                    };
                    fdy[s][k] = js;
                    q_enth += closure5::HF_SPECIES[s] * js;
                }
                gye[k] += q_enth;
                continue;
            }
            let (a, b) = (idx(i, f - 1), idx(i, f));
            let k = i * (ny + 1) + f;
            if is_solid(i, f - 1) || is_solid(i, f) {
                let (below_solid, above_solid) = (is_solid(i, f - 1), is_solid(i, f));
                if below_solid && above_solid {
                    continue;
                }
                let (txy, tyy, q_cond, q_cat, _j_wall) = wall_face_y(
                    below_solid,
                    i,
                    f,
                    &u,
                    &v,
                    ct,
                    species,
                    &state.rho,
                    nx,
                    g,
                    wall_opts,
                );
                gym_x[k] = txy;
                gym_y[k] = tyy;
                gye[k] = q_cond + q_cat;
                continue;
            }
            let mu_f = 0.5 * (ct.mu[a] + ct.mu[b]);
            let k_f = 0.5 * (ct.k[a] + ct.k[b]);
            let d_f = 0.5 * (ct.d[a] + ct.d[b]);
            let (ca, cb) = (g.centers_y[a], g.centers_y[b]);
            let dy = cb - ca;
            let dudy = (u[b] - u[a]) / dy;
            let dvdy = (v[b] - v[a]) / dy;
            let dtdy = (ct.t[b] - ct.t[a]) / dy;
            let div = 0.0 + dvdy;
            let tyy = mu_f * (2.0 * dvdy - 2.0 / 3.0 * div);
            let txy = mu_f * dudy;
            // rhs-flux convention: +k grad(T) along y (see the x sweep).
            let qy = k_f * dtdy;
            let rho_f = 0.5 * (state.rho[a] + state.rho[b]);
            let mut q_enth = 0.0f64;
            for (s, ys) in species.y.iter().enumerate() {
                let dydy = (ys[b] - ys[a]) / dy;
                let js = -rho_f * d_f * dydy;
                fdy[s][k] = js;
                q_enth += closure5::HF_SPECIES[s] * js;
            }
            gye[k] = 0.5 * (u[a] + u[b]) * txy + 0.5 * (v[a] + v[b]) * tyy + qy - q_enth;
            gym_x[k] = txy;
            gym_y[k] = tyy;
        }
    }
    // per-cell species-diffusion dt: the explicit diffusion number
    // per cell is dt < dx^2/(4 D/rho) locally. a global min over the
    // field lets the near-vacuum freestream (huge D/rho, uniform Y)
    // crush the march's step to nothing; the face-local pairing only
    // slows the cells the diffusion actually stiffens (the shock
    // layer), the same per-cell pairing the momentum bound uses.
    let cell_dt = |k: usize| -> f64 {
        let i = k % nx;
        let j = k / nx;
        let rho_c = state.rho[k].max(1e-9);
        let ddx = g.dxs[i] * g.dxs[i] / (4.0 * ct.d[k].max(1e-12) / rho_c);
        let ddy = g.dys[j] * g.dys[j] / (4.0 * ct.d[k].max(1e-12) / rho_c);
        ddx.min(ddy).min(dt)
    };
    // divergence update: flow state and species together.
    for j in 0..ny {
        for i in 0..nx {
            let c = idx(i, j);
            if is_solid(i, j) {
                continue;
            }
            // x faces: the face-local dt gates the whole viscous
            // update (stress, conduction, diffusion) so the near-
            // vacuum freestream cannot destabilize the shock layer
            // and vice versa; the pairing matches the species dt.
            let (fa, fb) = (i * ny + j, (i + 1) * ny + j);
            let dtx = cell_dt(c).min(if i + 1 < nx {
                cell_dt(idx(i + 1, j))
            } else {
                cell_dt(c)
            });
            state.mx[c] += dtx * (fxm[fb] - fxm[fa]) / g.dxs[i];
            state.my[c] += dtx * (fym[fb] - fym[fa]) / g.dxs[i];
            state.e[c] += dtx * (fe[fb] - fe[fa]) / g.dxs[i];
            for (s, ys) in species.y.iter_mut().enumerate() {
                // d(rho Y_s)/dt = -d(J_s)/dx: J already carries the
                // sign so the update subtracts the face difference of
                // J the same way the stress update does.
                ys[c] += -dtx * (fdx[s][fb] - fdx[s][fa]) / g.dxs[i] / state.rho[c].max(1e-9);
            }
            // y faces: same face-local pairing along y.
            let (da, db) = (i * (ny + 1) + j, i * (ny + 1) + j + 1);
            let dty = cell_dt(c).min(if j + 1 < ny {
                cell_dt(idx(i, j + 1))
            } else {
                cell_dt(c)
            });
            state.mx[c] += dty * (gym_x[db] - gym_x[da]) / g.dys[j];
            state.my[c] += dty * (gym_y[db] - gym_y[da]) / g.dys[j];
            state.e[c] += dty * (gye[db] - gye[da]) / g.dys[j];
            for (s, ys) in species.y.iter_mut().enumerate() {
                ys[c] += -dty * (fdy[s][db] - fdy[s][da]) / g.dys[j] / state.rho[c].max(1e-9);
            }
        }
    }
    // renormalize after diffusion (the single-D fluxes conserve sum
    // Y analytically; the discretization leaves O(eps) drift and the
    // overshoot clamp is cheap).
    for k in 0..nx * ny {
        if !is_solid(k % nx, k / nx) {
            species.sanitize(k);
        }
    }
    Ok(())
}

/// the no-slip isothermal + catalytic wall flux on one body x face.
/// mirrors body_wall_x's stencil: the one-sided quadratic through
/// the wall value and the first two fluid cells, exact for
/// quadratics.
#[allow(clippy::too_many_arguments)]
fn wall_face_x(
    left_solid: bool,
    f: usize,
    j: usize,
    u: &[f64],
    v: &[f64],
    ct: &CellTransport,
    species: &SpeciesState,
    rho: &[f64],
    nx: usize,
    g: &Grid2d,
    wall_opts: WallOpts,
) -> (f64, f64, f64, f64, f64) {
    let idx = |i: usize, jj: usize| jj * nx + i;
    // the fluid side is (f, j) when the left is solid, else (f-1, j)
    let (a, b, ok) = if left_solid {
        (idx(f, j), idx(f + 1, j), f + 1 < nx)
    } else {
        (idx(f - 1, j), idx(f - 2, j), f >= 2)
    };
    if !ok {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }
    let face_x = g.faces_x[f];
    let (c0, c1) = (g.centers_x[a], g.centers_x[b]);
    let mu_f = ct.mu[a];
    // no-slip: u = v = 0 at the wall
    let dudx = quad_deriv_at_wall(0.0, u[a], u[b], face_x, c0, c1);
    let dvdx = quad_deriv_at_wall(0.0, v[a], v[b], face_x, c0, c1);
    let txx = mu_f * dudx;
    let txy = mu_f * dvdx;
    let mut q_cond = 0.0;
    if wall_opts.temperature > 0.0 {
        let dtdn = quad_deriv_at_wall(wall_opts.temperature, ct.t[a], ct.t[b], face_x, c0, c1);
        // rhs-flux convention: +k grad(T) along the interior-ascending
        // axis (the conductor's q_wall), so the hot near-wall cell
        // loses energy to the cold wall.
        q_cond = ct.k[a] * dtdn;
    }
    // catalytic part: the wallward formation-enthalpy flux is energy
    // LEAVING the fluid, so the rhs flux stores its negative.
    let (q_cat, j_wall) = match wall_opts.catalytic {
        CatalyticWall::NonCatalytic => (0.0, 0.0),
        CatalyticWall::FullyCatalytic => {
            let mut y0 = [0.0f64; NS];
            let mut y1 = [0.0f64; NS];
            for (ys0, ys_src) in y0.iter_mut().zip(species.y.iter()) {
                *ys0 = ys_src[a];
            }
            for (ys1, ys_src) in y1.iter_mut().zip(species.y.iter()) {
                *ys1 = ys_src[b];
            }
            let y_wall =
                equilibrium_wall_composition(wall_opts.temperature.max(200.0), rho[a], &y0)
                    .unwrap_or(y0);
            // positions ascend away from the wall into the fluid
            let d0 = (c0 - face_x).abs();
            let d1 = (c1 - face_x).abs();
            let q = if d0 <= d1 {
                catalytic_surface_flux(&y_wall, &y0, &y1, rho[a], ct.d[a], [0.0, d0, d1])
            } else {
                catalytic_surface_flux(&y_wall, &y1, &y0, rho[a], ct.d[a], [0.0, d1, d0])
            };
            (-q, 0.0)
        }
    };
    (txx, txy, q_cond, q_cat, j_wall)
}

/// mirror of wall_face_x along y.
#[allow(clippy::too_many_arguments)]
fn wall_face_y(
    below_solid: bool,
    i: usize,
    f: usize,
    u: &[f64],
    v: &[f64],
    ct: &CellTransport,
    species: &SpeciesState,
    rho: &[f64],
    nx: usize,
    g: &Grid2d,
    wall_opts: WallOpts,
) -> (f64, f64, f64, f64, f64) {
    let idx = |ii: usize, jj: usize| jj * nx + ii;
    let (a, b, ok) = if below_solid {
        (idx(i, f), idx(i, f + 1), f + 1 < g.ny)
    } else {
        (idx(i, f - 1), idx(i, f - 2), f >= 2)
    };
    if !ok {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }
    let face_y = g.faces_y[f];
    let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
    let mu_f = ct.mu[a];
    let dudy = quad_deriv_at_wall(0.0, u[a], u[b], face_y, c0, c1);
    let dvdy = quad_deriv_at_wall(0.0, v[a], v[b], face_y, c0, c1);
    let txy = mu_f * dudy;
    let tyy = mu_f * dvdy;
    let mut q_cond = 0.0;
    if wall_opts.temperature > 0.0 {
        let dtdn = quad_deriv_at_wall(wall_opts.temperature, ct.t[a], ct.t[b], face_y, c0, c1);
        // rhs-flux convention: +k grad(T) along the interior-ascending
        // axis (the conductor's q_wall).
        q_cond = ct.k[a] * dtdn;
    }
    let (q_cat, j_wall) = match wall_opts.catalytic {
        CatalyticWall::NonCatalytic => (0.0, 0.0),
        CatalyticWall::FullyCatalytic => {
            let mut y0 = [0.0f64; NS];
            let mut y1 = [0.0f64; NS];
            for (ys0, ys_src) in y0.iter_mut().zip(species.y.iter()) {
                *ys0 = ys_src[a];
            }
            for (ys1, ys_src) in y1.iter_mut().zip(species.y.iter()) {
                *ys1 = ys_src[b];
            }
            let y_wall =
                equilibrium_wall_composition(wall_opts.temperature.max(200.0), rho[a], &y0)
                    .unwrap_or(y0);
            // positions ascend away from the wall into the fluid
            let d0 = (c0 - face_y).abs();
            let d1 = (c1 - face_y).abs();
            let q = if d0 <= d1 {
                catalytic_surface_flux(&y_wall, &y0, &y1, rho[a], ct.d[a], [0.0, d0, d1])
            } else {
                catalytic_surface_flux(&y_wall, &y1, &y0, rho[a], ct.d[a], [0.0, d1, d0])
            };
            // the wallward formation-enthalpy flux leaves the fluid:
            // the rhs flux stores its negative.
            (-q, 0.0)
        }
    };
    (txy, tyy, q_cond, q_cat, j_wall)
}

/// the one-sided quadratic derivative at the wall through (wall
/// value, first fluid, second fluid): same form as the conductor's
/// quad_deriv_at, kept local so this module carries no private deps
/// into the shared lane.
fn quad_deriv_at_wall(v_wall: f64, v0: f64, v1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    let w0 = (2.0 * y_wall - c0 - c1) / ((y_wall - c0) * (y_wall - c1));
    let w1 = (2.0 * y_wall - y_wall - c1) / ((c0 - y_wall) * (c0 - c1));
    let w2 = (2.0 * y_wall - y_wall - c0) / ((c1 - y_wall) * (c1 - c0));
    v_wall * w0 + v0 * w1 + v1 * w2
}

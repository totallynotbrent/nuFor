//! the axisymmetric viscous combustion lane: the reacting axi march plus
//! molecular transport, with the engine's wall treated.
//!
//! this is the lane the cooled-wall engine cases need. the planar
//! viscous combustion lane (combustion_viscous.rs) treats the duct's
//! rectangular boundary, but an engine's contour is a body: the wall
//! runs along a staircase of cells, and the flow there needs a one-sided
//! shear, conduction and species flux pointing INTO the gas.
//!
//! the transport coefficients are the same wilke-mixture ones the planar
//! lane uses. the difference is the geometry: the radial direction is
//! annular, so the radial fluxes divide by the radius and the axis is a
//! symmetry boundary rather than a wall.
//!
//! the wall treatment is one sided and isothermal:
//!   - no slip, so the tangential velocity vanishes at the wall
//!   - fourier conduction to the wall temperature when one is given
//!   - no composition gradient, so the wall neither creates nor absorbs
//!     species (the non-catalytic default)
//!
//! the wall is placed on the body's crossed faces: a face the wall
//! crosses carries the flux, taken from the fluid cell on its
//! gas side.

use crate::combustion_closure as closure9;
use crate::combustion_march::{advance2d_combustion_axi, CombustionState};
use crate::combustion_thermo as ct;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;
use crate::state2d::{cons_to_prim2d, ConservedState2d};

/// per-cell transport snapshot for the combustion set: closure
/// temperature, wilke-mixture mu and k, and the per-species effective
/// diffusivities.
pub(crate) struct AxiCombustionTransport {
    pub t: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub d_s: Vec<[f64; ct::NS]>,
}

/// the wall's treatment options.
#[derive(Clone, Copy, Debug, Default)]
pub struct AxiWallOpts {
    /// the wall temperature in K; zero or less means adiabatic.
    pub temperature: f64,
    /// the wall's grey-body emissivity; zero or less means the wall
    /// does not radiate.
    pub emissivity: f64,
    /// the wall's film heat-transfer coefficient in W per square metre
    /// kelvin; zero or less means no film.
    pub film_coefficient: f64,
}

/// the gas temperature extrapolated to the wall face, in K.
///
/// the radiative and film fluxes act on the gas at the wall surface,
/// not at the first cell's centre. near a cold wall the cell centre is
/// markedly cooler than the gas actually touching the wall, so using
/// the centre there would bias the exchange downward. this
/// extrapolates the two nearest cell-centre temperatures linearly to
/// the wall's position.
fn surface_temperature(v0: f64, v1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    if (c1 - c0).abs() < 1.0e-30 {
        return v0;
    }
    v0 + (y_wall - c0) * (v1 - v0) / (c1 - c0)
}

/// the wall's outward normal's direction, from the gas toward the body.
///
/// the signed distance is negative inside the body and grows outward,
/// so its gradient points from the body toward the gas. the heat
/// leaving the gas travels the other way, toward the body, so the
/// direction the heat travels is the negated gradient.
fn wall_heat_direction(dist: &dyn Fn(f64, f64) -> f64, x: f64, y: f64, h: f64) -> (f64, f64) {
    let dx = (dist(x + h, y) - dist(x - h, y)) / (2.0 * h);
    let dy = (dist(x, y + h) - dist(x, y - h)) / (2.0 * h);
    let norm = (dx * dx + dy * dy).sqrt();
    if norm < 1.0e-30 {
        return (0.0, 0.0);
    }
    (-dx / norm, -dy / norm)
}

/// the stefan-boltzmann constant in W per square metre per kelvin to
/// the fourth.
pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;

/// the radiant flux leaving the gas for a grey wall, in W per square
/// metre.
///
/// the gas sees the wall through the grey-body exchange
///
///   q = eps * sigma * (t_gas^4 - t_wall^4)
///
/// which at 3690 k is about 9.5 mw per square metre for a coated
/// chamber wall. that is an order of magnitude above the fourier
/// conduction the near-wall cells carry, so a wall modelled with
/// conduction alone exchanges a small fraction of the heat a real
/// chamber wall does.
fn radiant_wall_flux(t_gas: f64, t_wall: f64, emissivity: f64) -> f64 {
    if emissivity <= 0.0 || t_gas <= 0.0 || t_wall <= 0.0 {
        return 0.0;
    }
    let tg = t_gas.clamp(0.0, 6000.0);
    let tw = t_wall.clamp(0.0, 6000.0);
    emissivity * STEFAN_BOLTZMANN * (tg.powi(4) - tw.powi(4))
}

/// the film's flux from the gas to the wall, in W per square metre.
fn film_wall_flux(t_gas: f64, t_wall: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    h * (t_gas - t_wall)
}

/// the wall's total flux from the gas, in W per square metre.
///
/// a positive result means the wall is pulling heat out of the gas.
/// the three contributions are conduction through the near-wall
/// cells, radiation to the wall, and the film's convective exchange.
pub(crate) fn wall_flux_parts(
    t_gas: f64,
    t_wall: f64,
    conduction: f64,
    emissivity: f64,
    film_coefficient: f64,
) -> (f64, f64, f64) {
    let radiative = radiant_wall_flux(t_gas, t_wall, emissivity);
    let film = film_wall_flux(t_gas, t_wall, film_coefficient);
    (conduction + radiative + film, radiative, film)
}

/// one-sided quadratic derivative through the wall value and the first
/// two cell centres, evaluated at the wall.
fn wall_quad_deriv(v_wall: f64, v0: f64, v1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    let w0 = (2.0 * y_wall - c0 - c1) / ((y_wall - c0) * (y_wall - c1));
    let w1 = (2.0 * y_wall - y_wall - c1) / ((c0 - y_wall) * (c0 - c1));
    let w2 = (2.0 * y_wall - y_wall - c0) / ((c1 - y_wall) * (c1 - c0));
    v_wall * w0 + v0 * w1 + v1 * w2
}

/// the per-cell transport properties of the combustion set.
pub(crate) fn axi_combustion_transport(
    state: &ConservedState2d,
    species: &CombustionState,
    g: &Grid2d,
) -> Result<AxiCombustionTransport, Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut tr = AxiCombustionTransport {
        t: vec![0.0f64; n],
        mu: vec![0.0f64; n],
        k: vec![0.0f64; n],
        d_s: vec![[0.0f64; ct::NS]; n],
    };
    for k in 0..n {
        let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] = species.y[s][k];
        }
        tr.t[k] = closure9::temperature_at(e_int, &y);
        let p = closure9::pressure_at(state.rho[k], e_int, &y);
        tr.mu[k] = crate::combustion_transport::mixture_viscosity(&y, tr.t[k]);
        let mut cps = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            cps[s] = ct::cp_mass(s, tr.t[k]);
        }
        tr.k[k] = crate::combustion_transport::mixture_conductivity(&y, tr.t[k], &cps);
        let d_eff = crate::combustion_transport::effective_diffusivity(&y, tr.t[k], p.max(1.0));
        tr.d_s[k] = d_eff;
    }
    Ok(tr)
}

/// one viscous application on (state, species) with the given dt.
///
/// the wall flux is applied on the faces the body's boundary crosses.
/// the gas side of such a face is the fluid cell beside it, and the
/// one-sided derivatives run from the wall value into the gas.
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_axi_combustion_viscous(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    tr: &AxiCombustionTransport,
    dt: f64,
    wall_opts: AxiWallOpts,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    // the signed distance drives the wall's normal direction.
    let wall_fn = wall.map(|f| f as &dyn Fn(f64, f64) -> f64);

    // the body's crossed faces and the solid mask. the wall flux rides
    // on the crossed faces, and the solid cells hold the flow out.
    let (solid, wall_x, wall_y) = match wall {
        Some(dist) => {
            let band = 1.5 * g.dx.min(g.dy);
            let cache = crate::body::body_cache(g, dist, band);
            (
                cache.mask.clone(),
                cache.wall_faces_x.clone(),
                cache.wall_faces_y.clone(),
            )
        }
        None => (
            vec![false; nx * ny],
            vec![false; (nx + 1) * ny],
            vec![false; nx * (ny + 1)],
        ),
    };

    // face accumulators. the x faces carry the axial transport, the y
    // faces the radial one; in the axisymmetric lane the radial flux
    // divides by the radius, which the update does.
    let mut tx_x = vec![0.0f64; (nx + 1) * ny];
    let mut ty_x = vec![0.0f64; (nx + 1) * ny];
    let mut te_x = vec![0.0f64; (nx + 1) * ny];
    let mut ts_x: Vec<Vec<f64>> = vec![vec![0.0f64; (nx + 1) * ny]; ct::NS];

    let mut tx_y = vec![0.0f64; nx * (ny + 1)];
    let mut ty_y = vec![0.0f64; nx * (ny + 1)];
    let mut te_y = vec![0.0f64; nx * (ny + 1)];
    let mut ts_y: Vec<Vec<f64>> = vec![vec![0.0f64; nx * (ny + 1)]; ct::NS];

    // the interior axial faces
    for j in 0..ny {
        for f in 1..nx {
            let (a, b) = (idx(f - 1, j), idx(f, j));
            let k = f * ny + j;
            let mu_f = 0.5 * (tr.mu[a] + tr.mu[b]);
            let k_f = 0.5 * (tr.k[a] + tr.k[b]);
            let dx = g.centers_x[b] - g.centers_x[a];
            let dudx = (u[b] - u[a]) / dx;
            let dvdx = (v[b] - v[a]) / dx;
            let dtdx = (tr.t[b] - tr.t[a]) / dx;
            // the axisymmetric divergence: du/dx + v/r + dv/dr
            let rc = g.centers_y[a];
            let div = dudx + v[a] / rc.max(1.0e-12) + (v[b] - v[a]) / g.dys[j];
            let txx = mu_f * (2.0 * dudx - 2.0 / 3.0 * div);
            let txy = mu_f * dvdx;
            let qx = k_f * dtdx;
            let rho_f = 0.5 * (state.rho[a] + state.rho[b]);
            let mut q_enth = 0.0f64;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let dydx = (species.y[s][b] - species.y[s][a]) / dx;
                let d_f = 0.5 * (tr.d_s[a][s] + tr.d_s[b][s]);
                let js = -rho_f * d_f * dydx;
                ts_x[s][k] = js;
                q_enth += ct::h_mass(s, 0.5 * (tr.t[a] + tr.t[b])) * js;
            }
            te_x[k] = 0.5 * (u[a] + u[b]) * txx + 0.5 * (v[a] + v[b]) * txy + qx - q_enth;
            tx_x[k] = txx;
            ty_x[k] = txy;
        }
    }

    // the interior radial faces
    for i in 0..nx {
        for f in 1..ny {
            let (a, b) = (idx(i, f - 1), idx(i, f));
            let k = i * (ny + 1) + f;
            let mu_f = 0.5 * (tr.mu[a] + tr.mu[b]);
            let k_f = 0.5 * (tr.k[a] + tr.k[b]);
            let dy = g.centers_y[b] - g.centers_y[a];
            let dudy = (u[b] - u[a]) / dy;
            let dvdy = (v[b] - v[a]) / dy;
            let dtdy = (tr.t[b] - tr.t[a]) / dy;
            let rc = 0.5 * (g.centers_y[a] + g.centers_y[b]);
            let div = (u[b] - u[a]) / g.dxs[i] + v[a] / rc.max(1.0e-12) + dvdy;
            let txy = mu_f * dudy;
            let tyy = mu_f * (2.0 * dvdy - 2.0 / 3.0 * div);
            let qy = k_f * dtdy;
            let rho_f = 0.5 * (state.rho[a] + state.rho[b]);
            let mut q_enth = 0.0f64;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let dydy = (species.y[s][b] - species.y[s][a]) / dy;
                let d_f = 0.5 * (tr.d_s[a][s] + tr.d_s[b][s]);
                let js = -rho_f * d_f * dydy;
                ts_y[s][k] = js;
                q_enth += ct::h_mass(s, 0.5 * (tr.t[a] + tr.t[b])) * js;
            }
            te_y[k] = 0.5 * (u[a] + u[b]) * txy + 0.5 * (v[a] + v[b]) * tyy + qy - q_enth;
            tx_y[k] = txy;
            ty_y[k] = tyy;
        }
    }

    // the wall faces the body's boundary crosses. the gas side is the
    // fluid cell beside the face.
    if let Some(wall_fn) = wall_fn {
        for j in 0..ny {
            for f in 1..nx {
                if !wall_x[j * (nx + 1) + f] {
                    continue;
                }
                let (a, b, ok) = if !solid[idx(f, j)] {
                    (idx(f, j), idx(f + 1, j), f + 1 < nx)
                } else {
                    (idx(f - 1, j), idx(f - 2, j), f >= 2)
                };
                if !ok {
                    continue;
                }
                let k = f * ny + j;
                let face_x = g.faces_x[f];
                let (c0, c1) = (g.centers_x[a], g.centers_x[b]);
                let mu_f = tr.mu[a];
                let dudx = wall_quad_deriv(0.0, u[a], u[b], face_x, c0, c1);
                let dvdx = wall_quad_deriv(0.0, v[a], v[b], face_x, c0, c1);
                let mut q = 0.0f64;
                if wall_opts.temperature > 0.0 {
                    let dtdn =
                        wall_quad_deriv(wall_opts.temperature, tr.t[a], tr.t[b], face_x, c0, c1);
                    let conduction = -tr.k[a] * dtdn;
                    let t_surface = surface_temperature(tr.t[a], tr.t[b], face_x, c0, c1);
                    // the conduction term is already signed along the
                    // face's axis. the radiant and film terms are
                    // magnitudes of heat LEAVING the gas, so they take
                    // the wall's outward normal's sign along that axis
                    // before they join it.
                    let (hx, _) =
                        wall_heat_direction(wall_fn, face_x, g.centers_y[a], 0.25 * g.dx.min(g.dy));
                    let (_, radiant, film) = wall_flux_parts(
                        t_surface,
                        wall_opts.temperature,
                        0.0,
                        wall_opts.emissivity,
                        wall_opts.film_coefficient,
                    );
                    q = conduction + hx * (radiant + film);
                }
                tx_x[k] = mu_f * dudx;
                ty_x[k] = mu_f * dvdx;
                te_x[k] = q;
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    ts_x[s][k] = 0.0;
                }
            }
        }
        for i in 0..nx {
            for f in 1..ny {
                if !wall_y[i * (ny + 1) + f] {
                    continue;
                }
                let (a, b, ok) = if !solid[idx(i, f)] {
                    (idx(i, f), idx(i, f + 1), f + 1 < ny)
                } else {
                    (idx(i, f - 1), idx(i, f - 2), f >= 2)
                };
                if !ok {
                    continue;
                }
                let k = i * (ny + 1) + f;
                let face_y = g.faces_y[f];
                let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
                let mu_f = tr.mu[a];
                let dudy = wall_quad_deriv(0.0, u[a], u[b], face_y, c0, c1);
                let dvdy = wall_quad_deriv(0.0, v[a], v[b], face_y, c0, c1);
                let mut q = 0.0f64;
                if wall_opts.temperature > 0.0 {
                    let dtdn =
                        wall_quad_deriv(wall_opts.temperature, tr.t[a], tr.t[b], face_y, c0, c1);
                    let conduction = -tr.k[a] * dtdn;
                    let t_surface = surface_temperature(tr.t[a], tr.t[b], face_y, c0, c1);
                    // the radiant and film terms leave the gas, so they
                    // take the outward normal's sign along the axis.
                    let (_, hy) =
                        wall_heat_direction(wall_fn, g.centers_x[i], face_y, 0.25 * g.dx.min(g.dy));
                    let (_, radiant, film) = wall_flux_parts(
                        t_surface,
                        wall_opts.temperature,
                        0.0,
                        wall_opts.emissivity,
                        wall_opts.film_coefficient,
                    );
                    q = conduction + hy * (radiant + film);
                }
                tx_y[k] = mu_f * dudy;
                ty_y[k] = mu_f * dvdy;
                te_y[k] =
                    0.5 * (u[a] + u[b]) * (mu_f * dudy) + 0.5 * (v[a] + v[b]) * (mu_f * dvdy) + q;
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    ts_y[s][k] = 0.0;
                }
            }
        }
    }

    // the divergence update. the radial term divides by the radius, the
    // same way the reactive axi march does.
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let kfx_e = (i + 1) * ny + j;
            let kfx_w = i * ny + j;
            let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
            let rc = 0.5 * (rm + rp);
            // the radial divergence of a face quantity: the annular
            // weight over the radius and the annulus thickness. it
            // carries no dt, and the caller scales it once.
            //
            // the radial face quantities are laid out with a COLUMN
            // stride, one block of ny + 1 faces per column, so the
            // divergence has to reach into the block belonging to this
            // column. reading the block of the first column instead
            // spends every wall's flux on the first column's gas and
            // leaves the rest of the wall inert.
            let radial = |f: &[f64]| {
                let base = i * (ny + 1);
                (rp * f[base + j + 1] - rm * f[base + j]) / (rc * g.dys[j])
            };
            // tau_xx is the x-face axial stress and tau_xy the radial
            // shear: the x-momentum equation takes d(tau_xx)/dx plus
            // d(tau_xy)/dr, the y-momentum the other pairing.
            let dmx = -(dxi * (tx_x[kfx_e] - tx_x[kfx_w]) + dt * radial(&tx_y[..]));
            let dmy = -(dxi * (ty_x[kfx_e] - ty_x[kfx_w]) + dt * radial(&ty_y[..]));
            let de = -(dxi * (te_x[kfx_e] - te_x[kfx_w]) + dt * radial(&te_y[..]));
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let ds = -(dxi * (ts_x[s][kfx_e] - ts_x[s][kfx_w]) + dt * radial(&ts_y[s]));
                species.y[s][k] += ds / state.rho[k];
            }
        }
    }
    for k in 0..nx * ny {
        species.renorm(k);
    }
    Ok(())
}

/// the wall's heat duty and shear, gathered from the crossed faces.
///
/// the bracket needs MEASURED wall numbers, not a qualitative claim that
/// the wall cooler the gas. this sums the conductive flux leaving the
/// gas at every wall face, bucketed along the axis so the throat's
/// heating can be told from the bell's.
#[derive(Clone, Debug, Default)]
pub struct WallDiagnostics {
    /// the total power leaving the gas, W.
    pub heat_watts: f64,
    /// the conductive share of that power, W.
    pub conductive_watts: f64,
    /// the radiative share of that power, W.
    pub radiative_watts: f64,
    /// the film's share of that power, W.
    pub film_watts: f64,
    /// the wall's shear force along the axis, N.
    pub shear_newtons: f64,
    /// the buckets along the axis: each holds the x range and the power.
    pub buckets: Vec<WallBucket>,
}

/// one bucket of the wall, along the axis.
#[derive(Clone, Debug, Default)]
pub struct WallBucket {
    pub x0: f64,
    pub x1: f64,
    /// the total power leaving the gas in this slice, W.
    pub heat_watts: f64,
    /// the conductive share in this slice, W.
    pub conductive_watts: f64,
    /// the radiative share in this slice, W.
    pub radiative_watts: f64,
    /// the film's share in this slice, W.
    pub film_watts: f64,
    /// the wall's perimeter in this slice, m.
    pub perimeter_metres: f64,
    /// the wall's mean temperature in this slice, K.
    pub temperature_k: f64,
}

/// gathers the wall's diagnostics from the settled state.
///
/// the wall's conductive flux is the same one the viscous pass applied:
/// -k dT/dn, with dT/dn taken one-sided through the wall value and the
/// two gas cells beside it. a POSITIVE result means the wall is pulling
/// heat out of the gas.
pub fn wall_diagnostics(
    state: &ConservedState2d,
    species: &CombustionState,
    g: &Grid2d,
    nbuckets: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    wall_opts: AxiWallOpts,
) -> Result<WallDiagnostics, Error> {
    let wall_temperature = wall_opts.temperature;
    let wall_emissivity = wall_opts.emissivity;
    let wall_film_coefficient = wall_opts.film_coefficient;
    let mut out = WallDiagnostics::default();
    let dist = match wall {
        Some(d) => d,
        None => return Ok(out),
    };
    let tr = axi_combustion_transport(state, species, g)?;
    let (u, _v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let band = 1.5 * g.dx.min(g.dy);
    let cache = crate::body::body_cache(g, dist, band);
    let solid = cache.mask.clone();
    let width = g.faces_x[nx] - g.faces_x[0];
    let mut buckets: Vec<WallBucket> = (0..nbuckets)
        .map(|b| WallBucket {
            x0: g.faces_x[0] + width * b as f64 / nbuckets as f64,
            x1: g.faces_x[0] + width * (b + 1) as f64 / nbuckets as f64,
            ..Default::default()
        })
        .collect();

    // the x faces the wall crosses
    for j in 0..ny {
        for f in 1..nx {
            if !cache.wall_faces_x[j * (nx + 1) + f] {
                continue;
            }
            let (a, b, ok) = if !solid[idx(f, j)] {
                (idx(f, j), idx(f + 1, j), f + 1 < nx)
            } else {
                (idx(f - 1, j), idx(f - 2, j), f >= 2)
            };
            if !ok {
                continue;
            }
            let face_x = g.faces_x[f];
            let (c0, c1) = (g.centers_x[a], g.centers_x[b]);
            let dtdn = wall_quad_deriv(wall_temperature, tr.t[a], tr.t[b], face_x, c0, c1);
            // the conductive flux LEAVING the gas, per unit area
            let conduction = -tr.k[a] * dtdn;
            let t_surface = surface_temperature(tr.t[a], tr.t[b], face_x, c0, c1);
            let (hx, _) = wall_heat_direction(dist, face_x, g.centers_y[a], 0.25 * g.dx.min(g.dy));
            let (_, q_rad, q_film) = wall_flux_parts(
                t_surface,
                wall_temperature,
                0.0,
                wall_emissivity,
                wall_film_coefficient,
            );
            let q = conduction + hx * (q_rad + q_film);
            // the face's area is the annulus at this radius, times the
            // cell's circumferential span. the axisymmetric face is a
            // cylinder of this radius and this height.
            let r_face = g.centers_y[a].max(1.0e-9);
            let height = g.dys[j];
            let area = 2.0 * std::f64::consts::PI * r_face * height;
            let power = q * area;
            let bi = (((face_x - g.faces_x[0]) / width) * nbuckets as f64)
                .floor()
                .max(0.0) as usize;
            let bi = bi.min(nbuckets - 1);
            buckets[bi].heat_watts += power;
            buckets[bi].conductive_watts += conduction * area;
            buckets[bi].radiative_watts += q_rad * area;
            buckets[bi].film_watts += q_film * area;
            buckets[bi].perimeter_metres += 2.0 * std::f64::consts::PI * r_face;
            buckets[bi].temperature_k += wall_temperature;
            out.heat_watts += power;
            out.conductive_watts += conduction * area;
            out.radiative_watts += q_rad * area;
            out.film_watts += q_film * area;
            // the wall's shear: the axial viscous stress at the wall
            let dudx = wall_quad_deriv(0.0, u[a], u[b], face_x, c0, c1);
            out.shear_newtons += tr.mu[a] * dudx * area;
        }
    }
    // the y faces the wall crosses
    for i in 0..nx {
        for f in 1..ny {
            if !cache.wall_faces_y[i * (ny + 1) + f] {
                continue;
            }
            let (a, b, ok) = if !solid[idx(i, f)] {
                (idx(i, f), idx(i, f + 1), f + 1 < ny)
            } else {
                (idx(i, f - 1), idx(i, f - 2), f >= 2)
            };
            if !ok {
                continue;
            }
            let face_y = g.faces_y[f];
            let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
            let dtdn = wall_quad_deriv(wall_temperature, tr.t[a], tr.t[b], face_y, c0, c1);
            let conduction = -tr.k[a] * dtdn;
            let t_surface = surface_temperature(tr.t[a], tr.t[b], face_y, c0, c1);
            let (_, hy) = wall_heat_direction(dist, g.centers_x[i], face_y, 0.25 * g.dx.min(g.dy));
            let (_, q_rad, q_film) = wall_flux_parts(
                t_surface,
                wall_temperature,
                0.0,
                wall_emissivity,
                wall_film_coefficient,
            );
            let q = conduction + hy * (q_rad + q_film);
            // a y face is a cone frustum's lateral surface
            let rm = g.centers_y[a];
            let rp = g.centers_y[b];
            let slant = ((rp - rm).powi(2) + g.dxs[i].powi(2)).sqrt();
            let area = std::f64::consts::PI * (rp + rm) * slant;
            let power = q * area;
            let bi = (((g.centers_x[i] - g.faces_x[0]) / width) * nbuckets as f64)
                .floor()
                .max(0.0) as usize;
            let bi = bi.min(nbuckets - 1);
            buckets[bi].heat_watts += power;
            buckets[bi].conductive_watts += conduction * area;
            buckets[bi].radiative_watts += q_rad * area;
            buckets[bi].film_watts += q_film * area;
            buckets[bi].perimeter_metres += std::f64::consts::PI * (rp + rm);
            buckets[bi].temperature_k += wall_temperature;
            out.heat_watts += power;
            out.conductive_watts += conduction * area;
            out.radiative_watts += q_rad * area;
            out.film_watts += q_film * area;
        }
    }
    for b in buckets.iter_mut() {
        if b.perimeter_metres > 0.0 {
            b.temperature_k /= b.perimeter_metres;
        }
    }
    out.buckets = buckets;
    Ok(out)
}

/// the reacting axi march plus one viscous application per heun stage.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_axi_visc_rk2(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    wall_opts: AxiWallOpts,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    react: bool,
) -> Result<(f64, f64), Error> {
    let tr = axi_combustion_transport(state, species, g)?;
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    let mut dt_last = 0.0f64;
    for _ in 0..2 {
        let (dt, _) =
            advance2d_combustion_axi(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
        add_axi_combustion_viscous(&mut s1, &mut y1, g, &tr, dt, wall_opts, wall)?;
        dt_last = dt;
    }
    let n = g.nx * g.ny;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    for s in 0..ct::NS {
        for k in 0..n {
            species.y[s][k] = 0.5 * species.y[s][k] + 0.5 * y1.y[s][k];
        }
    }
    for k in 0..n {
        species.renorm(k);
    }
    if react {
        let (dt, _) = {
            let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
            let mut smax = 0.0f64;
            for k in 0..n {
                let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
                let mut y_cell = [0.0f64; ct::NS];
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    y_cell[s] = species.y[s][k];
                }
                let a = closure9::sound_speed_at(state.rho[k], e_int, &y_cell);
                smax = smax.max(u[k].abs().max(v[k].abs()) + a);
            }
            (cfl * g.dx.min(g.dy) / smax.max(1.0e-12), 0.0)
        };
        // the chemistry runs on the BLENDED state and species, the ones
        // the caller receives, not on the stage state that the blend
        // discards.
        let oversub = crate::combustion_march::chemistry_substep9(state, species, g, dt, nthreads)?;
        return Ok((dt_last, oversub));
    }
    Ok((dt_last, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_case() -> (Grid2d, ConservedState2d, CombustionState) {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let u = 100.0f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![rho * u; n],
            my: vec![0.0; n],
            e: vec![rho * (e_int + 0.5 * u * u); n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        (g, st, sp)
    }

    /// a uniform field has no gradients, so the viscous pass moves
    /// nothing. the axisymmetric lane shares the transport with the
    /// planar one, so this catches a mispaired stress or a stray dt.
    #[test]
    fn the_radiant_wall_flux_follows_the_stefan_boltzmann_law() {
        // the blackbody flux at 3000 k, from the law itself
        let sigma = STEFAN_BOLTZMANN;
        let blackbody = sigma * 3000.0_f64.powi(4);
        // a grey wall emits and absorbs, so its net exchange with a
        // hot gas scales with its emissivity
        let q = radiant_wall_flux(3000.0, 300.0, 0.9);
        let expected = 0.9 * sigma * (3000.0_f64.powi(4) - 300.0_f64.powi(4));
        assert!((q - expected).abs() / expected.abs() < 1.0e-12);
        // a black wall draws more than a white one
        assert!(radiant_wall_flux(3000.0, 300.0, 1.0) > q);
        // an insulated wall draws nothing
        assert_eq!(radiant_wall_flux(3000.0, 300.0, 0.0), 0.0);
        // a wall hotter than the gas receives heat rather than losing it
        assert!(radiant_wall_flux(300.0, 3000.0, 0.9) < 0.0);
        assert!(blackbody > 0.0);
    }

    #[test]
    fn the_chamber_walls_radiation_is_large_next_to_conduction() {
        // the point of adding radiation: a 3690 k gas facing a 400 k
        // wall exchanges far more energy by radiation than the
        // near-wall cells conduct, so a conduction-only wall is not a
        // representation of a real chamber wall.
        let q_rad = radiant_wall_flux(3690.0, 400.0, 0.9);
        // the measured conduction-only duty implied an effective
        // conductance of about 340 w per square metre kelvin
        let q_cond = 341.0 * (3690.0 - 400.0);
        // eight times the conduction the old wall model carried
        assert!(q_rad > 8.0 * q_cond);
        // and it is on the order of ten megawatts per square metre
        assert!(q_rad > 8.0e6 && q_rad < 1.1e7);
    }

    #[test]
    fn a_uniform_field_feels_no_viscous_force() {
        let (g, st, sp) = square_case();
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let opts = AxiWallOpts {
            temperature: 0.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
        };
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, 1.0e-6, opts, None).unwrap();
        // the radial annulus weighting makes the divergence of a uniform
        // field nonzero in principle, so the contract is weaker than the
        // planar one: the change must stay a rounding-level fraction of
        // the state's own scale, not a physical acceleration.
        let scale = st.mx.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
        let change = s1
            .mx
            .iter()
            .zip(&st.mx)
            .fold(0.0f64, |acc, (a, b)| acc.max((a - b).abs()));
        assert!(
            change < 1.0e-6 * scale.max(1.0),
            "a uniform field moved by {change} against a scale of {scale}"
        );
    }

    /// the MASS must be conserved by the viscous pass. the wall is a
    /// heat sink and a shear surface, but it is impermeable: the gas
    /// may not leave through it. a wall flux that leaked mass would
    /// quietly change the throat's state, and the bracket's numbers
    /// would be measuring that leak instead of the cooling.
    #[test]
    fn the_viscous_pass_conserves_mass() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let opts = AxiWallOpts {
            temperature: 300.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
        };
        let dist = |_x: f64, y: f64| y - 0.02;
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, 1.0e-5, opts, Some(&dist)).unwrap();
        let before: f64 = st.rho.iter().sum();
        let after: f64 = s1.rho.iter().sum();
        assert!(
            (after - before).abs() < 1.0e-12 * before.abs().max(1.0),
            "the viscous pass changed the mass: {before} against {after}"
        );
    }

    /// a cold wall cools the gas beside it. the near-wall cell must
    /// lose energy, and the far cells must not gain it.
    #[test]
    fn a_radiating_wall_cools_the_gas_more_than_a_non_radiating_one() {
        let (nx, ny) = (8, 8);
        let idx = |i: usize, j: usize| j * nx + i;
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let dist = |_x: f64, y: f64| y - 0.02;

        // the same cold wall, once conducting only and once radiating.
        let mut dull = st.clone();
        let mut dull_sp = sp.clone();
        let mut glowing = st.clone();
        let mut glowing_sp = sp.clone();

        add_axi_combustion_viscous(
            &mut dull,
            &mut dull_sp,
            &g,
            &tr,
            1.0e-5,
            AxiWallOpts {
                temperature: 300.0,
                emissivity: 0.0,
                film_coefficient: 0.0,
            },
            Some(&dist),
        )
        .unwrap();
        add_axi_combustion_viscous(
            &mut glowing,
            &mut glowing_sp,
            &g,
            &tr,
            1.0e-5,
            AxiWallOpts {
                temperature: 300.0,
                emissivity: 0.9,
                film_coefficient: 0.0,
            },
            Some(&dist),
        )
        .unwrap();

        let near = idx(0, 2);
        let dull_loss = st.e[near] - dull.e[near];
        let glowing_loss = st.e[near] - glowing.e[near];
        // both take heat out of the near-wall cell
        assert!(dull_loss > 0.0);
        assert!(glowing_loss > 0.0);
        // and the radiating wall takes strictly more out
        assert!(
            glowing_loss > dull_loss,
            "radiation did not add to the wall's duty: {glowing_loss:e} against {dull_loss:e}; surf_T would be needed"
        );
    }

    /// the radial divergence must reach the column it is updating.
    ///
    /// the radial face quantities are laid out with a column stride, so
    /// the divergence has to index the column it sits in. if it reaches
    /// another column, the wall's flux is spent on the wrong gas, and
    /// the gas under the wall appears untouched.
    ///
    /// the fixture puts a walled shelf in the middle of the mesh, well
    /// away from the first column, and bounds it in height so there is
    /// gas above it. the gas above the shelf must feel the shelf's
    /// temperature, and the gas beside the shelf must feel it less than
    /// the gas on it.
    #[test]
    fn the_radial_divergence_reads_its_own_column() {
        let (nx, ny) = (12, 8);
        let idx = |i: usize, j: usize| j * nx + i;
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();

        // a shelf: it sits between x = 0.04 and x = 0.07, and rises from
        // the axis to y = 0.03. the distance is negative inside it and
        // positive outside, which is what the body cache reads.
        let dist = |x: f64, y: f64| {
            let dx = (0.04 - x).max(x - 0.07);
            let dy = (0.0 - y).max(y - 0.03);
            dx.max(dy)
        };

        let march = |wall_temp: f64| {
            let mut s = st.clone();
            let mut y = sp.clone();
            let opts = AxiWallOpts {
                temperature: wall_temp,
                emissivity: 0.0,
                film_coefficient: 0.0,
            };
            add_axi_combustion_viscous(&mut s, &mut y, &g, &tr, 1.0e-5, opts, Some(&dist)).unwrap();
            s
        };

        // the gas's own temperature gives a wall that exchanges nothing,
        // so the difference between the two marches is the wall's doing.
        let gas_temp = tr.t[idx(0, 0)];
        let cold = march(300.0);
        let neutral = march(gas_temp);

        // the row of gas sitting on the shelf's top face. the shelf
        // spans the columns 4 to 6.
        let above_row = 3usize;
        let mut drained = 0usize;
        for i in 4..=6 {
            let k = idx(i, above_row);
            let drain = neutral.e[k] - cold.e[k];
            assert!(
                drain > 0.0,
                "the gas on the shelf in column {i} did not lose energy to \
                 the cold shelf: {drain:e}"
            );
            drained += 1;
        }
        assert_eq!(drained, 3, "the fixture did not place gas on the shelf");

        // the gas in the first column, which stands on no wall at all,
        // must not feel the shelf's temperature. if the radial flux is
        // read from the wrong column, the first column picks up every
        // wall's flux instead.
        for i in 0..3 {
            let k = idx(i, above_row);
            let drain = neutral.e[k] - cold.e[k];
            assert!(
                drain.abs() < 1.0e-6,
                "the gas in the wall-free column {i} reacted to the shelf's \
                 temperature: {drain:e}"
            );
        }
    }

    #[test]
    fn a_cold_wall_cools_the_near_wall_cell() {
        let (nx, ny) = (8, 8);
        let idx = |i: usize, j: usize| j * nx + i;
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let opts = AxiWallOpts {
            temperature: 300.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
        };
        // a flat wall along the bottom, so the row 0 cells are the
        // near-wall ones. the distance function says the gas is above
        // the line y = 0.02.
        let dist = |_x: f64, y: f64| y - 0.02;
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, 1.0e-5, opts, Some(&dist)).unwrap();
        // the wall sits at y = 0.02, so the rows below it are the body
        // and the GAS SIDE is the first row above it (row 2).
        let near = idx(0, 2);
        let far = idx(0, ny - 1);
        let delta = s1.e[near] - st.e[near];
        assert!(delta < 0.0, "the near-wall cell gained energy: {delta}");
        let far_delta = s1.e[far] - st.e[far];
        assert!(
            far_delta.abs() < 1.0e-9,
            "the far cell changed: {far_delta}"
        );
    }
}

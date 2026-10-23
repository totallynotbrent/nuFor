//! the viscous combustion lane: the reacting march plus molecular
//! transport. momentum and energy get shear and conduction through
//! the wilke-mixture transport of the combustion set, the species
//! get fick diffusion with the effective binary diffusivity, and a
//! no-slip isothermal wall carries one-sided shear, conduction, and
//! the wallward formation-enthalpy flux. this is the machinery the
//! cooled-wall engine cases ride.

use crate::combustion_closure as closure9;
use crate::combustion_march::{advance2d_combustion, CombustionState};
use crate::combustion_thermo as ct;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;
use crate::state2d::{cons_to_prim2d, ConservedState2d};

/// per-cell transport snapshot for the combustion set: closure
/// temperature and wilke-mixture mu, k, plus the effective
/// diffusivity per species.
pub(crate) struct CombustionTransport {
    pub t: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub d: Vec<f64>,
    pub d_s: Vec<[f64; ct::NS]>,
}

/// one-sided quadratic derivative through the wall value and the first
/// two cell centers, evaluated at the wall (the park lane's
/// quad_deriv_at_wall, restated locally).
fn wall_quad_deriv(v_wall: f64, v0: f64, v1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    let w0 = (2.0 * y_wall - c0 - c1) / ((y_wall - c0) * (y_wall - c1));
    let w1 = (2.0 * y_wall - y_wall - c1) / ((c0 - y_wall) * (c0 - c1));
    let w2 = (2.0 * y_wall - y_wall - c0) / ((c1 - y_wall) * (c1 - c0));
    v_wall * w0 + v0 * w1 + v1 * w2
}

/// the wall treatment options (mirror of the park lane's).
#[derive(Clone, Copy)]
pub(crate) struct CombustionWallOpts {
    pub temperature: f64,
}

pub(crate) fn combustion_transport(
    state: &ConservedState2d,
    species: &CombustionState,
    g: &Grid2d,
) -> Result<CombustionTransport, Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut tr = CombustionTransport {
        t: vec![0.0f64; n],
        mu: vec![0.0f64; n],
        k: vec![0.0f64; n],
        d: vec![0.0f64; n],
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
        // the scalar face diffusivity: the mass-weighted mean of the
        // per-species values, the same single-d field the park lane
        // feeds its face loops.
        let mut num = 0.0;
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            num += y[s] * d_eff[s];
        }
        tr.d[k] = num;
    }
    Ok(tr)
}

/// one viscous application on (state, species) with the given dt.
/// body-free for now (the engine cases are duct walls); the domain
/// no-slip faces carry the wall fluxes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_combustion_viscous(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    tr: &CombustionTransport,
    dt: f64,
    bc: &Boundaries2d,
    wall_opts: CombustionWallOpts,
) -> Result<(), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut fxm = vec![0.0f64; (nx + 1) * ny];
    let mut fym = vec![0.0f64; (nx + 1) * ny];
    let mut fe = vec![0.0f64; (nx + 1) * ny];
    let mut fdx: Vec<Vec<f64>> = vec![vec![0.0f64; (nx + 1) * ny]; ct::NS];
    let mut gym_x = vec![0.0f64; nx * (ny + 1)];
    let mut gym_y = vec![0.0f64; nx * (ny + 1)];
    let mut gye = vec![0.0f64; nx * (ny + 1)];
    let mut fdy: Vec<Vec<f64>> = vec![vec![0.0f64; nx * (ny + 1)]; ct::NS];

    // interior x faces
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
            let div = dudx;
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
                fdx[s][k] = js;
                q_enth += ct::h_mass(s, 0.5 * (tr.t[a] + tr.t[b])) * js;
            }
            fe[k] = 0.5 * (u[a] + u[b]) * txx + 0.5 * (v[a] + v[b]) * txy + qx - q_enth;
            fxm[k] = txx;
            fym[k] = txy;
        }
    }
    // domain west/east no-slip faces
    for j in 0..ny {
        let (face, a, b, into) = if matches!(bc.west, crate::solver2d::Bc2d::NoSlipWall) {
            (0usize, idx(0, j), idx(1, j), 1.0f64)
        } else if matches!(bc.east, crate::solver2d::Bc2d::NoSlipWall) {
            (nx, idx(nx - 1, j), idx(nx - 2, j), -1.0f64)
        } else {
            continue;
        };
        let k = face * ny + j;
        let (c0, c1) = (g.centers_x[a], g.centers_x[b]);
        let face_x = g.faces_x[face];
        let mu_f = tr.mu[a];
        let dudx = into * wall_quad_deriv(0.0, u[a], u[b], face_x, c0, c1);
        let dvdx = into * wall_quad_deriv(0.0, v[a], v[b], face_x, c0, c1);
        let txx = mu_f * dudx;
        let txy = mu_f * dvdx;
        let mut q = 0.0f64;
        if wall_opts.temperature > 0.0 {
            // the +x fourier component at the wall: -k dT/dx, so a
            // hot cell against a cold wall loses energy either side.
            q = -tr.k[a] * wall_quad_deriv(wall_opts.temperature, tr.t[a], tr.t[b], face_x, c0, c1);
        }
        fxm[k] = txx;
        fym[k] = txy;
        fe[k] = q;
    }
    // interior y faces
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
            let div = dvdy;
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
                fdy[s][k] = js;
                q_enth += ct::h_mass(s, 0.5 * (tr.t[a] + tr.t[b])) * js;
            }
            gye[k] = 0.5 * (u[a] + u[b]) * txy + 0.5 * (v[a] + v[b]) * tyy + qy - q_enth;
            gym_x[k] = txy;
            gym_y[k] = tyy;
        }
    }
    // domain south/north no-slip faces
    for i in 0..nx {
        let (face, a, b, flip) = if matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall) {
            (0usize, idx(i, 0), idx(i, 1), 1.0f64)
        } else if matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall) {
            (ny, idx(i, ny - 1), idx(i, ny - 2), -1.0f64)
        } else {
            continue;
        };
        let k = i * (ny + 1) + face;
        let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
        let face_y = g.faces_y[face];
        let mu_f = tr.mu[a];
        let dudy = flip * wall_quad_deriv(0.0, u[a], u[b], face_y, c0, c1);
        let dvdy = flip * wall_quad_deriv(0.0, v[a], v[b], face_y, c0, c1);
        let txy = mu_f * dudy;
        let tyy = mu_f * dvdy;
        let mut q = 0.0f64;
        if wall_opts.temperature > 0. {
            // the +y fourier component at the wall: -k dT/dy, so a
            // hot cell against a cold wall loses energy either side.
            q = -tr.k[a] * wall_quad_deriv(wall_opts.temperature, tr.t[a], tr.t[b], face_y, c0, c1);
        }
        gym_x[k] = txy;
        gym_y[k] = tyy;
        gye[k] = u[a] * txy + v[a] * tyy + q;
        // species at a no-slip wall: zero gradient (the wall neither
        // creates nor absorbs composition), the non-catalytic default.
        #[allow(clippy::needless_range_loop)]
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            fdy[s][k] = 0.0;
        }
    }

    // the divergence update, face-major layouts matching the park lane
    let dxs = &g.dxs;
    let dys = &g.dys;
    #[allow(clippy::needless_range_loop)]
    for j in 0..ny {
        #[allow(clippy::needless_range_loop)]
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = 1.0 / dxs[i];
            let dyi = 1.0 / dys[j];
            let kfx_e = (i + 1) * ny + j;
            let kfx_w = i * ny + j;
            let kfy_n = i * (ny + 1) + j + 1;
            let kfy_s = i * (ny + 1) + j;
            let dmx = -(dxi * (fxm[kfx_e] - fxm[kfx_w]) + dyi * (gym_x[kfy_n] - gym_x[kfy_s]));
            let dmy = -(dxi * (fym[kfx_e] - fym[kfx_w]) + dyi * (gym_y[kfy_n] - gym_y[kfy_s]));
            let de = -(dxi * (fe[kfx_e] - fe[kfx_w]) + dyi * (gye[kfy_n] - gye[kfy_s]));
            state.mx[k] += dt * dmx;
            state.my[k] += dt * dmy;
            state.e[k] += dt * de;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let ds = -(dxi * (fdx[s][kfx_e] - fdx[s][kfx_w])
                    + dyi * (fdy[s][kfy_n] - fdy[s][kfy_s]));
                species.y[s][k] += dt * ds / state.rho[k];
            }
        }
    }
    for k in 0..nx * ny {
        species.renorm(k);
    }
    Ok(())
}

/// the reacting march plus one viscous application per heun stage,
/// mirroring the park lane's visc_rk2 shape.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_visc_rk2(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    wall_temperature: f64,
    nthreads: usize,
    react: bool,
) -> Result<(f64, f64), Error> {
    let tr = combustion_transport(state, species, g)?;
    let wall_opts = CombustionWallOpts {
        temperature: wall_temperature,
    };
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    let mut dt_last = 0.0f64;
    for _ in 0..2 {
        let (dt, _) = advance2d_combustion(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads)?;
        add_combustion_viscous(&mut s1, &mut y1, g, &tr, dt, bc, wall_opts)?;
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
            (cfl * g.dx.min(g.dy) / smax.max(1e-12), 0.0)
        };
        let oversub = crate::combustion_march::chemistry_substep9(state, species, g, dt, nthreads)?;
        return Ok((dt_last, oversub));
    }
    Ok((dt_last, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the viscous pass must hold a uniform field: no gradients means
    /// no flux, so the state is unchanged modulo the inviscid part of
    /// the uniform duct.
    #[test]
    fn uniform_viscous_pass_is_a_no_op() {
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
        let tr = combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        add_combustion_viscous(
            &mut s1,
            &mut y1,
            &g,
            &tr,
            1.0e-6,
            &bc,
            CombustionWallOpts { temperature: 0.0 },
        )
        .unwrap();
        for k in 0..n {
            assert!((s1.e[k] - st.e[k]).abs() / st.e[k] < 1e-12, "energy moved");
            assert!((s1.mx[k] - st.mx[k]).abs() < 1e-9, "momentum moved");
        }
    }

    /// the isothermal wall pulls the near-wall energy down: a hot field
    /// against a cold wall must lose energy through the wall face.
    #[test]
    fn cold_wall_cools_the_near_wall_cell() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 1.0f64;
        let n = nx * ny;
        // uniform 3000 k field at rest
        let mut y = [0.0f64; ct::NS];
        y[ct::IH2O] = 0.9;
        y[ct::IOH] = 0.1;
        let e_int = closure9::internal_energy_at(3000.0, &y);
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::uniform(n, &y);
        let tr = combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::NoSlipWall,
        };
        add_combustion_viscous(
            &mut s1,
            &mut y1,
            &g,
            &tr,
            1.0e-5,
            &bc,
            CombustionWallOpts { temperature: 300.0 },
        )
        .unwrap();
        let k_wall = (ny - 1) * nx;
        assert!(
            s1.e[k_wall] < st.e[k_wall],
            "the cold wall must cool the near-wall cell"
        );
    }
}

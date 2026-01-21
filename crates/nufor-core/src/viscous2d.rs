//! 2d viscous navier-stokes terms layered on the euler update.
//!
//! the viscous flux is the diffusive counterpart of the euler flux: shear
//! stresses from velocity gradients and heat conduction from temperature
//! gradients. mu is the dynamic viscosity, pr the prandtl number linking
//! conduction to viscosity, and gamma sets the specific-heat factor in the
//! conductivity. a sutherland helper provides the temperature-dependent law.

use crate::eos2d::eos_pressure2d;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;
use crate::state2d::{cons_to_prim2d, ConservedState2d};

/// sutherland's law: viscosity rising with temperature, referenced at (mu0, t0).
pub fn sutherland_mu(mu0: f64, t: f64, t0: f64, s: f64) -> f64 {
    mu0 * (t / t0).powf(1.5) * (t0 + s) / (t + s)
}

/// mirror the interior field with one ghost layer so centred differences work.
fn pad2d(f: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let w = nx + 2;
    let mut out = vec![0.0; w * (ny + 2)];
    for j in 0..ny {
        for i in 0..nx {
            out[(j + 1) * w + (i + 1)] = f[j * nx + i];
        }
    }
    for j in 0..ny {
        out[(j + 1) * w] = out[(j + 1) * w + 1];
        out[(j + 1) * w + nx + 1] = out[(j + 1) * w + nx];
    }
    for i in 0..w {
        out[i] = out[i + w];
        out[(ny + 1) * w + i] = out[ny * w + i];
    }
    out
}

/// face-averaged primitives and their gradients on the face between two cells.
struct FaceGrad {
    u: f64,
    v: f64,
    gux: f64,
    guy: f64,
    gvx: f64,
    gvy: f64,
    gtx: f64,
    gty: f64,
}

/// average the cell-centred gradients of a and b onto their shared face.
fn face_grad(ga: &Grads, gb: &Grads, ua: f64, ub: f64, va: f64, vb: f64) -> FaceGrad {
    FaceGrad {
        u: 0.5 * (ua + ub),
        v: 0.5 * (va + vb),
        gux: 0.5 * (ga.ux + gb.ux),
        guy: 0.5 * (ga.uy + gb.uy),
        gvx: 0.5 * (ga.vx + gb.vx),
        gvy: 0.5 * (ga.vy + gb.vy),
        gtx: 0.5 * (ga.tx + gb.tx),
        gty: 0.5 * (ga.ty + gb.ty),
    }
}

/// cell-centred velocity and temperature gradients (central differencing).
#[derive(Clone, Copy)]
struct Grads {
    ux: f64,
    uy: f64,
    vx: f64,
    vy: f64,
    tx: f64,
    ty: f64,
}

/// stresses and heat flux from face gradients; kappa is the conductivity/mu.
fn viscous_flux(f: &FaceGrad, mu: f64, kappa: f64, gamma: f64) -> (f64, f64, f64, f64, f64) {
    let div = f.gux + f.gvy;
    let txx = f.gux * 2.0 - div * 2.0 / 3.0;
    let tyy = f.gvy * 2.0 - div * 2.0 / 3.0;
    let txy = f.guy + f.gvx;
    // heat conduction enters the energy flux as +k grad(T) (so a hot peak cools).
    let qx = f.gtx;
    let qy = f.gty;
    let _ = gamma;
    (
        mu * txx,
        mu * txy,
        mu * tyy,
        mu * kappa * qx,
        mu * kappa * qy,
    )
}

fn cell_grads(
    pu: &[f64],
    pv: &[f64],
    pt: &[f64],
    nx: usize,
    ny: usize,
    dx2: f64,
    dy2: f64,
) -> Vec<Grads> {
    let w = nx + 2;
    let mut out = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let (ip, jp) = (i + 1, j + 1);
            // x derivative uses the padded columns i and i+2 about cell column ip.
            let ux = (pu[jp * w + (i + 2)] - pu[jp * w + i]) / dx2;
            let vx = (pv[jp * w + (i + 2)] - pv[jp * w + i]) / dx2;
            let tx = (pt[jp * w + (i + 2)] - pt[jp * w + i]) / dx2;
            // y derivative uses the padded rows j and j+2 about cell row jp.
            let uy = (pu[(j + 2) * w + ip] - pu[j * w + ip]) / dy2;
            let vy = (pv[(j + 2) * w + ip] - pv[j * w + ip]) / dy2;
            let ty = (pt[(j + 2) * w + ip] - pt[j * w + ip]) / dy2;
            out.push(Grads {
                ux,
                uy,
                vx,
                vy,
                tx,
                ty,
            });
        }
    }
    out
}

/// add one viscous (diffusive) half-step to the state using central gradients.
pub fn add_viscous(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    mu: f64,
    pr: f64,
    dt: f64,
) -> Result<(), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    if mu <= 0.0 || pr <= 0.0 || dt <= 0.0 {
        return Ok(());
    }
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let t: Vec<f64> = p.iter().zip(&state.rho).map(|(pp, r)| pp / r).collect();
    let (pu, pv, pt) = (pad2d(&u, nx, ny), pad2d(&v, nx, ny), pad2d(&t, nx, ny));
    let grads = cell_grads(&pu, &pv, &pt, nx, ny, 2.0 * g.dx, 2.0 * g.dy);
    let kappa = gamma / ((gamma - 1.0) * pr);
    let (dtdx, dtdy) = (dt / g.dx, dt / g.dy);
    // x-face viscous flux (nx+1 faces per row); domain-edge faces carry zero flux.
    let mut fxm = vec![0.0; (nx + 1) * ny];
    let mut fym = vec![0.0; (nx + 1) * ny];
    let mut fe = vec![0.0; (nx + 1) * ny];
    for j in 0..ny {
        for f in 1..nx {
            let (a, b) = (idx(f - 1, j), idx(f, j));
            let g = face_grad(&grads[a], &grads[b], u[a], u[b], v[a], v[b]);
            let (txx, txy, _, qx, _) = viscous_flux(&g, mu, kappa, gamma);
            let k = f * ny + j;
            fxm[k] = txx;
            fym[k] = txy;
            fe[k] = g.u * txx + g.v * txy + qx;
        }
    }
    // y-face viscous flux (ny+1 faces per column); domain-edge faces carry zero flux.
    let mut gym_x = vec![0.0; nx * (ny + 1)];
    let mut gym_y = vec![0.0; nx * (ny + 1)];
    let mut gye = vec![0.0; nx * (ny + 1)];
    for i in 0..nx {
        for f in 1..ny {
            let (a, b) = (idx(i, f - 1), idx(i, f));
            let g = face_grad(&grads[a], &grads[b], u[a], u[b], v[a], v[b]);
            let (_, txy, tyy, _, qy) = viscous_flux(&g, mu, kappa, gamma);
            let k = i * (ny + 1) + f;
            gym_x[k] = txy;
            gym_y[k] = tyy;
            gye[k] = g.u * txy + g.v * tyy + qy;
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            let c = idx(i, j);
            // x-face divergence: faces i and i+1 bound cell i.
            let (a, b) = (i * ny + j, (i + 1) * ny + j);
            state.mx[c] += dtdx * (fxm[b] - fxm[a]);
            state.my[c] += dtdx * (fym[b] - fym[a]);
            state.e[c] += dtdx * (fe[b] - fe[a]);
            // y-face divergence: faces j and j+1 bound cell j.
            let (d, e) = (i * (ny + 1) + j, i * (ny + 1) + j + 1);
            state.mx[c] += dtdy * (gym_x[e] - gym_x[d]);
            state.my[c] += dtdy * (gym_y[e] - gym_y[d]);
            state.e[c] += dtdy * (gye[e] - gye[d]);
        }
    }
    Ok(())
}

/// the physical parameters of the viscous fluid.
#[derive(Debug, Clone, Copy)]
pub struct ViscParams {
    /// dynamic viscosity.
    pub mu: f64,
    /// prandtl number linking heat conduction to viscosity.
    pub pr: f64,
}

/// second-order heun step with viscous terms added at each stage.
pub fn advance2d_visc_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    v: ViscParams,
) -> Result<(f64, f64), Error> {
    use crate::solver2d::advance2d_capped;
    // explicit diffusion bounds the step; the binding coefficient is the larger
    // of momentum and heat diffusion (kappa = conductivity / mu).
    let mu = v.mu;
    let pr = v.pr;
    let kappa = gamma / ((gamma - 1.0) * pr);
    let rho_min = state
        .rho
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(1e-6);
    let dt_visc = 0.25 * g.dx.min(g.dy).powi(2) * rho_min / (mu * kappa.max(1.0));
    let mut s1 = state.clone();
    let (dt1, _) = advance2d_capped(&mut s1, g, gamma, cfl, muscl, bc, dt_visc)?;
    add_viscous(&mut s1, g, gamma, mu, pr, dt1)?;
    let (dt2, _) = advance2d_capped(&mut s1, g, gamma, cfl, muscl, bc, dt_visc)?;
    add_viscous(&mut s1, g, gamma, mu, pr, dt2)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

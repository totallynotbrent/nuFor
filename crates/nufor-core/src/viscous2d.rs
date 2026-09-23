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

/// the velocity pad: at a no-slip wall the ghost mirrors with a sign flip
/// (u_ghost = -u_cell), the reflection that makes the centred gradient at
/// the wall cell see the true wall shear; open sides copy.
fn pad2d_vel(f: &[f64], nx: usize, ny: usize, bc: &Boundaries2d) -> Vec<f64> {
    let w = nx + 2;
    let mut out = pad2d(f, nx, ny);
    let flip_s = matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall);
    let flip_n = matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall);
    if flip_s {
        for i in 0..w {
            out[i] = -out[i + w];
        }
    }
    if flip_n {
        for i in 0..w {
            out[(ny + 1) * w + i] = -out[ny * w + i];
        }
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

/// the face state on true positions: values and gradients are linearly
/// interpolated from the two cell centers to the actual face position; on
/// uniform spacing the face sits midway, so this is the plain average.
#[allow(clippy::too_many_arguments)]
fn face_grad_interp(
    ga: &Grads,
    gb: &Grads,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
    ca: f64,
    cb: f64,
    face: f64,
) -> FaceGrad {
    if (cb - ca).abs() < 1e-14 {
        return face_grad(ga, gb, ua, ub, va, vb);
    }
    let w = (face - ca) / (cb - ca);
    let mix = |a: f64, b: f64| a + w * (b - a);
    FaceGrad {
        u: mix(ua, ub),
        v: mix(va, vb),
        gux: mix(ga.ux, gb.ux),
        guy: mix(ga.uy, gb.uy),
        gvx: mix(ga.vx, gb.vx),
        gvy: mix(ga.vy, gb.vy),
        gtx: mix(ga.tx, gb.tx),
        gty: mix(ga.ty, gb.ty),
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

/// the viscosity/heat coefficients a face flux needs: kappa for the
/// molecular part, kappa_t for the eddy part (zero in the laminar limit).
#[derive(Clone, Copy)]
struct ViscCoeffs {
    kappa: f64,
    kappa_t: f64,
    gamma: f64,
}

/// the face flux with molecular and eddy viscosities contributing separately:
/// stresses use mu_lam + mu_t (face-averaged), while the heat flux splits so
/// each viscosity is paired with its own prandtl number (kappa for molecular,
/// kappa_t for turbulent; a zero kappa_t is the laminar limit).
fn face_flux_split(
    gf: &FaceGrad,
    mu_lam: f64,
    mu_t: &[f64],
    a: usize,
    b: usize,
    c: &ViscCoeffs,
) -> (f64, f64, f64, f64, f64) {
    if mu_t.is_empty() {
        return viscous_flux(gf, mu_lam, c.kappa, c.gamma);
    }
    let mu_f = mu_lam + 0.5 * (mu_t[a] + mu_t[b]);
    let (txx, txy, tyy, _qx, _qy) = viscous_flux(gf, mu_f, c.kappa, c.gamma);
    // the turbulent contribution replaces the molecular kappa pairing for the
    // eddy part: q = mu_lam*kappa*grad(T) + mu_t_face*kappa_t*grad(T).
    let mu_t_f = 0.5 * (mu_t[a] + mu_t[b]);
    let qx_total = (mu_lam * c.kappa + mu_t_f * c.kappa_t) * gf.gtx;
    let qy_total = (mu_lam * c.kappa + mu_t_f * c.kappa_t) * gf.gty;
    (txx, txy, tyy, qx_total, qy_total)
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

/// cell gradients on true positions: each cell takes the exact 3-point
/// lagrange derivative at its center over its two neighbors (ghost centers
/// reflected across the domain faces), quadratic-exact on any spacing and
/// the central difference when the axis is uniform.
fn cell_grads_metric(
    pu: &[f64],
    pv: &[f64],
    pt: &[f64],
    nx: usize,
    ny: usize,
    m: &MeshMetrics,
) -> Vec<Grads> {
    let w = nx + 2;
    let mut out = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let (ip, jp) = (i + 1, j + 1);
            let (xa, xb, xc) = m.x_nodes(i, nx);
            let (ya, yb, yc) = m.y_nodes(j, ny);
            let ux = quad_deriv_at(
                pu[jp * w + i],
                pu[jp * w + ip],
                pu[jp * w + (i + 2)],
                xa,
                xb,
                xc,
                xb,
            );
            let vx = quad_deriv_at(
                pv[jp * w + i],
                pv[jp * w + ip],
                pv[jp * w + (i + 2)],
                xa,
                xb,
                xc,
                xb,
            );
            let tx = quad_deriv_at(
                pt[jp * w + i],
                pt[jp * w + ip],
                pt[jp * w + (i + 2)],
                xa,
                xb,
                xc,
                xb,
            );
            let uy = quad_deriv_at(
                pu[j * w + ip],
                pu[jp * w + ip],
                pu[(j + 2) * w + ip],
                ya,
                yb,
                yc,
                yb,
            );
            let vy = quad_deriv_at(
                pv[j * w + ip],
                pv[jp * w + ip],
                pv[(j + 2) * w + ip],
                ya,
                yb,
                yc,
                yb,
            );
            let ty = quad_deriv_at(
                pt[j * w + ip],
                pt[jp * w + ip],
                pt[(j + 2) * w + ip],
                ya,
                yb,
                yc,
                yb,
            );
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

/// the per-axis coordinates the stretched stencils need: interior centers
/// plus the ghost centers reflected beyond the domain faces.
pub(crate) struct MeshMetrics {
    x_c: Vec<f64>,
    y_c: Vec<f64>,
    x_ghost: [f64; 2],
    y_ghost: [f64; 2],
}

impl MeshMetrics {
    /// build from a grid; only used when an axis is non-uniform.
    pub(crate) fn new(g: &Grid2d) -> Self {
        let x_c: Vec<f64> = (0..g.nx).map(|i| g.centers_x[i]).collect();
        let y_c: Vec<f64> = (0..g.ny).map(|j| g.centers_y[j * g.nx]).collect();
        let x_ghost = [
            2.0 * g.faces_x[0] - g.centers_x[0],
            2.0 * g.faces_x[g.nx] - g.centers_x[g.nx - 1],
        ];
        let y_ghost = [
            2.0 * g.faces_y[0] - g.centers_y[0],
            2.0 * g.faces_y[g.ny] - g.centers_y[(g.ny - 1) * g.nx],
        ];
        MeshMetrics {
            x_c,
            y_c,
            x_ghost,
            y_ghost,
        }
    }

    /// the three node positions for the x derivative of cell i.
    fn x_nodes(&self, i: usize, nx: usize) -> (f64, f64, f64) {
        let c = self.x_c[i];
        let l = if i == 0 {
            self.x_ghost[0]
        } else {
            self.x_c[i - 1]
        };
        let r = if i == nx - 1 {
            self.x_ghost[1]
        } else {
            self.x_c[i + 1]
        };
        (l, c, r)
    }

    /// the three node positions for the y derivative of cell j.
    fn y_nodes(&self, j: usize, ny: usize) -> (f64, f64, f64) {
        let c = self.y_c[j];
        let l = if j == 0 {
            self.y_ghost[0]
        } else {
            self.y_c[j - 1]
        };
        let r = if j == ny - 1 {
            self.y_ghost[1]
        } else {
            self.y_c[j + 1]
        };
        (l, c, r)
    }
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
    add_viscous_cells(state, g, gamma, mu, pr, None, dt)
}

/// the boundary-aware laminar form: no-slip wall faces carry the one-sided
/// wall shear so a boundary layer actually develops against a solid wall.
pub fn add_viscous_bc(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    mu: f64,
    pr: f64,
    bc: &Boundaries2d,
    dt: f64,
) -> Result<(), Error> {
    let turb = TurbCtx {
        mu_t: &[],
        pr_t: pr,
        bc,
    };
    add_viscous_cells(state, g, gamma, mu, pr, Some(turb), dt)
}

/// the derivative at y_eval of the quadratic through (y0, v0), (y1, v1),
/// (y2, v2): the exact 3-point lagrange derivative on true positions.
fn quad_deriv_at(v0: f64, v1: f64, v2: f64, y0: f64, y1: f64, y2: f64, y_eval: f64) -> f64 {
    let w0 = (2.0 * y_eval - y1 - y2) / ((y0 - y1) * (y0 - y2));
    let w1 = (2.0 * y_eval - y0 - y2) / ((y1 - y0) * (y1 - y2));
    let w2 = (2.0 * y_eval - y0 - y1) / ((y2 - y0) * (y2 - y1));
    v0 * w0 + v1 * w1 + v2 * w2
}

/// the wall-shear derivative on true positions: the quadratic through the
/// wall value and the first two cell centers, evaluated at the wall.
fn wall_shear_deriv(u0: f64, u1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    quad_deriv_at(0.0, u0, u1, y_wall, c0, c1, y_wall)
}

/// the wall-shear flux on one no-slip face, exact for any profile that is
/// quadratic near the wall (so both the linear sublayer and the laminar
/// parabola): the one-sided lagrange derivative through the wall value and
/// the first two cell centers; on a uniform grid it reduces to
/// du/dy|_wall = (9 u_0 - u_1) / (3 dy).
fn wall_shear_flux(u0: f64, u1: f64, v0: f64, v1: f64, mu_f: f64, dy: f64) -> (f64, f64) {
    let k = 1.0 / (3.0 * dy);
    (mu_f * (9.0 * u0 - u1) * k, mu_f * (9.0 * v0 - v1) * k)
}

/// the eddy-viscosity context for a turbulence-aware viscous pass: the
/// per-cell eddy viscosity, its prandtl number, and the boundaries (the
/// no-slip sides carry the one-sided wall shear).
#[derive(Clone, Copy)]
pub struct TurbCtx<'a> {
    /// per-cell eddy viscosity mu_t (empty means laminar).
    pub mu_t: &'a [f64],
    /// turbulent prandtl number for the eddy heat flux.
    pub pr_t: f64,
    /// boundary sides; no-slip walls get the wall shear flux.
    pub bc: &'a Boundaries2d,
}

/// the wall-cell y-gradient override with the same quadratic-consistent
/// stencil: du/dy at the wall cell = (3 u_0 + u_1) / (3 dy).
fn wall_cell_grad(g: &mut Grads, u0: f64, u1: f64, v0: f64, v1: f64, dy: f64, flip: f64) {
    let k = 1.0 / (3.0 * dy);
    g.uy = flip * (3.0 * u0 + u1) * k;
    g.vy = flip * (3.0 * v0 + v1) * k;
}

/// the stretched wall-cell y-gradient: the true derivative of the same
/// quadratic at the wall cell's center (the sign is physical, so unlike the
/// uniform shorthand there is no flip).
fn wall_cell_grad_st(
    g: &mut Grads,
    uv: (f64, f64),
    uv1: (f64, f64),
    y_wall: f64,
    c0: f64,
    c1: f64,
) {
    g.uy = quad_deriv_at(0.0, uv.0, uv1.0, y_wall, c0, c1, c0);
    g.vy = quad_deriv_at(0.0, uv.1, uv1.1, y_wall, c0, c1, c0);
}

/// the all-transmissive boundary set a laminar (wall-agnostic) pass uses.
static LAMINAR_BC: Boundaries2d = Boundaries2d {
    west: crate::solver2d::Bc2d::Transmissive,
    east: crate::solver2d::Bc2d::Transmissive,
    south: crate::solver2d::Bc2d::Transmissive,
    north: crate::solver2d::Bc2d::Transmissive,
};

/// the split-viscosity form of the viscous half-step: mu_lam is the constant
/// molecular dynamic viscosity and turb the eddy-viscosity context (None for
/// a laminar pass), each viscosity contributing to the heat flux through its
/// own prandtl number. a laminar call is byte-identical to add_viscous, and
/// a turbulent call applies the one-sided wall shear on no-slip sides.
pub fn add_viscous_cells(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    mu_lam: f64,
    pr: f64,
    turb: Option<TurbCtx>,
    dt: f64,
) -> Result<(), Error> {
    let bc = turb.map(|t| t.bc).unwrap_or(&LAMINAR_BC);
    let mu_t: &[f64] = turb.map(|t| t.mu_t).unwrap_or(&[]);
    let pr_t = turb.map(|t| t.pr_t).unwrap_or(pr);
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let turb = !mu_t.is_empty();
    if turb && mu_t.len() != nx * ny {
        return Err(Error::InvalidArgs);
    }
    let mu_max = if turb {
        mu_t.iter().cloned().fold(0.0f64, f64::max)
    } else {
        0.0
    };
    if (mu_lam <= 0.0 && mu_max <= 0.0) || pr <= 0.0 || pr <= 0.0 && !turb || pr_t <= 0.0 && turb {
        return Ok(());
    }
    if mu_lam <= 0.0 && !turb && dt <= 0.0 {
        return Ok(());
    }
    if dt <= 0.0 {
        return Ok(());
    }
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let t: Vec<f64> = p.iter().zip(&state.rho).map(|(pp, r)| pp / r).collect();
    let (pu, pv, pt) = (
        pad2d_vel(&u, nx, ny, bc),
        pad2d_vel(&v, nx, ny, bc),
        pad2d(&t, nx, ny),
    );
    // stretched axes take the exact position-aware stencils; uniform keeps
    // the scalar fast path bit-for-bit.
    let unif_x = (0..nx).all(|i| (g.dxs[i] - g.dxs[0]).abs() <= 1e-9 * g.dxs[0].abs().max(1e-12));
    let unif_y = (0..ny).all(|j| (g.dys[j] - g.dys[0]).abs() <= 1e-9 * g.dys[0].abs().max(1e-12));
    let m = if unif_x && unif_y {
        None
    } else {
        Some(MeshMetrics::new(g))
    };
    let grads = {
        let mut grads = match &m {
            None => cell_grads(&pu, &pv, &pt, nx, ny, 2.0 * g.dx, 2.0 * g.dy),
            Some(mm) => cell_grads_metric(&pu, &pv, &pt, nx, ny, mm),
        };
        // near-wall y-gradients use the quadratic-consistent one-sided
        // stencil so the wall cells balance exactly with the wall flux.
        let south = matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall);
        let north = matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall);
        if south {
            for i in 0..nx {
                let (a, b) = (idx(i, 0), idx(i, 1));
                match &m {
                    None => wall_cell_grad(&mut grads[a], u[a], u[b], v[a], v[b], g.dy, 1.0),
                    Some(mm) => wall_cell_grad_st(
                        &mut grads[a],
                        (u[a], v[a]),
                        (u[b], v[b]),
                        g.faces_y[0],
                        mm.y_c[0],
                        mm.y_c[1],
                    ),
                }
            }
        }
        if north {
            for i in 0..nx {
                let (a, b) = (idx(i, ny - 1), idx(i, ny - 2));
                match &m {
                    None => wall_cell_grad(&mut grads[a], u[a], u[b], v[a], v[b], g.dy, -1.0),
                    Some(mm) => wall_cell_grad_st(
                        &mut grads[a],
                        (u[a], v[a]),
                        (u[b], v[b]),
                        g.faces_y[ny],
                        mm.y_c[ny - 1],
                        mm.y_c[ny - 2],
                    ),
                }
            }
        }
        grads
    };
    let kappa = gamma / ((gamma - 1.0) * pr);
    let kappa_t = if turb {
        gamma / ((gamma - 1.0) * pr_t)
    } else {
        0.0
    };
    let coeffs = ViscCoeffs {
        kappa,
        kappa_t,
        gamma,
    };
    let (dtdx, dtdy) = (dt / g.dx, dt / g.dy);
    let (dxs_w, dys_w) = (&g.dxs, &g.dys);
    // x-face viscous flux (nx+1 faces per row); domain-edge faces carry zero flux.
    let mut fxm = vec![0.0; (nx + 1) * ny];
    let mut fym = vec![0.0; (nx + 1) * ny];
    let mut fe = vec![0.0; (nx + 1) * ny];
    for j in 0..ny {
        for f in 1..nx {
            let (a, b) = (idx(f - 1, j), idx(f, j));
            let gf = face_grad_interp(
                &grads[a],
                &grads[b],
                u[a],
                u[b],
                v[a],
                v[b],
                g.centers_x[a],
                g.centers_x[b],
                g.faces_x[f],
            );
            let (txx, txy, _, qx, _) = face_flux_split(&gf, mu_lam, mu_t, a, b, &coeffs);
            let k = f * ny + j;
            fxm[k] = txx;
            fym[k] = txy;
            fe[k] = gf.u * txx + gf.v * txy + qx;
        }
    }
    // y-face viscous flux (ny+1 faces per column); interior faces use the
    // centered pair, no-slip wall faces carry the one-sided wall shear, and
    // open domain-edge faces carry zero flux.
    let mut gym_x = vec![0.0; nx * (ny + 1)];
    let mut gym_y = vec![0.0; nx * (ny + 1)];
    let mut gye = vec![0.0; nx * (ny + 1)];
    let y_c0 = g.centers_y[0];
    let y_c1 = g.centers_y[g.nx];
    let y_c_n = g.centers_y[(g.ny - 1) * g.nx];
    let y_c_nm1 = g.centers_y[(g.ny - 2) * g.nx];
    for i in 0..nx {
        for f in 1..ny {
            let (a, b) = (idx(i, f - 1), idx(i, f));
            let cy_a = g.centers_y[a];
            let cy_b = g.centers_y[b];
            let gf = face_grad_interp(
                &grads[a],
                &grads[b],
                u[a],
                u[b],
                v[a],
                v[b],
                cy_a,
                cy_b,
                g.faces_y[f],
            );
            let (_, txy, tyy, _, qy) = face_flux_split(&gf, mu_lam, mu_t, a, b, &coeffs);
            let k = i * (ny + 1) + f;
            gym_x[k] = txy;
            gym_y[k] = tyy;
            gye[k] = gf.u * txy + gf.v * tyy + qy;
        }
        // bottom face f=0: no-slip walls shear the flow; other bcs free-slip.
        if matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall) {
            let (a, b) = (idx(i, 0), idx(i, 1));
            let mu_f = if mu_t.is_empty() {
                mu_lam
            } else {
                mu_lam + mu_t[a]
            };
            let (txy, tyy) = if unif_y {
                wall_shear_flux(u[a], u[b], v[a], v[b], mu_f, g.dy)
            } else {
                (
                    mu_f * wall_shear_deriv(u[a], u[b], g.faces_y[0], y_c0, y_c1),
                    mu_f * wall_shear_deriv(v[a], v[b], g.faces_y[0], y_c0, y_c1),
                )
            };
            let k = i * (ny + 1);
            gym_x[k] = txy;
            gym_y[k] = tyy;
            gye[k] = u[a] * txy + v[a] * tyy;
        }
        // top face f=ny. the outward normal is -y, so the shear the wall
        // exerts enters the flux difference with the opposite sign.
        if matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall) {
            let (a, b) = (idx(i, ny - 1), idx(i, ny - 2));
            let mu_f = if mu_t.is_empty() {
                mu_lam
            } else {
                mu_lam + mu_t[a]
            };
            // the uniform shorthand returns the south-oriented magnitude, so
            // the stretched branch mirrors it: negate the sign-true derivative.
            let (txy, tyy) = if unif_y {
                wall_shear_flux(u[a], u[b], v[a], v[b], mu_f, g.dy)
            } else {
                (
                    -mu_f * wall_shear_deriv(u[a], u[b], g.faces_y[ny], y_c_n, y_c_nm1),
                    -mu_f * wall_shear_deriv(v[a], v[b], g.faces_y[ny], y_c_n, y_c_nm1),
                )
            };
            let k = i * (ny + 1) + ny;
            gym_x[k] = -txy;
            gym_y[k] = -tyy;
            gye[k] = -(u[a] * txy + v[a] * tyy);
        }
    }
    for (j, &dyw) in dys_w.iter().enumerate() {
        for (i, &dxw) in dxs_w.iter().enumerate() {
            let c = idx(i, j);
            // x-face divergence: faces i and i+1 bound cell i.
            let (a, b) = (i * ny + j, (i + 1) * ny + j);
            let dxi = if unif_x { dtdx } else { dt / dxw };
            state.mx[c] += dxi * (fxm[b] - fxm[a]);
            state.my[c] += dxi * (fym[b] - fym[a]);
            state.e[c] += dxi * (fe[b] - fe[a]);
            // y-face divergence: faces j and j+1 bound cell j.
            let (d, e) = (i * (ny + 1) + j, i * (ny + 1) + j + 1);
            let dyj = if unif_y { dtdy } else { dt / dyw };
            state.mx[c] += dyj * (gym_x[e] - gym_x[d]);
            state.my[c] += dyj * (gym_y[e] - gym_y[d]);
            state.e[c] += dyj * (gye[e] - gye[d]);
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

/// the time-step cap from the sa diffusion: explicit stability wants
/// dt < dx^2 / (2 * (nu + nu_tilde_max)/sigma) per dimension, matching the
/// shape of the viscous bound so the same cfl-style safety factor applies.
pub fn sa_dt_cap(g: &Grid2d, nu_max: f64) -> f64 {
    let dcoef = nu_max / crate::sa::SIGMA;
    if dcoef <= 0.0 {
        return f64::INFINITY;
    }
    0.25 * g.dx.min(g.dy).powi(2) / dcoef
}

/// one coupled heun step: mean flow (euler + viscous with eddy viscosity)
/// and the sa transport, each stage using the other's updated state.
///
/// stage 1 advances both from the current state; stage 2 re-evaluates the
/// fluxes at the stage-1 state; the heun average closes the step. the
/// eddy-viscosity coupling uses stage-local nu_tilde so the momentum
/// diffusion and the transported field stay consistent.
pub fn advance2d_sa_rk2(
    state: &mut ConservedState2d,
    turb: &mut crate::turb2d::TurbState,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
) -> Result<f64, Error> {
    use crate::sa::eddy_viscosity;
    use crate::solver2d::advance2d_capped;
    use crate::turb2d::advance_turb;
    let n = g.nx * g.ny;
    let sa = turb.params;
    let d = &turb.d;
    if turb.nu_tilde.len() != n || d.len() != n {
        return Err(Error::InvalidArgs);
    }
    let nu = sa.mu; // molecular dynamic viscosity
                    // the sa cap uses the largest diffusivity in the field.
    let nt_max = turb
        .nu_tilde
        .iter()
        .cloned()
        .fold(sa.nu_tilde_inf, f64::max);
    let nu_max = nu / state.rho.iter().cloned().fold(f64::INFINITY, f64::min) + nt_max;
    let dt_cap = sa_dt_cap(g, nu_max);
    // eddy viscosity per cell from the current nu_tilde.
    let mu_t_of = |t: &crate::turb2d::TurbState, rho: &[f64]| -> Vec<f64> {
        t.nu_tilde
            .iter()
            .zip(rho)
            .map(|(&nt, &r)| eddy_viscosity(r, nt, nu / r.max(1e-12)))
            .collect()
    };
    // stage 1: euler + viscous + turb from (state, turb).
    let mut s1 = state.clone();
    let mut t1 = turb.clone();
    let (dt1, _) = advance2d_capped(&mut s1, g, gamma, cfl, muscl, bc, dt_cap)?;
    let mu_t1 = mu_t_of(&t1, &s1.rho);
    let ctx1 = TurbCtx {
        mu_t: &mu_t1,
        pr_t: sa.pr_t,
        bc,
    };
    add_viscous_cells(&mut s1, g, gamma, nu, sa.pr, Some(ctx1), dt1)?;
    advance_turb(&mut t1, &s1, g, bc, dt1)?;
    // stage 2 from the stage-1 state (heun's second evaluation).
    let (dt2, _) = advance2d_capped(&mut s1, g, gamma, cfl, muscl, bc, dt_cap)?;
    let mu_t1b = mu_t_of(&t1, &s1.rho);
    let ctx2 = TurbCtx {
        mu_t: &mu_t1b,
        pr_t: sa.pr_t,
        bc,
    };
    add_viscous_cells(&mut s1, g, gamma, nu, sa.pr, Some(ctx2), dt2)?;
    advance_turb(&mut t1, &s1, g, bc, dt2)?;
    // heun average of the original and the twice-advanced state.
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
        turb.nu_tilde[k] = 0.5 * turb.nu_tilde[k] + 0.5 * t1.nu_tilde[k];
    }
    Ok(dt1)
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
    add_viscous_bc(&mut s1, g, gamma, mu, pr, bc, dt1)?;
    let (dt2, _) = advance2d_capped(&mut s1, g, gamma, cfl, muscl, bc, dt_visc)?;
    add_viscous_bc(&mut s1, g, gamma, mu, pr, bc, dt2)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

//! 2d euler finite-volume step: muscl reconstruction, hllc flux, conservative update.

use crate::eos2d::eos_pressure2d;
use crate::grid2d::Grid2d;
use crate::hllc2d::{hllc_flux, FacePrim};
use crate::state2d::{cons_to_prim2d, ConservedState2d};
use crate::Error;

/// the van leer limiter: the harmonic mean, zero when the slopes disagree.
fn van_leer(dm: f64, dp: f64) -> f64 {
    if dm * dp <= 0.0 {
        0.0
    } else {
        2.0 * dm * dp / (dm + dp)
    }
}

/// muscl face states (left, right; each length n+1) from a transmissive-padded slice.
///
/// the padded slice p has length n+2 with p[0] and p[n+1] the ghost values; a
/// zero slope budget (muscl=false) reproduces a piecewise-constant, first-order
/// scheme and is used to measure the order the limiter buys.
fn face_states(p: &[f64], n: usize, muscl: bool) -> (Vec<f64>, Vec<f64>) {
    let mut d = vec![0.0; n];
    if muscl {
        for i in 0..n {
            let dm = p[i + 1] - p[i];
            let dp = p[i + 2] - p[i + 1];
            d[i] = van_leer(dm, dp);
        }
    }
    let mut fl = vec![0.0; n + 1];
    let mut fr = vec![0.0; n + 1];
    fl[0] = p[0];
    fr[0] = p[1] - 0.5 * d[0];
    fl[n] = p[n] + 0.5 * d[n - 1];
    fr[n] = p[n + 1];
    for f in 1..n {
        fl[f] = p[f] + 0.5 * d[f - 1];
        fr[f] = p[f + 1] - 0.5 * d[f];
    }
    (fl, fr)
}

/// transpose copies (rho,u,v,p) of a flat row-major slab into per-axis padded strips.
///
/// returns two closures avoided: this is just a helper collecting the face flux
/// arrays for the x (vertical) and y (horizontal) sweeps of the update.
#[derive(Clone)]
struct Sweep {
    mass: Vec<f64>,
    mx: Vec<f64>,
    my: Vec<f64>,
    e: Vec<f64>,
}

/// one conservative step for the 2d state over the grid.
pub fn advance2d(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
) -> Result<(f64, f64), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let dt = cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, nx * ny);

    // vertical (x) faces: one flux per row at each of the nx+1 vertical faces.
    let mut fx = vec![
        Sweep {
            mass: vec![0.0; nx + 1],
            mx: vec![0.0; nx + 1],
            my: vec![0.0; nx + 1],
            e: vec![0.0; nx + 1],
        };
        ny
    ];
    for j in 0..ny {
        let mut pr = vec![0.0; nx];
        let mut pu = vec![0.0; nx];
        let mut pv = vec![0.0; nx];
        let mut pp = vec![0.0; nx];
        for i in 0..nx {
            pr[i] = state.rho[idx(i, j)];
            pu[i] = u[idx(i, j)];
            pv[i] = v[idx(i, j)];
            pp[i] = p[idx(i, j)];
        }
        let pad = |c: &[f64]| -> Vec<f64> {
            let mut out = vec![0.0; nx + 2];
            out[0] = c[0];
            out[nx + 1] = c[nx - 1];
            out[1..=nx].copy_from_slice(c);
            out
        };
        let (rl, rr) = face_states(&pad(&pr), nx, muscl);
        let (ul, ur) = face_states(&pad(&pu), nx, muscl);
        let (vl, vr) = face_states(&pad(&pv), nx, muscl);
        let (pl, prr) = face_states(&pad(&pp), nx, muscl);
        for f in 0..=nx {
            let q = hllc_flux(
                gamma,
                FacePrim {
                    rho: rl[f],
                    u: ul[f],
                    v: vl[f],
                    p: pl[f],
                },
                FacePrim {
                    rho: rr[f],
                    u: ur[f],
                    v: vr[f],
                    p: prr[f],
                },
                0,
            );
            fx[j].mass[f] = q.mass;
            fx[j].mx[f] = q.mx;
            fx[j].my[f] = q.my;
            fx[j].e[f] = q.e;
        }
    }

    // horizontal (y) faces: one flux per column at each of the ny+1 faces.
    let mut fy = vec![
        Sweep {
            mass: vec![0.0; ny + 1],
            mx: vec![0.0; ny + 1],
            my: vec![0.0; ny + 1],
            e: vec![0.0; ny + 1],
        };
        nx
    ];
    for i in 0..nx {
        let mut cr = vec![0.0; ny];
        let mut cu = vec![0.0; ny];
        let mut cv = vec![0.0; ny];
        let mut cp = vec![0.0; ny];
        for j in 0..ny {
            cr[j] = state.rho[idx(i, j)];
            cu[j] = u[idx(i, j)];
            cv[j] = v[idx(i, j)];
            cp[j] = p[idx(i, j)];
        }
        let pad = |c: &[f64]| -> Vec<f64> {
            let mut out = vec![0.0; ny + 2];
            out[0] = c[0];
            out[ny + 1] = c[ny - 1];
            out[1..=ny].copy_from_slice(c);
            out
        };
        let (rl, rr) = face_states(&pad(&cr), ny, muscl);
        let (ul, ur) = face_states(&pad(&cu), ny, muscl);
        let (vl, vr) = face_states(&pad(&cv), ny, muscl);
        let (pl, prr) = face_states(&pad(&cp), ny, muscl);
        for f in 0..=ny {
            let q = hllc_flux(
                gamma,
                FacePrim {
                    rho: rl[f],
                    u: ul[f],
                    v: vl[f],
                    p: pl[f],
                },
                FacePrim {
                    rho: rr[f],
                    u: ur[f],
                    v: vr[f],
                    p: prr[f],
                },
                1,
            );
            fy[i].mass[f] = q.mass;
            fy[i].mx[f] = q.mx;
            fy[i].my[f] = q.my;
            fy[i].e[f] = q.e;
        }
    }

    let dtdx = dt / g.dx;
    let dtdy = dt / g.dy;
    let mut resid = 0.0f64;
    for (j, fxrow) in fx.iter().enumerate() {
        for (i, fycol) in fy.iter().enumerate() {
            let dr = dtdx * (fxrow.mass[i + 1] - fxrow.mass[i])
                + dtdy * (fycol.mass[j + 1] - fycol.mass[j]);
            let dm =
                dtdx * (fxrow.mx[i + 1] - fxrow.mx[i]) + dtdy * (fycol.mx[j + 1] - fycol.mx[j]);
            let dn =
                dtdx * (fxrow.my[i + 1] - fxrow.my[i]) + dtdy * (fycol.my[j + 1] - fycol.my[j]);
            let de = dtdx * (fxrow.e[i + 1] - fxrow.e[i]) + dtdy * (fycol.e[j + 1] - fycol.e[j]);
            state.rho[idx(i, j)] -= dr;
            state.mx[idx(i, j)] -= dm;
            state.my[idx(i, j)] -= dn;
            state.e[idx(i, j)] -= de;
            resid = resid
                .max(dr.abs())
                .max(dm.abs())
                .max(dn.abs())
                .max(de.abs());
        }
    }
    Ok((dt, resid))
}

/// the largest fast-characteristic speed over the cells, bounding the wave speeds.
fn smax_of(u: &[f64], v: &[f64], p: &[f64], rho: &[f64], gamma: f64, n: usize) -> f64 {
    let mut s = 0.0f64;
    for k in 0..n {
        let a = (gamma * p[k].max(1e-12) / rho[k]).sqrt();
        s = s.max(u[k].abs().max(v[k].abs()) + a);
    }
    s
}

/// second-order two-stage (heun) time step: explicit predictor-corrector so
/// muscl's spatial accuracy is not masked by first-order time integration.
pub fn advance2d_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
) -> Result<(f64, f64), Error> {
    // heun: s1 = st + dt l(st); s2 = s1 + dt l(s1); st = 0.5(st + s2).
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let dt = cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, g.nx * g.ny);
    let mut s1 = state.clone();
    advance2d(&mut s1, g, gamma, cfl, muscl)?;
    advance2d(&mut s1, g, gamma, cfl, muscl)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt, 0.0))
}

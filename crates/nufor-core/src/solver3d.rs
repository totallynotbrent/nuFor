//! 3d euler finite-volume step: muscl reconstruction, hllc flux in three axes.

use crate::eos3d::eos_pressure3d;
use crate::grid3d::Grid3d;
use crate::hllc3d::{hllc_flux3, FacePrim3};
use crate::state3d::{cons_to_prim3d, ConservedState3d};
use crate::Error;

/// the van leer limiter: harmonic mean, zero when the slopes disagree.
fn van_leer(dm: f64, dp: f64) -> f64 {
    if dm * dp <= 0.0 {
        0.0
    } else {
        2.0 * dm * dp / (dm + dp)
    }
}

/// the largest fast-characteristic speed, bounding the 3d wave speeds.
fn smax3d(u: &[f64], v: &[f64], w: &[f64], p: &[f64], rho: &[f64], gamma: f64) -> f64 {
    let mut s = 0.0f64;
    for k in 0..rho.len() {
        let a = (gamma * p[k].max(1e-12) / rho[k]).sqrt();
        s = s.max(u[k].abs().max(v[k].abs()).max(w[k].abs()) + a);
    }
    s
}

/// a primitive strip padded with one transmissive ghost at each end.
fn padded(strip: &[f64]) -> Vec<f64> {
    let n = strip.len();
    let mut out = vec![0.0; n + 2];
    out[0] = strip[0];
    out[n + 1] = strip[n - 1];
    out[1..=n].copy_from_slice(strip);
    out
}

/// muscl face states (left, right; each length n+1) from a padded strip.
fn face_states(p: &[f64], n: usize, muscl: bool) -> (Vec<f64>, Vec<f64>) {
    let mut d = vec![0.0; n];
    if muscl {
        for i in 0..n {
            d[i] = van_leer(p[i + 1] - p[i], p[i + 2] - p[i + 1]);
        }
    }
    let mut fl = vec![0.0; n + 1];
    let mut fr = vec![0.0; n + 1];
    for f in 0..=n {
        let l = if f == 0 { p[0] } else { p[f] + 0.5 * d[f - 1] };
        let r = if f == n {
            p[n + 1]
        } else {
            p[f + 1] - 0.5 * d[f]
        };
        fl[f] = l;
        fr[f] = r;
    }
    (fl, fr)
}

/// five parallel face-flux arrays returned from a 1d strip sweep.
type StripFlux = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

/// the five primitive strips of one axis-aligned line of cells.
struct Strips<'a> {
    rho: &'a [f64],
    u: &'a [f64],
    v: &'a [f64],
    w: &'a [f64],
    p: &'a [f64],
}

/// the hllc flux vector across the n+1 faces of a 1d primitive strip.
fn strip_flux(s: &Strips, gamma: f64, axis: usize, muscl: bool) -> StripFlux {
    let n = s.rho.len();
    let (pr, pu, pv, pw, pp) = (
        padded(s.rho),
        padded(s.u),
        padded(s.v),
        padded(s.w),
        padded(s.p),
    );
    let (rl, rr) = face_states(&pr, n, muscl);
    let (ul, ur) = face_states(&pu, n, muscl);
    let (vl, vr) = face_states(&pv, n, muscl);
    let (wl, wr) = face_states(&pw, n, muscl);
    let (pl, prr) = face_states(&pp, n, muscl);
    let (mut m, mut mx, mut my, mut mz, mut me) = (
        vec![0.0; n + 1],
        vec![0.0; n + 1],
        vec![0.0; n + 1],
        vec![0.0; n + 1],
        vec![0.0; n + 1],
    );
    for f in 0..=n {
        let q = hllc_flux3(
            gamma,
            FacePrim3 {
                rho: rl[f],
                u: ul[f],
                v: vl[f],
                w: wl[f],
                p: pl[f],
            },
            FacePrim3 {
                rho: rr[f],
                u: ur[f],
                v: vr[f],
                w: wr[f],
                p: prr[f],
            },
            axis,
        );
        m[f] = q.mass;
        mx[f] = q.mx;
        my[f] = q.my;
        mz[f] = q.mz;
        me[f] = q.e;
    }
    (m, mx, my, mz, me)
}

/// one conservative step for the 3d state over the grid.
pub fn advance3d(
    state: &mut ConservedState3d,
    g: &Grid3d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
) -> Result<(f64, f64), Error> {
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let idx = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
    let (u, v, w, et) = cons_to_prim3d(&state.rho, &state.mx, &state.my, &state.mz, &state.e)?;
    let p = eos_pressure3d(gamma, &state.rho, &et, &u, &v, &w)?;
    let dt = cfl * g.dx.min(g.dy).min(g.dz) / smax3d(&u, &v, &w, &p, &state.rho, gamma);

    // x-face fluxes: (nx+1) per (j,k) slab.
    let (mut fx_m, mut fx_x, mut fx_y, mut fx_z, mut fx_e) = (
        vec![0.0; (nx + 1) * ny * nz],
        vec![0.0; (nx + 1) * ny * nz],
        vec![0.0; (nx + 1) * ny * nz],
        vec![0.0; (nx + 1) * ny * nz],
        vec![0.0; (nx + 1) * ny * nz],
    );
    let (rm, rx, ry, rz, re) = (
        state.rho.as_slice(),
        u.as_slice(),
        v.as_slice(),
        w.as_slice(),
        p.as_slice(),
    );
    for k in 0..nz {
        for j in 0..ny {
            let (mut sr, mut su, mut sv, mut sw, mut sp) = (
                vec![0.0; nx],
                vec![0.0; nx],
                vec![0.0; nx],
                vec![0.0; nx],
                vec![0.0; nx],
            );
            for i in 0..nx {
                let c = idx(i, j, k);
                sr[i] = rm[c];
                su[i] = rx[c];
                sv[i] = ry[c];
                sw[i] = rz[c];
                sp[i] = re[c];
            }
            let (m, mx, my, mz, me) = strip_flux(
                &Strips {
                    rho: &sr,
                    u: &su,
                    v: &sv,
                    w: &sw,
                    p: &sp,
                },
                gamma,
                0,
                muscl,
            );
            let base = (k * ny + j) * (nx + 1);
            fx_m[base..base + nx + 1].copy_from_slice(&m);
            fx_x[base..base + nx + 1].copy_from_slice(&mx);
            fx_y[base..base + nx + 1].copy_from_slice(&my);
            fx_z[base..base + nx + 1].copy_from_slice(&mz);
            fx_e[base..base + nx + 1].copy_from_slice(&me);
        }
    }
    // y-face fluxes: (ny+1) per (i,k) slab.
    let (mut fy_m, mut fy_x, mut fy_y, mut fy_z, mut fy_e) = (
        vec![0.0; (ny + 1) * nx * nz],
        vec![0.0; (ny + 1) * nx * nz],
        vec![0.0; (ny + 1) * nx * nz],
        vec![0.0; (ny + 1) * nx * nz],
        vec![0.0; (ny + 1) * nx * nz],
    );
    for k in 0..nz {
        for i in 0..nx {
            let (mut sr, mut su, mut sv, mut sw, mut sp) = (
                vec![0.0; ny],
                vec![0.0; ny],
                vec![0.0; ny],
                vec![0.0; ny],
                vec![0.0; ny],
            );
            for j in 0..ny {
                let c = idx(i, j, k);
                sr[j] = rm[c];
                su[j] = rx[c];
                sv[j] = ry[c];
                sw[j] = rz[c];
                sp[j] = re[c];
            }
            let (m, mx, my, mz, me) = strip_flux(
                &Strips {
                    rho: &sr,
                    u: &su,
                    v: &sv,
                    w: &sw,
                    p: &sp,
                },
                gamma,
                1,
                muscl,
            );
            let base = (k * nx + i) * (ny + 1);
            fy_m[base..base + ny + 1].copy_from_slice(&m);
            fy_x[base..base + ny + 1].copy_from_slice(&mx);
            fy_y[base..base + ny + 1].copy_from_slice(&my);
            fy_z[base..base + ny + 1].copy_from_slice(&mz);
            fy_e[base..base + ny + 1].copy_from_slice(&me);
        }
    }
    // z-face fluxes: (nz+1) per (i,j) slab.
    let (mut fz_m, mut fz_x, mut fz_y, mut fz_z, mut fz_e) = (
        vec![0.0; (nz + 1) * nx * ny],
        vec![0.0; (nz + 1) * nx * ny],
        vec![0.0; (nz + 1) * nx * ny],
        vec![0.0; (nz + 1) * nx * ny],
        vec![0.0; (nz + 1) * nx * ny],
    );
    for j in 0..ny {
        for i in 0..nx {
            let (mut sr, mut su, mut sv, mut sw, mut sp) = (
                vec![0.0; nz],
                vec![0.0; nz],
                vec![0.0; nz],
                vec![0.0; nz],
                vec![0.0; nz],
            );
            for k in 0..nz {
                let c = idx(i, j, k);
                sr[k] = rm[c];
                su[k] = rx[c];
                sv[k] = ry[c];
                sw[k] = rz[c];
                sp[k] = re[c];
            }
            let (m, mx, my, mz, me) = strip_flux(
                &Strips {
                    rho: &sr,
                    u: &su,
                    v: &sv,
                    w: &sw,
                    p: &sp,
                },
                gamma,
                2,
                muscl,
            );
            let base = (j * nx + i) * (nz + 1);
            fz_m[base..base + nz + 1].copy_from_slice(&m);
            fz_x[base..base + nz + 1].copy_from_slice(&mx);
            fz_y[base..base + nz + 1].copy_from_slice(&my);
            fz_z[base..base + nz + 1].copy_from_slice(&mz);
            fz_e[base..base + nz + 1].copy_from_slice(&me);
        }
    }

    let (dtdx, dtdy, dtdz) = (dt / g.dx, dt / g.dy, dt / g.dz);
    let mut resid = 0.0f64;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let c = idx(i, j, k);
                // x-face pair i and i+1, y-face pair j and j+1, z-face pair k and k+1.
                let (a, b) = ((k * ny + j) * (nx + 1) + i, (k * ny + j) * (nx + 1) + i + 1);
                let (d, e) = ((k * nx + i) * (ny + 1) + j, (k * nx + i) * (ny + 1) + j + 1);
                let (f, h) = ((j * nx + i) * (nz + 1) + k, (j * nx + i) * (nz + 1) + k + 1);
                let dri = dtdx * (fx_m[b] - fx_m[a])
                    + dtdy * (fy_m[e] - fy_m[d])
                    + dtdz * (fz_m[h] - fz_m[f]);
                let dri_x = dtdx * (fx_x[b] - fx_x[a])
                    + dtdy * (fy_x[e] - fy_x[d])
                    + dtdz * (fz_x[h] - fz_x[f]);
                let dri_y = dtdx * (fx_y[b] - fx_y[a])
                    + dtdy * (fy_y[e] - fy_y[d])
                    + dtdz * (fz_y[h] - fz_y[f]);
                let dri_z = dtdx * (fx_z[b] - fx_z[a])
                    + dtdy * (fy_z[e] - fy_z[d])
                    + dtdz * (fz_z[h] - fz_z[f]);
                let dri_e = dtdx * (fx_e[b] - fx_e[a])
                    + dtdy * (fy_e[e] - fy_e[d])
                    + dtdz * (fz_e[h] - fz_e[f]);
                state.rho[c] -= dri;
                state.mx[c] -= dri_x;
                state.my[c] -= dri_y;
                state.mz[c] -= dri_z;
                state.e[c] -= dri_e;
                resid = resid.max(dri.abs()).max(dri_x.abs()).max(dri_y.abs());
            }
        }
    }
    Ok((dt, resid))
}

/// second-order two-stage (heun) 3d step so muscl's accuracy is not masked.
pub fn advance3d_rk2(
    state: &mut ConservedState3d,
    g: &Grid3d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
) -> Result<(f64, f64), Error> {
    let mut s1 = state.clone();
    let (dt1, _) = advance3d(&mut s1, g, gamma, cfl, muscl)?;
    advance3d(&mut s1, g, gamma, cfl, muscl)?;
    for k in 0..g.nx * g.ny * g.nz {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.mz[k] = 0.5 * state.mz[k] + 0.5 * s1.mz[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

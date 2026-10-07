//! the axisymmetric euler march: the 2d planar flux machinery with the
//! geometric source terms of a body of revolution.
//!
//! in cylindrical coordinates (x axial, r radial) the axisymmetric euler
//! equations keep the same face fluxes as the planar 2d case per unit
//! area, but the finite volume is an annulus, so the radial flux
//! difference carries the face radii and the radial momentum feels a
//! pressure source:
//!
//!   d/dt (rho) = -[ d(f_x)/dx + (1/r) d(r f_r)/dr ]
//!   d/dt (m_x) = -[ d(f_x^mx)/dx + (1/r) d(r f_r^mx)/dr ]
//!   d/dt (m_r) = -[ d(f_x^my)/dx + (1/r) d(r f_r^my)/dr ] + p / r
//!   d/dt (e)   = -[ d(f_x^e)/dx + (1/r) d(r f_r^e)/dr ]
//!
//! discretized on cell (i,j) with center radius r_c and radial face
//! radii r_- = r_j, r_+ = r_{j+1}:
//!
//!   u_t = -( f_x difference + (r_+ F_r[j+1] - r_- F_r[j]) / (r_c dr) )
//!         + (0, 0, p / r_c, 0)
//!
//! two consequences worth writing down: the south face of the domain is
//! the symmetry axis, and a slip-wall ghost (v -> -v, others copied)
//! makes the axis flux pure pressure transport with zero mass flow,
//! exactly the axis condition; and at the axis face r_- = 0, so the
//! first radial row only sees its outer face.
//!
//! the march reuses the planar flux machinery unchanged (fluxes per
//! unit area are identical in both formulations) and rewrites only the
//! update.

use crate::eos2d::eos_pressure2d;
use crate::grid2d::Grid2d;
use crate::hllc2d::{hllc_flux, FacePrim};
use crate::state2d::{cons_to_prim2d, ConservedState2d};
use crate::Error;

use crate::solver2d::{ghost_value, Boundaries2d};

/// the van leer limiter: the harmonic mean, zero when the slopes disagree.
fn van_leer(dm: f64, dp: f64) -> f64 {
    if dm * dp <= 0.0 {
        0.0
    } else {
        2.0 * dm * dp / (dm + dp)
    }
}

/// muscl face states (left, right; each length n+1) from a padded slice
/// of length n+2, identical to the planar march's version.
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

/// the largest fast-characteristic speed over the cells.
fn smax_of(u: &[f64], v: &[f64], p: &[f64], rho: &[f64], gamma: f64, n: usize) -> f64 {
    let mut s = 0.0f64;
    for k in 0..n {
        let a = (gamma * p[k].max(1e-12) / rho[k]).sqrt();
        s = s.max(u[k].abs().max(v[k].abs()) + a);
    }
    s
}

/// face flux arrays for one sweep direction.
struct Sweep {
    mass: Vec<f64>,
    mx: Vec<f64>,
    my: Vec<f64>,
    e: Vec<f64>,
}

impl Sweep {
    fn new(n: usize) -> Self {
        Sweep {
            mass: vec![0.0; n + 1],
            mx: vec![0.0; n + 1],
            my: vec![0.0; n + 1],
            e: vec![0.0; n + 1],
        }
    }
}

/// one axisymmetric euler step; compose with two half-steps for heun the
/// same way `advance2d_rk2` composes `advance2d`.
///
/// the south boundary is the symmetry axis: pass `Bc2d::SlipWall` there
/// (the ghost's negated v is exactly the axis condition).
pub fn advance2d_axi(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let dt = cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, nx * ny);

    let st: &ConservedState2d = state;
    let mut fx: Vec<Sweep> = (0..ny).map(|_| Sweep::new(nx)).collect();
    let mut fy: Vec<Sweep> = (0..nx).map(|_| Sweep::new(ny)).collect();

    for j in 0..ny {
        let (mut pr, mut pu, mut pv, mut pp) =
            (vec![0.0; nx], vec![0.0; nx], vec![0.0; nx], vec![0.0; nx]);
        for i in 0..nx {
            pr[i] = st.rho[idx(i, j)];
            pu[i] = u[idx(i, j)];
            pv[i] = v[idx(i, j)];
            pp[i] = p[idx(i, j)];
        }
        let y_row = g.centers_y[j * nx];
        let pad = |c: &[f64], var: usize| -> Vec<f64> {
            let mut out = vec![0.0; nx + 2];
            out[0] = ghost_value(&bc.west, c[0], var, 1, y_row);
            out[nx + 1] = ghost_value(&bc.east, c[nx - 1], var, 1, y_row);
            out[1..=nx].copy_from_slice(c);
            out
        };
        let flux_row = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            face_states(&padded, nx, muscl)
        };
        let (rl, rr) = flux_row(&pr, 0);
        let (ul, ur) = flux_row(&pu, 1);
        let (vl, vr) = flux_row(&pv, 2);
        let (pl, prr) = flux_row(&pp, 3);
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

    for i in 0..nx {
        let (mut cr, mut cu, mut cv, mut cp) =
            (vec![0.0; ny], vec![0.0; ny], vec![0.0; ny], vec![0.0; ny]);
        for j in 0..ny {
            cr[j] = st.rho[idx(i, j)];
            cu[j] = u[idx(i, j)];
            cv[j] = v[idx(i, j)];
            cp[j] = p[idx(i, j)];
        }
        let x_col = g.centers_x[i];
        let x_ghost_s = 2.0 * g.faces_x[0] - x_col;
        let x_ghost_n = 2.0 * g.faces_x[nx] - x_col;
        let pad = |c: &[f64], var: usize| -> Vec<f64> {
            let mut out = vec![0.0; ny + 2];
            out[0] = ghost_value(&bc.south, c[0], var, 2, x_ghost_s);
            out[ny + 1] = ghost_value(&bc.north, c[ny - 1], var, 2, x_ghost_n);
            out[1..=ny].copy_from_slice(c);
            out
        };
        let flux_col = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            face_states(&padded, ny, muscl)
        };
        let (rl, rr) = flux_col(&cr, 0);
        let (ul, ur) = flux_col(&cu, 1);
        let (vl, vr) = flux_col(&cv, 2);
        let (pl, prr) = flux_col(&cp, 3);
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

    // the annular update: r-weighted radial flux difference plus the
    // radial pressure source.
    let mut resid = 0.0f64;
    for j in 0..ny {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        for (i, fycol) in fy.iter().enumerate() {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let radial = |f: &Vec<f64>| (rp * f[j + 1] - rm * f[j]) / (rc * dr);
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i]) + dt * radial(&fycol.mass));
            let dmx = -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dt * radial(&fycol.mx));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dt * radial(&fycol.my)) + dt * p[k] / rc;
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dt * radial(&fycol.e));
            state.rho[k] += drho;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            resid = resid
                .max(drho.abs())
                .max(dmx.abs())
                .max(dmy.abs())
                .max(de.abs());
        }
    }
    Ok((dt, resid))
}

/// heun two-stage step over the axisymmetric step, mirroring
/// `advance2d_rk2`'s composition over `advance2d`: the increment comes
/// from the pre-step state (the first advance call returns it).
pub fn advance2d_axi_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
) -> Result<(f64, f64), Error> {
    let (dt, _) = {
        let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
        (
            cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, g.nx * g.ny),
            0.0,
        )
    };
    let mut s1 = state.clone();
    advance2d_axi(&mut s1, g, gamma, cfl, muscl, bc)?;
    advance2d_axi(&mut s1, g, gamma, cfl, muscl, bc)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt, 0.0))
}

//! 2d euler finite-volume step: muscl reconstruction, hllc flux, conservative update.

use crate::blasius::BlasiusProfile;
use crate::eos2d::eos_pressure2d;
use crate::grid2d::Grid2d;
use crate::hllc2d::{hllc_flux, FacePrim};
use crate::state2d::{cons_to_prim2d, ConservedState2d};
use crate::vectorize::apply_divergence;
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

/// muscl face states on a non-uniform axis: the limiter acts on physical
/// one-sided gradients and the reconstruction extrapolates to the true face
/// positions, with ghost centers reflected across the boundary faces.
fn face_states_axis(
    p: &[f64],
    n: usize,
    muscl: bool,
    faces: &[f64],
    centers: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let xg0 = 2.0 * faces[0] - centers[0];
    let xgn = 2.0 * faces[n] - centers[n - 1];
    let mut d = vec![0.0; n];
    if muscl {
        for i in 0..n {
            let (xl, pl) = if i == 0 {
                (xg0, p[0])
            } else {
                (centers[i - 1], p[i])
            };
            let (xr, pr) = if i == n - 1 {
                (xgn, p[n + 1])
            } else {
                (centers[i + 1], p[i + 2])
            };
            let gm = (p[i + 1] - pl) / (centers[i] - xl);
            let gp = (pr - p[i + 1]) / (xr - centers[i]);
            d[i] = van_leer(gm, gp);
        }
    }
    let mut fl = vec![0.0; n + 1];
    let mut fr = vec![0.0; n + 1];
    fl[0] = p[0];
    fr[0] = p[1] - d[0] * (centers[0] - faces[0]);
    fl[n] = p[n] + d[n - 1] * (faces[n] - centers[n - 1]);
    fr[n] = p[n + 1];
    for f in 1..n {
        fl[f] = p[f] + d[f - 1] * (faces[f] - centers[f - 1]);
        fr[f] = p[f + 1] - d[f] * (centers[f] - faces[f]);
    }
    (fl, fr)
}

/// true when every cell width on the axis matches the first to roundoff, so
/// the uniform-grid fast paths (and their exact float behavior) apply.
fn axis_uniform(ws: &[f64]) -> bool {
    let w0 = ws[0].abs().max(1e-12);
    ws.iter().all(|w| (w - ws[0]).abs() <= 1e-9 * w0)
}

/// an inflow profile prescribing u(y)/v(y) along an inflow plane instead of
/// one uniform state.
#[derive(Debug, Clone, PartialEq)]
pub enum InflowProfile {
    /// the blasius laminar layer, evaluated at the inflow plane.
    Blasius(BlasiusProfile),
    /// a tabulated profile: (y, u, v) rows, linearly interpolated.
    Table {
        /// the station y-coordinates, strictly increasing.
        ys: Vec<f64>,
        /// the streamwise velocity per row.
        us: Vec<f64>,
        /// the wall-normal velocity per row.
        vs: Vec<f64>,
    },
}

impl InflowProfile {
    /// the (u, v) prescribed at the inflow plane at height y; the plane's x
    /// station is where the profile is anchored (for blasius, the leading
    /// edge offset is folded in when the profile is built).
    pub fn at(&self, y: f64) -> (f64, f64) {
        match self {
            InflowProfile::Blasius(b) => (b.u(b.x_in, y), b.v(b.x_in, y)),
            InflowProfile::Table { ys, us, vs } => {
                if ys.is_empty() {
                    return (0.0, 0.0);
                }
                if y <= ys[0] {
                    return (us[0], vs[0]);
                }
                if y >= *ys.last().unwrap() {
                    let k = ys.len() - 1;
                    return (us[k], vs[k]);
                }
                let i = ys.partition_point(|&v| v < y);
                let (y0, y1) = (ys[i - 1], ys[i]);
                let w = if y1 > y0 { (y - y0) / (y1 - y0) } else { 0.0 };
                (
                    us[i - 1] + w * (us[i] - us[i - 1]),
                    vs[i - 1] + w * (vs[i] - vs[i - 1]),
                )
            }
        }
    }
}

/// a boundary condition applied to one side of the 2d domain.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Bc2d {
    /// open end; the ghost copies the interior state so waves leave freely.
    #[default]
    Transmissive,
    /// fixed freestream primitive state, used at a supersonic inflow face.
    SupersonicInflow { rho: f64, u: f64, v: f64, p: f64 },
    /// inflow carrying a boundary-layer profile: u(y)/v(y) from the profile,
    /// rho/p from the freestream values.
    ProfileInflow {
        profile: std::sync::Arc<InflowProfile>,
        rho: f64,
        p: f64,
    },
    /// the ghost copies the interior state, the usual supersonic outflow.
    SupersonicOutflow,
    /// solid wall; the ghost mirrors the normal velocity and copies the rest.
    SlipWall,
    /// viscous solid wall: both velocity components are reversed so the
    /// wall-average velocity vanishes (no slip), density and pressure copied.
    NoSlipWall,
}

/// the boundary condition on each of the four sides (west,east,south,north).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Boundaries2d {
    pub west: Bc2d,
    pub east: Bc2d,
    pub south: Bc2d,
    pub north: Bc2d,
}

/// the ghost value for one primitive at one side of a strip.
///
/// `var` is the variable index (0=rho,1=u,2=v,3=p) and `normal_var` the index of
/// the velocity normal to the face (1 for vertical walls, 2 for horizontal), so
/// a slip wall knows which component to reflect. a profile inflow additionally
/// needs `face_pos`, the coordinate along the inflow plane where the ghost
/// center sits (pass 0.0 when the side cannot be a profile).
fn ghost_value(bc: &Bc2d, interior: f64, var: usize, normal_var: usize, face_pos: f64) -> f64 {
    match bc {
        Bc2d::Transmissive | Bc2d::SupersonicOutflow => interior,
        Bc2d::SlipWall => {
            if var == normal_var {
                -interior
            } else {
                interior
            }
        }
        Bc2d::NoSlipWall => {
            if var == 1 || var == 2 {
                -interior
            } else {
                interior
            }
        }
        Bc2d::SupersonicInflow { rho, u, v, p } => match var {
            0 => *rho,
            2 => *v,
            3 => *p,
            _ => *u,
        },
        Bc2d::ProfileInflow { profile, rho, p } => match var {
            0 => *rho,
            3 => *p,
            // the ghost center sits one half-cell beyond the face; the
            // profile is continuous, so evaluating at the true position.
            1 => profile.at(face_pos).0,
            _ => profile.at(face_pos).1,
        },
    }
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
/// one conservative step, with the time increment additionally capped at max_dt.
///
/// cap is infinite for pure euler; the viscous wrapper passes the diffusive
/// stability bound so the explicit diffusion never outruns the time step.
pub(crate) fn advance2d_capped(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    max_dt: f64,
) -> Result<(f64, f64), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let dti = cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, nx * ny);
    let dt = dti.min(max_dt);

    // stretched-axis handling: uniform axes keep the exact scalar fast path,
    // non-uniform ones reconstruct and diverge with the true cell metrics.
    let unif_x = axis_uniform(&g.dxs);
    let unif_y = axis_uniform(&g.dys);
    let dtdxs: Vec<f64> = g.dxs.iter().map(|w| dt / w).collect();
    let dtdys: Vec<f64> = g.dys.iter().map(|w| dt / w).collect();

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
        // the row's y coordinate: a west/east ghost cell shares it, so a
        // profile inflow there evaluates at the row's own height.
        let y_row = g.centers_y[j * nx];
        let pad = |c: &[f64], var: usize| -> Vec<f64> {
            let mut out = vec![0.0; nx + 2];
            out[0] = ghost_value(&bc.west, c[0], var, 1, y_row);
            out[nx + 1] = ghost_value(&bc.east, c[nx - 1], var, 1, y_row);
            out[1..=nx].copy_from_slice(c);
            out
        };
        let x_centers: Vec<f64> = g.faces_x[..nx]
            .iter()
            .zip(g.faces_x[1..].iter())
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let row_faces_states = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            if unif_x {
                face_states(&padded, nx, muscl)
            } else {
                face_states_axis(&padded, nx, muscl, &g.faces_x, &x_centers)
            }
        };
        let (rl, rr) = row_faces_states(&pr, 0);
        let (ul, ur) = row_faces_states(&pu, 1);
        let (vl, vr) = row_faces_states(&pv, 2);
        let (pl, prr) = row_faces_states(&pp, 3);
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
        // the column's x coordinate and the reflected ghost-center x the
        // profile inflow would evaluate at beyond a horizontal face.
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
        let y_centers: Vec<f64> = g.faces_y[..ny]
            .iter()
            .zip(g.faces_y[1..].iter())
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let col_faces_states = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            if unif_y {
                face_states(&padded, ny, muscl)
            } else {
                face_states_axis(&padded, ny, muscl, &g.faces_y, &y_centers)
            }
        };
        let (rl, rr) = col_faces_states(&cr, 0);
        let (ul, ur) = col_faces_states(&cu, 1);
        let (vl, vr) = col_faces_states(&cv, 2);
        let (pl, prr) = col_faces_states(&cp, 3);
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
    for (j, fa) in fx.iter().enumerate() {
        // per-cell y factors; the uniform value when the y axis is uniform.
        let dyf = if unif_y { dtdy } else { dtdys[j] };
        if unif_x {
            // transverse (y) contribution per cell in this row, all four fields.
            let mut yr = vec![0.0; nx];
            let mut ym = vec![0.0; nx];
            let mut yn = vec![0.0; nx];
            let mut ye = vec![0.0; nx];
            for i in 0..nx {
                let fycol = &fy[i];
                yr[i] = dyf * (fycol.mass[j + 1] - fycol.mass[j]);
                ym[i] = dyf * (fycol.mx[j + 1] - fycol.mx[j]);
                yn[i] = dyf * (fycol.my[j + 1] - fycol.my[j]);
                ye[i] = dyf * (fycol.e[j + 1] - fycol.e[j]);
                resid = resid
                    .max((dtdx * (fa.mass[i + 1] - fa.mass[i]) + yr[i]).abs())
                    .max((dtdx * (fa.mx[i + 1] - fa.mx[i]) + ym[i]).abs())
                    .max((dtdx * (fa.my[i + 1] - fa.my[i]) + yn[i]).abs())
                    .max((dtdx * (fa.e[i + 1] - fa.e[i]) + ye[i]).abs());
            }
            let off = j * nx;
            apply_divergence(&mut state.rho[off..off + nx], &fa.mass, &yr, dtdx);
            apply_divergence(&mut state.mx[off..off + nx], &fa.mx, &ym, dtdx);
            apply_divergence(&mut state.my[off..off + nx], &fa.my, &yn, dtdx);
            apply_divergence(&mut state.e[off..off + nx], &fa.e, &ye, dtdx);
        } else {
            // stretched x: per-cell width ratios replace the constant dtdx.
            for i in 0..nx {
                let fycol = &fy[i];
                let yr = dyf * (fycol.mass[j + 1] - fycol.mass[j]);
                let ym = dyf * (fycol.mx[j + 1] - fycol.mx[j]);
                let yn = dyf * (fycol.my[j + 1] - fycol.my[j]);
                let ye = dyf * (fycol.e[j + 1] - fycol.e[j]);
                resid = resid
                    .max((dtdxs[i] * (fa.mass[i + 1] - fa.mass[i]) + yr).abs())
                    .max((dtdxs[i] * (fa.mx[i + 1] - fa.mx[i]) + ym).abs())
                    .max((dtdxs[i] * (fa.my[i + 1] - fa.my[i]) + yn).abs())
                    .max((dtdxs[i] * (fa.e[i + 1] - fa.e[i]) + ye).abs());
                let k = j * nx + i;
                state.rho[k] -= dtdxs[i] * (fa.mass[i + 1] - fa.mass[i]) + yr;
                state.mx[k] -= dtdxs[i] * (fa.mx[i + 1] - fa.mx[i]) + ym;
                state.my[k] -= dtdxs[i] * (fa.my[i + 1] - fa.my[i]) + yn;
                state.e[k] -= dtdxs[i] * (fa.e[i + 1] - fa.e[i]) + ye;
            }
        }
    }
    Ok((dt, resid))
}

/// the plain euler step (no time-step cap beyond the cfl bound).
pub fn advance2d(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
) -> Result<(f64, f64), Error> {
    advance2d_capped(state, g, gamma, cfl, muscl, bc, f64::INFINITY)
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

/// threaded 2d euler step: shards the x-face rows and y-face columns and the
/// update across a thread pool, then applies the deltas serially so the result
/// is bit-identical to the serial advance2d regardless of the thread count.
pub fn advance2d_par(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
) -> Result<(f64, f64), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let dt = cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, nx * ny);
    let unif_x = axis_uniform(&g.dxs);
    let unif_y = axis_uniform(&g.dys);
    // a shared reference for the read-only flux closures so they can be Sync.
    let st: &ConservedState2d = state;
    let x_centers: Vec<f64> = g.faces_x[..nx]
        .iter()
        .zip(g.faces_x[1..].iter())
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let y_centers: Vec<f64> = g.faces_y[..ny]
        .iter()
        .zip(g.faces_y[1..].iter())
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let row_flux = |j: usize| -> Sweep {
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
        let row_faces_states = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            if unif_x {
                face_states(&padded, nx, muscl)
            } else {
                face_states_axis(&padded, nx, muscl, &g.faces_x, &x_centers)
            }
        };
        let (rl, rr) = row_faces_states(&pr, 0);
        let (ul, ur) = row_faces_states(&pu, 1);
        let (vl, vr) = row_faces_states(&pv, 2);
        let (pl, prr) = row_faces_states(&pp, 3);
        let mut sw = Sweep {
            mass: vec![0.0; nx + 1],
            mx: vec![0.0; nx + 1],
            my: vec![0.0; nx + 1],
            e: vec![0.0; nx + 1],
        };
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
            sw.mass[f] = q.mass;
            sw.mx[f] = q.mx;
            sw.my[f] = q.my;
            sw.e[f] = q.e;
        }
        sw
    };
    let col_flux = |i: usize| -> Sweep {
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
        let col_faces_states = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            if unif_y {
                face_states(&padded, ny, muscl)
            } else {
                face_states_axis(&padded, ny, muscl, &g.faces_y, &y_centers)
            }
        };
        let (rl, rr) = col_faces_states(&cr, 0);
        let (ul, ur) = col_faces_states(&cu, 1);
        let (vl, vr) = col_faces_states(&cv, 2);
        let (pl, prr) = col_faces_states(&cp, 3);
        let mut sw = Sweep {
            mass: vec![0.0; ny + 1],
            mx: vec![0.0; ny + 1],
            my: vec![0.0; ny + 1],
            e: vec![0.0; ny + 1],
        };
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
            sw.mass[f] = q.mass;
            sw.mx[f] = q.mx;
            sw.my[f] = q.my;
            sw.e[f] = q.e;
        }
        sw
    };
    let build = |n: usize, flux: &(dyn Fn(usize) -> Sweep + Sync)| -> Vec<Sweep> {
        let threads = nthreads.min(n).max(1);
        let chunk = n.div_ceil(threads);
        let mut out: Vec<Option<Sweep>> = (0..n).map(|_| None).collect();
        std::thread::scope(|s| {
            let mut handles = Vec::new();
            for start in (0..n).step_by(chunk) {
                let end = (start + chunk).min(n);
                let flux = &flux;
                handles.push(s.spawn(move || (start, (start..end).map(flux).collect::<Vec<_>>())));
            }
            for h in handles {
                let (start, v) = h.join().unwrap();
                for (k, sw) in v.into_iter().enumerate() {
                    out[start + k] = Some(sw);
                }
            }
        });
        out.into_iter().map(|s| s.unwrap()).collect()
    };
    let fx = build(ny, &row_flux);
    let fy = build(nx, &col_flux);
    // per-cell deltas are independent; compute them in parallel, apply serially.
    let dtdx = dt / g.dx;
    let dtdy = dt / g.dy;
    let dtdxs: Vec<f64> = g.dxs.iter().map(|w| dt / w).collect();
    let dtdys: Vec<f64> = g.dys.iter().map(|w| dt / w).collect();
    let cells = g.nx * g.ny;
    let mut resid = 0.0f64;
    let threads = nthreads.min(cells).max(1);
    let chunk = cells.div_ceil(threads);
    let mut dr = vec![0.0; cells];
    let mut dm = vec![0.0; cells];
    let mut dn = vec![0.0; cells];
    let mut de = vec![0.0; cells];
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for start in (0..cells).step_by(chunk) {
            let end = (start + chunk).min(cells);
            let fx = &fx;
            let fy = &fy;
            let (mut drl, mut dml, mut dnl, mut del) = (
                vec![0.0; end - start],
                vec![0.0; end - start],
                vec![0.0; end - start],
                vec![0.0; end - start],
            );
            let mut res = 0.0f64;
            let (dtdx_i, dtdy_j) = if unif_x && unif_y {
                (dtdx, dtdy)
            } else {
                (0.0f64, 0.0f64)
            };
            for c in start..end {
                let i = c % nx;
                let j = c / nx;
                let k = c - start;
                let fxrow = &fx[j];
                let fycol = &fy[i];
                let (dxi, dyj) = if unif_x && unif_y {
                    (dtdx_i, dtdy_j)
                } else {
                    (
                        if unif_x { dtdx } else { dtdxs[i] },
                        if unif_y { dtdy } else { dtdys[j] },
                    )
                };
                drl[k] = dxi * (fxrow.mass[i + 1] - fxrow.mass[i])
                    + dyj * (fycol.mass[j + 1] - fycol.mass[j]);
                dml[k] =
                    dxi * (fxrow.mx[i + 1] - fxrow.mx[i]) + dyj * (fycol.mx[j + 1] - fycol.mx[j]);
                dnl[k] =
                    dxi * (fxrow.my[i + 1] - fxrow.my[i]) + dyj * (fycol.my[j + 1] - fycol.my[j]);
                del[k] = dxi * (fxrow.e[i + 1] - fxrow.e[i]) + dyj * (fycol.e[j + 1] - fycol.e[j]);
                res = res
                    .max(drl[k].abs())
                    .max(dml[k].abs())
                    .max(dnl[k].abs())
                    .max(del[k].abs());
            }
            handles.push(s.spawn(move || (start, drl, dml, dnl, del, res)));
        }
        for h in handles {
            let (start, drl, dml, dnl, del, res) = h.join().unwrap();
            let len = drl.len();
            dr[start..start + len].copy_from_slice(&drl);
            dm[start..start + len].copy_from_slice(&dml);
            dn[start..start + len].copy_from_slice(&dnl);
            de[start..start + len].copy_from_slice(&del);
            resid = resid.max(res);
        }
    });
    for c in 0..cells {
        state.rho[c] -= dr[c];
        state.mx[c] -= dm[c];
        state.my[c] -= dn[c];
        state.e[c] -= de[c];
    }
    Ok((dt, resid))
}

/// second-order two-stage (heun) time step: explicit predictor-corrector so
/// muscl's spatial accuracy is not masked by first-order time integration.
pub fn advance2d_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
) -> Result<(f64, f64), Error> {
    // heun: s1 = st + dt l(st); s2 = s1 + dt l(s1); st = 0.5(st + s2).
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let dt = cfl * g.dx.min(g.dy) / smax_of(&u, &v, &p, &state.rho, gamma, g.nx * g.ny);
    let mut s1 = state.clone();
    advance2d(&mut s1, g, gamma, cfl, muscl, bc)?;
    advance2d(&mut s1, g, gamma, cfl, muscl, bc)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt, 0.0))
}

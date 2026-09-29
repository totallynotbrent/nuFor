//! the flat-plate gallery: laminar blasius and turbulent sa plates solved
//! on the clustered grid the validations use, cached for the web ui.
//!
//! the laminar plate is cheap enough to solve on demand; the turbulent one
//! marches in a background thread at server start and reports progress so
//! the ui can show it arriving.

use nufor_core::{
    advance2d_sa_rk2, advance2d_visc_rk2, cons_to_prim2d, prim_to_cons2d, render_png_rect,
    stretched_grid2d, wall_distance2d, Bc2d, BlasiusProfile, Boundaries2d, Clustering,
    ConservedState2d, Grid2d, InflowProfile, SaParams, TurbState, ViscParams,
};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Mutex;

const GAMMA: f64 = 1.4;
const U_INF: f64 = 0.2;

/// the laminar plate: domain [0.2, 1.2] with the edge at x=0, nu = 1e-3.
pub fn laminar_grid(nx: usize, ny: usize) -> Grid2d {
    stretched_grid2d(
        nx,
        0.2,
        1.2,
        ny,
        0.0,
        1.0,
        Clustering {
            first_cell: 2.0e-3,
            growth: 1.2,
        },
    )
    .unwrap()
}

/// the turbulent plate: domain [0, 1.2], edge 0.2 upstream, mu = 4.8e-7.
pub fn turbulent_grid(nx: usize, ny: usize) -> Grid2d {
    stretched_grid2d(
        nx,
        0.0,
        1.2,
        ny,
        0.0,
        0.35,
        Clustering {
            first_cell: 5.18e-5,
            growth: 1.2,
        },
    )
    .unwrap()
}

/// the discrete wall skin friction per station: the quadratic-consistent
/// derivative the solver and the cli cf report use.
pub fn wall_cf(g: &Grid2d, st: &ConservedState2d, mu: f64) -> Vec<f64> {
    let (u, _, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let mut cf = vec![0.0; g.nx];
    for i in 0..g.nx {
        let (y0, y1) = (g.centers_y[i], g.centers_y[g.nx + i]);
        let (u0, u1) = (u[i], u[g.nx + i]);
        let w0 = (g.ymin - y1) / ((y0 - g.ymin) * (y0 - y1));
        let w1 = (g.ymin - y0) / ((y1 - g.ymin) * (y1 - y0));
        cf[i] = 2.0 * mu * (u0 * w0 + u1 * w1) / U_INF.powi(2);
    }
    cf
}

/// one solved plate case, ready for rendering and extraction.
pub struct Plate {
    pub grid: Grid2d,
    pub state: ConservedState2d,
    pub mu: f64,
    pub cf: Vec<f64>,
}

/// solve the laminar blasius plate to steady hold (several flow-throughs).
pub fn solve_laminar(nx: usize, ny: usize) -> Plate {
    let g = laminar_grid(nx, ny);
    let nu = 1.0e-3;
    let p = BlasiusProfile::new(U_INF, nu, 0.0)
        .unwrap()
        .anchored_at(0.2)
        .unwrap();
    let n = nx * ny;
    let mut u = vec![0.0; n];
    let mut v = vec![0.0; n];
    let rho = vec![1.0; n];
    let pr = vec![1.0; n];
    for k in 0..n {
        u[k] = p.u(g.centers_x[k], g.centers_y[k]);
        v[k] = p.v(g.centers_x[k], g.centers_y[k]);
    }
    let et: Vec<f64> = pr
        .iter()
        .zip(&u)
        .zip(&v)
        .map(|((pp, uu), vv)| pp / (GAMMA - 1.0) + 0.5 * (uu * uu + vv * vv))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: std::sync::Arc::new(InflowProfile::Blasius(p)),
            rho: 1.0,
            p: 1.0,
        },
        east: Bc2d::Transmissive,
        south: Bc2d::NoSlipWall,
        north: Bc2d::Transmissive,
    };
    let mut t = 0.0;
    while t < 0.4 {
        let (dt, _) = advance2d_visc_rk2(
            &mut st,
            &g,
            GAMMA,
            0.4,
            true,
            &bc,
            ViscParams { mu: nu, pr: 0.72 },
        )
        .unwrap();
        t += dt;
    }
    let cf = wall_cf(&g, &st, nu);
    Plate {
        grid: g,
        state: st,
        mu: nu,
        cf,
    }
}

/// march the turbulent sa plate to the given time, returning the case.
pub fn solve_turbulent(nx: usize, ny: usize, t_end: f64, progress: &AtomicU32) -> Plate {
    let g = turbulent_grid(nx, ny);
    let n = nx * ny;
    let mu = 4.8e-7;
    let nu = mu;
    let nu_tilde_inf = 3.0 * nu;
    let et = vec![1.0 / (GAMMA - 1.0) + 0.5 * U_INF * U_INF; n];
    let (mx, my, e) = prim_to_cons2d(&vec![1.0; n], &vec![U_INF; n], &vec![0.0; n], &et).unwrap();
    let mut st = ConservedState2d {
        rho: vec![1.0; n],
        mx,
        my,
        e,
    };
    let p = BlasiusProfile::new(U_INF, nu, -0.2)
        .unwrap()
        .anchored_at(0.0)
        .unwrap();
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: std::sync::Arc::new(InflowProfile::Blasius(p)),
            rho: 1.0,
            p: 1.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::NoSlipWall,
        north: Bc2d::SupersonicOutflow,
    };
    let mut turb = TurbState {
        nu_tilde: vec![nu_tilde_inf; n],
        d: wall_distance2d(&g, &bc, 0.35),
        params: SaParams {
            mu,
            pr: 0.72,
            pr_t: 0.9,
            nu_tilde_inf,
        },
    };
    let mut t = 0.0f64;
    while t < t_end {
        t += advance2d_sa_rk2(&mut st, &mut turb, &g, GAMMA, 0.3, true, &bc).unwrap();
        let pct = ((t / t_end) * 100.0) as u32;
        progress.store(pct.min(99), Ordering::Relaxed);
    }
    progress.store(100, Ordering::Relaxed);
    let cf = wall_cf(&g, &st, mu);
    Plate {
        grid: g,
        state: st,
        mu,
        cf,
    }
}

/// the shared turbulent-plate cache: solved once in the background, read by
/// the api handlers. a smaller grid keeps the warmup march under a minute.
/// the laminar plate solves on demand (about two seconds), so it needs no
/// cache slot.
pub struct PlateCache {
    pub turbulent: Mutex<Option<Plate>>,
    pub turb_progress: AtomicU32,
}

impl PlateCache {
    pub fn new() -> Self {
        Self {
            turbulent: Mutex::new(None),
            turb_progress: AtomicU32::new(0),
        }
    }

    /// kick off the background turbulent solve at server start.
    pub fn spawn_warmup(self: &std::sync::Arc<Self>) {
        let cache = std::sync::Arc::clone(self);
        std::thread::spawn(move || {
            let plate = solve_turbulent(96, 48, 5.0, &cache.turb_progress);
            *cache.turbulent.lock().unwrap() = Some(plate);
        });
    }
}

/// render one field of a plate to a png (rectangular: the plate is wide).
pub fn plate_image(plate: &Plate, field: &str) -> Result<Vec<u8>, String> {
    let g = &plate.grid;
    let (u, v, _) = cons_to_prim2d(
        &plate.state.rho,
        &plate.state.mx,
        &plate.state.my,
        &plate.state.e,
    )
    .map_err(|e| e.to_string())?;
    let (data, lo, hi) = match field {
        "mach" => {
            let m: Vec<f64> = (0..g.nx * g.ny)
                .map(|k| {
                    let pp = (GAMMA - 1.0)
                        * (plate.state.e[k]
                            - 0.5
                                * (plate.state.mx[k] * plate.state.mx[k]
                                    + plate.state.my[k] * plate.state.my[k])
                                / plate.state.rho[k]);
                    (u[k] * u[k] + v[k] * v[k]).sqrt()
                        / (GAMMA * pp / plate.state.rho[k]).sqrt().max(1e-12)
                })
                .collect();
            (m, 0.0, 0.25)
        }
        "p" => {
            let pr: Vec<f64> = (0..g.nx * g.ny)
                .map(|k| {
                    (GAMMA - 1.0)
                        * (plate.state.e[k]
                            - 0.5
                                * (plate.state.mx[k] * plate.state.mx[k]
                                    + plate.state.my[k] * plate.state.my[k])
                                / plate.state.rho[k])
                })
                .collect();
            (pr, 0.985, 1.015)
        }
        _ => (plate.state.rho.clone(), 0.985, 1.015),
    };
    render_png_rect(&data, g.nx, g.ny, lo, hi)
}

/// the json the validation dock consumes: our cf(x) plus the published
/// references evaluated at the same stations.
pub fn validation_json(plate: &Plate, turbulent: bool) -> String {
    let g = &plate.grid;
    let mut rows = String::new();
    for i in 0..g.nx {
        let x = g.centers_x[i];
        let re_x = U_INF * x / plate.mu;
        let cf = plate.cf[i];
        let (blasius, corr) = if turbulent {
            (String::new(), format!("{:.5e}", 0.0592 * re_x.powf(-0.2)))
        } else {
            (
                if re_x > 0.0 {
                    format!("{:.5e}", 0.664 / re_x.sqrt())
                } else {
                    String::new()
                },
                String::new(),
            )
        };
        rows.push_str(&format!(
            "{{\"x\":{x:.5},\"re_x\":{re_x:.1},\"cf\":{cf:.5e},\"blasius\":\"{blasius}\",\"corr\":\"{corr}\"}},"
        ));
    }
    rows.pop();
    format!("{{\"turbulent\":{turbulent},\"stations\":[{rows}]}}")
}

/// the wall-normal u(y) profile at the nearest station to x/xmax = frac,
/// against the exact blasius solution when the case is laminar.
pub fn profile_json(plate: &Plate, frac: f64, laminar: bool) -> String {
    let g = &plate.grid;
    let (u, _, _) = cons_to_prim2d(
        &plate.state.rho,
        &plate.state.mx,
        &plate.state.my,
        &plate.state.e,
    )
    .unwrap();
    let x_t = g.xmin + frac * (g.xmax - g.xmin);
    let i = (0..g.nx)
        .min_by(|a, b| {
            (g.centers_x[*a] - x_t)
                .abs()
                .partial_cmp(&(g.centers_x[*b] - x_t).abs())
                .unwrap()
        })
        .unwrap();
    let x = g.centers_x[i];
    let p = if laminar {
        Some(BlasiusProfile::new(U_INF, plate.mu, 0.0).unwrap())
    } else {
        None
    };
    let mut rows = String::new();
    for j in 0..g.ny {
        let k = j * g.nx + i;
        let y = g.centers_y[k];
        let exact = p
            .as_ref()
            .map(|pp| pp.u(x + 0.2, y))
            .map(|e| format!("{e:.5e}"))
            .unwrap_or_default();
        rows.push_str(&format!(
            "{{\"y\":{y:.5},\"u\":{:.5e},\"exact\":\"{exact}\"}},",
            u[k]
        ));
    }
    rows.pop();
    format!("{{\"x\":{x:.5},\"rows\":[{rows}]}}")
}

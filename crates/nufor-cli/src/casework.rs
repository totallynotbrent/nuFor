//! the case workflow: load, edit, run real case.toml files from the
//! browser, stream frames with the body silhouette, and read the
//! surface table (Cp along the body plus the integrated C_A).
//!
//! a case run occupies the same single-job slot as the generic march:
//! one live job at a time, polled by /api/cases/run-status, frames
//! baked into the same frame cache the playback control reads. the
//! server holds no case-specific code: everything comes from the
//! loaded case file.

use nufor_config::{load_case_config, CaseConfig};
use nufor_core::{
    advance2d_axi_rk2, advance2d_model_rk2, cons_to_prim2d, grid2d, render_png_rect,
    write_vtk2d_model, ConservedState2d, Grid2d,
};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::{Arc, Mutex};

/// one case job: the loaded config, the march, and the results.
pub struct CaseJob {
    pub name: String,
    pct: AtomicU32,
    steps: AtomicU32,
    t_now: Mutex<f64>,
    t_end: f64,
    nx: usize,
    ny: usize,
    equations: String,
    done: AtomicBool,
    failed: Mutex<Option<String>>,
    /// the final conserved field + grid, for the surface table and probes.
    pub final_state: Mutex<Option<(Grid2d, ConservedState2d)>>,
    /// the case's own freestream and gamma, so the surface table
    /// normalizes with the run's actual conditions.
    pub gamma: f64,
    pub rho_inf: f64,
    pub u_inf: f64,
    pub p_inf: f64,
    /// the selected case's name (the directory name), which may differ
    /// from metadata.name inside a copied file.
    pub case_name: String,
    /// the thermodynamic closure the run marched with, so the surface
    /// table reads pressure the way the march wrote energy.
    pub model: nufor_core::ThermoModel,
    /// the reference area for the integrated coefficient: pi*rb^2 for
    /// a sphere-cone, pi*r^2 for a sphere, else 1.
    pub a_ref: f64,
}

impl CaseJob {
    fn new(cfg: &CaseConfig) -> Self {
        let (rho_inf, u_inf, p_inf) = match &cfg.initial_condition {
            nufor_config::InitialCondition::Uniform { rho, u, p } => (*rho, *u, *p),
            _ => (1.0, 0.0, 1.0),
        };
        let a_ref = match &cfg.body {
            Some(nufor_config::BodySection::SphereCone { rb, .. }) => {
                std::f64::consts::PI * rb * rb
            }
            Some(nufor_config::BodySection::Sphere { r, .. }) => std::f64::consts::PI * r * r,
            _ => 1.0,
        };
        Self {
            name: cfg.metadata.name.clone(),
            case_name: String::new(),
            model: nufor_core::model_from_config(cfg),
            a_ref,
            gamma: cfg.physics.gamma,
            rho_inf,
            u_inf,
            p_inf,
            equations: format!("{:?}", cfg.physics.equations),
            t_end: cfg.time.final_time,
            nx: cfg.mesh.nx as usize,
            ny: cfg.mesh.ny.unwrap_or(cfg.mesh.nx) as usize,
            pct: AtomicU32::new(0),
            steps: AtomicU32::new(0),
            t_now: Mutex::new(0.0),
            done: AtomicBool::new(false),
            failed: Mutex::new(None),
            final_state: Mutex::new(None),
        }
    }

    pub fn status_json(&self) -> String {
        let pct = self.pct.load(Ordering::Relaxed);
        let done = self.done.load(Ordering::Relaxed);
        let failed = self.failed.lock().unwrap().clone();
        let t = *self.t_now.lock().unwrap();
        let steps = self.steps.load(Ordering::Relaxed);
        match failed {
            Some(msg) => format!(
                "{{\"running\":false,\"done\":true,\"failed\":true,\"error\":\"{msg}\",\"name\":\"{}\",\"nx\":{}}}",
                self.name, self.nx
            ),
            None => format!(
                "{{\"running\":{running},\"done\":{done},\"failed\":false,\"pct\":{pct},\"steps\":{steps},\"t\":{t:.4},\"t_end\":{t_end},\"name\":\"{name}\",\"nx\":{nx},\"ny\":{ny},\"equations\":\"{eq}\"}}",
                running = !done,
                t_end = self.t_end,
                name = self.name,
                nx = self.nx,
                ny = self.ny,
                eq = self.equations,
            ),
        }
    }
}

/// the shared case state: the cases directory, one live job.
pub struct CaseState {
    pub dir: PathBuf,
    job: Mutex<Option<Arc<CaseJob>>>,
    pub frames: Arc<crate::runstate::FrameCache>,
}

impl CaseState {
    pub fn new(dir: PathBuf, frames: Arc<crate::runstate::FrameCache>) -> Self {
        Self {
            dir,
            job: Mutex::new(None),
            frames,
        }
    }

    /// the case files in the directory, one entry per case.toml found.
    pub fn list_json(&self) -> String {
        let mut names: Vec<String> = Vec::new();
        if let Ok(rd) = std::fs::read_dir(&self.dir) {
            for e in rd.flatten() {
                let p = e.path();
                let is_case = p.is_dir() && p.join("case.toml").is_file();
                if is_case {
                    if let Some(n) = p.file_name() {
                        names.push(n.to_string_lossy().to_string());
                    }
                }
            }
        }
        names.sort();
        let items: Vec<String> = names.iter().map(|n| format!("\"{n}\"")).collect();
        format!("{{\"cases\":[{}]}}", items.join(","))
    }

    /// one case's raw toml text, for the editor.
    pub fn load(&self, name: &str) -> Result<String, String> {
        let p = self.case_path(name)?;
        std::fs::read_to_string(&p).map_err(|e| format!("load {name}: {e}"))
    }

    /// save edited toml text back to the case file, only after it
    /// parses and validates.
    pub fn save(&self, name: &str, text: &str) -> Result<(), String> {
        let p = self.case_path(name)?;
        // validate before writing: a bad edit must never land.
        nufor_config::parse_case_toml(text).map_err(|e| format!("{e}"))?;
        std::fs::write(&p, text).map_err(|e| format!("save {name}: {e}"))
    }

    /// delete a case directory.
    pub fn delete(&self, name: &str) -> Result<(), String> {
        let p = self
            .case_path(name)?
            .parent()
            .map(|d| d.to_path_buf())
            .ok_or_else(|| format!("delete {name}: no parent"))?;
        std::fs::remove_dir_all(p).map_err(|e| format!("delete {name}: {e}"))
    }

    /// create a new case from toml text (a save-as from the editor).
    pub fn create(&self, name: &str, text: &str) -> Result<(), String> {
        let cfg = nufor_config::parse_case_toml(text).map_err(|e| format!("{e}"))?;
        if cfg.metadata.name.is_empty() {
            return Err("case name must not be empty".into());
        }
        let p = self.dir.join(safe_name(name)?);
        std::fs::create_dir_all(&p).map_err(|e| format!("create {name}: {e}"))?;
        std::fs::write(p.join("case.toml"), text).map_err(|e| format!("create {name}: {e}"))
    }

    /// start a run of the named case if none is live; returns the job.
    pub fn start(&self, name: &str) -> Result<Arc<CaseJob>, String> {
        let p = self.case_path(name)?;
        let cfg = load_case_config(&p).map_err(|e| format!("{e}"))?;
        let mut guard = self.job.lock().unwrap();
        if let Some(old) = guard.as_ref() {
            if !old.done.load(Ordering::Relaxed) {
                return Ok(Arc::clone(old));
            }
        }
        let mut job = CaseJob::new(&cfg);
        job.case_name = name.to_string();
        let job = Arc::new(job);
        *guard = Some(Arc::clone(&job));
        drop(guard);
        let j = Arc::clone(&job);
        let frames = Arc::clone(&self.frames);
        std::thread::spawn(move || {
            if let Err(e) = march_case(&j, &cfg, &frames) {
                *j.failed.lock().unwrap() = Some(e);
            }
            j.pct.store(100, Ordering::Relaxed);
            j.done.store(true, Ordering::Relaxed);
        });
        Ok(job)
    }

    /// restore a finished case from its on-disk vtk: the results
    /// panel works for cases the cli ran, not only ui-run ones, and
    /// survives a server restart. the vtk carries density, pressure,
    /// energy, and mach; (rho, p) recovers the internal energy
    /// exactly through the run's closure, and kinetic energy follows
    /// from the momentum the file's velocity implies, so the frames,
    /// surface table, and probes all read as if the run just ended.
    pub fn restore(&self, name: &str) -> Result<Arc<CaseJob>, String> {
        let p = self.case_path(name)?;
        let cfg = load_case_config(&p).map_err(|e| format!("{e}"))?;
        let vtk = p
            .parent()
            .ok_or_else(|| format!("restore {name}: no parent"))?
            .join(format!("{name}.vtk"));
        let model = nufor_core::model_from_config(&cfg);
        let field = case_vtk_read(&vtk, model)?;
        let mut job = CaseJob::new(&cfg);
        job.case_name = name.to_string();
        let g = grid2d(
            cfg.mesh.nx as usize,
            cfg.mesh.ny.unwrap_or(cfg.mesh.nx) as usize,
            cfg.mesh.x0,
            cfg.mesh.x1,
            cfg.mesh.y0.unwrap_or(cfg.mesh.x0),
            cfg.mesh.y1.unwrap_or(cfg.mesh.x1),
        )
        .map_err(|e| e.to_string())?;
        let body = crate::case_body_fns(&cfg.body);
        let frames = Arc::clone(&self.frames);
        frames.frames_clear();
        for f in ["rho", "mach", "p"] {
            let png = render_with_body(&g, &field, f, &body, job.model);
            frames.frames_push((f, job.t_end, png));
        }
        *job.final_state.lock().unwrap() = Some((g, field));
        job.pct.store(100, Ordering::Relaxed);
        job.done.store(true, Ordering::Relaxed);
        let job = Arc::new(job);
        *self.job.lock().unwrap() = Some(Arc::clone(&job));
        Ok(job)
    }

    pub fn status_json(&self) -> String {
        let guard = self.job.lock().unwrap();
        match guard.as_ref() {
            Some(j) => j.status_json(),
            None => "{\"running\":false,\"done\":false,\"failed\":false}".to_string(),
        }
    }

    /// the finished run's surface table: Cp along the body and the
    /// integrated C_A, plus the silhouette for the overlays.
    pub fn surface_json(&self) -> String {
        let guard = self.job.lock().unwrap();
        let Some(job) = guard.as_ref() else {
            return "{\"error\":\"no run yet\"}".to_string();
        };
        let fs = job.final_state.lock().unwrap().clone();
        let Some((g, st)) = fs else {
            return "{\"error\":\"run not finished\"}".to_string();
        };
        let gamma = job.gamma;
        let (u, v, et) = match cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e) {
            Ok(x) => x,
            Err(e) => return format!("{{\"error\":\"{e}\"}}"),
        };
        let e_int: Vec<f64> = (0..st.rho.len())
            .map(|k| et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]))
            .collect();
        let p = match nufor_core::thermo_pressure(job.model, &st.rho, &e_int, &u, &v) {
            Ok(x) => x,
            Err(e) => return format!("{{\"error\":\"{e}\"}}"),
        };
        let (rho_inf, u_inf, p_inf) = (job.rho_inf, job.u_inf, job.p_inf);
        let q_inf = 0.5 * rho_inf * u_inf * u_inf;
        // the wall-band read: fluid cells adjacent to solid ones. the
        // silhouette comes from the mask: the topmost solid cell per
        // column.
        let mut rows: Vec<String> = Vec::new();
        let mut sil: Vec<String> = Vec::new();
        let mut fx = 0.0f64;
        let a_ref = job.a_ref;
        let mut prev: Option<(f64, f64)> = None;
        for i in 0..g.nx {
            // the top solid cell in this column, scanning from the axis
            // up, gives the silhouette where the body is.
            let mut top: Option<usize> = None;
            let e_ref = 1.0 / (gamma - 1.0);
            for j in 0..g.ny {
                let k = j * g.nx + i;
                if st.rho[k] == 1.0 && (st.e[k] - e_ref).abs() < 1e-9 {
                    top = Some(j);
                }
            }
            let x = g.centers_x[i];
            if let Some(tj) = top {
                let r = g.centers_y[tj * g.nx + i];
                sil.push(format!("[{x:.4},{r:.4}]"));
                // the wall cell: the first fluid cell above the solid.
                if tj + 1 < g.ny {
                    let k = (tj + 1) * g.nx + i;
                    let cp = (p[k] - p_inf) / q_inf;
                    rows.push(format!("[{x:.4},{cp:.4}]"));
                    if let Some((r0, cp0)) = prev {
                        let dr = r - r0;
                        let rm = 0.5 * (r + r0);
                        let pm = 0.5 * (cp + cp0) * q_inf + p_inf;
                        fx += pm * 2.0 * std::f64::consts::PI * rm * dr;
                    }
                    prev = Some((r, cp));
                }
            }
        }
        let c_a = fx / (q_inf * a_ref);
        format!(
            "{{\"name\":\"{name}\",\"cp\":[{cp_rows}],\"silhouette\":[{sil}],\"c_a\":{c_a:.4}}}",
            name = job.case_name,
            cp_rows = rows.join(","),
            sil = sil.join(",")
        )
    }

    fn case_path(&self, name: &str) -> Result<PathBuf, String> {
        let n = safe_name(name)?;
        let p = self.dir.join(&n).join("case.toml");
        if !p.is_file() {
            return Err(format!("no case named {n}"));
        }
        Ok(p)
    }
}

/// read a finished case's vtk snapshot back into a conserved state.
/// the file carries density, pressure (closure-consistent), energy
/// (total per volume), and mach. velocity direction can't be
/// recovered from a scalar field, so restored momentum carries the
/// magnitude along the freestream axis: exact for the capsule
/// family's axisymmetric flow upstream and over the flank, the one
/// approximation being the subsonic layer's radial component. the
/// panels (frames, surface cp, c_a, mach) all read as-run.
pub fn case_vtk_read(
    path: &std::path::Path,
    model: nufor_core::ThermoModel,
) -> Result<ConservedState2d, String> {
    let f = std::fs::read_to_string(path).map_err(|e| format!("read vtk: {e}"))?;
    let (nx, ny) = f
        .lines()
        .find_map(|l| {
            let l = l.trim();
            if let Some(rest) = l.strip_prefix("DIMENSIONS") {
                let mut it = rest.split_whitespace();
                let nx: usize = it.next()?.parse().ok()?;
                let ny: usize = it.next()?.parse().ok()?;
                return Some((nx, ny));
            }
            None
        })
        .ok_or_else(|| "read vtk: no DIMENSIONS".to_string())?;
    let field = |name: &str| -> Option<Vec<f64>> {
        let i = f.find(&format!("SCALARS {name}"))?;
        let lt = f[i..].find("LOOKUP_TABLE")? + i;
        let j = f[lt..].find('\n')? + lt + 1;
        Some(
            f[j..]
                .split_whitespace()
                .take(nx * ny)
                .filter_map(|x| x.parse::<f64>().ok())
                .collect(),
        )
    };
    let n = nx * ny;
    let rho = field("density").ok_or("read vtk: no density")?;
    if rho.len() != n {
        return Err("read vtk: short density block".into());
    }
    let e = field("energy").ok_or("read vtk: no energy")?;
    let mach = field("mach").unwrap_or_else(|| vec![0.0; n]);
    // recover the speed per cell: |v| = mach * a, with a from the
    // closure at the running estimate of the internal energy.
    let mut speed = vec![0.0f64; n];
    for k in 0..n {
        let et = e[k] / rho[k].max(1e-12);
        let mut ke = 0.0;
        for _ in 0..3 {
            let e_int = (et - ke).max(1.0e4);
            let a = match nufor_core::thermo_sound(model, &[rho[k]], &[e_int], &[0.0], &[0.0]) {
                Ok(v) => v[0],
                Err(_) => 300.0,
            };
            let s_k = mach[k] * a;
            ke = 0.5 * s_k * s_k;
        }
        // the iteration can overshoot in near-vacuum wake cells where
        // the closure's sound speed estimate is unstable; never let
        // the recovered kinetic energy exceed half the total.
        ke = ke.min(et * 0.5);
        speed[k] = (2.0 * ke).sqrt();
    }
    let mx: Vec<f64> = (0..n).map(|k| rho[k] * speed[k]).collect();
    Ok(ConservedState2d {
        rho,
        mx,
        my: vec![0.0; n],
        e,
    })
}

/// names are single path segments, letters/digits/dash/underscore/dot
/// only: no traversal.
fn safe_name(name: &str) -> Result<String, String> {
    let ok = !name.is_empty()
        && name.len() <= 64
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.')
        && !name.starts_with('.');
    if !ok {
        return Err(format!("invalid case name {name:?}"));
    }
    Ok(name.to_string())
}

impl CaseState {
    /// a line probe through the finished run's final field: physical
    /// coordinates, bilinear samples along the segment.
    pub fn probe_json(
        &self,
        field: &str,
        x0: f64,
        y0: f64,
        x1: f64,
        y1: f64,
        samples: usize,
    ) -> String {
        let guard = self.job.lock().unwrap();
        let Some(job) = guard.as_ref() else {
            return "{\"error\":\"no run yet\"}".to_string();
        };
        let fs = job.final_state.lock().unwrap().clone();
        let Some((g, st)) = fs else {
            return "{\"error\":\"run not finished\"}".to_string();
        };
        let n = g.nx * g.ny;
        let u: Vec<f64> = (0..n).map(|k| st.mx[k] / st.rho[k]).collect();
        let v: Vec<f64> = (0..n).map(|k| st.my[k] / st.rho[k]).collect();
        let e_int: Vec<f64> = (0..n)
            .map(|k| st.e[k] / st.rho[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]))
            .collect();
        let p = nufor_core::thermo_pressure(job.model, &st.rho, &e_int, &u, &v)
            .unwrap_or_else(|_| vec![1e-12; n])
            .iter()
            .map(|x| x.max(1e-12))
            .collect::<Vec<f64>>();
        let a = nufor_core::thermo_sound(job.model, &st.rho, &e_int, &u, &v)
            .unwrap_or_else(|_| vec![1.0; n])
            .iter()
            .map(|x| x.max(1e-6))
            .collect::<Vec<f64>>();
        let vals: Vec<f64> = match field {
            "mach" => (0..n)
                .map(|k| (u[k] * u[k] + v[k] * v[k]).sqrt() / a[k])
                .collect(),
            "p" => p,
            _ => st.rho.clone(),
        };
        let pts =
            nufor_core::probe_line(&vals, &g, x0, y0, x1, y1, samples.max(2)).unwrap_or_default();
        let rows: Vec<String> = pts
            .iter()
            .map(|(s, v)| format!("{{\"s\":{:.4},\"v\":{:.5}}}", s, v))
            .collect();
        format!(
            "{{\"field\":\"{}\",\"samples\":[{}]}}",
            field,
            rows.join(",")
        )
    }

    /// the mesh pane data for a loaded case: face coordinates, extents,
    /// and the body silhouette sampled from the case's own sdf. the pane
    /// draws exactly what the march will march on.
    pub fn mesh_json(&self, name: &str) -> Result<String, String> {
        let text = self.load(name)?;
        let cfg = nufor_config::parse_case_toml(&text).map_err(|e| format!("{e}"))?;
        let ny = cfg.mesh.ny.unwrap_or(cfg.mesh.nx);
        let g = grid2d(
            cfg.mesh.nx as usize,
            ny as usize,
            cfg.mesh.x0,
            cfg.mesh.x1,
            cfg.mesh.y0.unwrap_or(cfg.mesh.x0),
            cfg.mesh.y1.unwrap_or(cfg.mesh.x1),
        )
        .map_err(|e| e.to_string())?;

        let xs: Vec<String> = g.faces_x.iter().map(|v| format!("{v:.4}")).collect();
        let ys: Vec<String> = g.faces_y.iter().map(|v| format!("{v:.4}")).collect();
        let body = match crate::case_body_fns(&cfg.body) {
            None => String::new(),
            Some((dist, _)) => {
                // walk the zero contour of the sdf column by column: one
                // up-crossing per column for these convex analytic bodies.
                let mut pts: Vec<(f64, f64)> = Vec::new();
                for i in 0..g.nx {
                    let x = g.centers_x[i];
                    let mut prev = dist(x, g.ymin);
                    let mut prev_y = g.ymin;
                    for j in 0..g.ny {
                        let y = g.centers_y[j * g.nx + i];
                        let d = dist(x, y);
                        if prev < 0.0 && d >= 0.0 {
                            let t = prev / (prev - d);
                            pts.push((x, prev_y + t * (y - prev_y)));
                        }
                        prev = d;
                        prev_y = y;
                    }
                }
                let inner: Vec<String> = pts
                    .iter()
                    .map(|(x, y)| format!("[{x:.4},{y:.4}]"))
                    .collect();
                format!(",\"body\":[{}]", inner.join(","))
            }
        };
        Ok(format!(
        "{{\"x\":[{}],\"y\":[{}],\"nx\":{},\"ny\":{},\"xmin\":{:.4},\"xmax\":{:.4},\"ymin\":{:.4},\"ymax\":{:.4}{}",
        xs.join(","),
        ys.join(","),
        g.nx,
        g.ny,
        g.xmin,
        g.xmax,
        g.ymin,
        g.ymax,
        body
    )
        + "}")
    }
}

/// the march for a loaded case: 2d planar or axisymmetric euler, with
/// the body mask applied each step and frames baked with the
/// silhouette.
fn march_case(
    job: &CaseJob,
    cfg: &CaseConfig,
    frames: &crate::runstate::FrameCache,
) -> Result<(), String> {
    if matches!(cfg.physics.equations, nufor_config::Equations::Euler1d) {
        return march_case_1d(job, cfg, frames);
    }
    let model = nufor_core::model_from_config(cfg);
    let gamma = cfg.physics.gamma;
    let g = grid2d(
        cfg.mesh.nx as usize,
        cfg.mesh.ny.unwrap_or(cfg.mesh.nx) as usize,
        cfg.mesh.x0,
        cfg.mesh.x1,
        cfg.mesh.y0.unwrap_or(cfg.mesh.x0),
        cfg.mesh.y1.unwrap_or(cfg.mesh.x1),
    )
    .map_err(|e| e.to_string())?;
    let n = g.nx * g.ny;
    let (rho, u, p) = match &cfg.initial_condition {
        nufor_config::InitialCondition::Uniform { rho, u, p } => (*rho, *u, *p),
        _ => return Err("case runs need a uniform initial condition".into()),
    };
    let e_int = match model {
        nufor_core::ThermoModel::Perfect { gamma: g_ } => p / ((g_ - 1.0) * rho),
        nufor_core::ThermoModel::EqAir => {
            nufor_core::eqair_energy(&[rho], &[p]).map_err(|e| e.to_string())?[0]
        }
    };
    let e_cell = e_int + 0.5 * u * u;
    let mut st = ConservedState2d {
        rho: vec![rho; n],
        mx: vec![rho * u; n],
        my: vec![0.0; n],
        e: vec![rho * e_cell; n],
    };
    // the body mask closures from the case section
    let body = crate::case_body_fns(&cfg.body);
    if let Some((dist, normal)) = &body {
        nufor_core::apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
    }
    let bc = crate::case_boundaries(cfg);
    let cfl = cfg.numerics.cfl;
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let axi = matches!(cfg.physics.equations, nufor_config::Equations::EulerAxi);
    frames.frames_clear();
    let frame_dt = t_end / 40.0;
    let mut next_frame = frame_dt;
    let _ = axi;
    while t < t_end && steps < max_steps {
        let (dt, _) = if axi {
            advance2d_axi_rk2(&mut st, &g, model, cfl, true, &bc, 1)
        } else {
            advance2d_model_rk2(&mut st, &g, model, cfl, true, &bc, 1)
        }
        .map_err(|e| e.to_string())?;
        if let Some((dist, normal)) = &body {
            nufor_core::apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
        }
        t += dt;
        steps += 1;
        *job.t_now.lock().unwrap() = t.min(t_end);
        job.steps.store(steps as u32, Ordering::Relaxed);
        let pct = ((t / t_end) * 100.0) as u32;
        job.pct.store(pct.min(99), Ordering::Relaxed);
        if t >= next_frame {
            next_frame += frame_dt;
            for (f, _) in [("rho", 0), ("mach", 1), ("p", 2)] {
                let png = render_with_body(&g, &st, f, &body, model);
                frames.frames_push((f, t, png));
            }
        }
    }
    for (f, _) in [("rho", 0), ("mach", 1), ("p", 2)] {
        let png = render_with_body(&g, &st, f, &body, model);
        frames.frames_push((f, t_end, png));
    }
    // persist the finished field next to the case so restore() can
    // bring it back after a restart
    let vtk_path = std::path::Path::new("cases")
        .join(&job.case_name)
        .join(format!("{}.vtk", job.case_name));
    let _ = write_vtk2d_model(&vtk_path, &g, &st, model);
    *job.final_state.lock().unwrap() = Some((g, st));
    Ok(())
}

/// the density render with the body silhouette drawn over it: solid
/// cells paint as the near-black mask color.
/// the 1d case march: the shared solver on the case's grid, frames
/// rendered as line plots (the standard 1d results view).
fn march_case_1d(
    job: &CaseJob,
    cfg: &CaseConfig,
    frames: &crate::runstate::FrameCache,
) -> Result<(), String> {
    let gamma = cfg.physics.gamma;
    let g = crate::build_1d_mesh(&cfg.mesh)?;
    let n = g.centers.len();
    let mut st = crate::ic_from_case(&cfg.initial_condition, n, gamma)?;
    let (left, right) = (
        crate::bc1d(cfg.boundaries.left),
        crate::bc1d(cfg.boundaries.right),
    );
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let frame_dt = t_end / 40.0;
    let mut next_frame = frame_dt;
    while t < t_end && steps < max_steps {
        let (dt, _) = nufor_core::advance(&mut st, gamma, cfl, g.dx, left, right)
            .map_err(|e| e.to_string())?;
        t += dt;
        steps += 1;
        *job.t_now.lock().unwrap() = t.min(t_end);
        job.steps.store(steps as u32, Ordering::Relaxed);
        let pct = ((t / t_end) * 100.0) as u32;
        job.pct.store(pct.min(99), Ordering::Relaxed);
        if t >= next_frame {
            next_frame += frame_dt;
            for (f, _) in [("rho", 0), ("u", 1), ("p", 2)] {
                let png = nufor_core::render_line1d(&g.centers, &st, gamma, f);
                frames.frames_push((f, t, png));
            }
        }
    }
    for (f, _) in [("rho", 0), ("u", 1), ("p", 2)] {
        let png = nufor_core::render_line1d(&g.centers, &st, gamma, f);
        frames.frames_push((f, t_end, png));
    }
    // the final state rides the 2d final_state slot as a 1-row field so
    // the probe and surface endpoints keep their shape.
    let st2 = ConservedState2d {
        rho: st.rho.clone(),
        mx: st.m.clone(),
        my: vec![0.0; n],
        e: st.e.clone(),
    };
    let x0 = g.faces[0];
    let x1 = g.faces[g.faces.len() - 1];
    let g2 = grid2d(n, 1, x0, x1, 0.0, 1.0).map_err(|e| e.to_string())?;
    *job.final_state.lock().unwrap() = Some((g2, st2));
    Ok(())
}

fn render_with_body(
    g: &Grid2d,
    st: &ConservedState2d,
    field: &str,
    body: &Option<crate::CaseBodyFns>,
    model: nufor_core::ThermoModel,
) -> Vec<u8> {
    let n = g.nx * g.ny;
    let u: Vec<f64> = (0..n).map(|k| st.mx[k] / st.rho[k]).collect();
    let v: Vec<f64> = (0..n).map(|k| st.my[k] / st.rho[k]).collect();
    let e_int: Vec<f64> = (0..n)
        .map(|k| st.e[k] / st.rho[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]))
        .collect();
    let p = nufor_core::thermo_pressure(model, &st.rho, &e_int, &u, &v)
        .unwrap_or_else(|_| vec![1e-12; n]);
    let p: Vec<f64> = p.iter().map(|x| x.max(1e-12)).collect();
    let a =
        nufor_core::thermo_sound(model, &st.rho, &e_int, &u, &v).unwrap_or_else(|_| vec![1.0; n]);
    let a: Vec<f64> = a.iter().map(|x| x.max(1e-6)).collect();
    let mut data: Vec<f64> = match field {
        "mach" => (0..n)
            .map(|k| (u[k] * u[k] + v[k] * v[k]).sqrt() / a[k])
            .collect(),
        "p" => p.clone(),
        _ => st.rho.clone(),
    };
    // the body silhouette: solid cells carry the quiescent reference
    // (rho=1, u=v=0), so paint them dark for the overlay read.
    if let Some((dist, _)) = body {
        let (lo, hi) = if field == "mach" {
            (0.0, 2.0)
        } else if field == "p" {
            (0.0, 6.0)
        } else {
            (0.85, 3.2)
        };
        // mark solid cells with the mask value below the scale
        for j in 0..g.ny {
            for i in 0..g.nx {
                let (x, y) = (g.centers_x[j * g.nx + i], g.centers_y[j * g.nx + i]);
                if dist(x, y) < 0.0 {
                    data[j * g.nx + i] = lo - 1.0;
                }
            }
        }
        return render_png_rect(&data, g.nx, g.ny, lo, hi).unwrap_or_default();
    }
    let (lo, hi) = if field == "mach" {
        (0.0, 2.0)
    } else if field == "p" {
        (0.0, 6.0)
    } else {
        (0.85, 3.2)
    };
    render_png_rect(&data, g.nx, g.ny, lo, hi).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state() -> CaseState {
        let dir = std::env::temp_dir().join(format!("nf-cases-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        CaseState::new(dir, Arc::new(crate::runstate::FrameCache::new()))
    }

    #[test]
    fn names_reject_traversal() {
        assert!(safe_name("../etc").is_err());
        assert!(safe_name("a/b").is_err());
        assert!(safe_name(".hidden").is_err());
        assert!(safe_name("capsule-peak-q").is_ok());
    }

    #[test]
    fn create_load_save_delete_roundtrip() {
        let s = state();
        let find_case = || {
            for prefix in ["", "../", "../../"] {
                let p = format!("{prefix}cases/capsule-peak-q/case.toml");
                if let Ok(t) = std::fs::read_to_string(&p) {
                    return t;
                }
            }
            panic!("capsule case not found from the test cwd");
        };
        let text = find_case();
        s.create("rt-test", &text).unwrap();
        let loaded = s.load("rt-test").unwrap();
        assert!(loaded.contains("capsule-peak-q"));
        let edited = loaded.replace("final_time = 0.5", "final_time = 0.05");
        s.save("rt-test", &edited).unwrap();
        assert!(s.load("rt-test").unwrap().contains("final_time = 0.05"));
        // a bad edit must be rejected, file unchanged
        let bad = edited.replace("gamma = 1.4", "gamma = 0.5");
        assert!(s.save("rt-test", &bad).is_err());
        assert!(s.load("rt-test").unwrap().contains("gamma = 1.4"));
        let listing = s.list_json();
        assert!(listing.contains("rt-test"));
        s.delete("rt-test").unwrap();
        assert!(s.load("rt-test").is_err());
    }
}

//! the run orchestration for the web ui: a background march for the blast
//! cases with live progress, plus the frame cache that makes playback
//! actually animate instead of re-solving every tick.

use nufor_core::{
    advance2d_rk2, advance3d_rk2, grid2d, grid3d, render_png_rect, write_vtk2d, Boundaries2d,
    Bounds3d, ConservedState2d, Grid2d,
};
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::Mutex;

const GAMMA: f64 = 1.4;

/// one background blast run: marches in a thread while the ui polls status.
pub struct RunJob {
    pub dim: String,
    pub n: usize,
    pub t_end: f64,
    pub cfl: f64,
    pub name: String,
    pct: AtomicU32,
    steps: AtomicU32,
    t_now: Mutex<f64>,
    done: AtomicBool,
    failed: Mutex<Option<String>>,
}

impl RunJob {
    fn new(dim: &str, n: usize, t_end: f64, cfl: f64, name: &str) -> Self {
        Self {
            dim: dim.to_string(),
            n,
            t_end,
            cfl,
            name: name.to_string(),
            pct: AtomicU32::new(0),
            steps: AtomicU32::new(0),
            t_now: Mutex::new(0.0),
            failed: Mutex::new(None),
            done: AtomicBool::new(false),
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
                "{{\"running\":false,\"done\":true,\"failed\":true,\"error\":\"{msg}\",\"dim\":\"{}\",\"n\":{}}}",
                self.dim, self.n
            ),
            None => format!(
                "{{\"running\":{running},\"done\":{done},\"failed\":false,\"pct\":{pct},\"steps\":{steps},\"t\":{t:.4},\"t_end\":{t_end},\"dim\":\"{dim}\",\"n\":{n}}}",
                running = !done,
                t_end = self.t_end,
                dim = self.dim,
                n = self.n
            ),
        }
    }

    /// the march itself; runs on the job's own thread.
    fn march2d(&self) -> Result<(), String> {
        let g = grid2d(self.n, self.n, 0.0, 1.0, 0.0, 1.0).map_err(|e| e.to_string())?;
        let mut st = crate::serve::blast_ic2d(&g, GAMMA);
        let bc = Boundaries2d::default();
        let mut t = 0.0f64;
        let mut steps = 0usize;
        while t < self.t_end {
            let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, self.cfl, true, &bc)
                .map_err(|e| e.to_string())?;
            t += dt;
            steps += 1;
            *self.t_now.lock().unwrap() = t.min(self.t_end);
            self.steps.store(steps as u32, Ordering::Relaxed);
            let pct = ((t / self.t_end) * 100.0) as u32;
            self.pct.store(pct.min(99), Ordering::Relaxed);
        }
        let _ = std::fs::create_dir_all("results");
        let dst = format!("results/{}.vtk", self.name);
        let _ = write_vtk2d(std::path::Path::new(&dst), &g, &st, GAMMA);
        Ok(())
    }

    fn march3d(&self) -> Result<(), String> {
        let b = Bounds3d {
            xmin: 0.0,
            xmax: 1.0,
            ymin: 0.0,
            ymax: 1.0,
            zmin: 0.0,
            zmax: 1.0,
        };
        let g = grid3d(self.n, self.n, self.n, &b).map_err(|e| e.to_string())?;
        let mut st = crate::serve::blast_ic3d(&g, GAMMA);
        let mut t = 0.0f64;
        let mut steps = 0usize;
        while t < self.t_end {
            let (dt, _) =
                advance3d_rk2(&mut st, &g, GAMMA, self.cfl, true).map_err(|e| e.to_string())?;
            t += dt;
            steps += 1;
            *self.t_now.lock().unwrap() = t.min(self.t_end);
            self.steps.store(steps as u32, Ordering::Relaxed);
            let pct = ((t / self.t_end) * 100.0) as u32;
            self.pct.store(pct.min(99), Ordering::Relaxed);
        }
        Ok(())
    }
}

/// the frame cache: bake n_steps 2d blast frames once in the background so
/// playback animates from memory instead of solving per request.
pub struct FrameCache {
    frames: Mutex<Vec<(f64, Vec<u8>)>>,
    pub baking: AtomicBool,
}

impl FrameCache {
    pub fn new() -> Self {
        Self {
            frames: Mutex::new(Vec::new()),
            baking: AtomicBool::new(false),
        }
    }

    /// bake the default blast animation (2d, t from 0.02 to 0.35) once.
    pub fn spawn_bake(self: &std::sync::Arc<Self>, n: usize, field: &str) {
        if self.baking.load(Ordering::Relaxed) {
            return;
        }
        self.baking.store(true, Ordering::Relaxed);
        let cache = std::sync::Arc::clone(self);
        let field = field.to_string();
        std::thread::spawn(move || {
            let g = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
            let mut st = crate::serve::blast_ic2d(&g, GAMMA);
            let bc = Boundaries2d::default();
            let mut t = 0.02f64;
            let mut frames: Vec<(f64, Vec<u8>)> = Vec::new();
            let t_stop = 0.35f64;
            while t < t_stop {
                // march a fixed dt so the frames land on round t values
                let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.4, true, &bc).unwrap();
                if t + dt > t_stop {
                    break;
                }
                t += dt;
                if frames.is_empty() || (t - frames.last().unwrap().0).abs() > 0.009 {
                    let png = render_field2d(&g, &st, &field);
                    frames.push((t, png));
                }
            }
            *cache.frames.lock().unwrap() = frames;
            cache.baking.store(false, Ordering::Relaxed);
        });
    }

    /// the frame nearest the requested t, if baked.
    pub fn frame_near(&self, t: f64) -> Option<(f64, Vec<u8>)> {
        let frames = self.frames.lock().unwrap();
        if frames.is_empty() {
            return None;
        }
        frames
            .iter()
            .min_by(|a, b| {
                (a.0 - t)
                    .abs()
                    .partial_cmp(&(b.0 - t).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .cloned()
    }

    /// the sorted t values, for the ui slider ticks.
    pub fn times_json(&self) -> String {
        let frames = self.frames.lock().unwrap();
        let ts: Vec<String> = frames.iter().map(|(t, _)| format!("{t:.3}")).collect();
        format!("[{}]", ts.join(","))
    }
}

/// render one 2d field (mach or density) from a conserved state.
pub fn render_field2d(g: &Grid2d, st: &ConservedState2d, field: &str) -> Vec<u8> {
    let (lo, hi) = if field == "mach" {
        (0.0, 1.4)
    } else {
        (0.85, 3.0)
    };
    let data: Vec<f64> = match field {
        "mach" => (0..g.nx * g.ny)
            .map(|k| {
                let rho = st.rho[k];
                let u = st.mx[k] / rho;
                let v = st.my[k] / rho;
                let p = (GAMMA - 1.0) * (st.e[k] - 0.5 * rho * (u * u + v * v));
                (u * u + v * v).sqrt() / (GAMMA * p / rho).sqrt().max(1e-12)
            })
            .collect(),
        _ => st.rho.clone(),
    };
    render_png_rect(&data, g.nx, g.ny, lo, hi).unwrap_or_default()
}

/// the mesh description the mesh pane draws: face coordinates for the active
/// case grid (uniform blast, clustered plate, or imported rectilinear).
pub fn mesh_faces_json_for(kind: &str) -> String {
    match kind {
        "laminar" => {
            let g = crate::plates::laminar_grid(80, 40);
            faces_of(&g)
        }
        "turbulent" => {
            let g = crate::plates::turbulent_grid(80, 40);
            faces_of(&g)
        }
        _ => {
            let g = grid2d(32, 32, 0.0, 1.0, 0.0, 1.0).unwrap();
            faces_of(&g)
        }
    }
}

fn faces_of(g: &Grid2d) -> String {
    let xs: Vec<String> = g.faces_x.iter().map(|v| format!("{v:.6}")).collect();
    let ys: Vec<String> = g.faces_y.iter().map(|v| format!("{v:.6}")).collect();
    format!(
        "{{\"x\":[{}],\"y\":[{}],\"nx\":{},\"ny\":{},\"xmin\":{:.4},\"xmax\":{:.4},\"ymin\":{:.4},\"ymax\":{:.4}}}",
        xs.join(","),
        ys.join(","),
        g.nx,
        g.ny,
        g.xmin,
        g.xmax,
        g.ymin,
        g.ymax
    )
}

/// the shared run state: one live job at a time plus the frame cache.
pub struct RunState {
    pub job: Mutex<Option<std::sync::Arc<RunJob>>>,
    pub frames: std::sync::Arc<FrameCache>,
    pub run_count: AtomicU32,
}

impl RunState {
    pub fn new() -> Self {
        Self {
            job: Mutex::new(None),
            frames: std::sync::Arc::new(FrameCache::new()),
            run_count: AtomicU32::new(0),
        }
    }

    /// start a run if none is live; returns the job for status polling.
    pub fn start(
        &self,
        dim: &str,
        n: usize,
        t_end: f64,
        cfl: f64,
        name: &str,
    ) -> std::sync::Arc<RunJob> {
        let job = std::sync::Arc::new(RunJob::new(dim, n, t_end, cfl, name));
        {
            let mut guard = self.job.lock().unwrap();
            if let Some(old) = guard.as_ref() {
                if !old.done.load(Ordering::Relaxed) {
                    return std::sync::Arc::clone(old);
                }
            }
            *guard = Some(std::sync::Arc::clone(&job));
        }
        self.run_count.fetch_add(1, Ordering::Relaxed);
        let j = std::sync::Arc::clone(&job);
        std::thread::spawn(move || {
            let r = if j.dim == "3d" {
                j.march3d()
            } else {
                j.march2d()
            };
            if let Err(e) = r {
                *j.failed.lock().unwrap() = Some(e);
            }
            j.pct.store(100, Ordering::Relaxed);
            j.done.store(true, Ordering::Relaxed);
        });
        job
    }

    pub fn status_json(&self) -> String {
        let guard = self.job.lock().unwrap();
        match guard.as_ref() {
            Some(j) => j.status_json(),
            None => {
                format!(
                    "{{\"running\":false,\"done\":false,\"failed\":false,\"runs\":{}}}",
                    self.run_count.load(Ordering::Relaxed)
                )
            }
        }
    }
}

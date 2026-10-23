//! progress.json writer: a tiny live-run heartbeat the case viewer
//! polls. every march loop calls `tick` on its emit cadence and `finish`
//! at the end; the file is written tmp-then-rename so a reader never
//! sees a torn json.

use std::io::Write;
use std::path::{Path, PathBuf};

/// how many history samples to keep in the file (the page downsamples
/// further for plotting).
const MAX_SAMPLES: usize = 600;

pub struct Progress {
    path: PathBuf,
    kind: String,
    samples: Vec<Sample>,
    last_emit_step: u64,
    emit_every: u64,
    t0: std::time::Instant,
}

struct Sample {
    step: u64,
    t: f64,
    dt: f64,
    rho_min: f64,
    rho_max: f64,
    p_max: f64,
    wall: f64,
    // the cae readout: drag/lift coefficients and the peak wall heat
    // flux when the case carries a body, plus a march residual and the
    // peak temperature through the run's own thermo closure.
    cd: Option<f64>,
    cl: Option<f64>,
    q_max: Option<f64>,
    residual: f64,
    t_max: Option<f64>,
}

/// where the heartbeat lives: next to the case file, `progress.json`.
/// a case run from a path without an extension writes nothing.
pub fn path_for(case_path: &str) -> Option<PathBuf> {
    let p = Path::new(case_path);
    p.extension()?;
    p.parent().map(|d| d.join("progress.json"))
}

/// the optional cae readout a march can attach to a tick. everything
/// defaults to none, so a plain euler case costs nothing.
pub struct SampleExtra {
    pub cd: Option<f64>,
    pub cl: Option<f64>,
    pub q_max: Option<f64>,
    pub residual: f64,
    pub t_max: Option<f64>,
    /// the peak cell pressure through the run's own closure. none
    /// means the tick computes the perfect-gas form, which is right
    /// for the perfect-gas lanes only.
    pub p_peak: Option<f64>,
}

impl Default for SampleExtra {
    fn default() -> Self {
        SampleExtra {
            cd: None,
            cl: None,
            q_max: None,
            residual: 0.0,
            t_max: None,
            p_peak: None,
        }
    }
}

impl Progress {
    /// kind names the march family for the page header.
    pub fn new(case_path: &str, kind: &str) -> Option<Self> {
        let path = path_for(case_path)?;
        let slf = Self {
            path,
            kind: kind.to_string(),
            samples: Vec::new(),
            last_emit_step: 0,
            emit_every: 50,
            t0: std::time::Instant::now(),
        };
        slf.write(&format!(
            "{{\"kind\":\"{kind}\",\"state\":\"starting\",\"step\":0,\"t\":0,\"samples\":[]}}"
        ));
        Some(slf)
    }

    /// record one sample; the first step always samples (a short run
    /// must leave at least one point), then the emit cadence applies.
    /// `extra` carries the cae readout the caller has at hand: cd/cl,
    /// peak wall heat flux, and the march residual.
    pub fn tick(
        &mut self,
        step: u64,
        t: f64,
        dt: f64,
        st: &nufor_core::ConservedState2d,
        extra: SampleExtra,
    ) {
        if step != 1 && step < self.last_emit_step + self.emit_every {
            return;
        }
        self.push_sample(step, t, dt, st, extra);
    }

    /// the per-step heartbeat: the extrema scan stays on the emit
    /// cadence (it walks the whole field), but the header fields the
    /// console's progress box reads (step, t, dt) refresh every step
    /// so the percent and eta move continuously.
    pub fn step_tick(&mut self, step: u64, t: f64, dt: f64) {
        self.write_state("running", step, t, dt);
    }

    /// the 3d march ticks with the same sample shape (rho, p extrema).
    pub fn tick_3d(
        &mut self,
        step: u64,
        t: f64,
        dt: f64,
        st: &nufor_core::ConservedState3d,
        extra: SampleExtra,
    ) {
        if step != 1 && step < self.last_emit_step + self.emit_every {
            return;
        }
        let fake = nufor_core::ConservedState2d {
            rho: st.rho.clone(),
            mx: st.mx.clone(),
            my: st.my.clone(),
            e: st.e.clone(),
        };
        self.push_sample(step, t, dt, &fake, extra);
    }

    /// the 1d euler solver exposes its per-step log; the viewer gets
    /// the whole history at finish time instead of a live tick.
    pub fn finish_1d(
        &mut self,
        ok: bool,
        log: &[nufor_core::EulerLog],
        converged: bool,
        residual: f64,
    ) {
        let mut body = String::from("{\"kind\":\"1d\",\"state\":\"");
        body.push_str(if ok { "finished" } else { "failed" });
        body.push_str("\",\"step\":");
        body.push_str(&log.len().to_string());
        body.push_str(",\"t\":");
        if let Some(l) = log.last() {
            body.push_str(&format!("{:.6}", l.time));
        } else {
            body.push('0');
        }
        body.push_str(",\"dt\":0,\"converged\":");
        body.push_str(if converged { "true" } else { "false" });
        body.push_str(",\"residual\":");
        body.push_str(&format!("{:.3e}", residual));
        body.push_str(",\"samples\":[");
        let take = log.len().min(MAX_SAMPLES);
        let skip = log.len() - take;
        let mut first = true;
        for l in &log[skip..] {
            if !first {
                body.push(',');
            }
            first = false;
            body.push_str(&format!(
                "{{\"step\":{},\"t\":{:.6},\"dt\":{:.3e},\"residual\":{:.3e}}}",
                l.step, l.time, l.dt, l.residual
            ));
        }
        body.push_str("]}");
        self.write(&body);
    }

    fn push_sample(
        &mut self,
        step: u64,
        t: f64,
        dt: f64,
        st: &nufor_core::ConservedState2d,
        extra: SampleExtra,
    ) {
        self.last_emit_step = step;
        let n = st.rho.len();
        let mut rho_min = f64::INFINITY;
        let mut rho_max = f64::NEG_INFINITY;
        // the peak pressure comes from the run's own closure when the
        // caller supplies it: the perfect-gas form (gamma-1) * (e -
        // rho*ke) misreads the eqair and combustion lanes, whose
        // pressure is not a function of the internal energy alone.
        let p_max = match extra.p_peak {
            Some(p) => p,
            None => {
                let mut p_max = 0.0f64;
                for i in 0..n {
                    let e_int = st.e[i]
                        - 0.5 * (st.mx[i] * st.mx[i] + st.my[i] * st.my[i]) / st.rho[i].max(1e-30);
                    let p_i = 0.4 * e_int;
                    if p_i > 0.0 {
                        p_max = p_max.max(p_i);
                    }
                }
                p_max
            }
        };
        for i in 0..n {
            rho_min = rho_min.min(st.rho[i]);
            rho_max = rho_max.max(st.rho[i]);
        }
        let wall = self.t0.elapsed().as_secs_f64();
        self.samples.push(Sample {
            step,
            t,
            dt,
            rho_min,
            rho_max,
            p_max,
            wall,
            cd: extra.cd,
            cl: extra.cl,
            q_max: extra.q_max,
            residual: extra.residual,
            t_max: extra.t_max,
        });
        if self.samples.len() > MAX_SAMPLES {
            self.samples.remove(0);
        }
        self.write_state("running", step, t, dt);
    }

    fn write_state(&self, state: &str, step: u64, t: f64, dt: f64) {
        let samples: Vec<String> = self
            .samples
            .iter()
            .map(|s| {
                let mut o = format!(
                    "{{\"step\":{},\"t\":{:.6},\"dt\":{:.3e},\"rho_min\":{:.5},\"rho_max\":{:.5},\"p_max\":{:.1},\"wall\":{:.1}",
                    s.step, s.t, s.dt, s.rho_min, s.rho_max, s.p_max, s.wall
                );
                if let Some(cd) = s.cd {
                    o.push_str(&format!(",\"cd\":{cd:.4},\"cl\":{:.4}", s.cl.unwrap_or(0.0)));
                }
                if let Some(q) = s.q_max {
                    o.push_str(&format!(",\"q_max\":{q:.3e}"));
                }
                if s.residual > 0.0 {
                    o.push_str(&format!(",\"residual\":{:.3e}", s.residual));
                }
                if let Some(t) = s.t_max {
                    o.push_str(&format!(",\"t_max\":{t:.1}"));
                }
                o.push('}');
                o
            })
            .collect();
        self.write(&format!(
            "{{\"kind\":\"{}\",\"state\":\"{}\",\"step\":{},\"t\":{:.6},\"dt\":{:.3e},\"samples\":[{}]}}",
            self.kind,
            state,
            step,
            t,
            dt,
            samples.join(",")
        ));
    }

    /// running -> finished or failed; writes once, no more ticking.
    pub fn finish(&mut self, ok: bool, step: u64, t: f64) {
        self.write_state(if ok { "finished" } else { "failed" }, step, t, 0.0);
    }

    fn write(&self, body: &str) {
        if self.path.parent().is_some() {
            let tmp = self.path.with_extension("json.tmp");
            if let Ok(mut f) = std::fs::File::create(&tmp) {
                if f.write_all(body.as_bytes()).is_ok() {
                    let _ = std::fs::rename(&tmp, &self.path);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_writes_a_consistent_heartbeat() {
        let dir = std::env::temp_dir().join("nufor-progress-test");
        std::fs::create_dir_all(&dir).unwrap();
        let case = dir.join("case.toml");
        std::fs::write(&case, "").unwrap();
        let cpath = case.to_str().unwrap().to_string();
        let mut p = Progress::new(&cpath, "2d").unwrap();
        let st = nufor_core::ConservedState2d {
            rho: vec![1.0, 2.0],
            mx: vec![0.0; 2],
            my: vec![0.0; 2],
            e: vec![2.5e5, 2.5e5],
        };
        p.tick(
            1,
            1e-3,
            1e-3,
            &st,
            SampleExtra {
                cd: Some(0.5),
                cl: Some(0.1),
                ..Default::default()
            },
        );
        p.tick(50, 5e-2, 1e-3, &st, SampleExtra::default());
        p.finish(true, 50, 5e-2);
        let raw = std::fs::read_to_string(dir.join("progress.json")).unwrap();
        assert!(raw.contains("\"state\":\"finished\""));
        assert!(raw.contains("\"step\":50"));
        assert!(raw.contains("\"rho_max\":2"));
        // the tmp file is gone after the rename
        assert!(!dir.join("progress.json.tmp").exists());
    }

    #[test]
    fn a_short_run_still_records_its_first_step() {
        let dir = std::env::temp_dir().join("nufor-progress-short");
        std::fs::create_dir_all(&dir).unwrap();
        let case = dir.join("case.toml");
        std::fs::write(&case, "").unwrap();
        let cpath = case.to_str().unwrap().to_string();
        let mut p = Progress::new(&cpath, "2d").unwrap();
        let st = nufor_core::ConservedState2d {
            rho: vec![1.0; 3],
            mx: vec![0.0; 3],
            my: vec![0.0; 3],
            e: vec![2.5e5; 3],
        };
        // 12 steps: below the 50-step cadence, the heartbeat must still
        // carry the first sample.
        for s in 1..=12u64 {
            p.tick(s, s as f64 * 1e-3, 1e-3, &st, SampleExtra::default());
        }
        p.finish(true, 12, 12e-3);
        let raw = std::fs::read_to_string(dir.join("progress.json")).unwrap();
        assert!(raw.contains("\"samples\":[{\"step\":1"));
        assert!(raw.contains("\"state\":\"finished\""));
    }

    #[test]
    fn path_for_rejects_extensionless_paths() {
        assert!(path_for("sod").is_none());
        assert_eq!(
            path_for("cases/foo/case.toml"),
            Some(PathBuf::from("cases/foo/progress.json"))
        );
    }
}

//! the case viewer: `nufor watch <case.toml>` serves one page that polls
//! the run's progress.json (the march heartbeat) plus the finished-run
//! artifacts (run-log.txt, the vtk's presence). read-only window on
//! files the CLI already owns; the case itself still runs in a
//! terminal with `nufor run`.

use std::io::{Read, Write};
use std::net::TcpListener;
use std::path::{Path, PathBuf};

/// finished-field renders keyed by (field name, vtk mtime); a small
/// cache since each vtk parse is the expensive part. the tuple is the
/// png plus the color window the colorbar reads.
type FieldCache = std::collections::HashMap<(String, u128), (Vec<u8>, f64, f64)>;
static FIELD_CACHE: std::sync::OnceLock<std::sync::Mutex<FieldCache>> = std::sync::OnceLock::new();

fn field_cache() -> &'static std::sync::Mutex<FieldCache> {
    FIELD_CACHE.get_or_init(|| std::sync::Mutex::new(std::collections::HashMap::new()))
}

/// the per-cell values the probe reads, parsed once per vtk and cached
/// by mtime: density, pressure, temperature, mach, and the grid dims.
#[derive(Clone)]
struct CellValues {
    nx: usize,
    ny: usize,
    rho: Vec<f64>,
    p: Vec<f64>,
    t: Vec<f64>,
    mach: Vec<f64>,
}

static CELL_CACHE: std::sync::OnceLock<std::sync::Mutex<Option<(u128, CellValues)>>> =
    std::sync::OnceLock::new();

fn cell_cache() -> &'static std::sync::Mutex<Option<(u128, CellValues)>> {
    CELL_CACHE.get_or_init(|| std::sync::Mutex::new(None))
}

pub struct Watch {
    case_path: PathBuf,
    progress_path: Option<PathBuf>,
    pub dir: PathBuf,
    pub name: String,
}

const PAGE: &str = include_str!("watch_page.html");

/// the fields the api serves; a bogus name is a 404, not a silent rho.
const FIELDS: &[&str] = &["rho", "mach", "p", "t", "schlieren"];

impl Watch {
    pub fn new(case_path: &str) -> Option<Self> {
        let p = Path::new(case_path);
        if !p.exists() {
            return None;
        }
        let dir = p.parent()?.to_path_buf();
        // the march names its vtk/log files after the case's metadata
        // name, not the toml filename; the stem is only a fallback.
        let name = nufor_config::load_case_config(p)
            .map(|c| c.metadata.name)
            .ok()
            .or_else(|| {
                p.file_stem()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string())
            })?;
        Some(Self {
            case_path: p.to_path_buf(),
            progress_path: crate::progress::path_for(case_path),
            dir,
            name,
        })
    }

    pub fn run(&self, port: u16) -> std::io::Result<()> {
        let listener = TcpListener::bind(("0.0.0.0", port))?;
        println!(
            "watching case {} on http://0.0.0.0:{port}/ (Ctrl-C to stop)",
            self.name
        );
        for mut stream in listener.incoming().flatten() {
            let _ = self.serve_one(&mut stream);
        }
        Ok(())
    }

    fn serve_one(&self, stream: &mut std::net::TcpStream) -> std::io::Result<()> {
        let mut buf = [0u8; 4096];
        let n = stream.read(&mut buf).unwrap_or(0);
        let head = String::from_utf8_lossy(&buf[..n]).to_string();
        let path = head
            .lines()
            .next()
            .and_then(|l| l.split_whitespace().nth(1))
            .unwrap_or("/")
            .to_string();
        let path_only = path.split('?').next().unwrap_or("/").to_string();
        if path_only == "/api/field" || path_only == "/api/range" {
            let name = get_query(&path, "name").unwrap_or_else(|| "rho".to_string());
            if !FIELDS.contains(&name.as_str()) {
                return Self::write_response(
                    stream,
                    "404 Not Found",
                    "text/plain; charset=utf-8",
                    format!("unknown field: {name} (one of {})", FIELDS.join(", ")).into_bytes(),
                );
            }
        }
        if path_only == "/api/field" {
            let name = get_query(&path, "name").unwrap_or_else(|| "rho".to_string());
            return match self.field_png(&name) {
                Some(png) => Self::write_response(stream, "200 OK", "image/png", png),
                None => Self::write_response(
                    stream,
                    "404 Not Found",
                    "text/plain; charset=utf-8",
                    b"no finished field for this case (run it first)".to_vec(),
                ),
            };
        }
        if path_only == "/api/range" {
            let name = get_query(&path, "name").unwrap_or_else(|| "rho".to_string());
            let body = match self.field_render(&name) {
                Some((_, lo, hi)) => format!("{{\"lo\":{lo:.6},\"hi\":{hi:.6}}}"),
                None => "{\"error\":\"no finished field\"}".to_string(),
            };
            return Self::write_response(stream, "200 OK", "application/json", body.into_bytes());
        }
        if path_only == "/api/cell" {
            let i: usize = get_query(&path, "i")
                .and_then(|v| v.parse().ok())
                .unwrap_or(usize::MAX);
            let j: usize = get_query(&path, "j")
                .and_then(|v| v.parse().ok())
                .unwrap_or(usize::MAX);
            return match self.cell_json(i, j) {
                Some(body) => {
                    Self::write_response(stream, "200 OK", "application/json", body.into_bytes())
                }
                None => Self::write_response(
                    stream,
                    "404 Not Found",
                    "text/plain; charset=utf-8",
                    b"no finished field for this case (run it first)".to_vec(),
                ),
            };
        }
        let (status, ct, body) = match path_only.as_str() {
            "/" | "/index.html" => ("200 OK", "text/html; charset=utf-8", PAGE.to_string()),
            "/api/watch" => ("200 OK", "application/json", self.status_json()),
            "/api/case_info" => ("200 OK", "application/json", self.case_info_json()),
            "/api/case" => ("200 OK", "text/plain; charset=utf-8", self.case_toml()),
            "/api/log" => ("200 OK", "text/plain; charset=utf-8", self.run_log()),
            "/api/mesh" => match self.mesh_json() {
                Ok(j) => ("200 OK", "application/json", j),
                Err(e) => (
                    "200 OK",
                    "application/json",
                    format!("{{\"error\":\"{e}\"}}"),
                ),
            },
            _ => ("404 Not Found", "text/plain; charset=utf-8", String::new()),
        };
        Self::write_response(stream, status, ct, body.into_bytes())
    }

    fn write_response(
        stream: &mut std::net::TcpStream,
        status: &str,
        ct: &str,
        body: Vec<u8>,
    ) -> std::io::Result<()> {
        let resp = format!(
            "HTTP/1.1 {status}\r\nContent-Type: {ct}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        );
        stream.write_all(resp.as_bytes())?;
        stream.write_all(&body)?;
        Ok(())
    }

    fn status_json(&self) -> String {
        let prog = self
            .progress_path
            .as_ref()
            .and_then(|p| std::fs::read_to_string(p).ok())
            .map(|s| s.trim().to_string())
            .unwrap_or_else(|| {
                // a case with no heartbeat yet: the page says so plainly.
                "{\"kind\":\"none\",\"state\":\"no-run\",\"step\":0,\"t\":0,\"samples\":[]}"
                    .to_string()
            });
        // inject the case name and the sibling-artifact flags the page
        // renders as tabs.
        let has_vtk = self.dir.join(format!("{}.vtk", self.name)).exists();
        let has_log = self.dir.join("run-log.txt").exists();
        let has_last = self.dir.join("run-log-last.txt").exists();
        let mut obj = prog.trim_end_matches('}').to_string();
        obj.push_str(&format!(
            ",\"case\":\"{}\",\"vtk\":{},\"log\":{},\"last\":{}}}",
            self.name, has_vtk, has_log, has_last
        ));
        obj
    }

    fn case_toml(&self) -> String {
        std::fs::read_to_string(&self.case_path).unwrap_or_else(|_| "# case file unreadable".into())
    }

    /// the left rail's case readout: freestream, closure, numerics, and
    /// stop conditions, straight from the config. one fetch, no polling.
    fn case_info_json(&self) -> String {
        let Some(cfg) = nufor_config::load_case_config(&self.case_path).ok() else {
            return "{\"error\":\"case unreadable\"}".to_string();
        };
        let esc = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"");
        let eos = match cfg.physics.eos {
            nufor_config::Eos::Perfect => "perfect",
            nufor_config::Eos::EqAir => "eqair",
            nufor_config::Eos::Park2t => "park2t",
            nufor_config::Eos::Park5 => "park5",
        };
        let (frho, fu, fp) = match &cfg.initial_condition {
            nufor_config::InitialCondition::Uniform { rho, u, p } => {
                (format!("{rho:e}"), format!("{u:e}"), format!("{p:e}"))
            }
            nufor_config::InitialCondition::TwoState { left, .. } => (
                format!("{:e}", left.rho),
                format!("{:e}", left.u),
                format!("{:e}", left.p),
            ),
            nufor_config::InitialCondition::Blast { ambient, .. } => (
                format!("{:e}", ambient.rho),
                format!("{:e}", ambient.u),
                format!("{:e}", ambient.p),
            ),
        };
        let mach = {
            // the freestream sound speed through the closure: perfect gas
            // is closed form, eqair reads the table at the fs state.
            let (rho, p) = match &cfg.initial_condition {
                nufor_config::InitialCondition::Uniform { rho, p, .. } => (*rho, *p),
                nufor_config::InitialCondition::TwoState { left, .. } => (left.rho, left.p),
                nufor_config::InitialCondition::Blast { ambient, .. } => (ambient.rho, ambient.p),
            };
            let u = match &cfg.initial_condition {
                nufor_config::InitialCondition::Uniform { u, .. } => *u,
                nufor_config::InitialCondition::TwoState { left, .. } => left.u,
                nufor_config::InitialCondition::Blast { ambient, .. } => ambient.u,
            };
            match nufor_core::model_from_config(&cfg) {
                nufor_core::ThermoModel::Perfect { gamma } => {
                    let a = (gamma * p / rho.max(1e-12)).sqrt();
                    (u / a.max(1e-12)).abs()
                }
                m => {
                    // the closure sound speed at the freestream internal
                    // energy the perfect-gas relation would assign: close
                    // enough for the readout, since the eqair table at the
                    // freestream corner is near the perfect-gas value.
                    let e_int = 1.5 * p / rho.max(1e-12);
                    match nufor_core::thermo_sound(m, &[rho], &[e_int], &[0.0], &[0.0]) {
                        Ok(a) => (u / a.first().copied().unwrap_or(1e-12).max(1e-12)).abs(),
                        Err(_) => f64::NAN,
                    }
                }
            }
        };
        let viscous = cfg.physics.mu.unwrap_or(0.0) > 0.0;
        let wall_t = cfg.physics.wall_temperature.unwrap_or(0.0);
        let flux = match cfg.numerics.flux {
            nufor_config::FluxScheme::Hll => "hll",
            nufor_config::FluxScheme::Hllc => "hllc",
        };
        let recon = match cfg.numerics.reconstruction {
            nufor_config::Reconstruction::FirstOrder => "1st order",
            nufor_config::Reconstruction::Muscl => "muscl",
        };
        format!(
            "{{\"name\":\"{}\",\"eos\":\"{eos}\",\"viscous\":{viscous},\"wall_t\":{wall_t},\"frho\":{frho},\"fu\":{fu},\"fp\":{fp},\"mach\":{mach},\"flux\":\"{flux}\",\"recon\":\"{recon}\",\"cfl\":{},\"final_time\":{},\"max_steps\":{}}}",
            esc(&self.name),
            cfg.numerics.cfl,
            cfg.time.final_time,
            cfg.time.max_steps,
        )
    }

    /// the field render for the visual tab: the finished run's vtk,
    /// recovered to a state and painted through the same renderer the
    /// results ui uses, body silhouette included. ranges are taken
    /// from the data itself (percentiles over live cells) so any case
    /// reads sensibly, not just the ui's canonical ones. cached by
    /// vtk mtime because the parse is the expensive part. returns the
    /// png plus the color window so the page can draw a colorbar.
    fn field_render(&self, name: &str) -> Option<(Vec<u8>, f64, f64)> {
        let vtk = self.dir.join(format!("{}.vtk", self.name));
        let meta = std::fs::metadata(&vtk).ok()?;
        let mtime = meta
            .modified()
            .ok()?
            .duration_since(std::time::UNIX_EPOCH)
            .ok()?
            .as_nanos();
        let key = (name.to_string(), mtime);
        if let Some(cache) = field_cache().lock().unwrap().get(&key) {
            return Some(cache.clone());
        }
        let cfg = nufor_config::load_case_config(&self.case_path).ok()?;
        let model = nufor_core::model_from_config(&cfg);
        let field = crate::casework::case_vtk_read(&vtk, model).ok()?;
        let g = crate::build_2d_mesh_in(&cfg.mesh, Some(&self.dir)).ok()?;
        let body = crate::case_body_fns(&cfg.body);
        let (png, lo, hi) =
            crate::casework::render_with_body_ranged(&g, &field, name, &body, model);
        field_cache()
            .lock()
            .unwrap()
            .insert(key, (png.clone(), lo, hi));
        Some((png, lo, hi))
    }

    fn field_png(&self, name: &str) -> Option<Vec<u8>> {
        self.field_render(name).map(|(png, _, _)| png)
    }

    /// the probe's per-cell table: the vtk's own scalar columns, parsed
    /// once and cached by mtime so each hover costs a hashmap lookup.
    fn cell_values(&self) -> Option<CellValues> {
        let vtk = self.dir.join(format!("{}.vtk", self.name));
        let meta = std::fs::metadata(&vtk).ok()?;
        let mtime = meta
            .modified()
            .ok()?
            .duration_since(std::time::UNIX_EPOCH)
            .ok()?
            .as_nanos();
        if let Some((mt, vals)) = cell_cache().lock().unwrap().as_ref() {
            if *mt == mtime {
                return Some(vals.clone());
            }
        }
        let f = std::fs::read_to_string(&vtk).ok()?;
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
            .filter(|&(nx, ny)| nx > 0 && ny > 0)?;
        let n = nx * ny;
        let col = |name: &str| -> Option<Vec<f64>> {
            let i = f.find(&format!("SCALARS {name}"))?;
            let lt = f[i..].find("LOOKUP_TABLE")? + i;
            let j = f[lt..].find('\n')? + lt + 1;
            let v: Vec<f64> = f[j..]
                .split_whitespace()
                .take(n)
                .filter_map(|x| x.parse::<f64>().ok())
                .collect();
            (v.len() == n).then_some(v)
        };
        let vals = CellValues {
            nx,
            ny,
            rho: col("density")?,
            p: col("pressure")?,
            t: col("temperature").unwrap_or_else(|| vec![f64::NAN; n]),
            mach: col("mach").unwrap_or_else(|| vec![f64::NAN; n]),
        };
        *cell_cache().lock().unwrap() = Some((mtime, vals.clone()));
        Some(vals)
    }

    fn cell_json(&self, i: usize, j: usize) -> Option<String> {
        let v = self.cell_values()?;
        if i >= v.nx || j >= v.ny {
            return None;
        }
        let k = j * v.nx + i;
        let num = |x: f64| -> String {
            if x.is_finite() {
                format!("{x:e}")
            } else {
                "null".to_string()
            }
        };
        Some(format!(
            "{{\"i\":{i},\"j\":{j},\"rho\":{},\"p\":{},\"t\":{},\"mach\":{}}}",
            num(v.rho[k]),
            num(v.p[k]),
            num(v.t[k]),
            num(v.mach[k])
        ))
    }

    /// the mesh tab: the case's grid edges (downsampled) plus the body
    /// polygon verts when the case carries one, as json the page draws.
    fn mesh_json(&self) -> Result<String, String> {
        let cfg = nufor_config::load_case_config(&self.case_path).map_err(|e| e.to_string())?;
        let dir = self.case_path.parent().unwrap_or(Path::new("."));
        // 1d cases have no y axis to draw; report the cell count only.
        if cfg.mesh.ny.is_none() {
            return Ok(format!(
                "{{\"nx\":{},\"ny\":0,\"error\":\"1d case: no mesh to draw\"}}",
                cfg.mesh.nx
            ));
        }
        let g = crate::build_2d_mesh_in(&cfg.mesh, Some(dir))?;
        // cell edges from the cumulative widths; downsample x when wide.
        let mut xs = vec![g.xmin];
        let mut acc = g.xmin;
        for w in &g.dxs {
            acc += w;
            xs.push(acc);
        }
        let mut ys = vec![g.ymin];
        acc = g.ymin;
        for w in &g.dys {
            acc += w;
            ys.push(acc);
        }
        let stride = ((xs.len() - 1) / 80).max(1);
        let xs: Vec<String> = xs
            .iter()
            .step_by(stride)
            .chain(std::iter::once(&xs[xs.len() - 1]))
            .map(|v| format!("{v:.6}"))
            .collect();
        let stride = ((ys.len() - 1) / 80).max(1);
        let ys: Vec<String> = ys
            .iter()
            .step_by(stride)
            .chain(std::iter::once(&ys[ys.len() - 1]))
            .map(|v| format!("{v:.6}"))
            .collect();
        // the body polygon, verbatim from the case.
        let verts: Vec<String> = match &cfg.body {
            Some(nufor_config::BodySection::Polygon { verts }) => verts
                .iter()
                .map(|(x, y)| format!("[{x:.6},{y:.6}]"))
                .collect(),
            _ => Vec::new(),
        };
        Ok(format!(
            "{{\"nx\":{},\"ny\":{},\"xmin\":{:.6},\"xmax\":{:.6},\"ymin\":{:.6},\"ymax\":{:.6},\"xs\":[{}],\"ys\":[{}],\"verts\":[{}]}}",
            g.nx,
            g.ny,
            g.xmin,
            g.xmax,
            g.ymin,
            g.ymax,
            xs.join(","),
            ys.join(","),
            verts.join(",")
        ))
    }

    fn run_log(&self) -> String {
        // the live log first, then the last completed one, then nothing.
        for f in ["run-log.txt", "run-log-last.txt"] {
            let p = self.dir.join(f);
            if p.exists() {
                if let Ok(s) = std::fs::read_to_string(&p) {
                    return s;
                }
            }
        }
        "no run log yet".to_string()
    }
}

/// pull one query parameter out of a request path (?name=value).
fn get_query(path: &str, key: &str) -> Option<String> {
    let q = path.split('?').nth(1)?;
    for pair in q.split('&') {
        let mut it = pair.splitn(2, '=');
        if it.next() == Some(key) {
            return Some(it.next().unwrap_or("").to_string());
        }
    }
    None
}

/// the CLI entry: watch <case.toml> [port].
pub fn watch(args: &[String]) -> i32 {
    let Some(case) = args.get(2) else {
        eprintln!("watch needs a case file: nufor watch <case.toml> [port]");
        return 2;
    };
    let port = args
        .get(3)
        .and_then(|s| s.parse::<u16>().ok())
        .unwrap_or(8061);
    match Watch::new(case) {
        Some(w) => match w.run(port) {
            Ok(()) => 0,
            Err(e) => {
                eprintln!("watch error: {e}");
                1
            }
        },
        None => {
            eprintln!("case not found: {case}");
            2
        }
    }
}

/// spawn the viewer on a background thread for a case that is about
/// to run, so the cli can print the link before the march starts.
/// returns the port it bound (caller prints the url).
pub fn spawn_for_run(case_path: &str) -> Option<u16> {
    let w = Watch::new(case_path)?;
    for port in 8062..8070 {
        if let Ok(listener) = TcpListener::bind(("0.0.0.0", port)) {
            std::thread::spawn(move || {
                for mut stream in listener.incoming().flatten() {
                    let _ = w.serve_one(&mut stream);
                }
            });
            return Some(port);
        }
    }
    None
}

/// write the standalone html report for a finished run: the viewer
/// page with everything inlined (progress history, case toml, run
/// log, mesh geometry, every field render as a data url) so a single
/// file opens offline and prints to pdf for submission.
pub fn write_report(watch: &Watch, out_path: &std::path::Path) -> std::io::Result<()> {
    let mut fields = String::from("{");
    let mut ranges = String::from("{");
    let mut first = true;
    for name in FIELDS {
        if let Some((png, lo, hi)) = watch.field_render(name) {
            let b64 = b64_encode(&png);
            if !first {
                fields.push(',');
                ranges.push(',');
            }
            first = false;
            fields.push_str(&format!("\"{name}\":\"data:image/png;base64,{b64}\""));
            ranges.push_str(&format!("\"{name}\":{{\"lo\":{lo:.6},\"hi\":{hi:.6}}}"));
        }
    }
    fields.push('}');
    ranges.push('}');
    let payload = format!(
        "window.STANDALONE = {{
  watch: {watch},
  log: {log},
  mesh: {mesh},
  fields: {fields},
  ranges: {ranges}
}};",
        watch = serde_json_embed(&watch.status_json()),
        log = json_str(&watch.run_log()),
        mesh = serde_json_embed(
            &watch
                .mesh_json()
                .unwrap_or_else(|e| format!("{{\"error\":\"{e}\"}}"))
        ),
        fields = serde_json_embed(&fields),
        ranges = serde_json_embed(&ranges),
    );
    let mut html = PAGE.to_string();
    html.push_str(&format!("<script>\n{}\n</script>\n</body>", payload));
    std::fs::write(out_path, html)
}

/// standard base64, hand-rolled to keep the dependency list as is.
fn b64_encode(data: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let b = [
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ];
        let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        out.push(T[(n >> 18) as usize & 63] as char);
        out.push(T[(n >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 {
            T[(n >> 6) as usize & 63] as char
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            T[n as usize & 63] as char
        } else {
            '='
        });
    }
    out
}

/// escape a json payload for embedding inside a script block.
fn serde_json_embed(s: &str) -> String {
    s.replace('\\', "\\\\").replace("</", "<\\/")
}

/// encode a raw text value as one json string literal, script-safe.
fn json_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '<' => out.push_str("\\u003c"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn field_png_chain_reports_real_errors() {
        // the capsule e2e dir, only present when the ad-hoc run left it
        let dir = std::path::Path::new("/tmp/nufor-watch-2d");
        let case = dir.join("case.toml");
        if !case.exists() {
            return;
        }
        let vtk = dir.join("aeroshell-capsule-trim.vtk");
        assert!(std::fs::metadata(&vtk).is_ok(), "vtk missing");
        let cfg = nufor_config::load_case_config(&case);
        assert!(cfg.is_ok(), "config: {:?}", cfg.err());
        let cfg = cfg.unwrap();
        let model = nufor_core::model_from_config(&cfg);
        let field = crate::casework::case_vtk_read(&vtk, model);
        assert!(field.is_ok(), "vtk_read: {:?}", field.err());
        let g = crate::build_2d_mesh_in(&cfg.mesh, Some(dir));
        assert!(g.is_ok(), "mesh: {:?}", g.err());
        let field = field.unwrap();
        let g = g.unwrap();
        let body = crate::case_body_fns(&cfg.body);
        let (png, _, _) = crate::casework::render_with_body_ranged(&g, &field, "rho", &body, model);
        assert!(!png.is_empty(), "render produced empty png");
    }
}

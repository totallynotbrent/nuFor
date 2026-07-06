//! a minimal http server for the web ui skeleton.
//!
//! exposes the whole cli as a small data api (config, run, snapshot, export,
//! benchmark, history) plus one app shell page. the front end is intentionally
//! plain so an external design tool can restyle it against these endpoints.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::time::Instant;

use crate::webviews::app_html;

use nufor_core::{
    advance2d_rk2, cons_to_prim2d, euler_solve, grid1d, grid2d, prim_to_cons, prim_to_cons2d,
    probe_line, render_png, riemann, write_csv, write_h5, write_vtk, Boundaries2d, Boundary,
    ConservedState, ConservedState2d, EulerConfig, Grid2d, OutputState, PrimState,
    TerminationReason,
};

pub const DEFAULT_PORT: u16 = 8060;

/// which initial condition the run starts from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaseKind {
    Sod,
    Lax,
}

impl CaseKind {
    fn as_str(self) -> &'static str {
        match self {
            CaseKind::Sod => "sod",
            CaseKind::Lax => "lax",
        }
    }
    fn from_str(s: &str) -> CaseKind {
        if s == "lax" {
            CaseKind::Lax
        } else {
            CaseKind::Sod
        }
    }
}

/// a run configuration, edited from the case panel and passed to the solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunConfig {
    pub kind: CaseKind,
    pub n: usize,
    pub t: f64,
    pub gamma: f64,
    pub cfl: f64,
    pub bc: Boundary,
}

impl Default for RunConfig {
    fn default() -> Self {
        RunConfig {
            kind: CaseKind::Sod,
            n: 200,
            t: 0.2,
            gamma: 1.4,
            cfl: 0.5,
            bc: Boundary::Transmissive,
        }
    }
}

/// a run snapshot plus derived primitives, served as json.
#[derive(Debug, Clone)]
pub struct Snapshot {
    pub n: usize,
    pub gamma: f64,
    pub time: f64,
    pub centers: Vec<f64>,
    pub rho: Vec<f64>,
    pub m: Vec<f64>,
    pub e: Vec<f64>,
    pub u: Vec<f64>,
    pub p: Vec<f64>,
    pub mach: Vec<f64>,
}

/// the outcome of a solve kept in the server as the result summary.
#[derive(Debug, Clone, Copy)]
pub struct RunInfo {
    pub steps: usize,
    pub residual: f64,
    pub time: f64,
    pub reason: TerminationReason,
    pub n: usize,
    pub kind: CaseKind,
}

/// the mutable server state: current config, current snapshot, and run history.
pub struct Server {
    pub config: RunConfig,
    pub snap: Snapshot,
    pub history: Vec<RunInfo>,
}

/// the left/right primitive states for a case.
fn data(kind: CaseKind) -> ((f64, f64, f64), (f64, f64, f64)) {
    match kind {
        CaseKind::Sod => ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
        CaseKind::Lax => ((0.445, 0.698, 3.528), (0.5, 0.0, 0.571)),
    }
}

/// run the configured case to time t and return both the snapshot and outcome.
fn solve(cfg: &RunConfig) -> (Snapshot, RunInfo) {
    let g = grid1d(cfg.n, 0.0, 1.0).expect("grid");
    let (l, r) = data(cfg.kind);
    let mut state = ConservedState {
        rho: vec![0.0; cfg.n],
        m: vec![0.0; cfg.n],
        e: vec![0.0; cfg.n],
    };
    for (i, &x) in g.centers.iter().enumerate() {
        let (rl, u, p) = if x < 0.5 { l } else { r };
        state.rho[i] = rl;
        let et = p / ((cfg.gamma - 1.0) * rl) + 0.5 * u * u;
        let (mi, ei) = prim_to_cons(&[rl], &[u], &[et]).expect("prim");
        state.m[i] = mi[0];
        state.e[i] = ei[0];
    }
    let solver = EulerConfig {
        gamma: cfg.gamma,
        cfl: cfg.cfl,
        dx: 1.0 / cfg.n as f64,
        left: cfg.bc,
        right: cfg.bc,
        max_steps: (16.0 * cfg.t * cfg.n as f64) as usize + 300,
        t_end: cfg.t,
        tol: 0.0,
    };
    let result = euler_solve(&mut state, &solver).expect("solve");
    let mut u = vec![0.0; cfg.n];
    let mut p = vec![0.0; cfg.n];
    let mut mach = vec![0.0; cfg.n];
    for i in 0..cfg.n {
        u[i] = state.m[i] / state.rho[i];
        let et = state.e[i] / state.rho[i];
        p[i] = (cfg.gamma - 1.0) * state.rho[i] * (et - 0.5 * u[i] * u[i]);
        let a = (cfg.gamma * p[i] / state.rho[i]).sqrt();
        mach[i] = u[i].abs() / a;
    }
    let snap = Snapshot {
        n: cfg.n,
        gamma: cfg.gamma,
        time: result.time,
        centers: g.centers,
        rho: state.rho,
        m: state.m,
        e: state.e,
        u,
        p,
        mach,
    };
    let info = RunInfo {
        steps: result.steps,
        residual: result.residual,
        time: result.time,
        reason: result.reason,
        n: cfg.n,
        kind: cfg.kind,
    };
    (snap, info)
}

fn join(arr: &[f64]) -> String {
    arr.iter()
        .map(|v| format!("{v:.6}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn arr_json(key: &str, v: &[f64]) -> String {
    format!("{key}:[{}]", join(v))
}

/// json for the exact solution (density, velocity, pressure) plus the rho L1 error.
fn exact_json(cfg: &RunConfig, centers: &[f64], rho: &[f64]) -> String {
    let (l, r) = data(cfg.kind);
    let left = PrimState {
        rho: l.0,
        u: l.1,
        p: l.2,
    };
    let right = PrimState {
        rho: r.0,
        u: r.1,
        p: r.2,
    };
    let mut er = Vec::with_capacity(centers.len());
    let mut eu = Vec::with_capacity(centers.len());
    let mut ep = Vec::with_capacity(centers.len());
    let mut l1 = 0.0f64;
    for (i, &x) in centers.iter().enumerate() {
        let sol = riemann(left, right, cfg.gamma, x, cfg.t);
        er.push(sol.state.rho);
        eu.push(sol.state.u);
        ep.push(sol.state.p);
        l1 += (sol.state.rho - rho[i]).abs();
    }
    l1 /= centers.len() as f64;
    format!(
        "{{\"l1\":{:.3e},{},{},{}}}",
        l1,
        arr_json("\"rho\"", &er),
        arr_json("\"u\"", &eu),
        arr_json("\"p\"", &ep),
    )
}

/// the full result envelope the front end renders from.
fn envelope(server: &Server) -> String {
    let s = &server.snap;
    let info = server.history.last().copied().unwrap_or(RunInfo {
        steps: 0,
        residual: 0.0,
        time: s.time,
        reason: TerminationReason::TimeEnd,
        n: s.n,
        kind: server.config.kind,
    });
    let exact = exact_json(&server.config, &s.centers, &s.rho);
    let mut body = format!(
        "{{\"steps\":{},\"residual\":{:.3e},\"reason\":\"{:?}\",\"snapshot\":{{\"n\":{},\"gamma\":{},\"time\":{:.4},{},{},{},{},{}",
        info.steps,
        info.residual,
        info.reason,
        s.n,
        s.gamma,
        s.time,
        arr_json("\"centers\"", &s.centers),
        arr_json("\"rho\"", &s.rho),
        arr_json("\"m\"", &s.m),
        arr_json("\"e\"", &s.e),
        arr_json("\"u\"", &s.u),
    );
    body.push_str(&format!(
        ",{},{}",
        arr_json("\"p\"", &s.p),
        arr_json("\"mach\"", &s.mach),
    ));
    body.push_str(&format!("}},\"exact\":{exact}}}"));
    body
}

fn config_json(cfg: &RunConfig) -> String {
    format!(
        "{{\"kind\":\"{}\",\"n\":{},\"t\":{},\"gamma\":{},\"cfl\":{},\"boundary\":\"{}\"}}",
        cfg.kind.as_str(),
        cfg.n,
        cfg.t,
        cfg.gamma,
        cfg.cfl,
        match cfg.bc {
            Boundary::Transmissive => "transmissive",
            Boundary::Reflective => "reflective",
        },
    )
}

fn history_json(server: &Server) -> String {
    let rows: Vec<String> = server
        .history
        .iter()
        .enumerate()
        .map(|(i, r)| {
            format!(
                "{{\"index\":{},\"kind\":\"{}\",\"n\":{},\"steps\":{},\"time\":{:.4},\"residual\":{:.3e},\"reason\":\"{:?}\"}}",
                i + 1, r.kind.as_str(), r.n, r.steps, r.time, r.residual, r.reason
            )
        })
        .collect();
    format!("[{}]", rows.join(","))
}

/// run the solver a fixed number of steps at a few mesh sizes and report throughput.
fn one_solve(n: usize, steps: usize) -> f64 {
    let g = grid1d(n, 0.0, 1.0).unwrap();
    let (l, r) = data(CaseKind::Sod);
    let mut state = ConservedState {
        rho: vec![0.0; n],
        m: vec![0.0; n],
        e: vec![0.0; n],
    };
    for (i, &x) in g.centers.iter().enumerate() {
        let (rl, u, p) = if x < 0.5 { l } else { r };
        state.rho[i] = rl;
        let et = p / (0.4 * rl) + 0.5 * u * u;
        let (mi, ei) = prim_to_cons(&[rl], &[u], &[et]).unwrap();
        state.m[i] = mi[0];
        state.e[i] = ei[0];
    }
    let solver = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0 / n as f64,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: steps,
        t_end: f64::INFINITY,
        tol: 0.0,
    };
    let start = Instant::now();
    let _ = euler_solve(&mut state, &solver);
    start.elapsed().as_secs_f64()
}

fn benchmark_json(steps: usize) -> String {
    let mut rows = Vec::new();
    for &n in &[100usize, 500, 1000] {
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let secs = one_solve(n, steps);
            best = best.min(secs * 1e6 / (steps as f64 * n as f64));
        }
        let rate = 1e6 / best;
        rows.push(format!(
            "{{\"cells\":{},\"us_per_step_per_cell\":{:.3},\"cell_steps_per_second\":{:.0}}}",
            n, best, rate
        ));
    }
    format!("[{}]", rows.join(","))
}

/// write the snapshot in the requested format to a temp file and return its bytes.
fn export_bytes(snap: &Snapshot, format: &str) -> Result<Vec<u8>, String> {
    let ext = if format == "h5" {
        "h5"
    } else if format == "vtk" {
        "vtk"
    } else {
        "csv"
    };
    let path = std::env::temp_dir().join(format!("nufor_export_{}.{}", std::process::id(), ext));
    let st = OutputState {
        centers: &snap.centers,
        rho: &snap.rho,
        m: &snap.m,
        e: &snap.e,
        gamma: snap.gamma,
    };
    let res = match format {
        "h5" => write_h5(&path, &st, snap.time),
        "vtk" => write_vtk(&path, &st),
        _ => write_csv(&path, &st),
    };
    res.map(|_| std::fs::read(&path).unwrap_or_default())
        .map_err(|e| format!("{e}"))
}

fn url_decode(s: &str) -> String {
    s.replace('+', " ")
}

/// convert a query string into a key/value map.
fn parse_query(q: &str) -> Vec<(String, String)> {
    q.split('&')
        .filter(|p| !p.is_empty())
        .map(|p| {
            let (k, v) = p.split_once('=').unwrap_or((p, ""));
            (url_decode(k).to_string(), url_decode(v).to_string())
        })
        .collect()
}

/// build a response tuple from a status, content type, and text body.
fn respond(status: &str, ct: &'static str, body: String) -> (String, &'static str, Vec<u8>) {
    (status.to_string(), ct, body.into_bytes())
}

/// build a resolved 2d blast wave and extract one scalar field from it.
///
/// returns (field, colormap-low, colormap-high); the field layout is row-major
/// (j*nx+i) so both the png renderer and the line probe can read it.
fn blast_field(n: usize, t_run: f64, field: &str) -> (Vec<f64>, f64, f64) {
    const GAMMA: f64 = 1.4;
    let g: Grid2d = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let (cx, cy, r0): (f64, f64, f64) = (0.35, 0.5, 0.2);
    let (rho, mut p) = (vec![1.0; n * n], vec![1.0; n * n]);
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (g.centers_x[j * n + i], g.centers_y[j * n + i]);
            if (x - cx).powi(2) + (y - cy).powi(2) < r0.powi(2) {
                p[j * n + i] = 5.0;
            }
        }
    }
    let u = vec![0.0; n * n];
    let v = vec![0.0; n * n];
    let et: Vec<f64> = p
        .iter()
        .zip(&rho)
        .map(|(pp, r)| pp / (r * (GAMMA - 1.0)))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d::default();
    let mut t = 0.0;
    while t < t_run {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    let (u, v, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    match field {
        "mach" => {
            let m: Vec<f64> = (0..n * n)
                .map(|k| {
                    let pp = (GAMMA - 1.0)
                        * (st.e[k] - 0.5 * (st.mx[k] * st.mx[k] + st.my[k] * st.my[k]) / st.rho[k]);
                    (u[k] * u[k] + v[k] * v[k]).sqrt() / (GAMMA * pp / st.rho[k]).sqrt().max(1e-12)
                })
                .collect();
            (m, 0.0, 3.0)
        }
        "p" => {
            let pr: Vec<f64> = (0..n * n)
                .map(|k| {
                    (GAMMA - 1.0)
                        * (st.e[k] - 0.5 * (st.mx[k] * st.mx[k] + st.my[k] * st.my[k]) / st.rho[k])
                })
                .collect();
            (pr, 1.0, 5.0)
        }
        _ => (st.rho.clone(), 1.0, 2.6),
    }
}

/// render a resolved 2d blast wave field as a png.
fn blast_image(n: usize, t_run: f64, field: &str) -> Vec<u8> {
    let (data, lo, hi) = blast_field(n, t_run, field);
    render_png(&data, n, lo, hi).unwrap()
}

/// sample a 2d blast field along a segment and serialise (distance, value) pairs.
fn probe_json(n: usize, field: &str, x0: f64, y0: f64, x1: f64, y1: f64, samples: usize) -> String {
    let g = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).expect("grid");
    let (data, _, _) = blast_field(n, 0.10, field);
    let pts = probe_line(&data, &g, x0, y0, x1, y1, samples).unwrap_or_default();
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

fn line_coords(line: &str) -> [f64; 4] {
    match line {
        "V" => [0.5, 0.0, 0.5, 1.0],
        "D" => [0.1, 0.1, 0.9, 0.9],
        _ => [0.0, 0.5, 1.0, 0.5],
    }
}

fn compare_json(fields: &[String], n_list: &[usize], line: &str, samples: usize) -> String {
    let c = line_coords(line);
    let mut series: Vec<String> = Vec::new();
    for &n in n_list {
        let g = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).expect("grid");
        for field in fields {
            let (data, _, _) = blast_field(n, 0.10, field);
            let pts = probe_line(&data, &g, c[0], c[1], c[2], c[3], samples).unwrap_or_default();
            let points: Vec<String> = pts
                .iter()
                .map(|(s, v)| format!(r#"{{"s":{:.4},"v":{:.5}}}"#, s, v))
                .collect();
            let obj = format!(
                r#"{{"label":"{} · n={}","field":"{}","n":{},"points":[{}]}}"#,
                field,
                n,
                field,
                n,
                points.join(",")
            );
            series.push(obj);
        }
    }
    format!(r#"{{"line":"{}","series":[{}]}}"#, line, series.join(","))
}

/// (status, content-type, body) for a request path with its query string.
pub fn handle_request(path: &str, server: &mut Server) -> (String, &'static str, Vec<u8>) {
    let (route, query) = match path.split_once('?') {
        Some((r, q)) => (r, q),
        None => (path, ""),
    };
    let q = parse_query(query);
    let get = |k: &str| q.iter().find(|(a, _)| a == k).map(|(_, v)| v);
    match route {
        "/" | "/index.html" => respond("200 OK", "text/html; charset=utf-8", app_html()),
        "/api/result" | "/api/snapshot" => respond("200 OK", "application/json", envelope(server)),
        "/api/config" => respond("200 OK", "application/json", config_json(&server.config)),
        "/api/image" => {
            let n = get("n").and_then(|s| s.parse().ok()).unwrap_or(128);
            let field = get("field").map(String::as_str).unwrap_or("rho");
            let body = blast_image(n, 0.10, field);
            ("200 OK".to_string(), "image/png", body)
        }
        "/api/probe" => {
            let n = get("n").and_then(|s| s.parse().ok()).unwrap_or(128);
            let field = get("field").map(String::as_str).unwrap_or("rho");
            let p: Vec<f64> = [("x0", 0.0), ("y0", 0.5), ("x1", 1.0), ("y1", 0.5)]
                .iter()
                .map(|(k, d)| get(k).and_then(|s| s.parse().ok()).unwrap_or(*d))
                .collect();
            let samples = get("samples").and_then(|s| s.parse().ok()).unwrap_or(40);
            let body = probe_json(n, field, p[0], p[1], p[2], p[3], samples);
            respond("200 OK", "application/json", body)
        }
        "/api/compare" => {
            let fields: Vec<String> = get("fields")
                .map(|s| s.split(',').map(|x| x.to_string()).collect())
                .unwrap_or_else(|| vec!["rho".to_string(), "mach".to_string()]);
            let n_list: Vec<usize> = get("n")
                .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
                .unwrap_or_else(|| vec![64, 128]);
            let line = get("line").map(String::as_str).unwrap_or("H");
            let samples = get("samples").and_then(|s| s.parse().ok()).unwrap_or(48);
            let body = compare_json(&fields, &n_list, line, samples);
            respond("200 OK", "application/json", body)
        }
        "/api/history" => respond("200 OK", "application/json", history_json(server)),
        "/api/run" => {
            if let Some(v) = get("n").and_then(|s| s.parse().ok()) {
                server.config.n = v;
            }
            if let Some(v) = get("t").and_then(|s| s.parse().ok()) {
                server.config.t = v;
            }
            if let Some(v) = get("gamma").and_then(|s| s.parse().ok()) {
                server.config.gamma = v;
            }
            if let Some(v) = get("cfl").and_then(|s| s.parse().ok()) {
                server.config.cfl = v;
            }
            if let Some(v) = get("kind") {
                server.config.kind = CaseKind::from_str(v);
            }
            if let Some(v) = get("bc") {
                server.config.bc = if v == "reflective" {
                    Boundary::Reflective
                } else {
                    Boundary::Transmissive
                };
            }
            let (snap, info) = solve(&server.config);
            server.snap = snap;
            server.history.push(info);
            respond("200 OK", "application/json", envelope(server))
        }
        "/api/exact" => respond(
            "200 OK",
            "application/json",
            exact_json(&server.config, &server.snap.centers, &server.snap.rho),
        ),
        "/api/export" => {
            let format = get("format").map(String::as_str).unwrap_or("csv");
            match format {
                "h5" | "vtk" | "csv" => {
                    let (ctype, body) = match export_bytes(&server.snap, format) {
                        Ok(b) => (
                            match format {
                                "h5" => "application/octet-stream",
                                "vtk" => "text/plain; charset=utf-8",
                                _ => "text/csv; charset=utf-8",
                            },
                            b,
                        ),
                        Err(e) => ("text/plain; charset=utf-8", e.into_bytes()),
                    };
                    ("200 OK".to_string(), ctype, body)
                }
                _ => respond(
                    "400 Bad Request",
                    "text/plain; charset=utf-8",
                    "bad format".to_string(),
                ),
            }
        }
        "/api/benchmark" => {
            let steps = get("steps").and_then(|s| s.parse().ok()).unwrap_or(1500);
            respond("200 OK", "application/json", benchmark_json(steps))
        }
        _ => respond(
            "404 Not Found",
            "text/plain; charset=utf-8",
            "not found".to_string(),
        ),
    }
}

/// serves requests on all interfaces (0.0.0.0:port) until the process is stopped.
pub fn run(port: u16) -> std::io::Result<()> {
    let (snap, info) = solve(&RunConfig::default());
    let mut server = Server {
        config: RunConfig::default(),
        snap,
        history: vec![info],
    };
    let listener = TcpListener::bind(("0.0.0.0", port))?;
    println!("serving on http://0.0.0.0:{port}/ (Ctrl-C to stop)");
    for mut stream in listener.incoming().flatten() {
        let _ = serve_one(&mut stream, &mut server);
    }
    Ok(())
}

/// reads one request and writes its response.
fn serve_one(stream: &mut TcpStream, server: &mut Server) -> std::io::Result<()> {
    let mut buf = [0u8; 2048];
    let n = stream.read(&mut buf).unwrap_or(0);
    let head = String::from_utf8_lossy(&buf[..n]);
    let line = head.lines().next().unwrap_or("");
    let path = line.split_whitespace().nth(1).unwrap_or("/").to_string();
    let (status, ctype, body) = handle_request(&path, server);
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    stream.write_all(head.as_bytes())?;
    stream.write_all(&body)?;
    stream.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn server() -> Server {
        let (snap, info) = solve(&RunConfig::default());
        Server {
            config: RunConfig::default(),
            snap,
            history: vec![info],
        }
    }

    fn body(b: &[u8]) -> String {
        String::from_utf8_lossy(b).into_owned()
    }

    #[test]
    fn probe_endpoint_returns_samples_across_a_horizontal_midplane() {
        let mut s = server();
        let (st, ct, b) = handle_request(
            "/api/probe?n=64&field=rho&x0=0.05&y0=0.5&x1=0.95&y1=0.5&samples=9",
            &mut s,
        );
        assert_eq!(st, "200 OK");
        assert_eq!(ct, "application/json");
        let t = body(&b);
        // a structural check that avoids brittle quote-escape matching.
        assert!(t.starts_with('{') && t.ends_with('}'));
        assert!(t.contains("samples") && t.contains("rho"));
        // nine samples each contribute a comma between the s and v pair.
        assert!(t.matches(',').count() >= 9, "sample rows present: {}", t);
        // every value is finite (no nan/inf escapes the wire).
        assert!(!t.contains("nan") && !t.contains("inf"));
    }
    #[test]
    fn compare_endpoint_overlays_fields_across_resolutions() {
        let mut s = server();
        let (st, ct, b) = handle_request(
            "/api/compare?fields=rho,mach&n=32,64&line=H&samples=12",
            &mut s,
        );
        assert_eq!(st, "200 OK");
        assert_eq!(ct, "application/json");
        let t = body(&b);
        assert!(t.starts_with('{') && t.ends_with('}'));
        assert!(t.contains("series") && t.contains("label"));
        assert!(t.contains("rho") && t.contains("mach"));
        assert!(t.contains("n=32") && t.contains("n=64"));
        // two fields x two resolutions = four series.
        assert_eq!(t.matches(r#""label""#).count(), 4);
        assert!(
            !t.contains("nan") && !t.contains("inf"),
            "no invalid values: {}",
            t
        );
    }

    #[test]
    fn snapshot_endpoint_returns_the_envelope() {
        let mut s = server();
        let (st, ct, b) = handle_request("/api/result", &mut s);
        assert_eq!(st, "200 OK");
        assert_eq!(ct, "application/json");
        let t = body(&b);
        assert!(t.contains("\"n\":200"));
        assert!(t.contains("\"mach\":["));
        assert!(t.contains("\"exact\":"));
    }

    #[test]
    fn run_endpoint_reconfigures_and_runs() {
        let mut s = server();
        let (st, _, b) = handle_request(
            "/api/run?kind=lax&n=100&t=0.15&gamma=1.4&cfl=0.5&bc=transmissive",
            &mut s,
        );
        assert_eq!(st, "200 OK");
        assert!(body(&b).contains("\"n\":100"));
        assert_eq!(s.history.len(), 2);
        assert_eq!(s.snap.n, 100);
    }

    #[test]
    fn history_tracks_runs() {
        let mut s = server();
        let _ = handle_request("/api/run?kind=sod&n=50&t=0.1", &mut s);
        let _ = handle_request("/api/run?kind=lax&n=80&t=0.2", &mut s);
        let (_, _, b) = handle_request("/api/history", &mut s);
        let t = body(&b);
        assert!(t.contains("\"n\":50"));
        assert!(t.contains("\"n\":80"));
    }

    #[test]
    fn export_writes_all_three_formats() {
        let mut s = server();
        let (_, _, b) = handle_request("/api/export?format=csv", &mut s);
        assert!(body(&b).contains(",rho"));
        let (_, _, b) = handle_request("/api/export?format=vtk", &mut s);
        assert!(body(&b).contains("SCALARS"));
        let (st, ct, b) = handle_request("/api/export?format=h5", &mut s);
        assert_eq!(st, "200 OK");
        assert_eq!(ct, "application/octet-stream");
        // the hdf5 magic signature 89 48 44 46 (\x89HDF)
        assert!(b.starts_with(&[0x89, b'H', b'D', b'F']));
    }

    #[test]
    fn benchmark_returns_a_json_table() {
        let mut s = server();
        let (_, _, b) = handle_request("/api/benchmark?steps=300", &mut s);
        let t = body(&b);
        assert!(t.contains("\"cells\":"));
        assert!(t.contains("cell_steps_per_second"));
    }

    #[test]
    fn index_is_served_and_unknown_is_404() {
        let mut s = server();
        let (st, ct, b) = handle_request("/", &mut s);
        assert_eq!(st, "200 OK");
        assert_eq!(ct, "text/html; charset=utf-8");
        assert!(body(&b).contains("Solve"));
        let (st, _, b) = handle_request("/nope", &mut s);
        assert_eq!(st, "404 Not Found");
        assert!(body(&b) == "not found");
    }

    #[test]
    fn serves_data_over_http() {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let addr = listener.local_addr().unwrap();
        let mut server = server();
        std::thread::spawn(move || {
            if let Ok((mut stream, _)) = listener.accept() {
                let _ = serve_one(&mut stream, &mut server);
            }
        });
        use std::io::{Read, Write};
        let mut sock = std::net::TcpStream::connect(addr).unwrap();
        sock.write_all(b"GET /api/result HTTP/1.1\r\nHost: localhost\r\n\r\n")
            .unwrap();
        let mut resp = String::new();
        let mut buf = [0u8; 8192];
        loop {
            let n = sock.read(&mut buf).unwrap();
            if n == 0 {
                break;
            }
            resp.push_str(&String::from_utf8_lossy(&buf[..n]));
        }
        assert!(resp.starts_with("HTTP/1.1 200 OK"));
        assert!(resp.contains("\"rho\":["));
        assert!(resp.contains("\"mach\":["));
    }
}

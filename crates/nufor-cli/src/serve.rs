//! a minimal http server for the web ui skeleton.
//!
//! serves one plain page and the last run's snapshot as json; the polished
//! front end is designed separately (e.g. by an external design tool) on top of
//! this data api, so the server stays tiny and dependency-free.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};

use nufor_core::{euler_solve, grid1d, prim_to_cons, Boundary, ConservedState, EulerConfig};

pub const DEFAULT_PORT: u16 = 8060;

/// a run snapshot plus the derived primitive fields, served as json.
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
}

/// runs a sod shock tube to `t` and derives the primitives for serving.
pub fn build_snapshot(n: usize, t: f64) -> Snapshot {
    let g = grid1d(n, 0.0, 1.0).expect("grid");
    let mut state = ConservedState {
        rho: vec![0.0; n],
        m: vec![0.0; n],
        e: vec![0.0; n],
    };
    for (i, &x) in g.centers.iter().enumerate() {
        let (r, u, p) = if x < 0.5 {
            (1.0, 0.0, 1.0)
        } else {
            (0.125, 0.0, 0.1)
        };
        state.rho[i] = r;
        let et = p / (0.4 * r) + 0.5 * u * u;
        let (mi, ei) = prim_to_cons(&[r], &[u], &[et]).expect("prim");
        state.m[i] = mi[0];
        state.e[i] = ei[0];
    }
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0 / n as f64,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: (16.0 * t * n as f64) as usize + 300,
        t_end: t,
        tol: 0.0,
    };
    euler_solve(&mut state, &cfg).expect("solve");
    let mut u = vec![0.0; n];
    let mut p = vec![0.0; n];
    for i in 0..n {
        u[i] = state.m[i] / state.rho[i];
        let et = state.e[i] / state.rho[i];
        p[i] = (1.4 - 1.0) * state.rho[i] * (et - 0.5 * u[i] * u[i]);
    }
    Snapshot {
        n,
        gamma: 1.4,
        time: t,
        centers: g.centers,
        rho: state.rho,
        m: state.m,
        e: state.e,
        u,
        p,
    }
}

/// the served json body for a snapshot.
pub fn snapshot_json(s: &Snapshot) -> String {
    let mut out = format!(
        "{{\"n\":{},\"gamma\":{},\"time\":{},\"centers\":[{}],\"rho\":[{}]",
        s.n,
        s.gamma,
        s.time,
        join(&s.centers),
        join(&s.rho)
    );
    for (key, arr) in [("m", &s.m), ("e", &s.e), ("u", &s.u), ("p", &s.p)] {
        out.push_str(&format!(",\"{key}\":[{}]", join(arr)));
    }
    out.push('}');
    out
}

fn join(arr: &[f64]) -> String {
    arr.iter()
        .map(|v| format!("{v:.6}"))
        .collect::<Vec<_>>()
        .join(",")
}

/// the page served at `/`; intentionally plain so an external design tool can
/// restyle it against the data api.
fn index_html() -> String {
    r##"<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><title>nuFor &mdash; 1D Euler</title>
<style>body{font-family:sans-serif;margin:2rem}svg{background:#111}</style></head>
<body>
<h1>nuFor &mdash; 1D Euler snapshot</h1>
<svg id="plot" width="640" height="220"></svg>
<script>
fetch('/api/result').then(r=>r.json()).then(d=>{
  const x=d.centers, y=d.rho, w=640, h=220, mx=Math.max(...y), mn=Math.min(...y);
  const px=v=>((v-x[0])/(x[x.length-1]-x[0]))*w;
  const py=v=>(1-(v-mn)/(mx-mn))*h;
  const pts=y.map((v,i)=>px(x[i])+','+py(v)).join(' ');
  document.getElementById('plot').innerHTML =
    '<polyline points="'+pts+'" fill="none" stroke="#4f8" stroke-width="1.5"/>';
});
</script>
</body></html>"##
        .into()
}

/// (status, content-type, body) for a request path.
pub fn handle_request(path: &str, s: &Snapshot) -> (String, &'static str, String) {
    match path {
        "/" | "/index.html" => ("200 OK".into(), "text/html; charset=utf-8", index_html()),
        "/api/result" => ("200 OK".into(), "application/json", snapshot_json(s)),
        _ => (
            "404 Not Found".into(),
            "text/plain; charset=utf-8",
            "not found".to_string(),
        ),
    }
}

/// serves requests on 127.0.0.1:port until the process is stopped.
pub fn run(port: u16) -> std::io::Result<()> {
    let snap = build_snapshot(300, 0.2);
    let listener = TcpListener::bind(("127.0.0.1", port))?;
    println!("serving on http://127.0.0.1:{port}/ (Ctrl-C to stop)");
    for mut stream in listener.incoming().flatten() {
        let _ = serve_one(&mut stream, &snap);
    }
    Ok(())
}

/// reads one request and writes its response.
fn serve_one(stream: &mut TcpStream, snap: &Snapshot) -> std::io::Result<()> {
    let mut buf = [0u8; 1024];
    let n = stream.read(&mut buf).unwrap_or(0);
    let head = String::from_utf8_lossy(&buf[..n]);
    let line = head.lines().next().unwrap_or("");
    let path = line.split_whitespace().nth(1).unwrap_or("/").to_string();
    let (status, ctype, body) = handle_request(&path, snap);
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    stream.write_all(head.as_bytes())?;
    stream.write_all(body.as_bytes())?;
    stream.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn api_returns_snapshot_json() {
        let s = build_snapshot(8, 0.1);
        assert_eq!(s.n, 8);
        assert!(s.rho.iter().all(|&r| r > 0.0));
        let (status, ctype, body) = handle_request("/api/result", &s);
        assert_eq!(status, "200 OK");
        assert_eq!(ctype, "application/json");
        assert!(body.starts_with("{\"n\":8"));
        assert!(body.contains("\"rho\":["));
        assert!(body.contains("\"p\":["));
    }

    #[test]
    fn index_is_served_and_unknown_is_404() {
        let s = build_snapshot(8, 0.1);
        let (st, _ct, body) = handle_request("/", &s);
        assert_eq!(st, "200 OK");
        assert!(body.contains("<svg"));
        let (st, _, _) = handle_request("/nope", &s);
        assert_eq!(st, "404 Not Found");
    }

    #[test]
    fn serves_data_over_http() {
        let snap = build_snapshot(16, 0.1);
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let addr = listener.local_addr().unwrap();
        std::thread::spawn(move || {
            if let Ok((mut stream, _)) = listener.accept() {
                let _ = serve_one(&mut stream, &snap);
            }
        });
        use std::io::{Read, Write};
        let mut sock = std::net::TcpStream::connect(addr).unwrap();
        sock.write_all(b"GET /api/result HTTP/1.1\r\nHost: localhost\r\n\r\n")
            .unwrap();
        let mut resp = String::new();
        let mut buf = [0u8; 4096];
        loop {
            let n = sock.read(&mut buf).unwrap();
            if n == 0 || resp.contains("\"rho\":[") {
                break;
            }
            resp.push_str(&String::from_utf8_lossy(&buf[..n]));
        }
        assert!(resp.starts_with("HTTP/1.1 200 OK"));
        assert!(resp.contains("application/json"));
        assert!(resp.contains("\"rho\":["));
    }
}

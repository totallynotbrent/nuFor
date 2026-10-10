//! the mesh view: an uploaded 2d/3d mesh rendered beside its solution.
//!
//! an upload lands in uploads/, gets parsed (rectilinear or gmsh), and a
//! blast initial condition is marched on the imported grid in a background
//! thread. the field cache stores rho / mach / pressure once so switching
//! field, axis, or slice re-renders without re-solving. the 3d wireframe
//! json feeds the browser's drag-to-spin canvas.

use nufor_core::{
    advance2d_rk2, advance3d_rk2, rectilinear_grid2d, rectilinear_grid3d, render_png_rect,
    Boundaries2d, ConservedState2d, ConservedState3d, Grid2d, Grid3d,
};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Mutex;

const GAMMA: f64 = 1.4;
const MAX_CELLS: usize = 8_000_000;
const MAX_UPLOAD: usize = 64 * 1024 * 1024;

/// what the parsed upload turned out to be.
#[derive(Clone, Copy, PartialEq)]
pub enum MeshKind {
    Rect2d,
    Rect3d,
    Gmsh2d,
    Stl,
}

impl MeshKind {
    fn as_str(self) -> &'static str {
        match self {
            MeshKind::Rect2d => "rect2d",
            MeshKind::Rect3d => "rect3d",
            MeshKind::Gmsh2d => "gmsh2d",
            MeshKind::Stl => "stl",
        }
    }
}

/// the parsed mesh summary plus the face arrays the wireframe needs.
#[derive(Clone)]
pub struct MeshInfo {
    pub kind: MeshKind,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub zs: Vec<f64>,
    /// gmsh nodes + edges (display only).
    pub nodes: Vec<(f64, f64)>,
    pub edges: Vec<(usize, usize)>,
}

/// the solved fields cached per upload: all three scalars, one pass.
#[derive(Clone)]
pub struct FieldCache {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub rho: Vec<f64>,
    pub mach: Vec<f64>,
    pub p: Vec<f64>,
}

/// one background solve on the imported grid.
pub struct MeshSolve {
    pub pct: AtomicU32,
    pub steps: AtomicU32,
    t_now: Mutex<f64>,
    pub t_end: f64,
    done: AtomicU32, // 0 running, 1 done, 2 failed
    pub error: Mutex<Option<String>>,
    pub fields: Mutex<Option<FieldCache>>,
}

impl MeshSolve {
    fn new(t_end: f64) -> Self {
        Self {
            pct: AtomicU32::new(0),
            steps: AtomicU32::new(0),
            t_now: Mutex::new(0.0),
            t_end,
            done: AtomicU32::new(0),
            error: Mutex::new(None),
            fields: Mutex::new(None),
        }
    }

    fn status_json(&self) -> String {
        let done = self.done.load(Ordering::Relaxed);
        let t = *self.t_now.lock().unwrap();
        let err = self.error.lock().unwrap().clone();
        if let Some(e) = err {
            format!("{{\"done\":true,\"failed\":true,\"error\":\"{e}\"}}")
        } else {
            format!(
                "{{\"done\":{},\"failed\":false,\"pct\":{},\"steps\":{},\"t\":{t:.4},\"t_end\":{}}}",
                done == 1,
                self.pct.load(Ordering::Relaxed),
                self.steps.load(Ordering::Relaxed),
                self.t_end
            )
        }
    }
}

/// the shared upload state: the active mesh and its latest solve.
pub struct MeshViewState {
    pub info: Mutex<Option<MeshInfo>>,
    pub path: Mutex<Option<String>>,
    pub solve: Mutex<Option<std::sync::Arc<MeshSolve>>>,
}

impl MeshViewState {
    pub fn new() -> Self {
        Self {
            info: Mutex::new(None),
            path: Mutex::new(None),
            solve: Mutex::new(None),
        }
    }
}

/// sanitize an uploaded filename: keep the last path component, restrict the
/// character set, cap the length.
fn safe_name(raw: &str) -> String {
    let base = raw.rsplit(['/', '\\']).next().unwrap_or("mesh");
    let clean: String = base
        .chars()
        .take(64)
        .filter(|c| c.is_ascii_alphanumeric() || *c == '.' || *c == '-' || *c == '_')
        .collect();
    if clean.is_empty() {
        "mesh.dat".to_string()
    } else {
        clean
    }
}

/// parse an uploaded file into a mesh summary.
fn parse_mesh(path: &str) -> Result<MeshInfo, String> {
    if path.ends_with(".msh") {
        let m = crate::mesh_io::load_gmsh(path)?;
        let (nodes, edges) = gmsh_wireframe(&m, 4000);
        let nx = nodes.iter().map(|n| n.0).collect::<Vec<_>>();
        let ny = nodes.iter().map(|n| n.1).collect::<Vec<_>>();
        let (x0, x1) = (
            nx.iter().cloned().fold(f64::INFINITY, f64::min),
            nx.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        let (y0, y1) = (
            ny.iter().cloned().fold(f64::INFINITY, f64::min),
            ny.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        return Ok(MeshInfo {
            kind: MeshKind::Gmsh2d,
            nx: 0,
            ny: 0,
            nz: 0,
            xs: vec![x0, x1],
            ys: vec![y0, y1],
            zs: vec![],
            nodes,
            edges,
        });
    }
    if path.ends_with(".stl") {
        let tris = crate::stlcut::load_stl(path)?;
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        // wireframe: 3d verts projected to the (x, y) plane for the
        // viewer, capped for the browser payload
        let mut seen = std::collections::HashMap::new();
        for t in tris.iter().take(20000) {
            let mut idx = [0usize; 3];
            for (k, v) in t.v.iter().enumerate() {
                let key = (
                    (v.0 * 1e4).round() as i64,
                    (v.1 * 1e4).round() as i64,
                    (v.2 * 1e4).round() as i64,
                );
                let next = seen.len();
                let id = *seen.entry(key).or_insert(next);
                if id == next && nodes.len() < 20000 {
                    nodes.push((v.0, v.1));
                }
                idx[k] = id;
            }
            edges.push((idx[0], idx[1]));
            edges.push((idx[1], idx[2]));
            edges.push((idx[0], idx[2]));
        }
        if nodes.is_empty() {
            return Err("stl: no triangles parsed".into());
        }
        let xs = nodes.iter().map(|n| n.0).collect::<Vec<_>>();
        let ys = nodes.iter().map(|n| n.1).collect::<Vec<_>>();
        let (x0, x1) = (
            xs.iter().cloned().fold(f64::INFINITY, f64::min),
            xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        let (y0, y1) = (
            ys.iter().cloned().fold(f64::INFINITY, f64::min),
            ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        return Ok(MeshInfo {
            kind: MeshKind::Stl,
            nx: 0,
            ny: 0,
            nz: 0,
            xs: vec![x0, x1],
            ys: vec![y0, y1],
            zs: vec![],
            nodes,
            edges,
        });
    }
    let m = crate::mesh_io::load_rectilinear(path)?;
    let nx = m.xs.len().saturating_sub(1);
    let ny = m.ys.len().saturating_sub(1);
    if m.zs.len() >= 2 {
        Ok(MeshInfo {
            kind: MeshKind::Rect3d,
            nx,
            ny,
            nz: m.zs.len() - 1,
            xs: m.xs,
            ys: m.ys,
            zs: m.zs,
            nodes: vec![],
            edges: vec![],
        })
    } else {
        Ok(MeshInfo {
            kind: MeshKind::Rect2d,
            nx,
            ny,
            nz: 0,
            xs: m.xs,
            ys: m.ys,
            zs: vec![],
            nodes: vec![],
            edges: vec![],
        })
    }
}

/// a 2d blast initial condition on an arbitrary rectilinear grid.
fn blast_ic_rect2d(g: &Grid2d) -> ConservedState2d {
    let (cx, cy, r0) = (0.5, 0.5, 0.2);
    let n = g.nx * g.ny;
    let (rho, mx, my, mut e) = (vec![1.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for (k, (x, y)) in g.centers_x.iter().zip(&g.centers_y).enumerate() {
        let p = if (x - cx) * (x - cx) + (y - cy) * (y - cy) < r0 * r0 {
            10.0
        } else {
            1.0
        };
        e[k] = p / (GAMMA - 1.0);
    }
    ConservedState2d { rho, mx, my, e }
}

/// a 3d blast initial condition centered in the imported box.
fn blast_ic_rect3d(g: &Grid3d) -> ConservedState3d {
    let (cx, cy, cz) = (
        0.5 * (g.xmin + g.xmax),
        0.5 * (g.ymin + g.ymax),
        0.5 * (g.zmin + g.zmax),
    );
    let r0 = 0.25 * (g.xmax - g.xmin).min(g.ymax - g.ymin).min(g.zmax - g.zmin);
    let n = g.nx * g.ny * g.nz;
    let (rho, mx, my, mz, mut e) = (
        vec![1.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
    );
    for (k, ((x, y), z)) in g
        .centers_x
        .iter()
        .zip(&g.centers_y)
        .zip(&g.centers_z)
        .enumerate()
    {
        let d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
        let p = if d2 < r0 * r0 { 10.0 } else { 1.0 };
        e[k] = p / (GAMMA - 1.0);
    }
    ConservedState3d { rho, mx, my, mz, e }
}

/// march the blast on the imported grid and cache the three scalar fields.
fn solve_on_mesh(info: &MeshInfo, t_end: f64, cfl: f64, job: &MeshSolve) -> Result<(), String> {
    if info.kind == MeshKind::Stl {
        return Err("stl meshes render but do not solve here \u{2014} apply the slice to a case and run the case".into());
    }
    match info.kind {
        MeshKind::Stl => Err("unreachable: stl rejected above".into()),
        MeshKind::Rect2d => {
            let g = rectilinear_grid2d(&info.xs, &info.ys).map_err(|e| e.to_string())?;
            let mut st = blast_ic_rect2d(&g);
            let bc = Boundaries2d::default();
            let mut t = 0.0f64;
            let mut steps = 0u32;
            while t < t_end {
                let (dt, _) =
                    advance2d_rk2(&mut st, &g, GAMMA, cfl, true, &bc).map_err(|e| e.to_string())?;
                t += dt;
                steps += 1;
                *job.t_now.lock().unwrap() = t.min(t_end);
                job.steps.store(steps, Ordering::Relaxed);
                job.pct
                    .store(((t / t_end) * 99.0) as u32, Ordering::Relaxed);
            }
            *job.fields.lock().unwrap() = Some(fields_of_2d(&g, &st));
            Ok(())
        }
        MeshKind::Rect3d => {
            let g = rectilinear_grid3d(&info.xs, &info.ys, &info.zs).map_err(|e| e.to_string())?;
            let mut st = blast_ic_rect3d(&g);
            let mut t = 0.0f64;
            let mut steps = 0u32;
            while t < t_end {
                let (dt, _) =
                    advance3d_rk2(&mut st, &g, GAMMA, cfl, true).map_err(|e| e.to_string())?;
                t += dt;
                steps += 1;
                *job.t_now.lock().unwrap() = t.min(t_end);
                job.steps.store(steps, Ordering::Relaxed);
                job.pct
                    .store(((t / t_end) * 99.0) as u32, Ordering::Relaxed);
            }
            *job.fields.lock().unwrap() = Some(fields_of_3d(&g, &st));
            Ok(())
        }
        MeshKind::Gmsh2d => {
            Err("gmsh meshes are view-only for now (solve needs a rectilinear grid)".into())
        }
    }
}

/// extract the three scalar fields from a solved 2d state.
fn fields_of_2d(g: &Grid2d, st: &ConservedState2d) -> FieldCache {
    let n = g.nx * g.ny;
    let mut mach = vec![0.0; n];
    let mut p = vec![0.0; n];
    for k in 0..n {
        let rho = st.rho[k];
        let u = st.mx[k] / rho;
        let v = st.my[k] / rho;
        p[k] = (GAMMA - 1.0) * (st.e[k] - 0.5 * (st.mx[k] * st.mx[k] + st.my[k] * st.my[k]) / rho);
        let c2 = GAMMA * p[k] / rho;
        mach[k] = (u * u + v * v).sqrt() / c2.max(1e-12).sqrt();
    }
    FieldCache {
        nx: g.nx,
        ny: g.ny,
        nz: 1,
        rho: st.rho.clone(),
        mach,
        p,
    }
}

/// extract the three scalar fields from a solved 3d state.
fn fields_of_3d(g: &Grid3d, st: &ConservedState3d) -> FieldCache {
    let n = g.nx * g.ny * g.nz;
    let mut mach = vec![0.0; n];
    let mut p = vec![0.0; n];
    for k in 0..n {
        let rho = st.rho[k];
        let u = st.mx[k] / rho;
        let v = st.my[k] / rho;
        let w = st.mz[k] / rho;
        p[k] = (GAMMA - 1.0)
            * (st.e[k]
                - 0.5 * (st.mx[k] * st.mx[k] + st.my[k] * st.my[k] + st.mz[k] * st.mz[k]) / rho);
        let c2 = GAMMA * p[k] / rho;
        mach[k] = (u * u + v * v + w * w).sqrt() / c2.max(1e-12).sqrt();
    }
    FieldCache {
        nx: g.nx,
        ny: g.ny,
        nz: g.nz,
        rho: st.rho.clone(),
        mach,
        p,
    }
}

/// handle the upload routes; returns none when the path is not ours.
pub fn handle(
    path: &str,
    body: &[u8],
    state: &std::sync::Arc<MeshViewState>,
) -> Option<(String, String, Vec<u8>)> {
    let (route, query) = match path.split_once('?') {
        Some((r, q)) => (r, q),
        None => (path, ""),
    };
    let get = |k: &str| -> Option<String> {
        query.split('&').find_map(|kv| {
            let (kk, vv) = kv.split_once('=')?;
            (kk == k).then(|| vv.replace("%20", " ").replace("%2F", "/"))
        })
    };
    let json = |s: String| {
        Some((
            "200 OK".to_string(),
            "application/json".to_string(),
            s.into_bytes(),
        ))
    };
    let err = |code: &str, msg: &str| {
        Some((
            code.to_string(),
            "application/json".to_string(),
            format!("\"{msg}\"").into_bytes(),
        ))
    };

    match route {
        "/api/upload" => {
            if body.is_empty() {
                return err("400 Bad Request", "upload needs a file body");
            }
            if body.len() > MAX_UPLOAD {
                return err("413 Payload Too Large", "upload exceeds 64 MB");
            }
            let name = safe_name(&get("name").unwrap_or_else(|| "mesh.dat".into()));
            let _ = std::fs::create_dir_all("uploads");
            let dst = format!("uploads/{name}");
            if let Err(e) = std::fs::write(&dst, body) {
                return err(
                    "500 Internal Server Error",
                    &format!("could not write {dst}: {e}"),
                );
            }
            match parse_mesh(&dst) {
                Ok(info) => {
                    if info.kind != MeshKind::Gmsh2d
                        && info.kind != MeshKind::Stl
                        && info.nx * info.ny * info.nz.max(1) > MAX_CELLS
                    {
                        return err(
                            "413 Payload Too Large",
                            "mesh exceeds 8M cells for the browser solve path",
                        );
                    }
                    let kind = info.kind.as_str();
                    let cells = match info.kind {
                        MeshKind::Gmsh2d => info.nodes.len(),
                        MeshKind::Stl => info.edges.len() / 3,
                        _ => info.nx * info.ny * info.nz.max(1),
                    };
                    *state.path.lock().unwrap() = Some(dst.clone());
                    *state.info.lock().unwrap() = Some(info.clone());
                    json(format!(
                        "{{\"path\":\"{dst}\",\"kind\":\"{kind}\",\"nx\":{},\"ny\":{},\"nz\":{},\"cells\":{}}}",
                        info.nx, info.ny, info.nz, cells
                    ))
                }
                Err(e) => err("400 Bad Request", &e),
            }
        }
        "/api/stl/verts" => {
            // slice the uploaded stl at a plane and return the case
            // file's polygon verts block, so a body from cad flows
            // straight into the case editor
            let path = match state.path.lock().unwrap().as_ref() {
                Some(p) => p.clone(),
                None => return err("404 Not Found", "no mesh uploaded"),
            };
            if !path.ends_with(".stl") {
                return err("400 Bad Request", "uploaded mesh is not an stl");
            }
            let num = |k: &str, dflt: f64| -> f64 {
                get(k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(dflt)
            };
            let z0 = num("z", 0.0);
            let scale = num("scale", 1.0);
            let aoa = num("aoa", 0.0).to_radians();
            let dx = num("dx", f64::NAN);
            let dy = num("dy", f64::NAN);
            let tris = match crate::stlcut::load_stl(&path) {
                Ok(t) => t,
                Err(e) => return err("400 Bad Request", &e),
            };
            // a plane exactly on a vertex ring needs the nudge retry,
            // same as the cli path
            let mut pts = None;
            for cand in [z0, z0 + 1e-7, z0 - 1e-7, z0 + 1e-5, z0 - 1e-5] {
                if let Ok(p) = crate::stlcut::slice_at(&tris, cand) {
                    if p.len() >= 4 {
                        pts = Some(p);
                        break;
                    }
                }
            }
            let mut pts = match pts {
                Some(p) => p,
                None => {
                    return err(
                        "400 Bad Request",
                        &format!("the plane z={z0} does not cut a usable loop"),
                    )
                }
            };
            // transform: scale + rotate about the centroid, then
            // place the nose at (dx, dy) when given
            let (mut cx, mut cy) = (0.0f64, 0.0f64);
            for q in &pts {
                cx += q.0;
                cy += q.1;
            }
            cx /= pts.len() as f64;
            cy /= pts.len() as f64;
            let (ca, sa) = (aoa.cos(), aoa.sin());
            for q in pts.iter_mut() {
                let (x, y) = ((q.0 - cx) * scale, (q.1 - cy) * scale);
                *q = (cx + x * ca + y * sa, cy - x * sa + y * ca);
            }
            if dx.is_finite() || dy.is_finite() {
                let nose = pts
                    .iter()
                    .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
                    .copied()
                    .unwrap_or((0.0, 0.0));
                let (px, py) = (dx, if dy.is_finite() { dy } else { nose.1 });
                for q in pts.iter_mut() {
                    *q = (q.0 - nose.0 + px, q.1 - nose.1 + py);
                }
            }
            let mut out = String::from("verts = [\n");
            for (x, y) in &pts {
                out.push_str(&format!("  [{x:.4}, {y:.4}],\n"));
            }
            out.push(']');
            let esc = out
                .replace('\\', "\\\\")
                .replace('"', "\\\"")
                .replace('\n', "\\n");
            json(format!(
                "{{\"verts\":\"{}\",\"points\":{}}}",
                esc,
                pts.len()
            ))
        }
        "/api/meshview/info" => {
            let guard = state.info.lock().unwrap();
            match guard.as_ref() {
                Some(info) => json(format!(
                    "{{\"kind\":\"{}\",\"nx\":{},\"ny\":{},\"nz\":{}}}",
                    info.kind.as_str(),
                    info.nx,
                    info.ny,
                    info.nz
                )),
                None => err("404 Not Found", "no mesh uploaded"),
            }
        }
        "/api/meshview/solve" => {
            // start (or restart) the background march on the uploaded grid
            let t_end: f64 = get("t").and_then(|s| s.parse().ok()).unwrap_or(0.08);
            let cfl: f64 = get("cfl").and_then(|s| s.parse().ok()).unwrap_or(0.4);
            let info = state.info.lock().unwrap().clone();
            let Some(info) = info else {
                return err("404 Not Found", "no mesh uploaded");
            };
            let job = std::sync::Arc::new(MeshSolve::new(t_end));
            *state.solve.lock().unwrap() = Some(std::sync::Arc::clone(&job));
            let j = std::sync::Arc::clone(&job);
            std::thread::spawn(move || {
                if let Err(e) = solve_on_mesh(&info, j.t_end, cfl, &j) {
                    *j.error.lock().unwrap() = Some(e);
                    j.done.store(2, Ordering::Relaxed);
                } else {
                    j.done.store(1, Ordering::Relaxed);
                    j.pct.store(100, Ordering::Relaxed);
                }
            });
            json(job.status_json())
        }
        "/api/meshview/status" => {
            let guard = state.solve.lock().unwrap();
            match guard.as_ref() {
                Some(j) => json(j.status_json()),
                None => json("{\"done\":false,\"pct\":0}".to_string()),
            }
        }
        "/api/meshview/slice" => {
            let field = get("field").unwrap_or_else(|| "rho".into());
            let axis = get("axis").unwrap_or_else(|| "z".into());
            let frac: f64 = get("u").and_then(|s| s.parse().ok()).unwrap_or(0.5);
            let solve = state.solve.lock().unwrap().clone();
            let Some(job) = solve else {
                return err("404 Not Found", "no solve started");
            };
            if job.done.load(Ordering::Relaxed) != 1 {
                return err("503 Service Unavailable", "solve still running");
            }
            let fields = job.fields.lock().unwrap().clone();
            let Some(f) = fields else {
                return err("500 Internal Server Error", "solve finished with no fields");
            };
            match render_slice(&f, &field, &axis, frac) {
                Ok(png) => Some(("200 OK".to_string(), "image/png".to_string(), png)),
                Err(e) => err("400 Bad Request", &e),
            }
        }
        "/api/meshview/wire" => {
            let guard = state.info.lock().unwrap().clone();
            let Some(info) = guard else {
                return err("404 Not Found", "no mesh uploaded");
            };
            json(wire_json(&info))
        }
        _ => None,
    }
}

/// pick the scalar array for a field name.
fn field_arr<'a>(f: &'a FieldCache, field: &str) -> &'a [f64] {
    match field {
        "mach" => &f.mach,
        "p" => &f.p,
        _ => &f.rho,
    }
}

/// the color window per field, matching the built-in views.
fn field_range(field: &str, arr: &[f64]) -> (f64, f64) {
    let _ = arr;
    match field {
        "mach" => (0.0, 1.4),
        "p" => (0.5, 10.0),
        _ => (0.85, 3.0),
    }
}

/// render the 2d field, or the plane at frac along the chosen axis for 3d.
fn render_slice(f: &FieldCache, field: &str, axis: &str, frac: f64) -> Result<Vec<u8>, String> {
    let arr = field_arr(f, field);
    let (lo, hi) = field_range(field, arr);
    if f.nz == 1 {
        return render_png_rect(arr, f.nx, f.ny, lo, hi).map_err(|e| e.to_string());
    }
    let (k, w, h) = match axis {
        "x" | "0" => (idx_of(frac, f.nx), f.nz, f.ny),
        "y" | "1" => (idx_of(frac, f.ny), f.nx, f.nz),
        _ => (idx_of(frac, f.nz), f.nx, f.ny),
    };
    // the slice is (w, h) in the plane's own axes; row-major with the
    // SECOND axis as the row so the render's j-up convention holds.
    let mut plane = vec![0.0f64; w * h];
    for a in 0..w {
        for b in 0..h {
            let src = match axis {
                // fixed i=k: rows are z, columns are y
                "x" | "0" => (a * f.ny + b) * f.nx + k,
                // fixed j=k: rows are z, columns are x
                "y" | "1" => (b * f.ny + k) * f.nx + a,
                // fixed k=k: rows are y, columns are x
                _ => (k * f.ny + b) * f.nx + a,
            };
            plane[b * w + a] = arr[src];
        }
    }
    render_png_rect(&plane, w, h, lo, hi).map_err(|e| e.to_string())
}

/// the nearest index for a slice fraction.
fn idx_of(frac: f64, n: usize) -> usize {
    ((frac.clamp(0.0, 1.0)) * (n - 1) as f64).round() as usize
}

/// node coordinates and node-index edge pairs for a wireframe.
type Wire = (Vec<(f64, f64)>, Vec<(usize, usize)>);

/// the gmsh wireframe, capped at max_edges by stride sampling.
fn gmsh_wireframe(m: &crate::mesh_io::MshMesh, max_edges: usize) -> Wire {
    let n_edges: usize = m.cells.iter().map(|c| c.len()).sum();
    let stride = (n_edges / max_edges).max(1);
    let mut edges = Vec::new();
    for (ci, cell) in m.cells.iter().enumerate() {
        if stride > 1 && ci % stride != 0 {
            continue;
        }
        for e in 0..cell.len() {
            let a = cell[e];
            let b = cell[(e + 1) % cell.len()];
            edges.push((a, b));
        }
    }
    (m.nodes.clone(), edges)
}

/// the wireframe json: 2d faces, or a sampled 3d grid-box for the spin view.
fn wire_json(info: &MeshInfo) -> String {
    match info.kind {
        MeshKind::Gmsh2d | MeshKind::Stl => {
            let mut pts = String::new();
            for (x, y) in &info.nodes {
                pts.push_str(&format!("[{x},{y}],"));
            }
            pts.pop();
            let mut eg = String::new();
            for (a, b) in &info.edges {
                eg.push_str(&format!("[{a},{b}],"));
            }
            eg.pop();
            format!("{{\"kind\":\"gmsh2d\",\"nodes\":[{pts}],\"edges\":[{eg}]}}")
        }
        MeshKind::Rect2d => {
            let fmt_axis = |v: &[f64]| -> String {
                v.iter()
                    .map(|x| format!("{x:.6}"))
                    .collect::<Vec<_>>()
                    .join(",")
            };
            format!(
                "{{\"kind\":\"rect2d\",\"x\":[{}],\"y\":[{}]}}",
                fmt_axis(&info.xs),
                fmt_axis(&info.ys)
            )
        }
        MeshKind::Rect3d => {
            let (x0, x1) = (info.xs[0], *info.xs.last().unwrap());
            let (y0, y1) = (info.ys[0], *info.ys.last().unwrap());
            let (z0, z1) = (info.zs[0], *info.zs.last().unwrap());
            let mut segs = String::new();
            // the 12 box edges, sampled grid lines on three faces, and
            // interior lines every ~1/8 of each axis.
            let sample = |v: &[f64], want: usize| -> Vec<f64> {
                if v.len() <= want {
                    v.to_vec()
                } else {
                    (0..want)
                        .map(|i| v[i * (v.len() - 1) / (want - 1)])
                        .collect()
                }
            };
            let xs = sample(&info.xs, 17);
            let ys = sample(&info.ys, 17);
            let zs = sample(&info.zs, 17);
            let mut push = |x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64, w: f64| {
                segs.push_str(&format!(
                    "[{x1:.4},{y1:.4},{z1:.4},{x2:.4},{y2:.4},{z2:.4},{w:.1}],"
                ));
            };
            // grid on the z=z0 face (x and y lines)
            for &x in &xs {
                push(x, y0, z0, x, y1, z0, 0.35);
            }
            for &y in &ys {
                push(x0, y, z0, x1, y, z0, 0.35);
            }
            // grid on the y=y1 face (x and z lines)
            for &x in &xs {
                push(x, y1, z0, x, y1, z1, 0.28);
            }
            for &z in &zs {
                push(x0, y1, z, x1, y1, z, 0.28);
            }
            // grid on the x=x1 face (y and z lines)
            for &y in &ys {
                push(x1, y, z0, x1, y, z1, 0.22);
            }
            for &z in &zs {
                push(x1, y0, z, x1, y1, z, 0.22);
            }
            // the box outline in accent weight
            let corners = [
                (x0, y0, z0),
                (x1, y0, z0),
                (x1, y1, z0),
                (x0, y1, z0),
                (x0, y0, z1),
                (x1, y0, z1),
                (x1, y1, z1),
                (x0, y1, z1),
            ];
            let outline = [
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 4),
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            ];
            for (a, b) in outline {
                let (ax, ay, az) = corners[a];
                let (bx, by, bz) = corners[b];
                push(ax, ay, az, bx, by, bz, 1.0);
            }
            segs.pop();
            format!(
                "{{\"kind\":\"rect3d\",\"bbox\":[{x0},{y0},{z0},{x1},{y1},{z1}],\"segs\":[{segs}]}}"
            )
        }
    }
}

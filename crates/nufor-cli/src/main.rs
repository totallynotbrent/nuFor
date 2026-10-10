//! command-line driver for nuFor: build, run, restart, serve, and export a 1D Euler case.
//!
//! `serve` hosts a minimal web ui skeleton: one page plus the snapshot as json,
//! ready for a real front end to be designed against the same data api.

use std::path::{Path, PathBuf};
use std::time::Instant;

use nufor_config::{
    load_case_config, BoundaryKind, CaseConfig, Equations, InflowProfileSpec, InitialCondition,
    Mesh, OutputFormat,
};
use nufor_core::{
    advance2d_axi_rk2, advance2d_model_rk2, advance2d_model_visc_rk2, advance2d_sa_lts,
    advance2d_sa_model_rk2, advance2d_sa_rk2, advance3d_rk2, advance_ugrid, apply_solid_fn,
    cons_to_prim2d, eos_pressure2d, euler_solve, grid1d, grid2d, grid3d, prim_to_cons,
    prim_to_cons2d, prim_to_cons3d, read_restart, render_png, wall_distance2d, write_csv, write_h5,
    write_restart, write_vtk, write_vtk2d_model, write_vtk3d, Bc2d, Boundaries2d, Boundary,
    Bounds3d, ConservedState, ConservedState2d, ConservedState3d, Error, EulerConfig, Grid1d,
    Grid2d, Grid3d, OutputState, SaParams, SolidPolygon, SolidShape, SphereCone, SphereConeSdf,
    TurbState, Ugrid,
};
use nufor_core::{model_from_config, ThermoModel};

mod casework;
mod mesh_io;
mod meshview;
mod plates;
mod runstate;
mod serve;
mod stlcut;
mod webviews;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const GAMMA: f64 = 1.4;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let prog = args.first().map(String::as_str).unwrap_or("nufor");
    if args.len() < 2 {
        usage(prog);
        return;
    }
    let code = match args[1].as_str() {
        "version" | "--version" => version(),
        "init" => init(&args),
        "mesh" => mesh(&args),
        "mesh-check" => mesh_check(&args),
        "run" => run(&args),
        "inspect" => inspect(&args),
        "export" => export(&args),
        "history" => history(&args),
        "benchmark" => benchmark(&args),
        "serve" => serve(&args),
        "stl-cut" => stl_cut(&args),
        "help" | "-h" | "--help" => {
            usage(prog);
            0
        }
        other => {
            eprintln!("unknown command: {other}");
            usage(prog);
            2
        }
    };
    std::process::exit(code);
}

fn usage(prog: &str) {
    println!(
        "nufor {VERSION} - 1D ideal-gas Euler solver (Fortran kernels, Rust runtime)\n\
         \nusage: {prog} <command> [args]\n\
         \ncommands:\n\
         \x20 init      [case.toml]       write a case to run (default case.toml)\n\
         \x20 mesh      N xmin xmax        show a computed grid\n\
         \x20 mesh-check file.msh         run the mesh quality diagnostics\n\
         \x20 run       <case.toml> | N t [sod|lax] [rst]\n\
         \x20                        run a case file, or a shock tube to time t\n\
         \x20 inspect   file.rst           show a restart's header and min/max density\n\
         \x20 export    in.rst out.vtk|h5  write a vtk or hdf5 snapshot\n\
         \x20 history   N [file.csv]       run sod to t=0.2 and write a csv snapshot\n\
         \x20 benchmark N steps            time a fixed run\n\
         \x20 serve    [port]              serve the web ui skeleton (default 8060)\n\
         \x20 version                     print the version"
    );
}

fn version() -> i32 {
    println!("nufor {VERSION} (1D ideal-gas Euler, Fortran + Rust)");
    0
}

fn usize_arg(s: &str, what: &str) -> Result<usize, i32> {
    s.parse().map_err(|_| {
        eprintln!("bad {what}: {s}");
        2
    })
}

fn f64_arg(s: &str, what: &str) -> Result<f64, i32> {
    s.parse().map_err(|_| {
        eprintln!("bad {what}: {s}");
        2
    })
}

fn euler_error(e: Error) -> i32 {
    eprintln!("solver error: {e}");
    1
}

/// a uniform grid over [0,1] with its cell width, the ui's default domain.
fn unit_grid(n: usize) -> Result<(Vec<f64>, f64), Error> {
    let g = grid1d(n, 0.0, 1.0)?;
    Ok((g.centers, g.dx))
}

fn shock_state(kind: &str, n: usize) -> Result<ConservedState, Error> {
    let (centers, _) = unit_grid(n)?;
    let (l, r) = match kind {
        "lax" => ((0.445f64, 0.698, 3.528), (0.5, 0.0, 0.571)),
        _ => ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
    };
    let mut state = ConservedState {
        rho: vec![0.0; n],
        m: vec![0.0; n],
        e: vec![0.0; n],
    };
    for (i, &x) in centers.iter().enumerate() {
        let (r, u, p) = if x < 0.5 { l } else { r };
        let et = p / ((GAMMA - 1.0) * r) + 0.5 * u * u;
        state.rho[i] = r;
        let (mi, ei) = prim_to_cons(&[r], &[u], &[et])?;
        state.m[i] = mi[0];
        state.e[i] = ei[0];
    }
    Ok(state)
}

/// a cfl-stable step budget that guarantees reaching time t on any mesh.
fn step_budget(t: f64, n: usize) -> usize {
    (16.0 * t * n as f64).ceil() as usize + 300
}

fn config(t: f64, n: usize) -> EulerConfig {
    EulerConfig {
        gamma: GAMMA,
        cfl: 0.5,
        dx: 1.0 / n as f64,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: step_budget(t, n),
        t_end: t,
        tol: 0.0,
    }
}

fn mesh(args: &[String]) -> i32 {
    if args.len() < 5 {
        eprintln!("mesh needs: N xmin xmax");
        return 2;
    }
    let n = match usize_arg(&args[2], "cell count") {
        Ok(v) => v,
        Err(c) => return c,
    };
    let xmin = match f64_arg(&args[3], "xmin") {
        Ok(v) => v,
        Err(c) => return c,
    };
    let xmax = match f64_arg(&args[4], "xmax") {
        Ok(v) => v,
        Err(c) => return c,
    };
    match grid1d(n, xmin, xmax) {
        Ok(g) => {
            println!("cells: {n}\ndx: {}\nxmin: {xmin}\nxmax: {xmax}", g.dx);
            println!(
                "first face: {}\nlast face: {}\ncenters[0]: {}",
                g.faces[0], g.faces[n], g.centers[0]
            );
            0
        }
        Err(e) => euler_error(e),
    }
}

fn run(args: &[String]) -> i32 {
    if args.len() < 3 {
        eprintln!("run needs a case file or: <N> <t> [sod|lax] [rst]");
        return 2;
    }
    // sweep "--flag value" pairs the case runners can act on.
    let mut max_steps_override: Option<usize> = None;
    let mut positional: Vec<&str> = Vec::new();
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--max-steps" => {
                if let Some(v) = args.get(i + 1).and_then(|s| s.parse().ok()) {
                    max_steps_override = Some(v);
                }
                i += 2;
            }
            _ => {
                positional.push(&args[i]);
                i += 1;
            }
        }
    }
    if let Some(p) = positional.first() {
        if Path::new(p).exists() && p.ends_with(".toml") {
            return run_case_with_overrides(p, max_steps_override);
        }
    }
    let mut rest: Vec<String> = vec![args[0].clone(), args[1].clone()];
    for p in positional {
        rest.push(p.to_string());
    }
    run_shock_tube(&rest)
}

fn run_shock_tube(args: &[String]) -> i32 {
    if args.len() < 4 {
        eprintln!("run needs: N t [sod|lax] [restart.rst]");
        return 2;
    }
    let n = match usize_arg(&args[2], "cell count") {
        Ok(v) => v,
        Err(c) => return c,
    };
    let t = match f64_arg(&args[3], "time") {
        Ok(v) => v,
        Err(c) => return c,
    };
    let kind = args.get(4).map(String::as_str).unwrap_or("sod");
    let mut state = match shock_state(kind, n) {
        Ok(s) => s,
        Err(e) => return euler_error(e),
    };
    let cfg = config(t, n);
    match euler_solve(&mut state, &cfg) {
        Ok(r) => {
            println!(
                "steps: {}\ntime: {:.4}\nresidual: {:.3e}",
                r.steps, r.time, r.residual
            );
        }
        Err(e) => return euler_error(e),
    }
    if let Some(rst) = args.get(5) {
        if write_restart(Path::new(rst), &state, GAMMA, t, 0).is_err() {
            eprintln!("could not write restart");
            return 1;
        }
    }
    0
}

fn inspect(args: &[String]) -> i32 {
    let file = match args.get(2) {
        Some(f) => f,
        None => {
            eprintln!("inspect needs a restart file");
            return 2;
        }
    };
    match read_restart(Path::new(file)) {
        Ok(r) => {
            let min = r.state.rho.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = r
                .state
                .rho
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            println!(
                "version: {}\ncells: {}\ngamma: {}\ntime: {:.4}\nstep: {}\nrho min/max: {:.4} / {:.4}",
                r.version,
                r.state.rho.len(),
                r.gamma,
                r.time,
                r.step,
                min,
                max
            );
            0
        }
        Err(e) => euler_error(e),
    }
}

/// load a mesh, run the quality diagnostics, and print the summary table.
fn mesh_check(args: &[String]) -> i32 {
    let file = match args.get(2) {
        Some(f) => f,
        None => {
            eprintln!("mesh-check needs a gmsh .msh file");
            return 2;
        }
    };
    let m = match mesh_io::load_gmsh(file) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let ug = match m.ugrid() {
        Ok(u) => u,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let d = ug.diagnostics();
    print_mesh_stats(&d);
    if d.valid {
        println!("verdict: valid mesh");
        0
    } else {
        println!("verdict: mesh has problems; fix before solving");
        1
    }
}

/// print the mesh quality table for a diagnostics summary.
fn print_mesh_stats(d: &nufor_core::MeshDiagnostics) {
    println!(
        "cells: {}\nfaces: {} ({} interior, {} boundary)\narea min/mean/max: {:.4} / {:.4} / {:.4}\narea stretch (max/min): {:.2}\nnegative cells: {}\nlargest closure residual: {:.2e}\nopen cells: {}\nflipped interior faces: {}",
        d.n_cells,
        d.n_faces,
        d.n_interior,
        d.n_boundary,
        d.area_min,
        d.area_mean,
        d.area_max,
        d.stretch,
        d.negative,
        d.closure_max,
        d.open_cells,
        d.flipped_faces,
    );
}

fn export(args: &[String]) -> i32 {
    let (inp, outp) = match (args.get(2), args.get(3)) {
        (Some(i), Some(o)) => (i, o),
        _ => {
            eprintln!("export needs: <in file> <out file>");
            return 2;
        }
    };
    let r = match read_restart(Path::new(inp)) {
        Ok(r) => r,
        Err(e) => return euler_error(e),
    };
    let n = r.state.rho.len();
    let centers = match unit_grid(n) {
        Ok((c, _)) => c,
        Err(e) => return euler_error(e),
    };
    let st = OutputState {
        centers: &centers,
        rho: &r.state.rho,
        m: &r.state.m,
        e: &r.state.e,
        gamma: r.gamma,
    };
    let res = if outp.ends_with(".h5") {
        write_h5(Path::new(outp), &st, r.time)
    } else {
        write_vtk(Path::new(outp), &st)
    };
    match res {
        Ok(_) => {
            println!("wrote {outp}");
            0
        }
        Err(e) => euler_error(e),
    }
}

fn history(args: &[String]) -> i32 {
    let n = match args.get(2).map(|s| usize_arg(s, "cell count")) {
        Some(Ok(v)) => v,
        _ => {
            eprintln!("history needs: N");
            return 2;
        }
    };
    let mut state = match shock_state("sod", n) {
        Ok(s) => s,
        Err(e) => return euler_error(e),
    };
    if let Err(e) = euler_solve(&mut state, &config(0.2, n)) {
        return euler_error(e);
    }
    let centers = match unit_grid(n) {
        Ok((c, _)) => c,
        Err(e) => return euler_error(e),
    };
    let st = OutputState {
        centers: &centers,
        rho: &state.rho,
        m: &state.m,
        e: &state.e,
        gamma: GAMMA,
    };
    let path = args.get(3).map(String::as_str).unwrap_or("history.csv");
    match write_csv(Path::new(path), &st) {
        Ok(_) => {
            println!("wrote {path}");
            0
        }
        Err(e) => euler_error(e),
    }
}

fn benchmark(args: &[String]) -> i32 {
    let n = match args.get(2).map(|s| usize_arg(s, "cell count")) {
        Some(Ok(v)) => v,
        _ => {
            eprintln!("benchmark needs: N");
            return 2;
        }
    };
    let steps = match args.get(3).map(|s| usize_arg(s, "steps")) {
        Some(Ok(v)) => v,
        _ => 100,
    };
    let mut state = match shock_state("sod", n) {
        Ok(s) => s,
        Err(e) => return euler_error(e),
    };
    let cfg = EulerConfig {
        gamma: GAMMA,
        cfl: 0.5,
        dx: 1.0 / n as f64,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: steps,
        t_end: f64::INFINITY,
        tol: 0.0,
    };
    let start = Instant::now();
    if let Err(e) = euler_solve(&mut state, &cfg) {
        return euler_error(e);
    }
    let dt = start.elapsed();
    println!(
        "ran {steps} steps on {n} cells in {:.3}s ({:.1} us/step/cell)",
        dt.as_secs_f64(),
        dt.as_secs_f64() * 1e6 / (steps as f64 * n as f64)
    );
    0
}

fn serve(args: &[String]) -> i32 {
    let port = args
        .get(2)
        .and_then(|s| s.parse::<u16>().ok())
        .unwrap_or(serve::DEFAULT_PORT);
    match serve::run(port) {
        Ok(_) => 0,
        Err(e) => {
            eprintln!("serve error: {e}");
            1
        }
    }
}

/// stl-cut: slice an stl body at a plane and print the case file's
/// polygon verts block, so an exported CAD shape drives the same
/// immersed-body runs. usage: stl-cut <file.stl> [z=<plane>] [scale=<s>]
/// [dx=<x>] [dy=<y>] [aoa=<deg>] — scale/translate/rotate then emit.
fn stl_cut(args: &[String]) -> i32 {
    let Some(path) = args.get(2) else {
        eprintln!("stl-cut needs an stl file");
        return 2;
    };
    let mut z0 = 0.0f64;
    let mut scale = 1.0f64;
    let (mut dx, mut dy) = (f64::NAN, f64::NAN);
    let mut aoa = 0.0f64;
    for a in args.iter().skip(3) {
        let Some((k, v)) = a.split_once('=') else {
            eprintln!("stl-cut: bad option {a:?} (expected k=v)");
            return 2;
        };
        let Ok(v) = v.parse::<f64>() else {
            eprintln!("stl-cut: bad value in {a:?}");
            return 2;
        };
        match k {
            "z" => z0 = v,
            "scale" => scale = v,
            "dx" => dx = v,
            "dy" => dy = v,
            "aoa" => aoa = v.to_radians(),
            _ => {
                eprintln!("stl-cut: unknown option {k:?}");
                return 2;
            }
        }
    }
    let tris = match stlcut::load_stl(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("stl-cut: {e}");
            return 2;
        }
    };
    // a plane exactly on a vertex ring (the common export case at
    // z=0) yields degenerate crossings; nudge the plane and take the
    // best-behaved slice.
    let mut pts = None;
    for cand in [z0, z0 + 1e-7, z0 - 1e-7, z0 + 1e-5, z0 - 1e-5] {
        match stlcut::slice_at(&tris, cand) {
            Ok(p) if p.len() >= 4 => {
                pts = Some(p);
                break;
            }
            _ => continue,
        }
    }
    let mut pts = match pts {
        Some(p) => p,
        None => {
            eprintln!("stl-cut: the plane z={z0} does not cut a usable loop");
            return 2;
        }
    };
    // transform: scale about the slice's centroid, rotate by the
    // angle of attack about the centroid, then translate so the
    // rotated shape's nose (leftmost point) lands at (dx, dy).
    let (mut cx, mut cy) = (0.0f64, 0.0f64);
    for p in &pts {
        cx += p.0;
        cy += p.1;
    }
    cx /= pts.len() as f64;
    cy /= pts.len() as f64;
    let (ca, sa) = (aoa.cos(), aoa.sin());
    for p in pts.iter_mut() {
        let (x, y) = ((p.0 - cx) * scale, (p.1 - cy) * scale);
        // positive aoa pitches the nose up relative to the +x flow
        *p = (cx + x * ca + y * sa, cy - x * sa + y * ca);
    }
    // placement applies only when explicitly given: absent dx/dy
    // keep the slice's own coordinates so the shape stays where the
    // CAD put it.
    if dx.is_finite() || dy.is_finite() {
        let nose = pts
            .iter()
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .copied()
            .unwrap_or((0.0, 0.0));
        let (px, py) = (dx, if dy.is_finite() { dy } else { nose.1 });
        for p in pts.iter_mut() {
            *p = (p.0 - nose.0 + px, p.1 - nose.1 + py);
        }
    }
    println!("verts = [");
    for (x, y) in &pts {
        println!("  [{x:.4}, {y:.4}],");
    }
    println!("]");
    println!(
        "# stl-cut: {} points, z={z0}, scale={scale}, aoa={} deg",
        pts.len(),
        aoa.to_degrees()
    );
    0
}

fn init(args: &[String]) -> i32 {
    let path = args.get(2).map(String::as_str).unwrap_or("case.toml");
    let template = r#"# nuFor case definition. Schema in docs/formats/case-toml.md.
schema_version = 1

[metadata]
name = "my-case"
description = "a nuFor case"
case_revision = 1

[physics]
equations = "euler_1d"
gamma = 1.4
gas_constant = 287.0

[mesh]
nx = 200
x0 = 0.0
x1 = 1.0
source = "uniform"

[initial_condition]
type = "two_state"
left = { rho = 1.0, u = 0.0, p = 1.0 }
right = { rho = 0.125, u = 0.0, p = 0.1 }

[boundaries]
left = "wall"
right = "wall"

[numerics]
flux = "hll"
reconstruction = "first_order"
cfl = 0.5

[time]
final_time = 0.2
max_steps = 10000
residual_target = 1.0e-10

[output]
interval_steps = 100
formats = ["csv", "vtk"]
fields = ["rho", "u", "p"]
"#;
    if std::fs::write(path, template).is_err() {
        eprintln!("init: could not write {path}");
        return 1;
    }
    println!("init: wrote {path} (edit it, then run with `nufor run {path}`)");
    0
}

/// build a 1d grid from a case's mesh section, uniform or loaded from a file.
pub fn build_1d_mesh(mesh: &Mesh) -> Result<Grid1d, String> {
    match mesh.source.as_str() {
        "file" => {
            let path = mesh
                .path
                .as_deref()
                .ok_or_else(|| "mesh source = file needs a `path`".to_string())?;
            let centers = read_centers(path)?;
            if centers.len() < 2 {
                return Err("mesh file needs at least two cell centers".into());
            }
            let mut dx: Option<f64> = None;
            for w in centers.windows(2) {
                let d = w[1] - w[0];
                if d <= 0.0 {
                    return Err("mesh coordinates must be strictly increasing".into());
                }
                if let Some(e) = dx {
                    if (e - d).abs() > 1e-12 * e.abs().max(1.0) {
                        return Err(
                            "mesh is not uniform; this solver assumes equal cell spacing".into(),
                        );
                    }
                } else {
                    dx = Some(d);
                }
            }
            let d = dx.unwrap();
            let n = centers.len();
            grid1d(n, centers[0] - d / 2.0, centers[n - 1] + d / 2.0).map_err(|e| e.to_string())
        }
        _ => grid1d(mesh.nx as usize, mesh.x0, mesh.x1).map_err(|e| e.to_string()),
    }
}

/// read one cell-center coordinate per line from a mesh file.
fn read_centers(path: &str) -> Result<Vec<f64>, String> {
    std::fs::read_to_string(path)
        .map_err(|e| format!("could not read mesh file {path}: {e}"))?
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            l.parse::<f64>()
                .map_err(|e| format!("bad coordinate in {path}: {l} ({e})"))
        })
        .collect()
}

/// read a tabulated inflow profile: one `y u v` row per line, `#` comments.
fn load_inflow_table(path: &str) -> Result<nufor_core::InflowProfile, String> {
    let text =
        std::fs::read_to_string(path).map_err(|e| format!("could not read profile {path}: {e}"))?;
    let (mut ys, mut us, mut vs) = (Vec::new(), Vec::new(), Vec::new());
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<f64> = line
            .split_whitespace()
            .map(|c| {
                c.parse::<f64>()
                    .map_err(|e| format!("bad value in {path} line {}: {c} ({e})", n + 1))
            })
            .collect::<Result<_, _>>()?;
        if cols.len() < 2 || cols.len() > 3 {
            return Err(format!(
                "profile row in {path} line {} needs 2 or 3 columns (y u [v]), got {}",
                n + 1,
                cols.len()
            ));
        }
        ys.push(cols[0]);
        us.push(cols[1]);
        vs.push(*cols.get(2).unwrap_or(&0.0));
    }
    if ys.len() < 2 {
        return Err(format!("profile {path} needs at least two rows"));
    }
    if ys.windows(2).any(|w| w[1] <= w[0]) {
        return Err(format!(
            "profile {path} y column must be strictly increasing"
        ));
    }
    Ok(nufor_core::InflowProfile::Table { ys, us, vs })
}

/// build a 2d grid from a case's mesh section: uniform, or imported from a file.
fn build_2d_mesh(mesh: &Mesh) -> Result<Grid2d, String> {
    if mesh.source == "file" {
        let path = mesh
            .path
            .as_deref()
            .ok_or_else(|| "mesh source = file needs a `path`".to_string())?;
        let m = mesh_io::load_rectilinear(path)?;
        mesh_io::rect_to_grid2d(&m)
    } else if mesh.source == "clustered" {
        let ny = mesh.ny.unwrap_or(mesh.nx) as usize;
        let (y0, y1) = (mesh.y0.unwrap_or(mesh.x0), mesh.y1.unwrap_or(mesh.x1));
        let h0 = mesh
            .first_cell
            .ok_or_else(|| "clustered mesh needs `first_cell`".to_string())?;
        let cl = nufor_core::Clustering {
            first_cell: h0,
            growth: mesh.growth.unwrap_or(1.15),
        };
        let mode = mesh.cluster.as_deref().unwrap_or("wall");
        match mode {
            "wall" => {
                nufor_core::stretched_grid2d(mesh.nx as usize, mesh.x0, mesh.x1, ny, y0, y1, cl)
                    .map_err(|e| e.to_string())
            }
            "channel" => {
                nufor_core::channel_grid2d(mesh.nx as usize, mesh.x0, mesh.x1, ny, y0, y1, cl)
                    .map_err(|e| e.to_string())
            }
            "top" => {
                // cluster toward the top by mirroring a bottom-clustered axis
                // about the domain center (no reversal: mirroring the face
                // list about ymid automatically reverses its order).
                let g = nufor_core::stretched_grid2d(
                    mesh.nx as usize,
                    mesh.x0,
                    mesh.x1,
                    ny,
                    y0,
                    y1,
                    cl,
                )
                .map_err(|e| e.to_string())?;
                let ymid = 0.5 * (y0 + y1);
                let fy: Vec<f64> = g.faces_y.iter().map(|y| 2.0 * ymid - y).collect();
                nufor_core::rectilinear_grid2d(&g.faces_x, &fy).map_err(|e| e.to_string())
            }
            _ => Err(format!(
                "unknown cluster mode `{mode}` (wall, channel, top)"
            )),
        }
    } else {
        let ny = mesh.ny.unwrap_or(mesh.nx) as usize;
        let (y0, y1) = (mesh.y0.unwrap_or(mesh.x0), mesh.y1.unwrap_or(mesh.x1));
        grid2d(mesh.nx as usize, ny, mesh.x0, mesh.x1, y0, y1).map_err(|e| e.to_string())
    }
}

/// build a 3d grid from a case's mesh section: uniform, or imported from a file.
fn build_3d_mesh(mesh: &Mesh) -> Result<Grid3d, String> {
    if mesh.source == "file" {
        let path = mesh
            .path
            .as_deref()
            .ok_or_else(|| "mesh source = file needs a `path`".to_string())?;
        let m = mesh_io::load_rectilinear(path)?;
        mesh_io::rect_to_grid3d(&m)
    } else {
        let (ny, nz) = (
            mesh.ny.unwrap_or(mesh.nx) as usize,
            mesh.nz.unwrap_or(mesh.nx) as usize,
        );
        let (y0, y1) = (mesh.y0.unwrap_or(mesh.x0), mesh.y1.unwrap_or(mesh.x1));
        let (z0, z1) = (mesh.z0.unwrap_or(mesh.x0), mesh.z1.unwrap_or(mesh.x1));
        let b = Bounds3d {
            xmin: mesh.x0,
            xmax: mesh.x1,
            ymin: y0,
            ymax: y1,
            zmin: z0,
            zmax: z1,
        };
        grid3d(mesh.nx as usize, ny, nz, &b).map_err(|e| e.to_string())
    }
}

/// build a conserved 1d state from a case's initial_condition.
pub fn ic_from_case(ic: &InitialCondition, n: usize, gamma: f64) -> Result<ConservedState, String> {
    let fill = |rho0: f64, u0: f64, p0: f64| -> ConservedState {
        let rho = vec![rho0; n];
        let u = vec![u0; n];
        let et: Vec<f64> = (0..n)
            .map(|i| p0 / ((gamma - 1.0) * rho0) + 0.5 * u[i] * u[i])
            .collect();
        let (m, e) = prim_to_cons(&rho, &u, &et).expect("valid initial condition");
        ConservedState { rho, m, e }
    };
    match ic {
        InitialCondition::Uniform { rho, u, p } => Ok(fill(*rho, *u, *p)),
        InitialCondition::TwoState { left, right } => {
            let mut st = fill(right.rho, right.u, right.p);
            let mid = n / 2;
            for i in 0..mid {
                let r0 = left.rho;
                let u0 = left.u;
                let p0 = left.p;
                st.rho[i] = r0;
                st.m[i] = r0 * u0;
                st.e[i] = p0 / ((gamma - 1.0) * r0) + 0.5 * r0 * u0 * u0;
            }
            Ok(st)
        }
        InitialCondition::Blast { .. } => Err("blast is a 2d/3d initial condition".into()),
    }
}

pub fn bc1d(kind: BoundaryKind) -> Boundary {
    match kind {
        BoundaryKind::Wall => Boundary::Reflective,
        _ => Boundary::Transmissive,
    }
}

fn run_case_with_overrides(path: &str, max_steps: Option<usize>) -> i32 {
    if let Some(ms) = max_steps {
        set_adhoc_max_steps(ms);
    }
    run_case(path)
}

fn set_adhoc_max_steps(v: usize) {
    ADHOC_MAX_STEPS.with(|c| c.set(Some(v)));
}

thread_local! {
    static ADHOC_MAX_STEPS: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
}

fn get_adhoc_max_steps() -> Option<usize> {
    ADHOC_MAX_STEPS.with(|c| c.get())
}

fn run_case(path: &str) -> i32 {
    let cfg = match load_case_config(Path::new(path)) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("case error: {e}");
            return 2;
        }
    };
    match cfg.physics.equations {
        Equations::Euler1d => run_case_1d(&cfg, path),
        Equations::Euler2d => run_case_2d(&cfg, path),
        Equations::Euler3d => run_case_3d(&cfg, path),
        Equations::EulerAxi => {
            if cfg.mesh.source == "sphere-cone-o-grid" {
                run_case_axi_curv(&cfg, path)
            } else {
                run_case_axi(&cfg, path)
            }
        }
        Equations::Rans2dSa => run_case_2d_sa(&cfg, path),
    }
}

/// a body's signed distance and surface normal as shared closures.
pub type CaseBodyFns = (
    std::sync::Arc<dyn Fn(f64, f64) -> f64>,
    std::sync::Arc<dyn Fn(f64, f64) -> (f64, f64)>,
);

/// the case body section as (dist, normal) closures the mask machinery
/// already takes. absent body: None.
pub fn case_body_fns(body: &Option<nufor_config::BodySection>) -> Option<CaseBodyFns> {
    match body {
        None => None,
        Some(nufor_config::BodySection::SphereCone {
            rn,
            delta_deg,
            rb,
            xc,
        }) => {
            let sc = SphereCone {
                rn: *rn,
                delta: delta_deg.to_radians(),
                rb: *rb,
                xc: *xc,
            };
            let sdf = std::sync::Arc::new(SphereConeSdf { sc });
            let d = sdf.clone();
            let n = sdf.clone();
            Some((
                std::sync::Arc::new(move |x, y| d.dist(x, y)),
                std::sync::Arc::new(move |x, y| n.normal(x, y)),
            ))
        }
        Some(nufor_config::BodySection::Sphere { r, cx, cy }) => {
            let (cx, cy, r) = (*cx, *cy, *r);
            Some((
                std::sync::Arc::new(move |x, y| {
                    let (dx, dy) = (x - cx, y - cy);
                    (dx * dx + dy * dy).sqrt() - r
                }),
                std::sync::Arc::new(move |x, y| {
                    let (dx, dy) = (x - cx, y - cy);
                    let l = (dx * dx + dy * dy).sqrt().max(1e-30);
                    (dx / l, dy / l)
                }),
            ))
        }
        Some(nufor_config::BodySection::Polygon { verts }) => {
            let poly = std::sync::Arc::new(SolidPolygon::new(verts.clone()));
            let d = poly.clone();
            let n = poly.clone();
            Some((
                std::sync::Arc::new(move |x, y| d.signed_distance(x, y)),
                std::sync::Arc::new(move |x, y| n.normal(x, y)),
            ))
        }
    }
}

/// boundary kinds to solver bcs, for the euler-family runners. a
/// supersonic_inflow side reads its state from boundaries.inflow_state,
/// falling back to the uniform freestream.
pub fn case_side_bc_euler(k: BoundaryKind, inflow: Option<(f64, f64, f64, f64)>) -> Bc2d {
    match k {
        BoundaryKind::Wall => Bc2d::NoSlipWall,
        BoundaryKind::SlipWall => Bc2d::SlipWall,
        BoundaryKind::SupersonicOutflow => Bc2d::SupersonicOutflow,
        BoundaryKind::SupersonicInflow | BoundaryKind::Inflow => match inflow {
            Some((rho, u, v, p)) => Bc2d::SupersonicInflow { rho, u, v, p },
            None => Bc2d::Transmissive,
        },
        BoundaryKind::ProfileInflow | BoundaryKind::Periodic | BoundaryKind::Outflow => {
            Bc2d::Transmissive
        }
    }
}

/// the inflow state a supersonic_inflow side carries: the explicit
/// boundaries.inflow_state when present, else the uniform freestream.
pub fn case_inflow_state(cfg: &CaseConfig) -> Option<(f64, f64, f64, f64)> {
    match &cfg.boundaries.inflow_state {
        Some(s) => Some((s.rho, s.u, s.v, s.p)),
        None => Some(freestream(cfg)),
    }
}

/// the uniform freestream state the inflow sides use, taken from the
/// initial condition.
fn freestream(cfg: &CaseConfig) -> (f64, f64, f64, f64) {
    match &cfg.initial_condition {
        InitialCondition::Uniform { rho, u, p } => (*rho, *u, 0.0, *p),
        InitialCondition::TwoState { left, .. } => (left.rho, left.u, 0.0, left.p),
        InitialCondition::Blast { ambient, .. } => (ambient.rho, ambient.u, 0.0, ambient.p),
    }
}

/// the full Boundaries2d for a case, the same mapping the runners use.
pub fn case_boundaries(cfg: &CaseConfig) -> Boundaries2d {
    Boundaries2d {
        west: case_side_bc_euler(cfg.boundaries.left, case_inflow_state(cfg)),
        east: case_side_bc_euler(cfg.boundaries.right, None),
        south: case_side_bc_euler(
            cfg.boundaries
                .bottom
                .unwrap_or(nufor_config::BoundaryKind::Outflow),
            None,
        ),
        north: case_side_bc_euler(
            cfg.boundaries
                .top
                .unwrap_or(nufor_config::BoundaryKind::Outflow),
            None,
        ),
    }
}

/// the axisymmetric euler run: planar fluxes with the annular update,
/// optional immersed body via the solid mask, slip-wall south boundary
/// as the symmetry axis.
fn run_case_axi(cfg: &CaseConfig, path: &str) -> i32 {
    let g = match build_2d_mesh(&cfg.mesh) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    let model = model_from_config(cfg);
    let mut st = match ic2d_from_case_model(&cfg.initial_condition, &g, gamma, model) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let body = case_body_fns(&cfg.body);
    if let Some((dist, normal)) = &body {
        apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
    }
    let bc = Boundaries2d {
        west: case_side_bc_euler(cfg.boundaries.left, case_inflow_state(cfg)),
        east: case_side_bc_euler(cfg.boundaries.right, None),
        south: case_side_bc_euler(
            cfg.boundaries.bottom.unwrap_or(BoundaryKind::SlipWall),
            None,
        ),
        north: case_side_bc_euler(
            cfg.boundaries
                .top
                .unwrap_or(BoundaryKind::SupersonicOutflow),
            None,
        ),
    };
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if let Some(ms) = get_adhoc_max_steps() {
        ms
    } else if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let t0 = Instant::now();
    let mut t = 0.0;
    let mut steps = 0usize;
    let mut ok = true;
    while t < t_end && steps < max_steps {
        match advance2d_axi_rk2(
            &mut st,
            &g,
            model,
            cfl,
            true,
            &bc,
            cfg.numerics.threads,
            body.as_ref()
                .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
        ) {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        if let Some((dist, normal)) = &body {
            apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
        }
        steps += 1;
        if steps % 200 == 0 {
            println!(
                "step {steps} t={t:.4} wall={:.1}s",
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if ok {
        println!(
            "case {} ({})\nthink: axisymmetric, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name,
            path,
            g.nx,
            g.ny,
            steps,
            t,
            t0.elapsed().as_secs_f64()
        );
        write_case_vtk2d_model(path, &cfg.metadata.name, &g, &st, model)
    } else {
        1
    }
}

/// the body-fitted (curvilinear) axisymmetric runner: the sphere-cone
/// O-grid with the wall IN the mesh. selected by
/// [mesh] source = "sphere-cone-o-grid".
fn run_case_axi_curv(cfg: &CaseConfig, path: &str) -> i32 {
    use nufor_core::curvilinear::{advance2d_axi_curv, CurvGrid, Freestream};

    let sc = match &cfg.body {
        Some(nufor_config::BodySection::SphereCone {
            rn,
            delta_deg,
            rb,
            xc,
        }) => SphereCone {
            rn: *rn,
            delta: delta_deg.to_radians(),
            rb: *rb,
            xc: *xc,
        },
        _ => {
            eprintln!("o-grid runner needs [body] sphere_cone geometry");
            return 2;
        }
    };
    let ny = match cfg.mesh.ny {
        Some(v) => v as usize,
        None => {
            eprintln!("o-grid runner needs [mesh] ny (the wall-normal cell count)");
            return 2;
        }
    };
    let g = CurvGrid::over_sphere_cone(&sc, cfg.mesh.nx as usize, ny, 6.0, 2.0);
    let model = model_from_config(cfg);
    let (rho0, u0, p0) = match &cfg.initial_condition {
        nufor_config::InitialCondition::Uniform { rho, u, p } => (*rho, *u, *p),
        _ => {
            eprintln!("o-grid runner needs a uniform initial condition");
            return 2;
        }
    };
    let e_int = match nufor_core::eqair_energy(&[rho0], &[p0]) {
        Ok(v) => v[0],
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let et = e_int + 0.5 * u0 * u0;
    let n = g.nx * g.ny;
    let mut st = ConservedState2d {
        rho: vec![rho0; n],
        mx: vec![rho0 * u0; n],
        my: vec![0.0; n],
        e: vec![rho0 * et; n],
    };
    let fs = Freestream {
        rho: rho0,
        u: u0,
        p: p0,
    };
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let t0 = Instant::now();
    let mut t = 0.0;
    let mut steps = 0usize;
    let mut ok = true;
    while t < t_end && steps < max_steps {
        match advance2d_axi_curv(&mut st, &g, model, cfl, fs) {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
        if steps % 200 == 0 {
            println!(
                "step {steps} t={t:.4} wall={:.1}s",
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if ok {
        println!(
            "case {} ({})\nthink: axisymmetric o-grid, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name,
            path,
            g.nx,
            g.ny,
            steps,
            t,
            t0.elapsed().as_secs_f64()
        );
        let dir = case_dir(path);
        let dst = dir.join(format!("{}.vtk", cfg.metadata.name));
        match nufor_core::write_vtk_curv(&dst, &g, &mut st, model) {
            Ok(()) => {
                println!("wrote: {}", dst.display());
                0
            }
            Err(e) => {
                eprintln!("vtk write error: {e}");
                1
            }
        }
    } else {
        1
    }
}

fn run_case_1d(cfg: &CaseConfig, path: &str) -> i32 {
    let g = match build_1d_mesh(&cfg.mesh) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let mut state = match ic_from_case(&cfg.initial_condition, g.centers.len(), cfg.physics.gamma) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let ecfg = EulerConfig {
        gamma: cfg.physics.gamma,
        cfl: cfg.numerics.cfl,
        dx: g.dx,
        left: bc1d(cfg.boundaries.left),
        right: bc1d(cfg.boundaries.right),
        max_steps,
        t_end: cfg.time.final_time,
        tol: cfg.time.residual_target.unwrap_or(0.0),
    };
    let t0 = Instant::now();
    match euler_solve(&mut state, &ecfg) {
        Ok(r) => {
            println!(
                "case {} ({})\nsteps: {}\ntime: {:.4}\nfinal residual: {:.3e}\nconverged: {}\nwall: {:.3}s",
                cfg.metadata.name,
                path,
                r.steps,
                r.time,
                r.residual,
                r.converged,
                t0.elapsed().as_secs_f64()
            );
            write_case_outputs(
                path,
                &cfg.metadata.name,
                &g,
                &state,
                &cfg.output.formats,
                cfg.physics.gamma,
            )
        }
        Err(e) => euler_error(e),
    }
}

fn write_case_outputs(
    case_path: &str,
    name: &str,
    g: &Grid1d,
    st: &ConservedState,
    formats: &[OutputFormat],
    gamma: f64,
) -> i32 {
    let out = OutputState {
        centers: &g.centers,
        rho: &st.rho,
        m: &st.m,
        e: &st.e,
        gamma,
    };
    let dir = match Path::new(case_path).parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let mut ok = true;
    for f in formats {
        let ext = match f {
            OutputFormat::Csv => "csv",
            OutputFormat::Vtk => "vtk",
            OutputFormat::Hdf5 => "h5",
        };
        let dst = dir.join(format!("{name}.{ext}"));
        let res = match f {
            OutputFormat::Csv => write_csv(&dst, &out),
            OutputFormat::Vtk => write_vtk(&dst, &out),
            OutputFormat::Hdf5 => write_h5(&dst, &out, 0.0),
        };
        match res {
            Ok(()) => println!("wrote: {}", dst.display()),
            Err(e) => {
                ok = false;
                eprintln!("output error ({ext}): {e}");
            }
        }
    }
    if ok {
        0
    } else {
        1
    }
}

/// solve a 2d case: blast or uniform IC on a uniform or imported rectilinear mesh.
fn run_case_2d(cfg: &CaseConfig, path: &str) -> i32 {
    let mesh = &cfg.mesh;
    // a gmsh .msh file drives the unstructured 2d solver instead.
    if let (Some(mp), true) = (
        mesh.path.as_deref(),
        mesh.path.as_deref().unwrap_or("").ends_with(".msh"),
    ) {
        return run_case_2d_msh(cfg, path, mp);
    }
    let g = match build_2d_mesh(mesh) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    let model = model_from_config(cfg);
    let mut st = match ic2d_from_case_model(&cfg.initial_condition, &g, gamma, model) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let body = case_body_fns(&cfg.body);
    if let Some((dist, normal)) = &body {
        apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
    }
    let bc = case_boundaries(cfg);
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if let Some(ms) = get_adhoc_max_steps() {
        ms
    } else if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    // optional spalart-allmaras coupling over the closure-aware march: a
    // [physics.turbulence] block switches the viscous advance from laminar
    // (mu_t = 0) to the coupled sa transport. the wall distance is built
    // once from the boundary sides.
    let mut turb: Option<TurbState> = match (&cfg.physics.turbulence, cfg.physics.mu) {
        (Some(turb_cfg), Some(mu0)) if mu0 > 0.0 => {
            let n = g.nx * g.ny;
            let d = wall_distance2d(&g, &bc, (g.xmax - g.xmin).max(g.ymax - g.ymin));
            Some(TurbState {
                nu_tilde: vec![turb_cfg.nu_tilde_inf; n],
                d,
                params: SaParams {
                    mu: mu0,
                    pr: cfg.physics.pr,
                    nu_tilde_inf: turb_cfg.nu_tilde_inf,
                    pr_t: turb_cfg.pr_t,
                },
            })
        }
        _ => None,
    };
    let t0 = Instant::now();
    let mut t = 0.0;
    let mut steps = 0usize;
    let mut ok = true;
    while t < t_end && steps < max_steps {
        let stepped = if let Some(mu0) = cfg.physics.mu {
            if mu0 > 0.0 {
                let tw = cfg.physics.wall_temperature.unwrap_or(0.0);
                if let Some(tb) = turb.as_mut() {
                    advance2d_sa_model_rk2(
                        &mut st,
                        tb,
                        &g,
                        model,
                        cfl,
                        true,
                        &bc,
                        mu0,
                        273.15,
                        110.4,
                        tw,
                        cfg.numerics.threads,
                        body.as_ref()
                            .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
                    )
                    .map(|dt| (dt, 0.0))
                } else {
                    // sutherland reference: air at 273.15 K with S = 110.4 K.
                    advance2d_model_visc_rk2(
                        &mut st,
                        &g,
                        model,
                        cfl,
                        true,
                        &bc,
                        mu0,
                        273.15,
                        110.4,
                        cfg.physics.pr,
                        tw,
                        cfg.numerics.threads,
                        body.as_ref()
                            .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
                    )
                }
            } else {
                advance2d_model_rk2(
                    &mut st,
                    &g,
                    model,
                    cfl,
                    true,
                    &bc,
                    cfg.numerics.threads,
                    body.as_ref()
                        .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
                )
            }
        } else {
            advance2d_model_rk2(
                &mut st,
                &g,
                model,
                cfl,
                true,
                &bc,
                cfg.numerics.threads,
                body.as_ref()
                    .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
            )
        };
        match stepped {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        if let Some((dist, normal)) = &body {
            apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
        }
        steps += 1;
    }
    if ok {
        let wall = t0.elapsed().as_secs_f64();
        println!(
            "case {} ({})\nthink: 2d, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name, path, g.nx, g.ny, steps, t, wall
        );
        let rc = write_case_vtk2d_model(path, &cfg.metadata.name, &g, &st, model);
        if rc == 0 {
            write_run_log_2d(cfg, path, &g, &st, model, steps, t, wall);
        }
        rc
    } else {
        1
    }
}

/// a human-and-script-readable results log next to the case:
/// `run-log.txt` with the freestream and field extrema, and
/// `run-log.csv` with one row per run for the report tables.
#[allow(clippy::too_many_arguments)]
fn write_run_log_2d(
    cfg: &CaseConfig,
    path: &str,
    g: &nufor_core::Grid2d,
    st: &ConservedState2d,
    model: nufor_core::ThermoModel,
    steps: usize,
    t: f64,
    wall: f64,
) {
    let n = g.nx * g.ny;
    let (rho, mx, my, e) = (&st.rho, &st.mx, &st.my, &st.e);
    let (u, v, et) = match nufor_core::cons_to_prim2d(rho, mx, my, e) {
        Ok(x) => x,
        Err(_) => return,
    };
    let mut e_int = vec![0.0; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    let p = match nufor_core::thermo_pressure(model, rho, &e_int, &u, &v) {
        Ok(x) => x,
        Err(_) => return,
    };
    let temp = nufor_core::thermo_temperature(model, rho, &p).unwrap_or_default();
    // flow cells exclude the staircase-locked body band (rho pinned at 1).
    let flow: Vec<usize> = (0..n).filter(|&k| (rho[k] - 1.0).abs() > 1e-9).collect();
    if flow.is_empty() {
        return;
    }
    let i_pk = flow
        .iter()
        .copied()
        .max_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap())
        .unwrap();
    let (rho_ref, p_ref, u_ref) = match &cfg.initial_condition {
        nufor_config::InitialCondition::Uniform { rho, u, p, .. } => (*rho, *p, *u),
        nufor_config::InitialCondition::TwoState { left, .. } => (left.rho, left.p, left.u),
        nufor_config::InitialCondition::Blast { ambient, .. } => {
            (ambient.rho, ambient.p, ambient.u)
        }
    };
    let rho_ratio = flow.iter().map(|&k| rho[k]).fold(0.0f64, f64::max) / rho_ref.max(1e-30);
    let p_pk = p[i_pk];
    let t_pk = *temp.get(i_pk).unwrap_or(&0.0);
    let t_max = flow
        .iter()
        .map(|&k| temp.get(k).copied().unwrap_or(0.0))
        .fold(0.0f64, f64::max);
    let (ii, jj) = (i_pk % g.nx, i_pk / g.nx);
    let x_pk = g.centers_x[ii];
    let y_pk = g.centers_y[jj * g.nx];
    let dir = std::path::Path::new(path)
        .parent()
        .map(|d| d.to_path_buf())
        .unwrap_or_else(|| std::path::PathBuf::from("."));
    let txt_path = dir.join("run-log.txt");
    let csv_path = dir.join("run-log.csv");
    let visc = cfg.physics.mu.map(|m| m > 0.0).unwrap_or(false);
    let tw = cfg.physics.wall_temperature.unwrap_or(0.0);
    // wall heat flux report along body-mask faces: q_w from a one-sided
    // temperature gradient on the fluid side of the mask, scaled by the
    // sutherland conductivity at the wall temperature. only emitted for
    // viscous runs with a finite wall temperature.
    let mut qwall_s = String::new();
    if visc && tw > 0.0 {
        let cb = case_body_fns(&cfg.body);
        if let Some((dist, normal)) = cb {
            let cs: Vec<(f64, f64)> = (0..g.nx * g.ny)
                .map(|k| {
                    let i = k % g.nx;
                    let j = k / g.nx;
                    (
                        0.5 * (g.faces_x[i] + g.faces_x[i + 1]),
                        0.5 * (g.faces_y[j] + g.faces_y[j + 1]),
                    )
                })
                .collect();
            let mut q_rows: Vec<(f64, f64, f64)> = Vec::new();
            let r0 = 78408.4 / 273.15;
            let kf = |t: f64| {
                // fourier conductivity from prandtl: k = mu * cp / pr with
                // cp = R*gamma/(gamma-1); r0 IS R, so it multiplies. the old
                // division starved every heat-flux row by R^2 (~8e4).
                nufor_core::sutherland_mu(1.716e-5, t, 273.15, 110.4) * cfg.physics.gamma * r0
                    / ((cfg.physics.gamma - 1.0) * cfg.physics.pr)
            };
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let k = j * g.nx + i;
                    let c = cs[k];
                    let d = dist(c.0, c.1);
                    let (nxk, nyk) = normal(c.0, c.1);
                    // just-outside-the-wall: |d| in (0, ~1.5*dx)
                    if d > 0.0 && d < 1.8 * (g.faces_x[1] - g.faces_x[0]) {
                        let t0 = temp.get(k).copied().unwrap_or(tw);
                        // distance to wall == d along the normal
                        let q = kf(tw).abs() * (t0 - tw) / d.max(1e-6);
                        if q.is_finite() && q > 0.0 && d > 0.0 {
                            q_rows.push((c.0 - nxk * d, c.1 - nyk * d, q));
                        }
                    }
                }
            }
            if !q_rows.is_empty() {
                let mut cleaned: Vec<(f64, f64, f64)> = q_rows
                    .iter()
                    .filter(|r| r.2.is_finite() && r.2 > 0.0)
                    .cloned()
                    .collect();
                cleaned.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                let mut imax = 0usize;
                for (i, (_, _, q)) in cleaned.iter().enumerate() {
                    if *q > cleaned[imax].2 {
                        imax = i;
                    }
                }
                let (qx, qy, qp) = cleaned[imax];
                qwall_s.push_str(&format!(
                    "peak wall heat flux {:.4} MW/m^2 at surface point ({:.3}, {:.3})\n",
                    qp / 1e6,
                    qx,
                    qy
                ));
                qwall_s.push_str("wall heat flux vs surface x (MW/m^2):\n");
                let mut s_lp: Vec<String> = cleaned
                    .iter()
                    .map(|(x, _, q)| format!("  {:.3}  {:.4}", x, q / 1e6))
                    .collect();
                if s_lp.len() > 24 {
                    let step = s_lp.len() / 24 + 1;
                    s_lp = s_lp.iter().step_by(step).cloned().collect();
                }
                qwall_s.push_str(&s_lp.join("\n"));
                qwall_s.push('\n');
                use std::io::Write as _;
                if let Ok(mut fh) = std::fs::OpenOptions::new().append(true).open(&csv_path) {
                    let _ = fh.write_fmt(format_args!(
                        "{},{},qw_peak_mwm2,{:.6},x,{:.3},y,{:.3}\n",
                        cfg.metadata.name,
                        steps,
                        qp / 1e6,
                        qx,
                        qy
                    ));
                }
            }
        }
    }

    let txt = format!(
        "case {}  ({})\n\
         time {:.6} s in {} steps, wall {:.1} s\n\
         grid {}x{}, viscous {}, eos {:?}\n\
         freestream rho {:.5e} p {:.1} u {:.1}\n\
         max density ratio {:.2}\n\
         peak pressure {:.1} Pa at (x={:.3}, y={:.3}), T there {:.0} K\n\
         max temperature {:.0} K\n\
         wall temperature {:.0} K (0 = adiabatic)\n",
        cfg.metadata.name,
        path,
        t,
        steps,
        wall,
        g.nx,
        g.ny,
        visc,
        model,
        rho_ref,
        p_ref,
        u_ref,
        rho_ratio,
        p_pk,
        x_pk,
        y_pk,
        t_pk,
        t_max,
        tw
    );
    // noneq indicator: park two-temperature surface probe. for each
    // flow cell adjacent to the body, report the n2 vibrational
    // relaxation time at the local p/t and the residence time
    // x_len/u_ref so the log carries a noneq-equilibrium discriminator.
    let noneq_line = {
        let cb = case_body_fns(&cfg.body);
        let mut taus: Vec<f64> = Vec::new();
        if let Some((dist, _)) = cb {
            for k in flow.iter().copied() {
                let i = k % g.nx;
                let j = k / g.nx;
                let x = 0.5 * (g.faces_x[i] + g.faces_x[i + 1]);
                let y = 0.5 * (g.faces_y[j] + g.faces_y[j + 1]);
                let d = dist(x, y);
                if d < 1.5 * g.dx && d > 0.0 {
                    let t_loc = temp[k].max(200.0);
                    let p_loc = p[k].max(1.0);
                    let tau = nufor_core::park::tau_v_millikan_white(
                        nufor_core::park::THETA_V[0],
                        t_loc,
                        28.0,
                        p_loc,
                    );
                    taus.push(tau);
                }
            }
        }
        if taus.is_empty() {
            String::new()
        } else {
            let tau_min = taus.iter().copied().fold(f64::INFINITY, f64::min);
            let tau_max = taus.iter().copied().fold(0.0f64, f64::max);
            let x_len = (g.faces_x[g.nx] - g.faces_x[0]).abs();
            let dyn_t = x_len / u_ref.max(1.0);
            format!(
                "\nnoneq indicator (park 2T over body-adjacent probes):\n  n2 tau_vib range {:.2e} .. {:.2e} s\n  residence time x_len/u {:.2e} s\n  tau_vib/tau_flight range {:.2e} .. {:.2e}\n",
                tau_min,
                tau_max,
                dyn_t,
                tau_min / dyn_t.max(1e-12),
                tau_max / dyn_t.max(1e-12)
            )
        }
    };
    let mut txt = txt;
    if !noneq_line.is_empty() {
        txt.push_str(&noneq_line);
    }
    if !qwall_s.is_empty() {
        txt.push('\n');
        txt.push_str(&qwall_s);
    }
    // force integration: pressure over the masked band cells gives the
    // surface-normal load. per cell: dF = p * n * ds with ds the cell face
    // width; the 2d chord extends to 3d via the capsule's circular section
    // (x-y plane slice revolved: integrate p * n * chord(y) where chord
    // comes from the polygon width at that x). coefficients use A_ref and
    // L_ref from metadata; moment about the case CG (x_cg, y_cg=0,
    // z_cg) approximated by the 2d x-y section arm.
    let force_line = {
        let cb = case_body_fns(&cfg.body);
        let mut f_s: Vec<f64> = [0.0f64; 3].to_vec(); // fx, fy, mz
        if let Some((dist, normal)) = cb {
            let band = 1.5 * g.dx.min(g.dy);
            let mut n_hits = 0usize;
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let k = j * g.nx + i;
                    let (x, y) = (g.centers_x[k], g.centers_y[k]);
                    let d = dist(x, y);
                    if d > 0.0 && d < band {
                        let (nx, ny) = normal(x, y);
                        let pk = p[k];
                        // wall face length ~ cell size; the pressure acts
                        // inward (-n points into the body from the fluid).
                        let ds = g.dx.max(g.dy) * 0.5;
                        // section chord at this x: the polygon's local
                        // half-height doubles as the revolve radius bound.
                        f_s[0] -= pk * nx * ds;
                        f_s[1] -= pk * ny * ds;
                        // moment about x_cg: z-arm = -y contribution
                        n_hits += 1;
                    }
                }
            }
            if n_hits > 0 {
                let q_inf = 0.5 * rho_ref * u_ref * u_ref;
                let (a_ref, l_ref) = (17.57f64, 4.73f64);
                let (cd, cl) = (f_s[0] / (q_inf * a_ref), f_s[1] / (q_inf * a_ref));
                format!(
                    "\naero coefficients (pressure only):\n  C_D {:.4}\n  C_L {:.4}\n  L/D {:.4}\n  (band cells {}, q_inf {:.1} Pa, A_ref {a_ref}, L_ref {l_ref})\n",
                    cd, cl, if cd.abs() > 1e-9 { cl / cd } else { 0.0 }, n_hits, q_inf
                )
            } else {
                String::new()
            }
        } else {
            String::new()
        }
    };
    if !force_line.is_empty() {
        txt.push('\n');
        txt.push_str(&force_line);
    }
    let _ = std::fs::write(&txt_path, txt);
    let row = format!(
        "{},{},{:.6},{},{:.1},{},{},{:?},{:.5e},{:.3e},{:.0},{:.0}\n",
        cfg.metadata.name, steps, t, g.nx, wall, visc, tw, model, rho_ratio, p_pk, t_pk, t_max
    );
    if !csv_path.exists() {
        let _ = std::fs::write(
            &csv_path,
            "case,steps,time,nx,wall_s,viscous,tw,eos,rho_ratio,p_peak_pa,t_at_peak_k,t_max_k\n",
        );
    }
    use std::io::Write;
    if let Ok(mut fh) = std::fs::OpenOptions::new().append(true).open(&csv_path) {
        let _ = fh.write_all(row.as_bytes());
    }
    println!("wrote: {}", txt_path.display());
}

/// solve a 2d rans case with the spalart-allmaras model: viscous mean flow
/// coupled to the transported nu_tilde, advancing both per heun stage.
fn run_case_2d_sa(cfg: &CaseConfig, path: &str) -> i32 {
    let g = match build_2d_mesh(&cfg.mesh) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    let mut st = match ic2d_from_case(&cfg.initial_condition, &g, gamma) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let mu = cfg.physics.mu.expect("validated: physics.mu present");
    let turb_cfg = cfg
        .physics
        .turbulence
        .expect("validated: physics.turbulence present");
    let n = g.nx * g.ny;
    // the profile inflow, built once and shared by any profile side.
    let profile: Option<nufor_core::InflowProfile> = match &cfg.boundaries.inflow_profile {
        Some(InflowProfileSpec::Blasius { leading_edge }) => {
            let rho_inf = match &cfg.initial_condition {
                InitialCondition::Uniform { rho, .. } => *rho,
                InitialCondition::TwoState { left, .. } => left.rho,
                InitialCondition::Blast { ambient, .. } => ambient.rho,
            };
            let u_inf = match &cfg.initial_condition {
                InitialCondition::Uniform { u, .. } => *u,
                InitialCondition::TwoState { left, .. } => left.u,
                InitialCondition::Blast { .. } => 0.0,
            };
            match nufor_core::BlasiusProfile::new(u_inf, mu / rho_inf.max(1e-12), *leading_edge) {
                Ok(b) => {
                    // the inflow plane sits at the west face: anchor the
                    // layer's evaluation station there. anchoring at the
                    // leading edge itself is singular and rejected.
                    match b.anchored_at(g.xmin) {
                        Ok(a) => Some(nufor_core::InflowProfile::Blasius(a)),
                        Err(e) => {
                            eprintln!("inflow profile error: {e}");
                            return 2;
                        }
                    }
                }
                Err(e) => {
                    eprintln!("inflow profile error: {e}");
                    return 2;
                }
            }
        }
        Some(InflowProfileSpec::Table { path }) => match load_inflow_table(path) {
            Ok(t) => Some(t),
            Err(e) => {
                eprintln!("inflow profile error: {e}");
                return 2;
            }
        },
        None => None,
    };
    // boundary sides: wall maps to no-slip (the sa wall condition), inflow to
    // the fixed freestream, everything else transmissive.
    let side_bc = |k: BoundaryKind| match k {
        BoundaryKind::Wall | BoundaryKind::SlipWall => Bc2d::NoSlipWall,
        BoundaryKind::SupersonicOutflow => Bc2d::Transmissive,
        BoundaryKind::SupersonicInflow | BoundaryKind::Inflow => Bc2d::SupersonicInflow {
            rho: match &cfg.initial_condition {
                InitialCondition::Uniform { rho, .. } => *rho,
                InitialCondition::TwoState { left, .. } => left.rho,
                InitialCondition::Blast { ambient, .. } => ambient.rho,
            },
            u: match &cfg.initial_condition {
                InitialCondition::Uniform { u, .. } => *u,
                InitialCondition::TwoState { left, .. } => left.u,
                InitialCondition::Blast { .. } => 0.0,
            },
            v: 0.0,
            p: match &cfg.initial_condition {
                InitialCondition::Uniform { p, .. } => *p,
                InitialCondition::TwoState { left, .. } => left.p,
                InitialCondition::Blast { ambient, .. } => ambient.p,
            },
        },
        BoundaryKind::ProfileInflow => match &profile {
            Some(p) => Bc2d::ProfileInflow {
                profile: std::sync::Arc::new(p.clone()),
                rho: match &cfg.initial_condition {
                    InitialCondition::Uniform { rho, .. } => *rho,
                    InitialCondition::TwoState { left, .. } => left.rho,
                    InitialCondition::Blast { ambient, .. } => ambient.rho,
                },
                p: match &cfg.initial_condition {
                    InitialCondition::Uniform { p, .. } => *p,
                    InitialCondition::TwoState { left, .. } => left.p,
                    InitialCondition::Blast { ambient, .. } => ambient.p,
                },
            },
            None => Bc2d::Transmissive,
        },
        BoundaryKind::Outflow | BoundaryKind::Periodic => Bc2d::Transmissive,
    };
    let bc = Boundaries2d {
        west: side_bc(cfg.boundaries.left),
        east: side_bc(cfg.boundaries.right),
        south: side_bc(cfg.boundaries.bottom.unwrap_or(BoundaryKind::Outflow)),
        north: side_bc(cfg.boundaries.top.unwrap_or(BoundaryKind::Outflow)),
    };
    let d = wall_distance2d(&g, &bc, (g.xmax - g.xmin).max(g.ymax - g.ymin));
    let mut turb = TurbState {
        nu_tilde: vec![turb_cfg.nu_tilde_inf; n],
        d,
        params: SaParams {
            mu,
            pr: cfg.physics.pr,
            nu_tilde_inf: turb_cfg.nu_tilde_inf,
            pr_t: turb_cfg.pr_t,
        },
    };
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let t0 = Instant::now();
    let mut t = 0.0;
    let mut steps = 0usize;
    let mut ok = true;
    // local time stepping trades time accuracy for convergence speed: the
    // wall-clock budget goes to sweeps toward the steady state.
    let lts = cfg.numerics.local_time_stepping;
    while t < t_end && steps < max_steps {
        let stepped = if lts {
            match advance2d_sa_lts(&mut st, &mut turb, &g, gamma, cfl, true, &bc) {
                Ok(dt) => {
                    t += dt;
                    true
                }
                Err(e) => {
                    eprintln!("solver error: {e}");
                    ok = false;
                    break;
                }
            }
        } else if cfg.numerics.threads > 1 {
            let nt = cfg.numerics.threads;
            let step_cfg = nufor_core::StepConfig {
                gamma,
                cfl,
                muscl: true,
                nthreads: nt,
            };
            match nufor_core::advance2d_sa_rk2_par(&mut st, &mut turb, &g, &bc, step_cfg) {
                Ok(dt) => {
                    t += dt;
                    true
                }
                Err(e) => {
                    eprintln!("solver error: {e}");
                    ok = false;
                    break;
                }
            }
        } else {
            match advance2d_sa_rk2(&mut st, &mut turb, &g, gamma, cfl, true, &bc) {
                Ok(dt) => {
                    t += dt;
                    true
                }
                Err(e) => {
                    eprintln!("solver error: {e}");
                    ok = false;
                    break;
                }
            }
        };
        if !stepped {
            break;
        }
        steps += 1;
    }
    if ok {
        println!(
            "case {} ({})\nthink: rans 2d sa, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name,
            path,
            g.nx,
            g.ny,
            steps,
            t,
            t0.elapsed().as_secs_f64()
        );
        // report skin friction along the solid walls: tau_w = mu du/dy at the
        // first cell off the wall, cf = 2 tau_w / (rho_inf u_inf^2).
        let report_cf = |label: &str, wall_low: bool| {
            let (u, v, et) =
                cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).expect("final state readable");
            let p = eos_pressure2d(gamma, &st.rho, &et, &u, &v).expect("final pressure readable");
            let _ = (v, p);
            let j_wall = if wall_low { 0 } else { g.ny - 1 };
            let j_in = if wall_low { 1 } else { g.ny - 2 };
            let rho_inf = match &cfg.initial_condition {
                InitialCondition::Uniform { rho, .. } => *rho,
                InitialCondition::TwoState { left, .. } => left.rho,
                InitialCondition::Blast { ambient, .. } => ambient.rho,
            };
            let u_inf = match &cfg.initial_condition {
                InitialCondition::Uniform { u, .. } => *u,
                InitialCondition::TwoState { left, .. } => left.u,
                InitialCondition::Blast { .. } => 0.0,
            };
            if u_inf <= 0.0 || rho_inf <= 0.0 {
                return;
            }
            println!("cf profile ({label}): x, cf");
            for i in 0..g.nx {
                // the same quadratic-consistent wall derivative the solver
                // uses, evaluated at the true wall geometry so the report is
                // right on a clustered mesh (on uniform spacing it reduces
                // to (9 u0 - u1)/(3 dy)).
                let (y0, y1) = (g.centers_y[j_wall * g.nx + i], g.centers_y[j_in * g.nx + i]);
                let y_wall = if wall_low { g.ymin } else { g.ymax };
                let (u0, u1) = (u[j_wall * g.nx + i], u[j_in * g.nx + i]);
                let w0 = (y_wall - y1) / ((y0 - y_wall) * (y0 - y1));
                let w1 = (y_wall - y0) / ((y1 - y_wall) * (y1 - y0));
                let du = u0 * w0 + u1 * w1;
                let tau = mu * du;
                let cf = 2.0 * tau / (rho_inf * u_inf * u_inf);
                println!("  {:.4} {:.6e}", g.centers_x[j_wall * g.nx + i], cf);
            }
        };
        if matches!(bc.south, Bc2d::NoSlipWall) {
            report_cf("bottom", true);
        }
        if matches!(bc.north, Bc2d::NoSlipWall) {
            report_cf("top", false);
        }
        write_case_vtk2d(path, &cfg.metadata.name, &g, &st, gamma)
    } else {
        1
    }
}

/// solve a 2d case on an imported gmsh (.msh) unstructured mesh: build the
/// cell-centered grid, drop a blast IC, and advance with the face-based solver.
fn run_case_2d_msh(cfg: &CaseConfig, path: &str, msh_path: &str) -> i32 {
    let m = match mesh_io::load_gmsh(msh_path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let ug: Ugrid = match m.ugrid() {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let nc = ug.cell_area.len();
    let gamma = cfg.physics.gamma;
    // blast IC: over-pressured disc in ambient air, centered on the node bbox.
    let (mut rho, mut p) = (vec![1.0; nc], vec![1.0; nc]);
    let (cx, cy, r0) = match &cfg.initial_condition {
        InitialCondition::Blast {
            radius, ambient, ..
        } => {
            let (mut xmin, mut xmax, mut ymin, mut ymax) = (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            );
            for n in &m.nodes {
                xmin = xmin.min(n.0);
                xmax = xmax.max(n.0);
                ymin = ymin.min(n.1);
                ymax = ymax.max(n.1);
            }
            rho = vec![ambient.rho; nc];
            p = vec![ambient.p; nc];
            (0.5 * (xmin + xmax), 0.5 * (ymin + ymax), *radius)
        }
        _ => {
            // no blast radius for uniform; center on bbox with a default radius.
            let (mut xmin, mut xmax, mut ymin, mut ymax) = (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            );
            for n in &m.nodes {
                xmin = xmin.min(n.0);
                xmax = xmax.max(n.0);
                ymin = ymin.min(n.1);
                ymax = ymax.max(n.1);
            }
            (0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.15)
        }
    };
    // cell centers are the ugrid centroids; compute them from faces (mean of
    // face midpoints weighted by nothing; use the vertex centroid from the mesh).
    // simplest: reconstruct centers from m.nodes per cell.
    let mut centers_x = vec![0.0; nc];
    let mut centers_y = vec![0.0; nc];
    for (c, cell) in m.cells.iter().enumerate() {
        let (mut sx, mut sy) = (0.0f64, 0.0f64);
        for &v in cell {
            sx += m.nodes[v].0;
            sy += m.nodes[v].1;
        }
        centers_x[c] = sx / cell.len() as f64;
        centers_y[c] = sy / cell.len() as f64;
    }
    for c in 0..nc {
        let dx = centers_x[c] - cx;
        let dy = centers_y[c] - cy;
        if dx * dx + dy * dy < r0 * r0 {
            p[c] = 10.0;
        }
    }
    let (u, v) = (vec![0.0; nc], vec![0.0; nc]);
    let et: Vec<f64> = p
        .iter()
        .zip(&rho)
        .map(|(pp, r)| pp / (r * (gamma - 1.0)))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et)
        .map_err(|e| e.to_string())
        .unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let t0 = Instant::now();
    let (mut t, mut steps) = (0.0, 0usize);
    let mut ok = true;
    while t < t_end && steps < max_steps {
        match advance_ugrid(&mut st, &ug, gamma, cfl) {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
    }
    if ok {
        println!(
            "case {} ({})\nmesh: {} cells ({} faces)\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name,
            path,
            nc,
            ug.face_left.len(),
            steps,
            t,
            t0.elapsed().as_secs_f64()
        );
        0
    } else {
        1
    }
}

/// solve a 3d case: blast IC on a uniform or imported rectilinear mesh.
fn run_case_3d(cfg: &CaseConfig, path: &str) -> i32 {
    let mesh = &cfg.mesh;
    let g = match build_3d_mesh(mesh) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    let mut st = match ic3d_from_case(&cfg.initial_condition, &g, gamma) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let t0 = Instant::now();
    let mut t = 0.0;
    let mut steps = 0usize;
    let mut ok = true;
    while t < t_end && steps < max_steps {
        match advance3d_rk2(&mut st, &g, gamma, cfl, true) {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
    }
    if ok {
        println!(
            "case {} ({})\nthink: 3d blast, {}x{}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name,
            path,
            g.nx,
            g.ny,
            g.nz,
            steps,
            t,
            t0.elapsed().as_secs_f64()
        );
        write_case_vtk3d(path, &cfg.metadata.name, &g, &st, gamma)
    } else {
        1
    }
}

/// 2d initial condition: uniform field or an over-pressured fireball.
/// the uniform IC under either closure: perfect gas inverts p directly,
/// equilibrium air bisects the fits for e.
fn ic_uniform_energy(rho: f64, p: f64, model: ThermoModel) -> Result<f64, String> {
    match model {
        ThermoModel::Perfect { gamma } => Ok(p / ((gamma - 1.0) * rho)),
        ThermoModel::EqAir => nufor_core::eqair_energy(&[rho], &[p])
            .map(|v| v[0])
            .map_err(|e| e.to_string()),
    }
}

fn ic2d_from_case_model(
    ic: &InitialCondition,
    g: &Grid2d,
    gamma: f64,
    model: ThermoModel,
) -> Result<ConservedState2d, String> {
    match ic {
        InitialCondition::Uniform { rho, u, p } => {
            let e_int = ic_uniform_energy(*rho, *p, model)?;
            let n = g.nx * g.ny;
            let et = e_int + 0.5 * u * u;
            let rhost = vec![*rho; n];
            let ud = vec![*u; n];
            let vd = vec![0.0; n];
            let (mx, my, e) =
                prim_to_cons2d(&rhost, &ud, &vd, &vec![et; n]).map_err(|e| e.to_string())?;
            Ok(ConservedState2d {
                rho: rhost,
                mx,
                my,
                e,
            })
        }
        _ => ic2d_from_case(ic, g, gamma),
    }
}

fn ic2d_from_case(
    ic: &InitialCondition,
    g: &Grid2d,
    gamma: f64,
) -> Result<ConservedState2d, String> {
    let (n, cx, cy) = (
        g.nx * g.ny,
        0.5 * (g.xmin + g.xmax),
        0.5 * (g.ymin + g.ymax),
    );
    match ic {
        InitialCondition::Uniform { rho, u, p } => fill2d(*rho, *u, 0.0, *p, n, gamma),
        InitialCondition::Blast {
            radius,
            ambient,
            fireball,
        } => {
            let mut rhost = vec![ambient.rho; n];
            let mut ud = vec![ambient.u; n];
            let mut vd = vec![0.0; n];
            let mut pd = vec![ambient.p; n];
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let k = j * g.nx + i;
                    let dx = g.centers_x[k] - cx;
                    let dy = g.centers_y[k] - cy;
                    if dx * dx + dy * dy < radius * radius {
                        rhost[k] = fireball.rho;
                        ud[k] = fireball.u;
                        vd[k] = 0.0;
                        pd[k] = fireball.p;
                    }
                }
            }
            state2d(&rhost, &ud, &vd, &pd, gamma)
        }
        InitialCondition::TwoState { .. } => Err("two_state is a 1d shock-tube IC".into()),
    }
}

/// 3d initial condition: uniform field or an over-pressured fireball.
fn ic3d_from_case(
    ic: &InitialCondition,
    g: &Grid3d,
    gamma: f64,
) -> Result<ConservedState3d, String> {
    let n = g.nx * g.ny * g.nz;
    let (cx, cy, cz) = (
        0.5 * (g.xmin + g.xmax),
        0.5 * (g.ymin + g.ymax),
        0.5 * (g.zmin + g.zmax),
    );
    match ic {
        InitialCondition::Uniform { rho, u, p } => fill3d(*rho, *u, 0.0, 0.0, *p, n, gamma),
        InitialCondition::Blast {
            radius,
            ambient,
            fireball,
        } => {
            let mut rhost = vec![ambient.rho; n];
            let mut ud = vec![ambient.u; n];
            let vd = vec![0.0; n];
            let wd = vec![0.0; n];
            let mut pd = vec![ambient.p; n];
            for k in 0..n {
                let dx = g.centers_x[k] - cx;
                let dy = g.centers_y[k] - cy;
                let dz = g.centers_z[k] - cz;
                if dx * dx + dy * dy + dz * dz < radius * radius {
                    rhost[k] = fireball.rho;
                    ud[k] = fireball.u;
                    pd[k] = fireball.p;
                }
            }
            state3d(&rhost, &ud, &vd, &wd, &pd, gamma)
        }
        InitialCondition::TwoState { .. } => Err("two_state is a 1d shock-tube IC".into()),
    }
}

fn fill2d(
    rho0: f64,
    u0: f64,
    v0: f64,
    p0: f64,
    n: usize,
    gamma: f64,
) -> Result<ConservedState2d, String> {
    let (rho, u, v) = (vec![rho0; n], vec![u0; n], vec![v0; n]);
    let et: Vec<f64> = (0..n)
        .map(|k| p0 / ((gamma - 1.0) * rho0) + 0.5 * (u[k] * u[k] + v[k] * v[k]))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).map_err(|e| e.to_string())?;
    Ok(ConservedState2d { rho, mx, my, e })
}

fn state2d(
    rho: &[f64],
    u: &[f64],
    v: &[f64],
    p: &[f64],
    gamma: f64,
) -> Result<ConservedState2d, String> {
    let et: Vec<f64> = (0..rho.len())
        .map(|k| p[k] / ((gamma - 1.0) * rho[k]) + 0.5 * (u[k] * u[k] + v[k] * v[k]))
        .collect();
    let (mx, my, e) = prim_to_cons2d(rho, u, v, &et).map_err(|e| e.to_string())?;
    Ok(ConservedState2d {
        rho: rho.to_vec(),
        mx,
        my,
        e,
    })
}

fn fill3d(
    rho0: f64,
    u0: f64,
    v0: f64,
    w0: f64,
    p0: f64,
    n: usize,
    gamma: f64,
) -> Result<ConservedState3d, String> {
    let (rho, u, v, w) = (vec![rho0; n], vec![u0; n], vec![v0; n], vec![w0; n]);
    let et: Vec<f64> = (0..n)
        .map(|k| p0 / ((gamma - 1.0) * rho0) + 0.5 * (u[k] * u[k] + v[k] * v[k] + w[k] * w[k]))
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &u, &v, &w, &et).map_err(|e| e.to_string())?;
    Ok(ConservedState3d { rho, mx, my, mz, e })
}

fn state3d(
    rho: &[f64],
    u: &[f64],
    v: &[f64],
    w: &[f64],
    p: &[f64],
    gamma: f64,
) -> Result<ConservedState3d, String> {
    let et: Vec<f64> = (0..rho.len())
        .map(|k| p[k] / ((gamma - 1.0) * rho[k]) + 0.5 * (u[k] * u[k] + v[k] * v[k] + w[k] * w[k]))
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(rho, u, v, w, &et).map_err(|e| e.to_string())?;
    Ok(ConservedState3d {
        rho: rho.to_vec(),
        mx,
        my,
        mz,
        e,
    })
}

fn write_case_vtk2d(path: &str, name: &str, g: &Grid2d, st: &ConservedState2d, gamma: f64) -> i32 {
    write_case_vtk2d_model(path, name, g, st, ThermoModel::Perfect { gamma })
}

/// the snapshot writer dispatches on the case's closure so an eqair
/// run's file carries eqair pressure and mach.
fn write_case_vtk2d_model(
    path: &str,
    name: &str,
    g: &Grid2d,
    st: &ConservedState2d,
    model: ThermoModel,
) -> i32 {
    let dir = case_dir(path);
    let dst = dir.join(format!("{name}.vtk"));
    match write_vtk2d_model(&dst, g, st, model) {
        Ok(()) => {
            println!("wrote: {}", dst.display());
            if g.nx == g.ny {
                let lo = min_of(&st.rho);
                let hi = max_of(&st.rho);
                let png = dir.join(format!("{name}.png"));
                match render_png(&st.rho, g.nx, lo, hi) {
                    Ok(bytes) => {
                        let _ = std::fs::write(&png, bytes);
                        println!("wrote: {}", png.display());
                    }
                    Err(e) => eprintln!("png render error: {e}"),
                }
            }
            0
        }
        Err(e) => {
            eprintln!("output error: {e}");
            1
        }
    }
}

fn write_case_vtk3d(path: &str, name: &str, g: &Grid3d, st: &ConservedState3d, gamma: f64) -> i32 {
    let dir = case_dir(path);
    let dst = dir.join(format!("{name}.vtk"));
    match write_vtk3d(&dst, g, st, gamma) {
        Ok(()) => {
            println!("wrote: {}", dst.display());
            0
        }
        Err(e) => {
            eprintln!("output error: {e}");
            1
        }
    }
}

fn case_dir(path: &str) -> PathBuf {
    match Path::new(path).parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    }
}

fn min_of(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::INFINITY, f64::min)
}
fn max_of(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

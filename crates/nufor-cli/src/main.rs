//! command-line driver for nuFor: build, run, restart, serve, and export a 1D Euler case.
//!
//! `serve` hosts a minimal web ui skeleton: one page plus the snapshot as json,
//! ready for a real front end to be designed against the same data api.

use std::path::{Path, PathBuf};
use std::time::Instant;

use nufor_config::{
    load_case_config, BoundaryKind, CaseConfig, Equations, InflowProfileSpec, InitialCondition,
    Mesh, OutputFormat, Reconstruction,
};
use nufor_core::ffi_mpi::{
    march_euler2d_mpi, march_rans2d_mpi, march_viscous2d_mpi, mpi_bcast_state, mpi_launched,
    mpi_rank, mpi_shutdown,
};
use nufor_core::{
    advance2d_axi_lts, advance2d_axi_rk2, advance2d_model_rk2, advance2d_model_visc_lts,
    advance_ugrid_recon,
    advance2d_model_visc_rk2, advance2d_sa_lts, advance2d_sa_model_rk2, advance3d_rk2,
    advance_ugrid, apply_solid_fn, apply_solid_ref, body_cache, cons_to_prim2d, eos_pressure2d,
    euler_solve, grid1d, grid2d, grid3d, prim_to_cons, prim_to_cons2d, prim_to_cons3d,
    read_restart, read_restart2d, render_png, wall_distance2d, write_csv, write_h5, write_restart,
    write_restart2d, write_vtk, write_vtk2d_combustion, write_vtk2d_model, write_vtk2d_park5, write_vtk3d, Bc2d,
    Boundaries2d, Boundary, Bounds3d, ConservedState, ConservedState2d, ConservedState3d, Error,
    EulerConfig, Grid1d, Grid2d, Grid3d, OutputState, SaParams, SolidPolygon, SolidShape,
    SphereCone, SphereConeSdf, TurbState, Ugrid,
};
use nufor_core::{model_from_config, park2t, ThermoModel};

mod casework;
mod mesh_io;
mod progress;
mod stlcut;
mod watch;

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
        "watch" => watch::watch(&args),
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
         \x20 watch    <case.toml> [port]  live case viewer (default 8061)\n\
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
    // boundary group census: how many faces each physical group claimed.
    let mut groups = std::collections::BTreeMap::new();
    let mut untagged = 0usize;
    for f in 0..ug.face_left.len() {
        if ug.face_right[f] >= 0 {
            continue;
        }
        if ug.face_group[f] == 0 {
            untagged += 1;
        } else {
            *groups.entry(ug.face_group[f]).or_insert(0usize) += 1;
        }
    }
    if !groups.is_empty() || untagged > 0 {
        let gstr = groups
            .iter()
            .map(|(g, n)| format!("group {g}: {n}"))
            .collect::<Vec<_>>()
            .join(", ");
        println!("boundary faces: {untagged} untagged, {gstr}");
    }
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

/// build_2d_mesh with a case-file directory: a relative mesh path
/// resolves against the case's own directory first (so `nufor run
/// cases/x/case.toml` works from anywhere), falling back to the raw
/// path for absolute paths and repo-root-invoked runs.
fn build_2d_mesh_in(mesh: &Mesh, case_dir: Option<&std::path::Path>) -> Result<Grid2d, String> {
    if mesh.source == "file" {
        let path = mesh
            .path
            .as_deref()
            .ok_or_else(|| "mesh source = file needs a `path`".to_string())?;
        let m = match case_dir {
            Some(dir) => {
                let joined = dir.join(path);
                match mesh_io::load_rectilinear(joined.to_str().unwrap_or(path)) {
                    Ok(m) => m,
                    Err(_) => mesh_io::load_rectilinear(path)?,
                }
            }
            None => mesh_io::load_rectilinear(path)?,
        };
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
    // the live viewer boots with the run: a background thread serves
    // the case console and the march prints its link up front.
    if let Some(port) = watch::spawn_for_run(path) {
        println!("viewer: http://127.0.0.1:{port}/  (live case console)");
    }
    let rc = match cfg.physics.equations {
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
    };
    // a finished run leaves a standalone report next to the case: one
    // html file with everything inlined, openable offline or printable
    // to pdf for submission.
    if rc == 0 {
        if let Some(w) = watch::Watch::new(path) {
            let report = w.dir.join(format!("{}-report.html", cfg.metadata.name));
            match watch::write_report(&w, &report) {
                Ok(()) => {
                    println!("report: {}", report.display());
                    println!("viewer closed with the run; reopen: nufor watch {}", path);
                }
                Err(e) => eprintln!("report write error: {e}"),
            }
        }
    }
    rc
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
pub fn case_side_bc_euler(
    k: BoundaryKind,
    inflow: Option<(f64, f64, f64, f64)>,
    gas_const: f64,
    gamma: f64,
    p_back: Option<f64>,
) -> Bc2d {
    match k {
        BoundaryKind::Wall => Bc2d::NoSlipWall,
        BoundaryKind::SlipWall => Bc2d::SlipWall,
        BoundaryKind::SupersonicOutflow => Bc2d::SupersonicOutflow,
        BoundaryKind::SupersonicInflow | BoundaryKind::Inflow => match inflow {
            Some((rho, u, v, p)) => Bc2d::SupersonicInflow { rho, u, v, p },
            None => Bc2d::Transmissive,
        },
        // the reservoir inflow reads its stagnation state from
        // inflow_state: pt = the state's p, tt from rho/p through the
        // case gas constant.
        BoundaryKind::ReservoirInflow => match inflow {
            Some((rho, _u, _v, p)) if rho > 0.0 && p > 0.0 => Bc2d::ReservoirInflow {
                pt: p,
                tt: p / (rho * gas_const),
                gas_const,
                gamma,
            },
            _ => Bc2d::Transmissive,
        },
        BoundaryKind::SubsonicOutflow => match p_back {
            Some(p) if p > 0.0 => Bc2d::SubsonicOutflow { p_back: p },
            _ => Bc2d::Transmissive,
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
    let pb = cfg.boundaries.outflow_state.as_ref().map(|s| s.p);
    let (rc, gm) = (cfg.physics.gas_constant, cfg.physics.gamma);
    Boundaries2d {
        west: case_side_bc_euler(cfg.boundaries.left, case_inflow_state(cfg), rc, gm, pb),
        east: case_side_bc_euler(cfg.boundaries.right, None, rc, gm, pb),
        south: case_side_bc_euler(
            cfg.boundaries
                .bottom
                .unwrap_or(nufor_config::BoundaryKind::Outflow),
            None,
            rc,
            gm,
            pb,
        ),
        north: case_side_bc_euler(
            cfg.boundaries
                .top
                .unwrap_or(nufor_config::BoundaryKind::Outflow),
            None,
            rc,
            gm,
            pb,
        ),
    }
}

/// the axisymmetric euler run: planar fluxes with the annular update,
/// optional immersed body via the solid mask, slip-wall south boundary
/// as the symmetry axis.
fn run_case_axi(cfg: &CaseConfig, path: &str) -> i32 {
    let g = match build_2d_mesh_in(&cfg.mesh, Some(&case_dir(path))) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    // restart support mirrors the planar runner: a [time] restart_in
    // seeds the march (mesh dims must match), a restart_out saves the
    // final state for chaining.
    let mut restart_step = 0u64;
    let mut restart_time = 0.0f64;
    let model = model_from_config(cfg);
    let mut st = match &cfg.time.restart_in {
        Some(rpath) => match read_restart2d(Path::new(rpath)) {
            Ok(rd) => {
                if rd.nx != g.nx || rd.ny != g.ny {
                    eprintln!(
                        "restart mesh mismatch: case {}x{} vs file {}x{}",
                        g.nx, g.ny, rd.nx, rd.ny
                    );
                    return 2;
                }
                restart_step = rd.step;
                restart_time = rd.time;
                rd.state
            }
            Err(e) => {
                eprintln!("restart read error: {e}");
                return 2;
            }
        },
        None => match ic2d_from_case_model(&cfg.initial_condition, &g, gamma, model) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("initial condition error: {e}");
                return 2;
            }
        },
    };
    let body = case_body_fns(&cfg.body);
    if let Some((dist, normal)) = &body {
        apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
    }
    let mut bc = case_boundaries(cfg);
    // the axisymmetric defaults: the south side is the symmetry axis
    // (slip wall) and the north side vents unless the case says more.
    if cfg.boundaries.bottom.is_none() {
        bc.south = Bc2d::SlipWall;
    }
    if cfg.boundaries.top.is_none() {
        bc.north = Bc2d::SupersonicOutflow;
    }
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
    let mut t = restart_time;
    let mut steps = 0usize;
    let mut ok = true;
    let mut prog = progress::Progress::new(path, "2d-axi");
    while t < t_end && steps < max_steps {
        let stepped_axi = if cfg.numerics.local_time_stepping {
            advance2d_axi_lts(
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
        } else {
            advance2d_axi_rk2(
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
        match stepped_axi {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        };
        if let Some((dist, normal)) = &body {
            apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
        }
        steps += 1;
        if let Some(p) = prog.as_mut() {
            let dt_last = match stepped_axi {
                Ok((dt, _)) => dt,
                Err(_) => 0.0,
            };
            p.tick(
                steps as u64,
                t,
                dt_last,
                &st,
                progress::SampleExtra::default(),
            );
        }
        if steps % 200 == 0 {
            println!(
                "step {steps} t={t:.4} wall={:.1}s",
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
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
        let rc = write_case_vtk2d_model(
            path,
            &cfg.metadata.name,
            &g,
            &st,
            model,
            None,
            None,
        );
        if rc == 0 {
            write_run_log_2d(
                cfg,
                path,
                &g,
                &st,
                model,
                None,
                None,
                steps,
                t,
                t0.elapsed().as_secs_f64(),
            );
        }
        if let Some(rpath) = &cfg.time.restart_out {
            if let Err(e) = write_restart2d(
                Path::new(rpath),
                &st,
                g.nx,
                g.ny,
                gamma,
                t,
                restart_step + steps as u64,
            ) {
                eprintln!("restart write error: {e}");
            }
        }
        rc
    } else {
        if let Some(p) = prog.as_mut() {
            p.finish(false, steps as u64, t);
        }
        1
    }
}

/// march a uniform-IC unstructured case to its final time and report.
fn march_msh_uniform(
    cfg: &CaseConfig,
    path: &str,
    ug: Ugrid,
    st: &mut ConservedState2d,
    gamma: f64,
    model: ThermoModel,
) -> i32 {
    let cfl = cfg.numerics.cfl;
    let t_end = cfg.time.final_time;
    let max_steps = if cfg.time.max_steps > 0 {
        cfg.time.max_steps as usize
    } else {
        1_000_000
    };
    let nc = ug.cell_area.len();
    let t0 = Instant::now();
    let (mut t, mut steps) = (0.0, 0usize);
    let mut ok = true;
    let mut prog = progress::Progress::new(path, "2d-msh");
    while t < t_end && steps < max_steps {
        let stepped = match cfg.numerics.reconstruction {
            Reconstruction::Muscl => advance_ugrid_recon(st, &ug, gamma, cfl, model, true),
            Reconstruction::FirstOrder => advance_ugrid_recon(st, &ug, gamma, cfl, model, false),
        };
        let dt_last = match stepped {
            Ok((dt, _)) => dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        };
        t += dt_last;
        steps += 1;
        if let Some(p) = prog.as_mut() {
            p.tick(
                steps as u64,
                t,
                dt_last,
                st,
                progress::SampleExtra::default(),
            );
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
        if let Some(rpath) = &cfg.time.restart_out {
            // the structured restart format requires nx*ny dims, which
            // an unstructured cell count cannot honor; dump the raw
            // state (n, then rho/mx/my/e) for offline inspection.
            let mut out = Vec::with_capacity(8 + 4 * nc * 8);
            out.extend_from_slice(&(nc as u64).to_le_bytes());
            for a in [&st.rho, &st.mx, &st.my, &st.e] {
                for v in a {
                    out.extend_from_slice(&v.to_le_bytes());
                }
            }
            if let Err(e) = std::fs::write(rpath, &out) {
                eprintln!("state dump failed: {e}");
            } else {
                println!("state dump: {rpath}");
            }
        }
        println!(
            "case {} ({})\nmesh: {} cells ({} faces)\nsteps: {}\ntime: {:.5}\nwall: {:.3}s",
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
        if let Some(p) = prog.as_mut() {
            p.finish(false, steps as u64, t);
        }
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
    let mut prog = progress::Progress::new(path, "2d-ogrid");
    while t < t_end && steps < max_steps {
        let stepped_curv = advance2d_axi_curv(&mut st, &g, model, cfl, fs);
        match stepped_curv {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
        if let Some(p) = prog.as_mut() {
            let dt_last = match stepped_curv {
                Ok((dt, _)) => dt,
                Err(_) => 0.0,
            };
            p.tick(
                steps as u64,
                t,
                dt_last,
                &st,
                progress::SampleExtra::default(),
            );
        }
        if steps % 200 == 0 {
            println!(
                "step {steps} t={t:.4} wall={:.1}s",
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
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
        if let Some(p) = prog.as_mut() {
            p.finish(false, steps as u64, t);
        }
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
    let mut prog = progress::Progress::new(path, "1d");
    match euler_solve(&mut state, &ecfg) {
        Ok(r) => {
            if let Some(p) = prog.as_mut() {
                p.finish_1d(true, &r.log, r.converged, r.residual);
            }
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
        Err(e) => {
            if let Some(p) = prog.as_mut() {
                p.finish_1d(false, &[], false, 0.0);
            }
            euler_error(e)
        }
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

#[allow(clippy::too_many_arguments)]
/// solve a 2d case: blast or uniform IC on a uniform or imported rectilinear mesh.
/// the cae readout the heartbeat samples on its emit cadence: surface
/// force coefficients integrated over the body band, the peak wall
/// heat flux when the case is viscous, and a density-change residual
/// as the march's convergence measure. everything optional so a plain
/// shock tube case pays nothing.
fn progress_readout(
    cfg: &CaseConfig,
    g: &Grid2d,
    st: &nufor_core::ConservedState2d,
    body: &Option<CaseBodyFns>,
    model: ThermoModel,
    species: Option<&nufor_core::species::SpeciesState>,
    combustion: Option<&nufor_core::combustion_march::CombustionState>,
    viscous: bool,
) -> progress::SampleExtra {
    let n = st.rho.len();
    // cell pressure through the run's own closure
    let u: Vec<f64> = (0..n).map(|k| st.mx[k] / st.rho[k].max(1e-12)).collect();
    let v: Vec<f64> = (0..n).map(|k| st.my[k] / st.rho[k].max(1e-12)).collect();
    let e_int: Vec<f64> = (0..n)
        .map(|k| st.e[k] / st.rho[k].max(1e-12) - 0.5 * (u[k] * u[k] + v[k] * v[k]))
        .collect();
    let p: Vec<f64> = match (species, combustion) {
        (_, Some(sp)) => (0..n)
            .map(|k| {
                let mut y = [0.0f64; nufor_core::combustion_thermo::NS];
                #[allow(clippy::needless_range_loop)]
                #[allow(clippy::needless_range_loop)]
            for s in 0..nufor_core::combustion_thermo::NS {
                    y[s] = sp.y[s][k];
                }
                nufor_core::combustion_closure::pressure_at(st.rho[k], e_int[k], &y)
            })
            .collect(),
        (Some(sp), None) => (0..n)
            .map(|k| {
                let mut y = [0.0f64; 5];
                for (ys, ysrc) in y.iter_mut().zip(sp.y.iter()) {
                    *ys = ysrc[k];
                }
                nufor_core::species_closure::pressure_at(st.rho[k], e_int[k], &y)
            })
            .collect(),
        (None, None) => nufor_core::thermo_pressure(model, &st.rho, &e_int, &u, &v)
            .unwrap_or_else(|_| vec![0.0; n]),
    };
    // residual proxy: the max density change this step would imply at
    // the last dt. the true per-step delta needs the previous field,
    // so this uses the freestream deviation, a coarse but honest
    // convergence signal: it decays as the field settles.
    let rho_ref = match &cfg.initial_condition {
        nufor_config::InitialCondition::Uniform { rho, .. } => *rho,
        nufor_config::InitialCondition::TwoState { left, .. } => left.rho,
        nufor_config::InitialCondition::Blast { ambient, .. } => ambient.rho,
    };
    let mut residual = 0.0f64;
    for k in 0..n {
        residual = residual.max((st.rho[k] - rho_ref).abs());
    }
    // peak temperature through the same closure the field render uses
    let t_max: Option<f64> = match (species, combustion) {
        (_, Some(sp)) => (0..n)
            .map(|k| {
                let mut y = [0.0f64; nufor_core::combustion_thermo::NS];
                #[allow(clippy::needless_range_loop)]
                #[allow(clippy::needless_range_loop)]
            for s in 0..nufor_core::combustion_thermo::NS {
                    y[s] = sp.y[s][k];
                }
                nufor_core::combustion_closure::temperature_at(e_int[k], &y)
            })
            .fold(None::<f64>, |a, t| {
                Some(match a {
                    Some(m) if m >= t => m,
                    _ => t,
                })
            }),
        (_, None) => nufor_core::thermo_temperature(model, &st.rho, &p)
            .ok()
            .and_then(|ts| {
                ts.iter().cloned().fold(None::<f64>, |a, t| {
                    Some(match a {
                        Some(m) if m >= t => m,
                        _ => t,
                    })
                })
            }),
    };
    let mut extra = progress::SampleExtra {
        residual,
        t_max,
        ..Default::default()
    };
    let Some((dist, normal)) = body else {
        return extra;
    };
    // freestream dynamic pressure for the coefficients
    let (rho_fs, u_fs) = match &cfg.initial_condition {
        nufor_config::InitialCondition::Uniform { rho, u, .. } => (*rho, *u),
        nufor_config::InitialCondition::TwoState { left, .. } => (left.rho, left.u),
        nufor_config::InitialCondition::Blast { ambient, .. } => (ambient.rho, ambient.u),
    };
    let q_inf = 0.5 * rho_fs * u_fs * u_fs;
    if q_inf <= 0.0 {
        return extra;
    }
    let band = 1.5 * g.dx.min(g.dy);
    let ds = g.dx.max(g.dy) * 0.5;
    let (mut fx, mut fy) = (0.0f64, 0.0f64);
    let mut q_max = 0.0f64;
    let mut hits = 0usize;
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let d = dist(g.centers_x[k], g.centers_y[k]);
            if d > 0.0 && d < band {
                let (nx, ny) = normal(g.centers_x[k], g.centers_y[k]);
                fx -= p[k] * nx * ds;
                fy -= p[k] * ny * ds;
                hits += 1;
                // wall heat flux when a wall temperature is set: the
                // band cell's radial gradient toward the surface gives
                // q = k * dT/dn with k from the sutherland law at the
                // near-wall state.
                if viscous {
                    let tw = cfg.physics.wall_temperature.unwrap_or(0.0);
                    if tw > 0.0 {
                        let mu = cfg.physics.mu.unwrap_or(0.0);
                        let pr = 0.71;
                        let k_th = mu * 1005.0 / pr;
                        let dtdn =
                            ((e_int[k] + p[k] / st.rho[k].max(1e-12)) / 718.0 - tw) / d.max(1e-9);
                        q_max = q_max.max(k_th * dtdn.abs());
                    }
                }
            }
        }
    }
    if hits > 0 {
        // 2d section force -> coefficient per unit depth; A_ref from
        // metadata when present, else the body's own frontal height.
        // reference area: the body's frontal height times unit depth;
        // the band cell count over ny approximates the wetted height.
        let a_ref = g.dy * hits as f64 / g.ny as f64;
        extra.cd = Some(fx / (q_inf * a_ref));
        extra.cl = Some(fy / (q_inf * a_ref));
        extra.q_max = if q_max > 0.0 { Some(q_max) } else { None };
    }
    extra
}

fn run_case_2d(cfg: &CaseConfig, path: &str) -> i32 {
    let mesh = &cfg.mesh;
    // a gmsh .msh file drives the unstructured 2d solver instead.
    if let (Some(mp), true) = (
        mesh.path.as_deref(),
        mesh.path.as_deref().unwrap_or("").ends_with(".msh"),
    ) {
        return run_case_2d_msh(cfg, path, mp);
    }
    let g = match build_2d_mesh_in(mesh, Some(&case_dir(path))) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    let model = model_from_config(cfg);
    // a [time] restart_in seeds the march from a prior run's state; the
    // mesh dims must match the case file's.
    let mut restart_step = 0u64;
    let mut restart_time = 0.0f64;
    let mut st = match &cfg.time.restart_in {
        Some(rpath) => {
            let rd = match read_restart2d(Path::new(rpath)) {
                Ok(rd) => rd,
                Err(e) => {
                    eprintln!("restart read error: {e}");
                    return 2;
                }
            };
            if rd.nx != g.nx || rd.ny != g.ny {
                eprintln!(
                    "restart mesh {}x{} does not match case {}x{}",
                    rd.nx, rd.ny, g.nx, g.ny
                );
                return 2;
            }
            restart_step = rd.step;
            restart_time = rd.time;
            rd.state
        }
        None => match ic2d_from_case_model(&cfg.initial_condition, &g, gamma, model) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("initial condition error: {e}");
                return 2;
            }
        },
    };
    let body = case_body_fns(&cfg.body);
    // the viscous lane owns its wall: the solid fill takes the freestream
    // reference (the eqair closure reads the scalar inert reference as a
    // cold dense slug) and the slip band projection stays off so the
    // no-slip wall fluxes do the braking.
    let visc_no_slip = cfg.physics.mu.is_some_and(|m| m > 0.0);
    let fs_ref = match (&cfg.initial_condition, model) {
        (nufor_config::InitialCondition::Uniform { rho, u, p }, ThermoModel::EqAir { .. }) => {
            Some((*rho, *p, *u))
        }
        _ => None,
    };
    if let Some((dist, normal)) = &body {
        if let Some((rho_fs, p_fs, u_fs)) = fs_ref {
            apply_solid_ref(
                &mut st,
                &g,
                dist.as_ref(),
                normal.as_ref(),
                gamma,
                Some((rho_fs, p_fs)),
                visc_no_slip,
            );
            let _ = u_fs;
        } else {
            apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
        }
    }
    let bc = case_boundaries(cfg);
    // park5 cases carry the composition beside the flow state: pure
    // nitrogen for the dissociation gates, standard air otherwise.
    let mut park5_species = if matches!(model, ThermoModel::Park5 { .. }) {
        let built = match cfg.physics.gas.as_str() {
            "n2" => {
                nufor_core::species::SpeciesState::uniform(g.nx * g.ny, &[1.0, 0.0, 0.0, 0.0, 0.0])
            }
            _ => Ok(nufor_core::species::SpeciesState::air(g.nx * g.ny)),
        };
        Some(match built {
            Ok(s) => s,
            Err(e) => {
                eprintln!("park5 initial composition error: {e}");
                return 2;
            }
        })
    } else {
        None
    };
    // the combustion lane carries its own 9-species state: the
    // stoichiometric unburned fill (a chamber case reseeds through
    // the inflow boundary). chemistry off freezes the composition.
    let mut combustion_state = if matches!(model, ThermoModel::Combustion9 { .. }) {
        Some(nufor_core::combustion_march::CombustionState::stoich_h2o2(g.nx * g.ny))
    } else {
        None
    };
    let combustion_react = match cfg.physics.chemistry.as_deref() {
        Some("off") => false,
        Some("on") | None => true,
        Some(other) => {
            eprintln!("unknown chemistry '{other}' (use on or off)");
            return 2;
        }
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
    // mpi path: under mpirun (or NFOR_MPI=1) the eqair euler march runs in
    // the fortran core across ranks. laminar viscous eqair cases route
    // through the fortran viscous step (two-stage heun euler + viscous per
    // call); turbulence stays serial until its port. body cases march in
    // single-step chunks so rank 0 can re-project the band cells between
    // steps and broadcast the refreshed state; body-free cases march in
    // one call. viscous body cases stay on the serial rust lane: the
    // fortran viscous kernel wall-treats domain faces only, so a masked
    // body under mpirun would silently skip the no-slip band.
    let visc_lam = cfg.physics.mu.is_some_and(|m| m > 0.0) && turb.is_none();
    let rans_sa = turb.is_some() && cfg.physics.mu.is_some_and(|m| m > 0.0);
    let visc_body = (visc_lam || rans_sa) && body.is_some();
    if mpi_launched()
        && matches!(model, ThermoModel::EqAir { .. })
        && !visc_body
        && (cfg.physics.mu.map_or(true, |m| m == 0.0) || visc_lam || rans_sa)
    {
        let is_root = mpi_rank().unwrap_or(0) == 0;
        let mut solid_mask = vec![false; g.nx * g.ny];
        if let Some((dist, _)) = &body {
            let band = 1.5 * g.dx.min(g.dy);
            solid_mask = body_cache(&g, dist.as_ref(), band).mask.clone();
        }
        let t0m = Instant::now();
        let mut t = restart_time;
        let mut steps = 0usize;
        let mut ok = true;
        // rank 0 owns the progress heartbeat for the chunked march.
        let mut prog = if is_root {
            progress::Progress::new(path, "2d-mpi")
        } else {
            None
        };
        // the viscous lane: same chunked loop, the fortran viscous heun
        // step per call. the mu/pr/wall params mirror the serial
        // advance2d_model_visc_rk2 call site.
        let mu0 = cfg.physics.mu.filter(|m| *m > 0.0);
        let pr = cfg.physics.pr;
        let tw = cfg.physics.wall_temperature.unwrap_or(0.0);
        while t < t_end && steps < max_steps {
            let pr_t = turb.as_ref().map_or(0.9, |tb| tb.params.pr_t);
            let r = match (&mu0, &turb) {
                (Some(m), Some(tb)) => march_rans2d_mpi(
                    &st,
                    tb,
                    &g,
                    model,
                    &bc,
                    cfl,
                    1,
                    &solid_mask,
                    *m,
                    pr,
                    pr_t,
                    tw,
                ),
                (Some(m), None) => march_viscous2d_mpi(
                    &st,
                    &g,
                    model,
                    &bc,
                    cfl,
                    f64::INFINITY,
                    1,
                    &solid_mask,
                    *m,
                    pr,
                    tw,
                )
                .map(|r| (r, None)),
                (None, _) => {
                    march_euler2d_mpi(&st, &g, model, &bc, cfl, f64::INFINITY, 1, &solid_mask)
                        .map(|r| (r, None))
                }
            };
            if let Some(nt) = r.as_ref().ok().and_then(|(_, nt)| nt.clone()) {
                // keep the local turb field in lockstep with the gathered
                // fortran transport output
                if let Some(tb) = turb.as_mut() {
                    tb.nu_tilde = nt;
                }
            }
            #[allow(unused_assignments)]
            let mut dt_chunk = 0.0f64;
            match r {
                Ok((mr, _)) => {
                    if let Some(gathered) = mr.state {
                        st = gathered;
                    }
                    // one step per call: the march's internal clock starts
                    // at zero each call, so accumulate the step's dt here.
                    t += mr.dt_last;
                    dt_chunk = mr.dt_last;
                    if is_root {
                        eprintln!(
                            "chunk {steps}: t={t:.6e} dt_last={:.6e} rho_min={:.6} rho_max={:.6}",
                            mr.dt_last,
                            st.rho.iter().cloned().fold(f64::INFINITY, f64::min),
                            st.rho.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                        );
                    }
                }
                Err(e) => {
                    eprintln!("solver error at step {steps}: {e}");
                    ok = false;
                    break;
                }
            }
            if let Some((dist, normal)) = &body {
                apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
            }
            // refresh every rank's full state from rank 0 so the next chunk
            // carves identical slabs (and band projection lands everywhere)
            let _ = mpi_bcast_state(&mut st);
            steps += 1;
            if let Some(p) = prog.as_mut() {
                p.tick(
                    steps as u64,
                    t,
                    dt_chunk,
                    &st,
                    progress::SampleExtra::default(),
                );
            }
        }
        if !ok {
            if let Some(p) = prog.as_mut() {
                p.finish(false, steps as u64, t);
            }
            mpi_shutdown();
            return 1;
        }
        let wall = t0m.elapsed().as_secs_f64();
        if !is_root {
            mpi_shutdown();
            return 0;
        }
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
        println!(
            "case {} ({}) [mpi]\nthink: 2d, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name, path, g.nx, g.ny, steps, t, wall
        );
        eprintln!(
            "dbg: rank0 state rho [{:.6},{:.6}] e [{:.6e},{:.6e}] steps={steps} t={t:.6e}",
            st.rho.iter().cloned().fold(f64::INFINITY, f64::min),
            st.rho.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            st.e.iter().cloned().fold(f64::INFINITY, f64::min),
            st.e.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        let rc = write_case_vtk2d_model(
            path,
            &cfg.metadata.name,
            &g,
            &st,
            model,
            None,
            None,
        );
        if rc == 0 {
            write_run_log_2d(
                cfg,
                path,
                &g,
                &st,
                model,
                park5_species.as_ref(),
                combustion_state.as_ref(),
                steps,
                t,
                wall,
            );
        }
        if let Some(rpath) = &cfg.time.restart_out {
            if let Err(e) = write_restart2d(
                Path::new(rpath),
                &st,
                g.nx,
                g.ny,
                gamma,
                t,
                restart_step + steps as u64,
            ) {
                eprintln!("restart write error: {e}");
            }
        }
        mpi_shutdown();
        return rc;
    }

    // under mpirun but NOT routed through the mpi branch (viscous,
    // non-eqair closure, or turbulence): only rank 0 marches and
    // writes. non-root ranks would just duplicate the serial run and
    // clobber the output files.
    if mpi_launched() && mpi_rank().unwrap_or(0) != 0 {
        mpi_shutdown();
        return 0;
    }

    // the park two-temperature field: vibrational energy per cell, J/kg,
    // seeded at the freestream equilibrium value for the initial
    // translational temperature. the lane advances with the euler march
    // (perfect-gas fluxes) plus the landau-teller substep per step.
    let park_lane = matches!(cfg.physics.eos, nufor_config::Eos::Park2t);
    let mut park = park2t::Park2tState::from_uniform(&g, 0.0);
    let mut park_solid = vec![false; g.nx * g.ny];
    if park_lane {
        // seed the vibrational pool at the initial translational
        // temperature: freestream equilibrium (casseau: vibrational
        // modes excited at 220 K already).
        let t_seed = match cfg.initial_condition {
            nufor_config::InitialCondition::Uniform { rho, p, .. } => {
                p / (rho * cfg.physics.gas_constant)
            }
            _ => 220.0,
        };
        park = park2t::Park2tState::from_uniform(&g, park2t::e_vib_of_t_v(t_seed));
        if let Some((dist, _)) = &body {
            let band = 1.5 * g.dx.min(g.dy);
            park_solid = body_cache(&g, dist.as_ref(), band).mask.clone();
        }
    }

    let t0 = Instant::now();
    let mut t = restart_time;
    let mut steps = 0usize;
    let mut ok = true;
    let mut prog = progress::Progress::new(path, "2d");
    while t < t_end && steps < max_steps {
        let stepped = if let Some(sp) = combustion_state.as_mut() {
            let mu_case = cfg.physics.mu.unwrap_or(0.0);
            if mu_case > 0.0 {
                let tw = cfg.physics.wall_temperature.unwrap_or(0.0);
                nufor_core::combustion_viscous::advance2d_combustion_visc_rk2(
                    &mut st,
                    sp,
                    &g,
                    cfl,
                    true,
                    &bc,
                    tw,
                    cfg.numerics.threads,
                    combustion_react,
                )
            } else {
                nufor_core::combustion_march::advance2d_combustion_rk2(
                    &mut st,
                    sp,
                    &g,
                    cfl,
                    true,
                    &bc,
                    cfg.numerics.threads,
                    combustion_react,
                )
            }
        } else if let Some(sp) = park5_species.as_mut() {
            // the park5 lane first: a park5 case with mu > 0 routes to
            // the species viscous march, and the generic mu branch below
            // must never see a park5 model (it rejects it).
            let mu_case = cfg.physics.mu.unwrap_or(0.0);
            if mu_case > 0.0 {
                let tw = cfg.physics.wall_temperature.unwrap_or(0.0);
                let cat = match cfg.physics.catalytic_wall.as_deref() {
                    Some("full") | Some("fully") | Some("fully_catalytic") => {
                        nufor_core::catalytic::CatalyticWall::FullyCatalytic
                    }
                    Some("none") | None => nufor_core::catalytic::CatalyticWall::NonCatalytic,
                    Some(other) => {
                        eprintln!("unknown catalytic_wall '{other}' (use none or full)");
                        ok = false;
                        break;
                    }
                };
                nufor_core::species_march::advance2d_park5_visc_rk2(
                    &mut st,
                    sp,
                    &g,
                    cfl,
                    true,
                    &bc,
                    tw,
                    cat,
                    cfg.numerics.threads,
                    body.as_ref()
                        .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
                )
            } else {
                nufor_core::species_march::advance2d_park5_rk2(
                    &mut st,
                    sp,
                    &g,
                    cfl,
                    true,
                    &bc,
                    cfg.numerics.threads,
                    body.as_ref()
                        .map(|b| b.0.as_ref() as &dyn Fn(f64, f64) -> f64),
                )
            }
        } else if let Some(mu0) = cfg.physics.mu {
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
                } else if cfg.numerics.local_time_stepping {
                    // steady-state acceleration: per-cell increments, so
                    // only the converged state is meaningful.
                    advance2d_model_visc_lts(
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
                    .map(|(dt, _)| (dt, 0.0))
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
                    .map(|(dt, _)| (dt, 0.0))
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
            if let Some((rho_fs, p_fs, _u_fs)) = fs_ref {
                apply_solid_ref(
                    &mut st,
                    &g,
                    dist.as_ref(),
                    normal.as_ref(),
                    gamma,
                    Some((rho_fs, p_fs)),
                    visc_no_slip,
                );
            } else {
                apply_solid_fn(&mut st, &g, dist.as_ref(), normal.as_ref(), gamma);
            }
        }
        if park_lane {
            let dt_last = match stepped {
                Ok((dt, _)) => dt,
                Err(_) => break,
            };
            park2t::park2t_substep(&mut st, &mut park, &g, gamma, dt_last, &park_solid, false)
                .map_err(|e| eprintln!("park2t substep error: {e}"))
                .ok();
        }
        steps += 1;
        if let Some(p) = prog.as_mut() {
            let dt_last = match stepped {
                Ok((dt, _)) => dt,
                Err(_) => 0.0,
            };
            // the cae readout rides the emit cadence only: on the steps
            // the heartbeat actually samples, integrate the surface
            // forces and the peak wall heat flux. off steps pay a
            // branch and nothing else. the cadence emits at step 1
            // then every 50th (1, 51, 101, ...), so the same gate
            // keeps the readout aligned with a kept sample.
            if steps % 50 == 1 {
                // the cadence sample: full extrema scan plus the cae
                // readout, then the heavy tick.
                let extra = progress_readout(
                    cfg,
                    &g,
                    &st,
                    &body,
                    model,
                    park5_species.as_ref(),
                    combustion_state.as_ref(),
                    visc_no_slip,
                );
                p.tick(steps as u64, t, dt_last, &st, extra);
            } else {
                // every other step: refresh the header the progress
                // box reads, without walking the field.
                p.step_tick(steps as u64, t, dt_last);
            }
        }
        if steps % 200 == 0 {
            eprintln!(
                "step {steps} t={t:.4} wall={:.1}s",
                t0.elapsed().as_secs_f64()
            );
        }
        // interval snapshots: a long march should not lose everything
        // when a late step goes bad; the interval vtk keeps the last
        // checkpoint on disk (same name as the final write).
        let iv = cfg.output.interval_steps as usize;
        if iv > 0 && steps % iv == 0 {
            let rc = write_case_vtk2d_model(
                path,
                &cfg.metadata.name,
                &g,
                &st,
                model,
                park5_species.as_ref(),
                combustion_state.as_ref(),
            );
            if rc != 0 {
                eprintln!("interval output failed at step {steps}");
            }
            // the run-log ride-along is the expensive half of the
            // interval write (it walks the field again); keep it on
            // the 200-step cadence even when the vtk writes faster.
            if steps % 200 == 0 {
                write_run_log_2d(
                    cfg,
                    path,
                    &g,
                    &st,
                    model,
                    park5_species.as_ref(),
                    combustion_state.as_ref(),
                    steps,
                    t,
                    t0.elapsed().as_secs_f64(),
                );
                let dst = case_dir(path).join("run-log-last.txt");
                let _ = std::fs::rename(case_dir(path).join("run-log.txt"), dst);
            }
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
        let wall = t0.elapsed().as_secs_f64();
        println!(
            "case {} ({})\nthink: 2d, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name, path, g.nx, g.ny, steps, t, wall
        );
        let rc = write_case_vtk2d_model(
            path,
            &cfg.metadata.name,
            &g,
            &st,
            model,
            park5_species.as_ref(),
            combustion_state.as_ref(),
        );
        if rc == 0 && park_lane {
            // the vibrational temperature field: the two-temperature
            // readout (t_v per cell from the stored vibrational energy).
            let tv: Vec<f64> = park.e_v.iter().map(|e| park2t::t_vib_of_e_v(*e)).collect();
            match write_scalar_vtk2d(
                &case_dir(path).join(format!("{}-tv.vtk", cfg.metadata.name)),
                &cfg.metadata.name,
                &g,
                &tv,
                "t_vib",
            ) {
                Ok(()) => println!("wrote: t_vib field appended"),
                Err(e) => eprintln!("t_vib write error: {e}"),
            }
        }
        if rc == 0 {
            write_run_log_2d(
                cfg,
                path,
                &g,
                &st,
                model,
                park5_species.as_ref(),
                combustion_state.as_ref(),
                steps,
                t,
                wall,
            );
        }
        if let Some(rpath) = &cfg.time.restart_out {
            if let Err(e) = write_restart2d(
                Path::new(rpath),
                &st,
                g.nx,
                g.ny,
                gamma,
                t,
                restart_step + steps as u64,
            ) {
                eprintln!("restart write error: {e}");
            }
        }
        rc
    } else {
        if let Some(p) = prog.as_mut() {
            p.finish(false, steps as u64, t);
        }
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
    species: Option<&nufor_core::species::SpeciesState>,
    combustion: Option<&nufor_core::combustion_march::CombustionState>,
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
    // park5 computes p and T through the species closure per cell;
    // the other models keep the shared table paths.
    let (p, temp) = if let Some(sp) = combustion {
        let mut pp = vec![0.0f64; n];
        let mut tt = vec![0.0f64; n];
        for k in 0..n {
            let mut y = [0.0f64; nufor_core::combustion_thermo::NS];
            #[allow(clippy::needless_range_loop)]
            for s in 0..nufor_core::combustion_thermo::NS {
                y[s] = sp.y[s][k];
            }
            pp[k] = nufor_core::combustion_closure::pressure_at(rho[k], e_int[k], &y);
            tt[k] = nufor_core::combustion_closure::temperature_at(e_int[k], &y);
        }
        (pp, tt)
    } else if let Some(sp) = species {
        let mut pp = vec![0.0f64; n];
        let mut tt = vec![0.0f64; n];
        for k in 0..n {
            let mut y = [0.0f64; 5];
            for (ys, ysrc) in y.iter_mut().zip(sp.y.iter()) {
                *ys = ysrc[k];
            }
            pp[k] = nufor_core::species_closure::pressure_at(rho[k], e_int[k], &y);
            tt[k] = nufor_core::species_closure::temperature_at(e_int[k], &y);
        }
        (pp, tt)
    } else {
        let p = match nufor_core::thermo_pressure(model, rho, &e_int, &u, &v) {
            Ok(x) => x,
            Err(_) => return,
        };
        let temp = nufor_core::thermo_temperature(model, rho, &p).unwrap_or_default();
        (p, temp)
    };
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
    let cat_wall = match cfg.physics.catalytic_wall.as_deref() {
        Some("full") | Some("fully") | Some("fully_catalytic") => {
            nufor_core::catalytic::CatalyticWall::FullyCatalytic
        }
        _ => nufor_core::catalytic::CatalyticWall::NonCatalytic,
    };
    // wall heat flux report along body-mask faces: q_w from a one-sided
    // temperature gradient on the fluid side of the mask, scaled by the
    // sutherland conductivity at the wall temperature. only emitted for
    // viscous runs with a finite wall temperature.
    let mut qwall_s = String::new();
    if visc && tw > 0.0 {
        let cb = case_body_fns(&cfg.body);
        if let Some((dist, _normal)) = cb {
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
            // fourier conductivity from prandtl: k = mu * cp / pr with
            // cp = R*gamma/(gamma-1) and the case gas constant. park5
            // cases use the wilke mixture conductivity at the near-wall
            // composition instead, so the log matches the march.
            let kf_suth = |t: f64| {
                nufor_core::sutherland_mu(1.716e-5, t, 273.15, 110.4)
                    * cfg.physics.gamma
                    * cfg.physics.gas_constant
                    / ((cfg.physics.gamma - 1.0) * cfg.physics.pr)
            };
            let kf = |t: f64, k: usize| {
                if let Some(sp) = species {
                    let mut y = [0.0f64; 5];
                    for (ys, ysrc) in y.iter_mut().zip(sp.y.iter()) {
                        *ys = ysrc[k];
                    }
                    nufor_core::species_props::mixture_conductivity(&y, t)
                } else {
                    kf_suth(t)
                }
            };
            // the solid mask the march itself used: a wall face is a fluid
            // cell with a solid neighbor. the flux is the same one-sided
            // quadratic through the pinned wall temperature and the first
            // two fluid cells the solver's conductor evaluates, so the log
            // matches the wall energy budget the march actually applied.
            let mask = nufor_core::solid_mask(g, dist.as_ref());
            let solid = |i: usize, j: usize| -> bool { i < g.nx && j < g.ny && mask[j * g.nx + i] };
            let temp_at = |k: usize| temp.get(k).copied().unwrap_or(tw);
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let k = j * g.nx + i;
                    if solid(i, j) {
                        continue;
                    }
                    // wall faces on x: fluid cell with an east or west solid
                    let (a, m, face_x) = if i + 1 < g.nx && solid(i + 1, j) {
                        if i < 2 {
                            continue;
                        }
                        (i, i - 1, g.faces_x[i + 1])
                    } else if i > 0 && solid(i - 1, j) {
                        if i + 2 >= g.nx {
                            continue;
                        }
                        (i, i + 1, g.faces_x[i])
                    } else {
                        continue;
                    };
                    let ka = j * g.nx + a;
                    let km = j * g.nx + m;
                    // lagrange derivative at the face through (tw@face,
                    // T_a, T_m) on true centers
                    let (y0, y1, y2) = (
                        face_x,
                        0.5 * (g.faces_x[a] + g.faces_x[a + 1]),
                        0.5 * (g.faces_x[m] + g.faces_x[m + 1]),
                    );
                    let w0 = (2.0 * y0 - y1 - y2) / ((y0 - y1) * (y0 - y2));
                    let w1 = (2.0 * y0 - y0 - y2) / ((y1 - y0) * (y1 - y2));
                    let w2 = (2.0 * y0 - y0 - y1) / ((y2 - y0) * (y2 - y1));
                    let dtdn = tw * w0 + temp_at(ka) * w1 + temp_at(km) * w2;
                    let mut q = (kf(tw, ka) * dtdn).abs();
                    // park5 fully-catalytic: add the formation-enthalpy
                    // flux the atoms carry into the wall. the same
                    // one-sided composition stencil the march applies;
                    // non-catalytic adds exactly zero.
                    if let Some(sp) = species {
                        if cat_wall == nufor_core::catalytic::CatalyticWall::FullyCatalytic {
                            let mut y0 = [0.0f64; 5];
                            let mut y1v = [0.0f64; 5];
                            for s in 0..5 {
                                y0[s] = sp.y[s][ka];
                                y1v[s] = sp.y[s][km];
                            }
                            let y_wall = nufor_core::catalytic::equilibrium_wall_composition(
                                tw.max(200.0),
                                rho[ka],
                                &y0,
                            )
                            .unwrap_or(y0);
                            // positions along x: wall face, then the two
                            // fluid-cell centers (the same ascending triple
                            // the march's wall stencil evaluates)
                            let d_eff = {
                                let mut y = [0.0f64; 5];
                                for (ys, ysrc) in y.iter_mut().zip(sp.y.iter()) {
                                    *ys = ysrc[ka];
                                }
                                let t_a = temp_at(ka);
                                let p_a = nufor_core::species_closure::pressure_at(
                                    rho[ka], e_int[ka], &y,
                                );
                                nufor_core::species_props::effective_diffusivity(
                                    &y,
                                    t_a,
                                    p_a.max(1.0),
                                )
                            };
                            // wallward = toward the solid: the derivative
                            // the surface flux evaluates is along +x, so an
                            // east-wall face (solid at i+1, wallward = -x)
                            // negates it; west-wall faces keep it.
                            let wallward_sign = if m == a - 1 { -1.0 } else { 1.0 };
                            q += wallward_sign
                                * nufor_core::catalytic::catalytic_surface_flux(
                                    &y_wall,
                                    &y0,
                                    &y1v,
                                    rho[ka],
                                    d_eff,
                                    [
                                        face_x,
                                        0.5 * (g.faces_x[a] + g.faces_x[a + 1]),
                                        0.5 * (g.faces_x[m] + g.faces_x[m + 1]),
                                    ],
                                )
                                .max(0.0);
                        }
                    }
                    if q.is_finite() && q > 0.0 {
                        let c = cs[k];
                        q_rows.push((c.0, c.1, q));
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
    // continuum-validity probe: bird's gradient-length local
    // knudsen over the flow cells, with the 0.05 breakdown threshold.
    // the molecular viscosity is sutherland at the local temperature
    // when the case gives a reference mu, else the perfect-gas N2 law.
    let continuum_line = {
        let gas_const = cfg.physics.gas_constant;
        let l_ref = cfg
            .physics
            .reference
            .as_ref()
            .map(|r| r.l)
            .filter(|l| *l > 0.0)
            .unwrap_or_else(|| (g.faces_x[g.nx] - g.faces_x[0]).abs());
        let mu_ref = cfg.physics.mu.filter(|m| *m > 0.0);
        let mu_cell: Vec<f64> = if let Some(mu0) = mu_ref {
            temp.iter()
                .map(|&t| nufor_core::sutherland_mu(mu0, t.max(1.0), 273.15, 110.4))
                .collect()
        } else {
            temp.iter()
                .map(|&t| nufor_core::fay_riddell::mu_n2(t.max(1.0)))
                .collect()
        };
        let mut solid = vec![false; g.nx * g.ny];
        if let Some((dist, _)) = case_body_fns(&cfg.body) {
            solid = nufor_core::solid_mask(g, dist.as_ref()).to_vec();
        }
        let pr = nufor_core::continuum_probe::continuum_probe(
            g,
            rho,
            &temp,
            &mu_cell,
            gas_const,
            &solid,
            rho_ref,
            temp.first().copied().unwrap_or(220.0),
            mu_cell.first().copied().unwrap_or(1.7e-5),
            l_ref,
        );
        let (wi, wj) = (pr.worst_cell % g.nx, pr.worst_cell / g.nx);
        let (wx, wy) = (g.centers_x[wi], g.centers_y[wj * g.nx]);
        format!(
            "\ncontinuum probe (bird gradient-length knudsen):\n  kn_inf (lambda_inf/l_ref) {:.3e}\n  max kn_gll(rho) {:.3e} at ({:.3}, {:.3})\n  max kn_gll(t) {:.3e}\n  cells above breakdown threshold {:.1} pct\n",
            pr.kn_inf,
            pr.kn_rho_max,
            wx,
            wy,
            pr.kn_t_max,
            100.0 * pr.frac_breakdown
        )
    };
    let mut txt = txt;
    if !noneq_line.is_empty() {
        txt.push_str(&noneq_line);
    }
    txt.push_str(&continuum_line);
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
    let model = nufor_core::model_from_config(cfg);
    // the fortran mpi rans lane covers the eqair closure without a body
    // (masked bodies stay serial: the fortran kernel wall-treats domain
    // faces only). run_case_2d builds the same turb state from the case
    // config and routes it through march_rans2d_mpi under mpirun.
    if mpi_launched() && matches!(model, ThermoModel::EqAir { .. }) && cfg.body.is_none() {
        return run_case_2d(cfg, path);
    }
    let g = match build_2d_mesh_in(&cfg.mesh, Some(&case_dir(path))) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let gamma = cfg.physics.gamma;
    let mut st = match ic2d_from_case_model(&cfg.initial_condition, &g, gamma, model) {
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
        // the reservoir inflow and subsonic back-pressure outflow are
        // euler-lane bcs; the sa lane treats them as transmissive until
        // its own nozzle story exists.
        BoundaryKind::ReservoirInflow | BoundaryKind::SubsonicOutflow => Bc2d::Transmissive,
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
    let mut prog = progress::Progress::new(path, "2d-sa");
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
        } else {
            // the closure-aware sa driver: per-cell sutherland viscosity and
            // the case's eos for pressure/temperature, threaded on
            // cfg.numerics.threads. matches the fortran mpi rans lane.
            match nufor_core::advance2d_sa_model_rk2(
                &mut st,
                &mut turb,
                &g,
                model,
                cfl,
                true,
                &bc,
                mu,
                273.15,
                110.4,
                cfg.physics.wall_temperature.unwrap_or(0.0),
                cfg.numerics.threads.max(1),
                None,
            ) {
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
        if let Some(p) = prog.as_mut() {
            p.tick(steps as u64, t, 0.0, &st, progress::SampleExtra::default());
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
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
        write_case_vtk2d_model(path, &cfg.metadata.name, &g, &st, model, None, None)
    } else {
        1
    }
}

/// solve a 2d case on an imported gmsh (.msh) unstructured mesh: build the
/// cell-centered grid, drop a blast IC, and advance with the face-based solver.
fn run_case_2d_msh(cfg: &CaseConfig, path: &str, msh_path: &str) -> i32 {
    let model = model_from_config(cfg);
    let m = match mesh_io::load_gmsh(msh_path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let mut ug: Ugrid = match m.ugrid() {
        Ok(g) => g,
        Err(e) => {
            eprintln!("mesh error: {e}");
            return 2;
        }
    };
    let nc = ug.cell_area.len();
    let gamma = cfg.physics.gamma;
    // a uniform IC seeds every cell with the freestream state and
    // matches the inflow boundary to it; the blast IC keeps its demo
    // path below.
    if let nufor_config::InitialCondition::Uniform { rho, u, p } = &cfg.initial_condition {
        ug.freestream = [*rho, *u, 0.0, *p];
        // total specific energy = internal + kinetic; prim_to_cons2d
        // stores e = rho * et verbatim.
        let et_i = p / (rho * (gamma - 1.0));
        let et = vec![et_i + 0.5 * u * u; nc];
        let (mx, my, e) = prim_to_cons2d(&vec![*rho; nc], &vec![*u; nc], &vec![0.0; nc], &et)
            .map_err(|e| e.to_string())
            .unwrap();
        let mut st = ConservedState2d {
            rho: vec![*rho; nc],
            mx,
            my,
            e,
        };
        return march_msh_uniform(cfg, path, ug, &mut st, gamma, model);
    }
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
    let mut prog = progress::Progress::new(path, "2d-msh");
    while t < t_end && steps < max_steps {
        let stepped_msh = match cfg.numerics.reconstruction {
            Reconstruction::Muscl => advance_ugrid_recon(&mut st, &ug, gamma, cfl, model, true),
            Reconstruction::FirstOrder => advance_ugrid(&mut st, &ug, gamma, cfl, model),
        };
        match stepped_msh {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
        if let Some(p) = prog.as_mut() {
            let dt_last = match stepped_msh {
                Ok((dt, _)) => dt,
                Err(_) => 0.0,
            };
            p.tick(
                steps as u64,
                t,
                dt_last,
                &st,
                progress::SampleExtra::default(),
            );
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
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
        if let Some(p) = prog.as_mut() {
            p.finish(false, steps as u64, t);
        }
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
    let mut prog = progress::Progress::new(path, "3d");
    while t < t_end && steps < max_steps {
        let stepped_3d = advance3d_rk2(&mut st, &g, gamma, cfl, true);
        match stepped_3d {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error at step {steps}: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
        if let Some(p) = prog.as_mut() {
            let dt_last = match stepped_3d {
                Ok((dt, _)) => dt,
                Err(_) => 0.0,
            };
            p.tick_3d(
                steps as u64,
                t,
                dt_last,
                &st,
                progress::SampleExtra::default(),
            );
        }
    }
    if ok {
        if let Some(p) = prog.as_mut() {
            p.finish(true, steps as u64, t);
        }
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
        if let Some(p) = prog.as_mut() {
            p.finish(false, steps as u64, t);
        }
        1
    }
}

/// 2d initial condition: uniform field or an over-pressured fireball.
/// the uniform IC under either closure: perfect gas inverts p directly,
/// equilibrium air bisects the fits for e.
fn ic_uniform_energy(rho: f64, p: f64, model: ThermoModel) -> Result<f64, String> {
    match model {
        ThermoModel::Perfect { gamma } => Ok(p / ((gamma - 1.0) * rho)),
        ThermoModel::EqAir { .. } => nufor_core::eqair_energy(&[rho], &[p])
            .map(|v| v[0])
            .map_err(|e| e.to_string()),
        ThermoModel::Park5 { .. } => {
            // uniform park5 ICs start as standard air; e_int carries the
            // formation offset through the mixture closure.
            Ok(nufor_core::species_march::air_ic_energy(rho, p))
        }
        ThermoModel::Combustion9 { .. } => {
            // uniform combustion ICs start as the stoichiometric
            // unburned fill; e_int carries the formation offset.
            Ok(nufor_core::combustion_march::h2o2_ic_energy(rho, p))
        }
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

/// the snapshot writer dispatches on the case's closure so an eqair
/// run's file carries eqair pressure and mach.
fn write_case_vtk2d_model(
    path: &str,
    name: &str,
    g: &Grid2d,
    st: &ConservedState2d,
    model: ThermoModel,
    species: Option<&nufor_core::species::SpeciesState>,
    combustion: Option<&nufor_core::combustion_march::CombustionState>,
) -> i32 {
    let dir = case_dir(path);
    let dst = dir.join(format!("{name}.vtk"));
    if let (Some(sp), true) = (combustion, matches!(model, ThermoModel::Combustion9 { .. })) {
        return match write_vtk2d_combustion(&dst, g, st, sp) {
            Ok(()) => {
                println!("wrote: {}", dst.display());
                0
            }
            Err(e) => {
                eprintln!("output error: {e}");
                1
            }
        };
    }
    if let (Some(sp), true) = (species, matches!(model, ThermoModel::Park5 { .. })) {
        return match write_vtk2d_park5(&dst, g, st, sp) {
            Ok(()) => {
                println!("wrote: {}", dst.display());
                0
            }
            Err(e) => {
                eprintln!("output error: {e}");
                1
            }
        };
    }
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

/// write one scalar field as a standalone structured-points vtk
/// (the two-temperature readout: t_vib per cell).
fn write_scalar_vtk2d(
    path: &Path,
    name: &str,
    g: &Grid2d,
    f: &[f64],
    label: &str,
) -> std::io::Result<()> {
    let n = g.nx * g.ny;
    let mut s = String::with_capacity(n * 24);
    s.push_str("# vtk DataFile Version 3.0\n");
    s.push_str(&format!("{name} {label} snapshot\n"));
    s.push_str("ASCII\nDATASET STRUCTURED_POINTS\n");
    s.push_str(&format!(
        "DIMENSIONS {} {} 1\nORIGIN {} {} 0\nSPACING {} {} 1\n",
        g.nx, g.ny, g.xmin, g.ymin, g.dx, g.dy
    ));
    s.push_str(&format!("POINT_DATA {n}\n"));
    s.push_str(&format!("SCALARS {label} double 1\nLOOKUP_TABLE default\n"));
    for v in f {
        s.push_str(&format!("{v:.6}\n"));
    }
    std::fs::write(path, s)
}

fn min_of(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::INFINITY, f64::min)
}
fn max_of(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

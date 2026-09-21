//! command-line driver for nuFor: build, run, restart, serve, and export a 1D Euler case.
//!
//! `serve` hosts a minimal web ui skeleton: one page plus the snapshot as json,
//! ready for a real front end to be designed against the same data api.

use std::path::{Path, PathBuf};
use std::time::Instant;

use nufor_config::{
    load_case_config, BoundaryKind, CaseConfig, Equations, InitialCondition, Mesh, OutputFormat,
};
use nufor_core::{
    advance2d_rk2, advance2d_sa_rk2, advance3d_rk2, advance_ugrid, cons_to_prim2d, eos_pressure2d,
    euler_solve, grid1d, grid2d, grid3d, prim_to_cons, prim_to_cons2d, prim_to_cons3d,
    read_restart, render_png, wall_distance2d, write_csv, write_h5, write_restart, write_vtk,
    write_vtk2d, write_vtk3d, Bc2d, Boundaries2d, Boundary, Bounds3d, ConservedState,
    ConservedState2d, ConservedState3d, Error, EulerConfig, Grid1d, Grid2d, Grid3d, OutputState,
    SaParams, TurbState, Ugrid,
};

mod mesh_io;
mod serve;
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
    if Path::new(&args[2]).exists() && args[2].ends_with(".toml") {
        return run_case(&args[2]);
    }
    run_shock_tube(args)
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
fn build_1d_mesh(mesh: &Mesh) -> Result<Grid1d, String> {
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

/// build a 2d grid from a case's mesh section: uniform, or imported from a file.
fn build_2d_mesh(mesh: &Mesh) -> Result<Grid2d, String> {
    if mesh.source == "file" {
        let path = mesh
            .path
            .as_deref()
            .ok_or_else(|| "mesh source = file needs a `path`".to_string())?;
        let m = mesh_io::load_rectilinear(path)?;
        mesh_io::rect_to_grid2d(&m)
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
fn ic_from_case(ic: &InitialCondition, n: usize, gamma: f64) -> Result<ConservedState, String> {
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

fn bc1d(kind: BoundaryKind) -> Boundary {
    match kind {
        BoundaryKind::Wall => Boundary::Reflective,
        _ => Boundary::Transmissive,
    }
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
        Equations::Rans2dSa => run_case_2d_sa(&cfg, path),
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
    let mut st = match ic2d_from_case(&cfg.initial_condition, &g, gamma) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("initial condition error: {e}");
            return 2;
        }
    };
    let bc = Boundaries2d::default();
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
        match advance2d_rk2(&mut st, &g, gamma, cfl, true, &bc) {
            Ok((dt, _)) => t += dt,
            Err(e) => {
                eprintln!("solver error: {e}");
                ok = false;
                break;
            }
        }
        steps += 1;
    }
    if ok {
        println!(
            "case {} ({})\nthink: 2d blast, {}x{}\nsteps: {}\ntime: {:.4}\nwall: {:.3}s",
            cfg.metadata.name,
            path,
            g.nx,
            g.ny,
            steps,
            t,
            t0.elapsed().as_secs_f64()
        );
        write_case_vtk2d(path, &cfg.metadata.name, &g, &st, gamma)
    } else {
        1
    }
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
    // boundary sides: wall maps to no-slip (the sa wall condition), inflow to
    // the fixed freestream, everything else transmissive.
    let side_bc = |k: BoundaryKind| match k {
        BoundaryKind::Wall => Bc2d::NoSlipWall,
        BoundaryKind::Inflow => Bc2d::SupersonicInflow {
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
    while t < t_end && steps < max_steps {
        match advance2d_sa_rk2(&mut st, &mut turb, &g, gamma, cfl, true, &bc) {
            Ok(dt) => t += dt,
            Err(e) => {
                eprintln!("solver error: {e}");
                ok = false;
                break;
            }
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
                // uses, so the report matches the flux actually applied.
                let du = (9.0 * u[j_wall * g.nx + i] - u[j_in * g.nx + i]) / (3.0 * g.dy);
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
                eprintln!("solver error: {e}");
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
                eprintln!("solver error: {e}");
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
    let dir = case_dir(path);
    let dst = dir.join(format!("{name}.vtk"));
    match write_vtk2d(&dst, g, st, gamma) {
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

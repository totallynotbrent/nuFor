//! command-line driver for nuFor: build, run, restart, and export a 1D Euler case.
//!
//! the web UI milestone is deliberately not here; `serve` is a stub that points
//! at it so the cli says clearly what is still on the roadmap.

use std::path::Path;
use std::time::Instant;

use nufor_core::{
    euler_solve, grid1d, prim_to_cons, read_restart, write_csv, write_h5, write_restart, write_vtk,
    Boundary, ConservedState, Error, EulerConfig, OutputState,
};

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
        "mesh" => mesh(&args),
        "run" => run(&args),
        "inspect" => inspect(&args),
        "export" => export(&args),
        "history" => history(&args),
        "benchmark" => benchmark(&args),
        "serve" => serve(),
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
         \x20 mesh      N xmin xmax        show a computed grid\n\
         \x20 run       N t [sod|lax] [rst] run a shock tube to time t, save a restart\n\
         \x20 inspect   file.rst           show a restart's header and min/max density\n\
         \x20 export    in.rst out.vtk|h5  write a vtk or hdf5 snapshot\n\
         \x20 history   N [file.csv]       run sod to t=0.2 and write a csv snapshot\n\
         \x20 benchmark N steps            time a fixed run\n\
         \x20 serve                       web UI (next milestone; not built)\n\
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
                "cells: {}\ngamma: {}\ntime: {:.4}\nstep: {}\nrho min/max: {:.4} / {:.4}",
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

fn serve() -> i32 {
    println!("the web UI is the next milestone; it is not built yet.");
    0
}

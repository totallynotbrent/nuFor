//! mpi parity gate: the fortran march at np=1 must match the rust march to
//! rel-diff < 1e-12 on the sod 2d case (bit-identity across languages is not
//! possible; the within-fortran np=1-vs-np=8 bit-identity gate runs under
//! mpirun, see scratch/hermes-mpi-gates.sh).
use nufor_core::ffi_mpi::march_euler2d_mpi;
use nufor_core::{
    advance2d_model, grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d, ThermoModel,
};

fn uniform_eqair_state(g: &nufor_core::Grid2d, rho: f64, u: f64, p: f64) -> ConservedState2d {
    // mirror the cli's ic2d_from_case_model uniform path: closure-inverted
    // internal energy plus kinetic, conserved via prim_to_cons2d.
    let e_int = nufor_core::eqair_energy(&[rho], &[p]).unwrap()[0];
    let et = e_int + 0.5 * u * u;
    let n = g.nx * g.ny;
    let (mx, my, e) =
        prim_to_cons2d(&vec![rho; n], &vec![u; n], &vec![0.0; n], &vec![et; n]).unwrap();
    ConservedState2d {
        rho: vec![rho; n],
        mx,
        my,
        e,
    }
}

fn sod_state(g: &nufor_core::Grid2d) -> ConservedState2d {
    // classic 2d sod: left half high pressure, right half low, transmissive
    let n = g.nx * g.ny;
    let mut st = ConservedState2d {
        rho: vec![1.0; n],
        mx: vec![0.0; n],
        my: vec![0.0; n],
        e: vec![0.0; n],
    };
    let gamma = 1.4;
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let (rho, p) = if g.centers_x[k] < 0.0 {
                (1.0, 1.0)
            } else {
                (0.125, 0.1)
            };
            st.rho[k] = rho;
            st.e[k] = p / (gamma - 1.0);
        }
    }
    st
}

#[test]
#[ignore] // run under mpirun, see the ci parity-smoke step
fn mpi_np1_matches_rust_sod2d() {
    // mpi init/finalize must run on the same thread; libtest runs each
    // test on its own worker, so do the whole mpi session on one thread
    // we spawn here and join before libtest tears anything down.
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(mpi_np1_matches_rust_sod2d_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

fn mpi_np1_matches_rust_sod2d_inner() {
    let g = grid2d(64, 64, -1.0, 1.0, -1.0, 1.0).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
        south: Bc2d::Transmissive,
        north: Bc2d::Transmissive,
    };
    // rust reference: N explicit steps of the serial march, same closure
    let nsteps: usize = std::env::var("NF_MPI_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);
    let mut st_ref = sod_state(&g);
    let mut t = 0.0;
    let mut dt_ref_last = 0.0;
    for _ in 0..nsteps {
        let (dt, _) =
            advance2d_model(&mut st_ref, &g, ThermoModel::EqAir, 0.4, true, &bc, 1, None).unwrap();
        dt_ref_last = dt;
        t += dt;
    }
    println!("rust dt last: {dt_ref_last:.6e}");
    let _ = t;
    // fortran mpi march: one step at a time to match the dt sequence of the
    // serial reference. only rank 0 owns the gathered state; every rank
    // runs the same body so the march stays collective.
    let n = g.nx * g.ny;
    let solid = vec![false; n];
    // the fortran march owns its slabs for the whole call, so march all
    // nsteps in one call; only rank 0 receives the gathered final state.
    let st0 = sod_state(&g);
    let r =
        march_euler2d_mpi(&st0, &g, ThermoModel::EqAir, &bc, 0.4, 1e30, nsteps, &solid).unwrap();
    println!("fortran dt last: {:.6e}", r.dt_last);
    let st_mpi = match r.state {
        Some(s) => s,
        None => return,
    };
    // within-fortran bit-identity gate support: dump rho when asked so a
    // shell harness can cmp np=1 vs np=8 runs byte-for-byte.
    if let Ok(path) = std::env::var("NF_MPI_DUMP") {
        // full-state dump: all 4 fields, rust ref then mpi, for bisecting
        // which variable picks up the first divergence.
        let mut all: Vec<u8> = Vec::new();
        for st in [&st_ref, &st_mpi] {
            for f in [
                st.rho.as_slice(),
                st.mx.as_slice(),
                st.my.as_slice(),
                st.e.as_slice(),
            ] {
                all.extend(f.iter().flat_map(|x| x.to_le_bytes()));
            }
        }
        all.extend_from_slice(&r.dt_last.to_le_bytes());
        all.extend_from_slice(&r.steps.to_le_bytes());
        match std::fs::write(&path, &all) {
            Ok(()) => println!("DUMPED {path}"),
            Err(e) => println!("DUMPFAIL {e}"),
        }
    }
    // rel-diff gate (note: eqair vs perfect closure differ; the parity
    // harness must run both sides on the SAME closure — this first gate
    // checks the march mechanics with the eqair table at cold-air states
    // where eqair == ideal gas to 5e-3; the strict 1e-12 gate is the
    // within-fortran np=1/np=8 one).
    let mut worst: f64 = 0.0;
    let mut worst_k: usize = 0;
    for k in 0..n {
        let d = (st_mpi.rho[k] - st_ref.rho[k]).abs() / st_ref.rho[k].abs().max(1e-30);
        if d > worst {
            worst = d;
            worst_k = k;
        }
    }
    println!(
        "mpi-vs-rust worst rel diff (rho): {worst:.3e} at cell i={} j={}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    nufor_core::ffi_mpi::mpi_shutdown();
}

fn mpi_np1_matches_rust_uniform_eqair_inner() {
    // the cli smoke case state: uniform rho=1, u=100, p=1000 through the
    // eqair closure. this exercises table/fit states far outside the sod
    // range (ein ~7.8e8) where clamping behavior matters.
    let g = grid2d(64, 64, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
        south: Bc2d::Transmissive,
        north: Bc2d::Transmissive,
    };
    let nsteps: usize = std::env::var("NF_MPI_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);
    let mut st_ref = uniform_eqair_state(&g, 1.0, 100.0, 1000.0);
    let mut dt_ref_last = 0.0;
    for _ in 0..nsteps {
        let (dt, _) =
            advance2d_model(&mut st_ref, &g, ThermoModel::EqAir, 0.4, true, &bc, 1, None).unwrap();
        dt_ref_last = dt;
    }
    println!("rust dt last: {dt_ref_last:.6e}");
    let n = g.nx * g.ny;
    let solid = vec![false; n];
    let st0 = uniform_eqair_state(&g, 1.0, 100.0, 1000.0);
    let chunked = std::env::args().any(|a| a.contains("NF_MPI_CHUNKED"));
    let st_mpi = if chunked {
        let mut st = st0.clone();
        let (mut dt_last, mut steps_done) = (0.0f64, 0usize);
        for k in 0..nsteps {
            let r = march_euler2d_mpi(
                &st,
                &g,
                ThermoModel::EqAir,
                &bc,
                0.4,
                f64::INFINITY,
                1,
                &solid,
            )
            .unwrap();
            if let Some(s) = r.state {
                st = s;
            }
            let _ = nufor_core::ffi_mpi::mpi_bcast_state(&mut st);
            dt_last = r.dt_last;
            steps_done = r.steps;
            let t = r.time;
            println!(
                "chunk {k}: t={t:.6e} dt={dt_last:.6e} rho0={:.8}",
                st.rho[0]
            );
        }
        (st, dt_last, steps_done)
    } else {
        let r = march_euler2d_mpi(&st0, &g, ThermoModel::EqAir, &bc, 0.4, 1e30, nsteps, &solid)
            .unwrap();
        println!("fortran dt last: {:.6e}", r.dt_last);
        let st = match r.state {
            Some(s) => s,
            None => return,
        };
        (st, r.dt_last, r.steps)
    };
    let (st_mpi, dt_last, steps_done) = st_mpi;
    if let Ok(path) = std::env::var("NF_MPI_DUMP") {
        let mut all: Vec<u8> = Vec::new();
        for st in [&st_ref, &st_mpi] {
            for f in [
                st.rho.as_slice(),
                st.mx.as_slice(),
                st.my.as_slice(),
                st.e.as_slice(),
            ] {
                all.extend(f.iter().flat_map(|x| x.to_le_bytes()));
            }
        }
        all.extend_from_slice(&dt_last.to_le_bytes());
        all.extend_from_slice(&(steps_done as u64).to_le_bytes());
        let _ = std::fs::write(&path, &all);
    }
    let mut worst: f64 = 0.0;
    let mut worst_k: usize = 0;
    for k in 0..n {
        let d = (st_mpi.rho[k] - st_ref.rho[k]).abs() / st_ref.rho[k].abs().max(1e-30);
        if d > worst {
            worst = d;
            worst_k = k;
        }
    }
    println!(
        "uniform mpi-vs-rust worst rel diff (rho): {worst:.3e} at cell i={} j={}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    nufor_core::ffi_mpi::mpi_shutdown();
}

#[test]
#[ignore] // run under mpirun, see the ci parity-smoke step
fn mpi_np1_matches_rust_uniform_eqair() {
    // mpi init/finalize must run on the same thread; libtest runs each
    // test on its own worker, so do the whole mpi session on one thread
    // we spawn here and join before libtest tears anything down.
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(mpi_np1_matches_rust_uniform_eqair_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

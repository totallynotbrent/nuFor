//! mpi parity gate: the fortran march at np=1 must match the rust march to
//! rel-diff < 1e-12 on the sod 2d case (bit-identity across languages is not
//! possible; the within-fortran np=1-vs-np=8 bit-identity gate runs under
//! mpirun, see scratch/mpi-gates.sh).
use nufor_core::ffi_mpi::{march_euler2d_mpi, march_viscous2d_mpi};
use nufor_core::{
    advance2d_model, advance2d_model_rk2, advance2d_model_visc_rk2, grid2d, prim_to_cons2d, Bc2d,
    Boundaries2d, ConservedState2d, ThermoModel,
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
        let (dt, _) = advance2d_model_rk2(
            &mut st_ref,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            0.4,
            true,
            &bc,
            1,
            None,
        )
        .unwrap();
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
    let r = march_euler2d_mpi(
        &st0,
        &g,
        ThermoModel::EqAir { r_gas: 287.0 },
        &bc,
        0.4,
        1e30,
        nsteps,
        &solid,
    )
    .unwrap();
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
        let (dt, _) = advance2d_model_rk2(
            &mut st_ref,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            0.4,
            true,
            &bc,
            1,
            None,
        )
        .unwrap();
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
                ThermoModel::EqAir { r_gas: 287.0 },
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
        let r = march_euler2d_mpi(
            &st0,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            &bc,
            0.4,
            1e30,
            nsteps,
            &solid,
        )
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

// the laminar viscous gate: the fortran viscous heun march vs the rust
// advance2d_model_visc_rk2 twin on a body-free eqair case with sutherland
// transport. this is the lane the viscous mpi port added; the euler gates
// above never touch the viscous kernels.
fn mpi_np1_visc_matches_rust_inner() {
    let g = grid2d(64, 64, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
        south: Bc2d::Transmissive,
        north: Bc2d::Transmissive,
    };
    let (mu0, t0_ref, s_par, pr, tw) = (1.716e-5, 273.15, 110.4, 0.71, 0.0);
    let nsteps: usize = std::env::var("NF_MPI_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);
    // capsule-regime seed (the eqair fits' intended envelope; the old
    // rho=1.0 seed saturated the energy bracket and made every viscous
    // flux identically zero) with a wall-adjacent perturbation so the
    // gradients, shear, and conduction are actually exercised.
    let (rho0, u0, p0) = (0.004f64, 6960.0f64, 287.0f64);
    let perturb = |st: &mut nufor_core::ConservedState2d| {
        let n = g.nx * g.ny;
        for k in 0..n {
            let (i, j) = (k % g.nx, k / g.nx);
            if j == 0 && i % 7 == 0 {
                st.rho[k] *= 1.05;
                st.mx[k] *= 1.10;
                st.e[k] *= 1.02;
            }
        }
    };
    let mut st_ref = uniform_eqair_state(&g, rho0, u0, p0);
    perturb(&mut st_ref);
    let mut dt_ref_last = 0.0;
    for _ in 0..nsteps {
        let (dt, _) = advance2d_model_visc_rk2(
            &mut st_ref,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            0.4,
            true,
            &bc,
            mu0,
            t0_ref,
            s_par,
            pr,
            tw,
            1,
            None,
        )
        .unwrap();
        dt_ref_last = dt;
    }
    println!("rust visc dt last: {dt_ref_last:.6e}");
    let n = g.nx * g.ny;
    let solid = vec![false; n];
    let mut st0 = uniform_eqair_state(&g, rho0, u0, p0);
    perturb(&mut st0);
    let r = march_viscous2d_mpi(
        &st0,
        &g,
        ThermoModel::EqAir { r_gas: 287.0 },
        &bc,
        0.4,
        1e30,
        nsteps,
        &solid,
        mu0,
        pr,
        tw,
    )
    .unwrap();
    println!("fortran visc dt last: {:.6e}", r.dt_last);
    let st_mpi = match r.state {
        Some(s) => s,
        None => return,
    };
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
        "viscous mpi-vs-rust worst rel diff (rho): {worst:.3e} at cell i={} j={}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    assert!(worst < 1e-6, "viscous mpi-vs-rust drift {worst:.3e}");
    nufor_core::ffi_mpi::mpi_shutdown();
}

#[test]
#[ignore] // run under mpirun, see the ci parity-smoke step
fn mpi_np1_visc_matches_rust() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(mpi_np1_visc_matches_rust_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

// the body-bearing gate: a supersonic freestream over a capsule polygon
// at angle of attack. this is the case class the body-free gates could
// never see: nonzero mx AND my, solid cells, band projection, and tgas1
// blocks beyond the cold region. the chunked mpi march must track the
// serial march (which closes walls through the dist closure) to 1e-6.
fn capsule_body_state(
    g: &nufor_core::Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    normal: &dyn Fn(f64, f64) -> (f64, f64),
) -> ConservedState2d {
    let n = g.nx * g.ny;
    let e_int = nufor_core::eqair_energy(&[0.01], &[685.0]).unwrap()[0];
    let et = e_int + 0.5 * 6691.0 * 6691.0;
    let (mx, my, e) = prim_to_cons2d(
        &vec![0.01; n],
        &vec![6691.0; n],
        &vec![0.0; n],
        &vec![et; n],
    )
    .unwrap();
    let mut st = ConservedState2d {
        rho: vec![0.01; n],
        mx,
        my,
        e,
    };
    nufor_core::apply_solid_fn(&mut st, g, dist, normal, 1.4);
    st
}

fn mpi_np1_matches_rust_capsule_inner() {
    let g = grid2d(96, 56, -1.0, 7.0, -3.5, 3.5).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho: 0.01,
            u: 6691.0,
            v: 0.0,
            p: 685.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SupersonicInflow {
            rho: 0.01,
            u: 6691.0,
            v: 0.0,
            p: 685.0,
        },
        north: Bc2d::SupersonicOutflow,
    };
    // trimmed sphere-cone aeroshell at 22.5 deg, sized to this grid
    let poly = std::sync::Arc::new(nufor_core::SolidPolygon::new(vec![
        (0.9, 0.0),
        (1.05, 0.44),
        (1.5, 0.86),
        (2.1, 1.24),
        (2.85, 1.57),
        (3.7, 1.83),
        (4.6, 1.99),
        (5.35, 2.0),
        (5.9, 1.82),
        (6.2, 1.5),
        (6.3, 1.1),
        (6.3, -1.1),
        (6.2, -1.5),
        (5.9, -1.82),
        (5.35, -2.0),
        (4.6, -1.99),
        (3.7, -1.83),
        (2.85, -1.57),
        (2.1, -1.24),
        (1.5, -0.86),
        (1.05, -0.44),
    ]));
    let pd = poly.clone();
    let pn = poly.clone();
    let dist = move |x: f64, y: f64| pd.signed_distance(x, y);
    let normal = move |x: f64, y: f64| nufor_core::SolidShape::normal(&*pn, x, y);
    let band = 1.5 * g.dx.min(g.dy);
    let mask = nufor_core::body_cache(&g, &dist, band).mask.clone();

    let nsteps: usize = std::env::var("NF_MPI_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    // serial reference: same chunking convention as the cli mpi path
    // (one step, then the band projection between steps)
    let mut st_ref = capsule_body_state(&g, &dist, &normal);
    let mut dt_ref_last = 0.0;
    for _ in 0..nsteps {
        let (dt, _) = advance2d_model(
            &mut st_ref,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            0.4,
            true,
            &bc,
            1,
            Some(&dist),
        )
        .unwrap();
        nufor_core::apply_solid_fn(&mut st_ref, &g, &dist, &normal, 1.4);
        dt_ref_last = dt;
    }
    println!("rust capsule dt last: {dt_ref_last:.6e}");

    let st0 = capsule_body_state(&g, &dist, &normal);
    let mut st = st0.clone();
    let mut dt_last = 0.0;
    for k in 0..nsteps {
        let r = march_euler2d_mpi(
            &st,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            &bc,
            0.4,
            f64::INFINITY,
            1,
            &mask,
        )
        .unwrap();
        if let Some(s) = r.state {
            st = s;
        }
        let _ = nufor_core::ffi_mpi::mpi_bcast_state(&mut st);
        nufor_core::apply_solid_fn(&mut st, &g, &dist, &normal, 1.4);
        dt_last = r.dt_last;
        println!("capsule chunk {k}: dt={:.6e}", r.dt_last);
    }
    println!("fortran capsule dt last: {dt_last:.6e}");

    let n = g.nx * g.ny;
    let mut worst: f64 = 0.0;
    let mut worst_k: usize = 0;
    for k in 0..n {
        let d = (st.rho[k] - st_ref.rho[k]).abs() / st_ref.rho[k].abs().max(1e-30);
        if d > worst {
            worst = d;
            worst_k = k;
        }
    }
    println!(
        "capsule mpi-vs-rust worst rel diff (rho): {worst:.3e} at cell i={} j={}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    assert!(
        worst < 1e-6,
        "capsule body parity broken: worst rel diff {worst:.3e} at cell {}",
        worst_k
    );
    nufor_core::ffi_mpi::mpi_shutdown();
}

#[test]
#[ignore] // run under mpirun, see the ci parity-smoke step
fn mpi_np1_matches_rust_capsule() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(mpi_np1_matches_rust_capsule_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

// the sa rans gate: body-free eqair case with a no-slip south wall so the
// sa transport sees a real wall (wall distance, production, damping).
// the mpi path (np=1) must track the serial advance2d_sa_model_rk2 to
// 1e-6 over n steps on the state AND nu_tilde. non-uniform seed so the
// viscous + transport fluxes are actually exercised.
fn mpi_np1_sa_matches_rust_inner() {
    use nufor_core::ffi_mpi::march_rans2d_mpi;
    use nufor_core::{advance2d_sa_model_rk2, SaParams, TurbState};

    let g = grid2d(48, 32, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
        south: Bc2d::NoSlipWall,
        north: Bc2d::Transmissive,
    };
    let (mu0, pr, pr_t, tw) = (1.716e-5, 0.71, 0.9, 0.0);
    let nsteps: usize = std::env::var("NF_MPI_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let n = g.nx * g.ny;
    // capsule-regime seed (the eqair fits' envelope; rho=1.0 saturates the
    // energy bracket) with a wall-adjacent perturbation so gradients,
    // shear, vorticity, and the transport are all nonzero.
    let (rho0, u0, p0) = (0.004f64, 6960.0f64, 287.0f64);
    let mut st_ref = uniform_eqair_state(&g, rho0, u0, p0);
    let nopert = std::env::var("NF_GATE_NOPERT").is_ok();
    for k in 0..n {
        let (i, j) = (k % g.nx, k / g.nx);
        if !nopert && j == 0 && i % 7 == 0 {
            st_ref.rho[k] *= 1.05;
            st_ref.mx[k] *= 1.10;
            st_ref.e[k] *= 1.02;
        }
    }
    let d = nufor_core::wall_distance2d(&g, &bc, 1.0);
    let nu = mu0 / rho0;
    let mut turb_ref = TurbState {
        nu_tilde: vec![3.0 * nu; n],
        d,
        params: SaParams {
            mu: mu0,
            pr,
            nu_tilde_inf: 3.0 * nu,
            pr_t,
        },
    };
    let mut dt_ref_last = 0.0;
    for _ in 0..nsteps {
        let dt = advance2d_sa_model_rk2(
            &mut st_ref,
            &mut turb_ref,
            &g,
            ThermoModel::EqAir { r_gas: 287.0 },
            0.4,
            true,
            &bc,
            mu0,
            273.15,
            110.4,
            tw,
            1,
            None,
        )
        .unwrap();
        dt_ref_last = dt;
    }
    println!("rust sa dt last: {dt_ref_last:.6e}");

    let mut st0 = uniform_eqair_state(&g, rho0, u0, p0);
    for k in 0..n {
        let (i, j) = (k % g.nx, k / g.nx);
        if !nopert && j == 0 && i % 7 == 0 {
            st0.rho[k] *= 1.05;
            st0.mx[k] *= 1.10;
            st0.e[k] *= 1.02;
        }
    }
    let d2 = nufor_core::wall_distance2d(&g, &bc, 1.0);
    let turb0 = TurbState {
        nu_tilde: vec![3.0 * nu; n],
        d: d2,
        params: SaParams {
            mu: mu0,
            pr,
            nu_tilde_inf: 3.0 * nu,
            pr_t,
        },
    };
    let solid = vec![false; n];
    let (r, nt_out) = march_rans2d_mpi(
        &st0,
        &turb0,
        &g,
        ThermoModel::EqAir { r_gas: 287.0 },
        &bc,
        0.4,
        nsteps,
        &solid,
        mu0,
        pr,
        pr_t,
        tw,
    )
    .unwrap();
    println!("fortran sa dt last: {:.6e}", r.dt_last);
    let st_mpi = match r.state {
        Some(s) => s,
        None => return,
    };
    let nt_mpi = match nt_out {
        Some(v) => v,
        None => return,
    };
    // per-component field-scale denominators: a per-cell ratio blows up on
    // dust values (my at the wall row sits ~1e-9 against a 1e-2 field scale)
    // and would gate noise, not parity.
    let scale = |v: &[f64]| v.iter().map(|x| x.abs()).fold(0.0f64, f64::max).max(1e-30);
    let (s_rho, s_mx, s_my, s_e, s_nt) = (
        scale(&st_ref.rho),
        scale(&st_ref.mx),
        scale(&st_ref.my),
        scale(&st_ref.e),
        scale(&turb_ref.nu_tilde),
    );
    let mut worst: f64 = 0.0;
    let mut worst_k: usize = 0;
    let mut worst_field = "";
    let per_cell: Vec<[f64; 5]> = st_mpi
        .rho
        .iter()
        .zip(&st_ref.rho)
        .zip(&st_mpi.mx)
        .zip(&st_ref.mx)
        .zip(&st_mpi.my)
        .zip(&st_ref.my)
        .zip(&st_mpi.e)
        .zip(&st_ref.e)
        .zip(&nt_mpi)
        .zip(&turb_ref.nu_tilde)
        .map(|(((((((((r0, r1), m0), m1), y0), y1), e0), e1), t0), t1)| {
            [
                (r0 - r1).abs() / s_rho,
                (m0 - m1).abs() / s_mx,
                (y0 - y1).abs() / s_my,
                (e0 - e1).abs() / s_e,
                (t0 - t1).abs() / s_nt,
            ]
        })
        .collect();
    for (k, fields) in per_cell.iter().enumerate() {
        for (name, d) in [
            ("rho", fields[0]),
            ("mx", fields[1]),
            ("my", fields[2]),
            ("e", fields[3]),
            ("nt", fields[4]),
        ] {
            if d > worst {
                worst = d;
                worst_k = k;
                worst_field = name;
            }
        }
    }
    println!(
        "sa mpi-vs-rust worst rel diff: {worst:.3e} on {worst_field} at cell i={} j={}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    {
        let fields: [(&str, Vec<f64>, f64); 5] = [
            (
                "rho",
                st_mpi
                    .rho
                    .iter()
                    .zip(&st_ref.rho)
                    .map(|(a, b)| (a - b).abs())
                    .collect(),
                s_rho,
            ),
            (
                "mx",
                st_mpi
                    .mx
                    .iter()
                    .zip(&st_ref.mx)
                    .map(|(a, b)| (a - b).abs())
                    .collect(),
                s_mx,
            ),
            (
                "my",
                st_mpi
                    .my
                    .iter()
                    .zip(&st_ref.my)
                    .map(|(a, b)| (a - b).abs())
                    .collect(),
                s_my,
            ),
            (
                "e",
                st_mpi
                    .e
                    .iter()
                    .zip(&st_ref.e)
                    .map(|(a, b)| (a - b).abs())
                    .collect(),
                s_e,
            ),
            (
                "nt",
                nt_mpi
                    .iter()
                    .zip(&turb_ref.nu_tilde)
                    .map(|(a, b)| (a - b).abs())
                    .collect(),
                s_nt,
            ),
        ];
        for (name, dv, s) in &fields {
            let (mut km, mut mx) = (0usize, 0.0f64);
            for (k, &x) in dv.iter().enumerate() {
                if x > mx {
                    mx = x;
                    km = k;
                }
            }
            println!(
                "  field {name}: worst abs {mx:.3e} rel {:.3e} at i={} j={}",
                mx / s,
                km % g.nx,
                km / g.nx
            );
        }
    }
    assert!(
        worst < 1e-6,
        "sa mpi-vs-rust drift {worst:.3e} on {worst_field} at cell {}",
        worst_k
    );
    nufor_core::ffi_mpi::mpi_shutdown();
}

#[test]
#[ignore] // run under mpirun, see the ci parity-smoke step
fn mpi_np1_sa_matches_rust() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(mpi_np1_sa_matches_rust_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

//! mpi parity gate: the fortran march at np=1 must match the rust march to
//! rel-diff < 1e-12 on the sod 2d case (bit-identity across languages is not
//! possible; the within-fortran np=1-vs-np=8 bit-identity gate runs under
//! mpirun).
//!
//! run each gate in its own process, the way the ci does. the gates share
//! one mpi communicator, so a process that runs them back to back has the
//! previous gate's collective state live when the next one starts.
//! `tools/regression.sh` runs them one per process.
use nufor_core::ffi_mpi::{march_euler2d_mpi, march_viscous2d_mpi};
use nufor_core::{
    advance2d_model, advance2d_model_rk2, advance2d_model_visc_rk2, grid2d, prim_to_cons2d, Bc2d,
    Boundaries2d, ConservedState2d, ThermoModel,
};

mod common;
use common::mpi_report;

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
        None => {
            // a rank holding no gathered state still has to reach the
            // collective finalize, or the launcher tears the job down
            // while rank 0 is still walking its report to disk.
            mpi_report::fence();
            nufor_core::ffi_mpi::mpi_shutdown();
            return;
        }
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
    mpi_report::number("euler_worst_rel_diff", worst);
    mpi_report::integer("euler_worst_cell_i", (worst_k % g.nx) as i64);
    mpi_report::integer("euler_worst_cell_j", (worst_k / g.nx) as i64);
    mpi_report::number("euler_dt_last_serial", dt_ref_last);
    mpi_report::number("euler_dt_last_mpi", r.dt_last);
    mpi_report::flush();
    mpi_report::fence();
    assert!(
        worst == 0.0,
        "the np=1 mpi march must reproduce the serial march exactly: {worst:.3e}"
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
            None => {
                // this rank has no gathered state, but it still has to
                // reach the fence or the launcher tears the job down
                // while rank 0 walks its report to disk.
                mpi_report::fence();
                nufor_core::ffi_mpi::mpi_shutdown();
                return;
            }
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
    mpi_report::number("uniform_worst_rel_diff", worst);
    mpi_report::integer("uniform_worst_cell_i", (worst_k % g.nx) as i64);
    mpi_report::integer("uniform_worst_cell_j", (worst_k / g.nx) as i64);
    mpi_report::flush();
    mpi_report::fence();
    assert!(
        worst == 0.0,
        "the np=1 mpi march must reproduce the serial march exactly: {worst:.3e}"
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
        None => {
            // a rank holding no gathered state still has to reach the
            // collective finalize, or the launcher tears the job down
            // while rank 0 is still walking its report to disk.
            mpi_report::fence();
            nufor_core::ffi_mpi::mpi_shutdown();
            return;
        }
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
    mpi_report::number("viscous_worst_rel_diff", worst);
    mpi_report::integer("viscous_worst_cell_i", (worst_k % g.nx) as i64);
    mpi_report::integer("viscous_worst_cell_j", (worst_k / g.nx) as i64);
    mpi_report::flush();
    mpi_report::fence();
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
            t: None,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SupersonicInflow {
            rho: 0.01,
            u: 6691.0,
            v: 0.0,
            p: 685.0,
            t: None,
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
    // per-column profile: which columns carry the disagreement, so a
    // seam-local fault and a global one look different.
    let mut col_profile = vec![0.0f64; g.nx];
    for k in 0..n {
        let d = (st.rho[k] - st_ref.rho[k]).abs() / st_ref.rho[k].abs().max(1e-30);
        let i_col = k % g.nx;
        if d > col_profile[i_col] {
            col_profile[i_col] = d;
        }
        if d > worst {
            worst = d;
            worst_k = k;
        }
    }
    let hot_cols: Vec<usize> = col_profile
        .iter()
        .enumerate()
        .filter(|(_, v)| **v > 1.0e-9)
        .map(|(i, _)| i)
        .collect();
    mpi_report::integer("body_columns_above_1e9", hot_cols.len() as i64);
    if let (Some(first), Some(last)) = (hot_cols.first(), hot_cols.last()) {
        mpi_report::integer("body_first_hot_column", *first as i64);
        mpi_report::integer("body_last_hot_column", *last as i64);
        mpi_report::text(
            "body_hot_columns",
            &hot_cols
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
    }
    if let Ok(path) = std::env::var("NF_CAPSULE_COLUMNS") {
        let mut csv = String::from("column,max_rel_diff\n");
        for (i, v) in col_profile.iter().enumerate() {
            csv.push_str(&format!("{i},{v:e}\n"));
        }
        std::fs::write(path, csv).ok();
    }
    println!(
        "capsule mpi-vs-rust worst rel diff (rho): {worst:.3e} at cell i={} j={}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    mpi_report::number("body_worst_rel_diff", worst);
    mpi_report::integer("body_worst_cell_i", (worst_k % g.nx) as i64);
    mpi_report::integer("body_worst_cell_j", (worst_k / g.nx) as i64);
    mpi_report::flush();
    mpi_report::fence();
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
        None => {
            // a rank holding no gathered state still has to reach the
            // collective finalize, or the launcher tears the job down
            // while rank 0 is still walking its report to disk.
            mpi_report::fence();
            nufor_core::ffi_mpi::mpi_shutdown();
            return;
        }
    };
    let nt_mpi = match nt_out {
        Some(v) => v,
        None => {
            // this rank has no gathered plane, but it still has to reach
            // the fence or the launcher tears the job down while rank 0
            // walks its report to disk.
            mpi_report::fence();
            nufor_core::ffi_mpi::mpi_shutdown();
            return;
        }
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
    mpi_report::number("sa_worst_rel_diff", worst);
    mpi_report::text("sa_worst_field", worst_field);
    mpi_report::integer("sa_worst_cell_i", (worst_k % g.nx) as i64);
    mpi_report::integer("sa_worst_cell_j", (worst_k / g.nx) as i64);
    mpi_report::flush();
    mpi_report::fence();
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

/// phase f3 transport gate: every rank carves its x-slab of a recognizable
/// multi-field payload, `gather_fields_mpi` reassembles it on rank 0, and the
/// assembled planes must reproduce the original full-domain fields exactly
/// (mpi doubles carry no arithmetic). non-zero rank returns after the
/// collective. exercises the nf-general gather at np=1 and np=2 in ci.
#[test]
#[ignore] // run under mpirun -np N, see the ci parity-smoke step
fn mpi_gather_nf_fields_reassembles_exactly() {
    use nufor_core::ffi_mpi::{gather_fields_mpi, mpi_shutdown, MpiSession};
    use nufor_core::partition::{assemble_slabs, carve_slab, rank_columns};
    let (nx, ny, nf) = (37usize, 11usize, 3usize);
    let n = nx * ny;
    MpiSession::new().unwrap();
    let (rank, _i0, nx_local) = MpiSession::case_init(nx, ny).unwrap();
    let nx_l = nx_local as usize;
    // full-domain field f: distinguishable values so a transposed slab
    // cannot accidentally pass.
    let full: Vec<Vec<f64>> = (0..nf)
        .map(|f| (0..n).map(|k| (f + 1) as f64 * 1000.0 + k as f64).collect())
        .collect();
    let slab: Vec<f64> = full
        .iter()
        .flat_map(|plane| carve_slab(plane, nx, ny, _i0 as usize, nx_l))
        .collect();
    let gathered = gather_fields_mpi(&slab, nf, nx_l * ny, nx, ny).unwrap();
    if rank == 0 {
        let mut equal = 0i64;
        let mut differing = 0i64;
        for f in 0..nf {
            let plane = &gathered[f * n..(f + 1) * n];
            for (got, want) in plane.iter().zip(full[f].iter()) {
                if got.to_bits() == want.to_bits() {
                    equal += 1;
                } else {
                    differing += 1;
                }
            }
            assert_eq!(plane, &full[f], "gathered field {f} diverged");
        }
        mpi_report::integer("gather_compared_entries", equal + differing);
        mpi_report::integer("gather_equal_entries", equal);
        mpi_report::integer("gather_differing_entries", differing);
        mpi_report::integer("gather_world_size", num_ranks() as i64);
        mpi_report::flush();
        mpi_report::fence();
        // and the same slabs reassemble through the pure-rust path
        let cols = rank_columns(nx, num_ranks());
        let slabs: Vec<Vec<f64>> = (0..num_ranks())
            .scan(0usize, |acc, r| {
                let s = carve_slab(&full[0], nx, ny, *acc, cols[r]);
                *acc += cols[r];
                Some(s)
            })
            .collect();
        assert_eq!(assemble_slabs(&slabs, &cols, nx, ny), full[0]);
    }
    mpi_shutdown();
}

// the two-temperature gate: the fortran 2t lane (heun march + the nufor2t
// plane substep, seam halo included) must reproduce the serial rust lane
// (advance2d_model_rk2 followed by park2t_substep) exactly. the e_v seed is
// a lopsided blob so donor advection pushes mass across the rank seam, and
// the wall-adjacent perturbations give the march real gradients.
#[test]
#[ignore] // run under mpirun -np N, see the ci parity-smoke step
fn mpi_np1_2t_matches_rust() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(mpi_np1_2t_matches_rust_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

fn twot_seed(g: &nufor_core::Grid2d) -> (ConservedState2d, Vec<f64>) {
    // the calm seed the exact-zero euler lanes march on: uniform flow keeps
    // the conserved-state comparison bit-exact, while the lopsided e_v blob
    // below makes the donor advection nonzero so the plane transport is what
    // the gate measures.
    let (rho0, u0, p0) = (1.0f64, 100.0f64, 1000.0f64);
    let st = uniform_eqair_state(g, rho0, u0, p0);
    let n = g.nx * g.ny;
    // vibrational seed: gaussian blobs at several x stations, placed so a
    // blob straddles a rank seam for every split factor the gate runs
    // (np2 cuts at 32, np3 near 21/43, np4 at 16, np5 near 13/26/39/51).
    // each blob is modulated along y so the donor stencil is exercised in
    // both directions.
    let t_v0 = 1500.0f64;
    let e_base = nufor_core::park2t::e_vib_of_t_v(t_v0);
    let centres = [13.0f64, 21.0, 32.0, 43.0, 51.0];
    let sig = 4.0f64;
    let mut ev = vec![0.0f64; n];
    for (k, evk) in ev.iter_mut().enumerate() {
        let i = (k % g.nx) as f64;
        let j = (k / g.nx) as f64;
        let mut amp = 0.0;
        for c in centres {
            let d = (i - c) / sig;
            amp += (-0.5 * d * d).exp();
        }
        *evk = e_base * amp * (0.6 + 0.4 * ((j / g.ny as f64) * std::f64::consts::TAU).cos());
    }
    (st, ev)
}

fn mpi_np1_2t_matches_rust_inner() {
    use nufor_core::ffi_mpi::{march_2t2d_mpi, mpi_shutdown};
    use nufor_core::park2t::{self, Park2tState};
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
        .unwrap_or(6);
    let n = g.nx * g.ny;
    let solid = vec![false; n];
    let (st0, ev0) = twot_seed(&g);

    // rust reference: rk2 march, then the serial 2t substep on the
    // post-march state with the march's dt, per step.
    let mut st_ref = st0.clone();
    let mut park = Park2tState { e_v: ev0.clone() };
    let mut dt_ref_last = 0.0;
    for step_i in 0..nsteps {
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
        .unwrap_or_else(|e| panic!("rk2 reference step {step_i} failed: {:?}", e));
        dt_ref_last = dt;
        park2t::park2t_substep(&mut st_ref, &mut park, &g, 1.4, dt, &solid, false)
            .unwrap_or_else(|e| panic!("park2t_substep step {step_i} failed: {:?}", e));
    }

    // fortran mpi lane: same seed, one call for all steps.
    let (r, ev_out) = march_2t2d_mpi(
        &st0,
        &ev0,
        &g,
        ThermoModel::EqAir { r_gas: 287.0 },
        &bc,
        0.4,
        1e30,
        nsteps,
        &solid,
    )
    .unwrap();
    let ev_mpi = ev_out.expect("the 2t lane returns the plane on every rank");
    let st_mpi = match r.state {
        Some(s) => s,
        None => {
            mpi_report::fence();
            mpi_shutdown();
            return;
        }
    };
    // conserved planes must match the rust march bit for bit (the 2t lane
    // never writes them, same guarantee as the euler lane). counted per
    // field with to_bits so a near-zero momentum cell cannot fake a large
    // relative diff.
    let mut state_worst: f64 = 0.0;
    let mut state_differing = [0i64; 4];
    let fields: [(&str, &Vec<f64>, &Vec<f64>); 4] = [
        ("rho", &st_mpi.rho, &st_ref.rho),
        ("mx", &st_mpi.mx, &st_ref.mx),
        ("my", &st_mpi.my, &st_ref.my),
        ("e", &st_mpi.e, &st_ref.e),
    ];
    for fidx in 0..4 {
        let (_, fa, fb) = &fields[fidx];
        for k in 0..n {
            let (a, b) = (fa[k], fb[k]);
            if a.to_bits() != b.to_bits() {
                state_differing[fidx] += 1;
                if fidx == 0 {
                    let d = (a - b).abs() / b.abs().max(1e-30);
                    if d > state_worst {
                        state_worst = d;
                    }
                }
            }
        }
    }
    println!(
        "2t state differing bits: rho {} mx {} my {} e {} of {n}",
        state_differing[0], state_differing[1], state_differing[2], state_differing[3]
    );
    // the plane: worst relative and absolute disagreement against the
    // serial substep.
    let mut worst_rel: f64 = 0.0;
    let mut worst_abs: f64 = 0.0;
    let mut worst_k = 0usize;
    let mut differing = 0i64;
    for (k, (a, b)) in ev_mpi.iter().zip(park.e_v.iter()).enumerate() {
        let (a, b) = (*a, *b);
        let ad = (a - b).abs();
        let rd = ad / b.abs().max(1e-30);
        if a.to_bits() != b.to_bits() {
            differing += 1;
        }
        if rd > worst_rel {
            worst_rel = rd;
            worst_abs = ad;
            worst_k = k;
        }
    }
    println!(
        "2t plane worst rel diff {worst_rel:.3e} (abs {worst_abs:.3e}) at cell i={} j={}; differing bits {differing}/{n}",
        worst_k % g.nx,
        worst_k / g.nx
    );
    // liveness: the substep has to move the plane, or the parity above
    // would hold for an inert field. the peak value falls as the vibronic
    // pool relaxes, so the largest departure from the seed proves the
    // donor advection and the landau-teller sink both ran.
    let mut liveness: f64 = 0.0;
    for k in 0..n {
        let d = (ev_mpi[k] - ev0[k]).abs();
        if d > liveness {
            liveness = d;
        }
    }
    let seed_peak = ev0.iter().fold(0.0f64, |a, b| a.max(*b));
    let mpi_peak = ev_mpi.iter().fold(0.0f64, |a, b| a.max(*b));
    let ref_peak = park.e_v.iter().fold(0.0f64, |a, b| a.max(*b));
    println!(
        "2t plane liveness: peak change {liveness:.3e} on a seed peak of {seed_peak:.3e} (final mpi peak {mpi_peak:.3e}, serial peak {ref_peak:.3e})"
    );
    mpi_report::number("twot_worst_rel_diff", worst_rel);
    mpi_report::number("twot_worst_abs_diff", worst_abs);
    mpi_report::number("twot_state_worst_rel_diff", state_worst);
    mpi_report::integer("twot_state_differing_rho", state_differing[0]);
    mpi_report::integer("twot_state_differing_mx", state_differing[1]);
    mpi_report::integer("twot_state_differing_my", state_differing[2]);
    mpi_report::integer("twot_state_differing_e", state_differing[3]);
    mpi_report::integer("twot_differing_entries", differing);
    mpi_report::integer("twot_compared_entries", n as i64);
    mpi_report::number("twot_plane_liveness", liveness);
    mpi_report::number("twot_seed_peak", seed_peak);
    mpi_report::integer("twot_world_size", mpi_report::world_size());
    mpi_report::number("twot_dt_last_serial", dt_ref_last);
    mpi_report::number("twot_dt_last_mpi", r.dt_last);
    mpi_report::flush();
    mpi_report::fence();
    assert_eq!(
        state_differing.iter().sum::<i64>(),
        0,
        "the 2t lane disturbed the conserved march: {:?}",
        state_differing
    );
    assert_eq!(
        differing, 0,
        "the mpi 2t plane must match the serial substep bit for bit: {worst_rel:.3e}"
    );
    assert!(
        liveness > 1.0e-9 * seed_peak.max(1.0),
        "the 2t lane did not move the plane (peak change {liveness:.3e}), so the \
         parity above would hold for an inert field"
    );
    mpi_shutdown();
}

/// five-species carrier gate: every rank pushes a distinguishable park5
/// composition slab through the fortran carrier, the seam-face
/// disagreement between neighboring ranks must be bitwise zero at any
/// np, and the gathered global composition must reproduce the seed
/// exactly (mpi doubles carry no arithmetic).
#[test]
#[ignore] // run under mpirun -np N, see the ci parity-smoke step
fn mpi_mpipark5_carrier_seam_and_gather() {
    use nufor_core::ffi_mpi::{
        mpi_shutdown, mpipark5_carrier, mpipark5_gather, mpipark5_init, MpiSession,
    };
    use nufor_core::partition::carve_slab;
    let (nx, ny, ns) = (32usize, 6usize, 5usize);
    let n = nx * ny;
    MpiSession::new().unwrap();
    let (_rank, _i0, nx_local) = MpiSession::case_init(nx, ny).unwrap();
    let (prank, pi0, pnl) = mpipark5_init(nx as i32, ny as i32).unwrap();
    assert_eq!(pnl as usize, nx_local.try_into().unwrap());
    let nx_l = pnl as usize;
    let nn = nx_l * ny;
    // full-domain five planes with distinguishable values per cell.
    let full: Vec<Vec<f64>> = (0..ns)
        .map(|s| {
            (0..n)
                .map(|k| (s + 1) as f64 * 1000.0 + k as f64 * 0.125)
                .collect()
        })
        .collect();
    // the fortran decomposition reports a one-based column origin; rust slices are zero-based.
    let i0_off = pi0 as usize - 1;
    let slab_y: Vec<f64> = full
        .iter()
        .flat_map(|p| carve_slab(p, nx, ny, i0_off, nx_l))
        .collect();
    // the carrier divides by density, so seed a real one; zero momentum
    // keeps the primitive conversion exact.
    let rho_l = vec![1.0f64; nn];
    let zeros = vec![0.0f64; nn];
    let energy = vec![360.0f64; nn];
    let bc_vals = [0.0f64; 20];
    let y_fill = [0.0f64; 5];
    let (rg, _ug, _vg, _eg, yg, seam) = mpipark5_carrier(
        nx_l, ny, &rho_l, &zeros, &zeros, &energy, &slab_y, 0, 0, &bc_vals, &y_fill,
    )
    .unwrap();
    assert_eq!(seam, 0.0, "seam faces disagree at rank {prank}: {seam:e}");
    // the carrier preserves the density it was handed on every interior cell.
    for j in 0..ny {
        for i in 0..nx_l {
            assert_eq!(
                rg[j * (nx_l + 4) + i + 2],
                1.0,
                "carrier altered density at rank {prank}"
            );
        }
    }
    // gather: the padded interior strip reassembles the seeded planes.
    let gathered_opt = mpipark5_gather(nx, ny, &yg).unwrap();
    if prank == 0 {
        let gathered = gathered_opt.expect("rank 0 receives the assembly");
        for s in 0..ns {
            let plane = &gathered[s * n..(s + 1) * n];
            assert_eq!(plane, &full[s], "gathered plane {s} diverged from seed");
        }
        mpi_report::integer("mpipark5_seam_zero", if seam == 0.0 { 1 } else { 0 });
        mpi_report::integer("mpipark5_world_size", num_ranks() as i64);
        mpi_report::flush();
        mpi_report::fence();
    }
    mpi_shutdown();
}

/// world size from the launcher's environment (openmpi/prte stamp every
/// rank's OMPI_COMM_WORLD_SIZE; fall back to 1 outside mpirun).
fn num_ranks() -> usize {
    std::env::var("OMPI_COMM_WORLD_SIZE")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&n: &usize| n >= 1)
        .unwrap_or(1)
}

//! distributed restart gate: a march that is stopped halfway and resumed
//! from a restart file must land on exactly the state the uninterrupted
//! march reaches.
//!
//! the gate covers the defect class that broke the scattered state: the
//! conserved planes travelled to every rank, but the sa turbulence plane did
//! not, so the non-root ranks re-seeded it from its initial value every
//! chunk. the state came out bit-identical while the turbulence field
//! drifted at the rank seam, and a restart therefore reproduced a march that
//! never existed.
//!
//! the march runs in chunks the way the cli drives it: each chunk advances,
//! gathers the state and the turbulence plane, and every rank re-carves the
//! next chunk's slab from the same numbers. the split point falls inside the
//! chunk sequence, so the resumed half starts on a rank seam.
//!
//! run under mpirun, see the ci parity-smoke step:
//!   mpirun -np 2 $(ls -t target/debug/deps/mpi_restart-*) \
//!     --exact mpi_split_restart_reproduces_continuous --nocapture --ignored

use nufor_core::ffi_mpi::{
    march_rans2d_mpi, mpi_bcast_field, mpi_bcast_state, mpi_rank, mpi_shutdown,
};
use nufor_core::{
    grid2d, wall_distance2d, Bc2d, Boundaries2d, ConservedState2d, Restart2dMeta, SaParams,
    ThermoModel, TurbState,
};

mod common;
use common::mpi_report;

/// the capsule-regime seed the eqair closure fits cover.
const RHO0: f64 = 0.004;
const U0: f64 = 6960.0;
const P0: f64 = 287.0;

fn uniform_eqair_state(g: &nufor_core::Grid2d) -> ConservedState2d {
    let n = g.nx * g.ny;
    let e_int = nufor_core::eqair_energy(&[RHO0], &[P0]).unwrap()[0];
    let et = e_int + 0.5 * U0 * U0;
    let (mx, my, e) =
        nufor_core::prim_to_cons2d(&vec![RHO0; n], &vec![U0; n], &vec![0.0; n], &vec![et; n])
            .unwrap();
    ConservedState2d {
        rho: vec![RHO0; n],
        mx,
        my,
        e,
    }
}

fn perturbed_state(g: &nufor_core::Grid2d) -> ConservedState2d {
    let mut st = uniform_eqair_state(g);
    for k in 0..g.nx * g.ny {
        let (i, j) = (k % g.nx, k / g.nx);
        if j == 0 && i % 7 == 0 {
            st.rho[k] *= 1.05;
            st.mx[k] *= 1.10;
            st.e[k] *= 1.02;
        }
    }
    st
}

fn turbidity_state(g: &nufor_core::Grid2d, bc: &Boundaries2d) -> TurbState {
    let n = g.nx * g.ny;
    let (mu0, pr, pr_t) = (1.716e-5, 0.71, 0.9);
    let nu = mu0 / RHO0;
    let d = wall_distance2d(g, bc, 1.0);
    // a non-uniform seed so the transport actually moves the field
    let mut nu_tilde = vec![3.0 * nu; n];
    for (k, value) in nu_tilde.iter_mut().enumerate().take(n) {
        let (i, j) = (k % g.nx, k / g.nx);
        *value *= 1.0 + 0.1 * ((i % 11) as f64) / 10.0 + 0.05 * ((j % 7) as f64) / 7.0;
    }
    TurbState {
        nu_tilde,
        d,
        params: SaParams {
            mu: mu0,
            pr,
            nu_tilde_inf: 3.0 * nu,
            pr_t,
        },
    }
}

/// marches `chunks` chunks of `chunk_len` steps and returns the final state
/// and turbulence plane on rank 0.
fn march_chunks(
    g: &nufor_core::Grid2d,
    bc: &Boundaries2d,
    st0: &ConservedState2d,
    turb0: &TurbState,
    chunk_len: usize,
    chunks: usize,
) -> (ConservedState2d, Vec<f64>) {
    let (mu0, pr, pr_t, tw) = (turb0.params.mu, turb0.params.pr, turb0.params.pr_t, 0.0);
    let solid = vec![false; g.nx * g.ny];
    let mut st = st0.clone();
    let mut turb = turb0.clone();
    let mut t = 0.0f64;
    let mut last_nt: Vec<f64> = turb.nu_tilde.clone();
    for _ in 0..chunks {
        let (_mr, nt_out) = march_rans2d_mpi(
            &st,
            &turb,
            g,
            ThermoModel::EqAir { r_gas: 287.0 },
            bc,
            0.4,
            chunk_len,
            &solid,
            mu0,
            pr,
            pr_t,
            tw,
        )
        .expect("the sa march");
        if let Some(nt) = nt_out {
            turb.nu_tilde.copy_from_slice(&nt);
        }
        let nt = turb.nu_tilde.clone();
        last_nt = nt;
        // rank 0 owns the assembled state; every rank re-carves the next
        // chunk's slab from the same numbers.
        if let Some(gathered) = _mr.state.clone() {
            st = gathered;
        }
        mpi_bcast_state(&mut st).expect("broadcast the state");
        t += _mr.dt_last;
    }
    let _ = t;
    (st, last_nt)
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

#[test]
#[ignore] // run under mpirun, see the ci parity-smoke step
fn mpi_split_restart_reproduces_continuous() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(split_restart_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

fn split_restart_inner() {
    let g = grid2d(48, 32, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
        south: Bc2d::NoSlipWall,
        north: Bc2d::Transmissive,
    };
    let st0 = perturbed_state(&g);
    let turb0 = turbidity_state(&g, &bc);

    // the uninterrupted march: eight chunks of two steps.
    let (cont_st, cont_nt) = march_chunks(&g, &bc, &st0, &turb0, 2, 8);

    // the split march: four chunks, then a restart file, then four chunks
    // resumed from that file.
    let (half_st, half_nt) = march_chunks(&g, &bc, &st0, &turb0, 2, 4);

    // every rank now holds the same full-domain state (the state broadcast
    // at the end of each chunk) and the same global turbulence plane.
    let rank = mpi_rank().unwrap_or(0);

    // the restart file is written once by rank 0 and read by every rank,
    // the way the cli does it. the directory name cannot carry the process
    // id: every rank is its own process, so the ranks would each look in a
    // different place.
    let dir = std::env::var_os("NUFOUR_MPI_RESTART_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("nufor_mpi_restart_gate"));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("half.rst");
    if rank == 0 {
        nufor_core::write_restart2d_aux(
            &path,
            &half_st,
            Restart2dMeta {
                nx: g.nx,
                ny: g.ny,
                gamma: 1.4,
                time: 0.0,
                step: 8,
            },
            &half_nt,
        )
        .expect("write the restart");
    }

    // contract 1: the turbulence plane has to be the SAME on every rank.
    // the conserved state travels to all of them after each chunk, and the
    // turbulence plane has to travel with it; when it does not, the
    // non-root ranks march on a stale field and the plane picked up at the
    // rank seam describes a march that never happened.
    //
    // the broadcast is collective, so it doubles as the fence that
    // publishes rank 0's restart file to the other ranks. it has to come
    // after that write.
    let mut root_nt = half_nt.clone();
    mpi_bcast_field(&mut root_nt).expect("broadcast the reference plane");
    let disagree = max_abs_diff(&root_nt, &half_nt);
    println!("rank {rank}: the turbulence plane differs from rank 0 by {disagree:.3e}",);
    mpi_report::number("restart_plane_disagreement", disagree);
    mpi_report::integer("restart_world_size", mpi_report::world_size());
    assert_eq!(
        disagree, 0.0,
        "rank {rank} marches on a different turbulence plane than rank 0"
    );

    // every rank reads the restart the way the cli does, then marches on.
    let rd = nufor_core::read_restart2d(&path).expect("read the restart");
    assert_eq!(rd.nx, g.nx);
    assert_eq!(rd.ny, g.ny);
    assert_eq!(rd.step, 8);
    let resumed_nt = rd.aux.expect("the restart carries the turbulence plane");
    assert_eq!(resumed_nt, half_nt, "the restart plane came back altered");
    let mut turb_resumed = turbidity_state(&g, &bc);
    turb_resumed.nu_tilde = resumed_nt;

    let (res_st, res_nt) = march_chunks(&g, &bc, &rd.state, &turb_resumed, 2, 4);

    if rank == 0 {
        for (name, a, b) in [
            ("rho", &cont_st.rho, &res_st.rho),
            ("mx", &cont_st.mx, &res_st.mx),
            ("my", &cont_st.my, &res_st.my),
            ("e", &cont_st.e, &res_st.e),
            ("nu_tilde", &cont_nt, &res_nt),
        ] {
            let d = max_abs_diff(a, b);
            println!("restart vs continuous: {name} max abs diff {d:.3e}");
            mpi_report::number(&format!("restart_{}_max_abs_diff", name), d);
            assert!(
                d <= 0.0,
                "the resumed march diverged from the continuous one on {name}: {d:.3e}"
            );
        }
        mpi_report::flush();
        std::fs::remove_dir_all(&dir).ok();
    }
    mpi_report::fence();
    mpi_shutdown();
}

// ---------------------------------------------------------------------------
// two-temperature restart gate: the e_v aux plane has to travel through the
// restart file exactly like the sa nu_tilde plane. a march stopped halfway
// must resume from the file and land bit-identical on the uninterrupted
// march, plane included.
//

fn twot_restart_seed(g: &nufor_core::Grid2d) -> (ConservedState2d, Vec<f64>) {
    // calm uniform seed (bit-exact march territory), blobs at x stations that
    // straddle a rank seam for every split the gate runs.
    let (rho0, u0, p0) = (1.0f64, 100.0f64, 1000.0f64);
    let st = uniform_eqair_state_grid(g, rho0, u0, p0);
    let n = g.nx * g.ny;
    let e_base = nufor_core::park2t::e_vib_of_t_v(1500.0);
    let centres = [13.0f64, 21.0, 32.0, 43.0];
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

fn uniform_eqair_state_grid(
    g: &nufor_core::Grid2d,
    rho0: f64,
    u0: f64,
    p0: f64,
) -> ConservedState2d {
    let n = g.nx * g.ny;
    let e_int = nufor_core::eqair_energy(&[rho0], &[p0]).unwrap()[0];
    let et = e_int + 0.5 * u0 * u0;
    let (mx, my, e) =
        nufor_core::prim_to_cons2d(&vec![rho0; n], &vec![u0; n], &vec![0.0; n], &vec![et; n])
            .unwrap();
    ConservedState2d {
        rho: vec![rho0; n],
        mx,
        my,
        e,
    }
}

/// marches `chunks` chunks of `chunk_len` steps on the 2t lane; every rank
/// ends each chunk holding the same full-domain state and the same global
/// e_v plane, the way the cli drives it.
fn march_2t_chunks(
    g: &nufor_core::Grid2d,
    bc: &Boundaries2d,
    st0: &ConservedState2d,
    ev0: &[f64],
    chunk_len: usize,
    chunks: usize,
) -> (ConservedState2d, Vec<f64>) {
    let solid = vec![false; g.nx * g.ny];
    let mut st = st0.clone();
    let mut ev = ev0.to_vec();
    for _ in 0..chunks {
        let (mr, ev_out) = nufor_core::ffi_mpi::march_2t2d_mpi(
            &st,
            &ev,
            g,
            ThermoModel::EqAir { r_gas: 287.0 },
            bc,
            0.4,
            1e30,
            chunk_len,
            &solid,
        )
        .expect("the 2t march");
        if let Some(s) = mr.state.clone() {
            st = s;
        }
        if let Some(v) = ev_out {
            ev.copy_from_slice(&v);
        }
        mpi_bcast_state(&mut st).expect("broadcast the state");
    }
    (st, ev)
}

#[test]
#[ignore] // run under mpirun -np N, see the ci parity-smoke step
fn mpi_split_restart_2t_reproduces_continuous() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(split_restart_2t_inner)
        .expect("spawn mpi test thread")
        .join()
        .expect("mpi test thread panicked");
}

fn split_restart_2t_inner() {
    use nufor_core::ffi_mpi::{mpi_bcast_field, mpi_rank, mpi_shutdown};
    let g = grid2d(48, 32, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d {
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
        south: Bc2d::Transmissive,
        north: Bc2d::Transmissive,
    };
    let (st0, ev0) = twot_restart_seed(&g);

    // uninterrupted: eight chunks of two steps.
    let (cont_st, cont_ev) = march_2t_chunks(&g, &bc, &st0, &ev0, 2, 8);

    // split: four chunks, a restart file carrying the e_v plane as the aux
    // plane, then four chunks resumed from that file.
    let (half_st, half_ev) = march_2t_chunks(&g, &bc, &st0, &ev0, 2, 4);

    let rank = mpi_rank().unwrap_or(0);
    let dir = std::env::var_os("NUFOUR_MPI_RESTART_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("nufor_mpi_restart_gate"));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("half2t.rst");
    if rank == 0 {
        nufor_core::write_restart2d_aux(
            &path,
            &half_st,
            Restart2dMeta {
                nx: g.nx,
                ny: g.ny,
                gamma: 1.4,
                time: 0.0,
                step: 8,
            },
            &half_ev,
        )
        .expect("write the 2t restart");
    }

    // the broadcast is collective: it fences rank 0's file write to the
    // other ranks and re-confirms every rank holds the same plane.
    let mut root_ev = half_ev.clone();
    mpi_bcast_field(&mut root_ev).expect("broadcast the reference plane");
    let disagree = max_abs_diff(&root_ev, &half_ev);
    mpi_report::number("restart2t_plane_disagreement", disagree);
    mpi_report::integer("restart2t_world_size", mpi_report::world_size());
    assert_eq!(
        disagree, 0.0,
        "rank {rank} marched on a different e_v plane than rank 0"
    );

    let rd = nufor_core::read_restart2d(&path).expect("read the 2t restart");
    assert_eq!(rd.nx, g.nx);
    assert_eq!(rd.ny, g.ny);
    assert_eq!(rd.step, 8);
    let resumed_ev = rd.aux.expect("the restart carries the e_v plane");
    assert_eq!(resumed_ev, half_ev, "the e_v plane came back altered");

    let (res_st, res_ev) = march_2t_chunks(&g, &bc, &rd.state, &resumed_ev, 2, 4);

    if rank == 0 {
        for (name, a, b) in [
            ("rho", &cont_st.rho, &res_st.rho),
            ("mx", &cont_st.mx, &res_st.mx),
            ("my", &cont_st.my, &res_st.my),
            ("e", &cont_st.e, &res_st.e),
            ("e_v", &cont_ev, &res_ev),
        ] {
            let d = max_abs_diff(a, b);
            println!("2t restart vs continuous: {name} max abs diff {d:.3e}");
            mpi_report::number(&format!("restart2t_{}_max_abs_diff", name), d);
            assert!(
                d <= 0.0,
                "the resumed 2t march diverged from the continuous one on {name}: {d:.3e}"
            );
        }
        mpi_report::flush();
        std::fs::remove_dir_all(&dir).ok();
    }
    mpi_report::fence();
    mpi_shutdown();
}

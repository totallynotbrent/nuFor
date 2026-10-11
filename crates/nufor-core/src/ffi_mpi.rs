//! the mpi march bridge: rust loads the case and ships flat arrays to the
//! fortran core, which does decomposition, halo exchange, and the march.
//! every call is one phase boundary; rust never calls mpi directly.

use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;
use crate::state2d::ConservedState2d;
use crate::thermo::ThermoModel;

// the fortran bridge surface (nuforbridge.f90, bind(c))
extern "C" {
    fn nfor_mpi_init(ierr: *mut i32);
    fn nfor_mpi_finalize(ierr: *mut i32);
    fn nfor_mpi_case_init(
        nx: i32,
        ny: i32,
        out_rank: *mut i32,
        out_i0: *mut i32,
        out_nx_local: *mut i32,
        ierr: *mut i32,
    );
    fn nfor_set_tgas1(
        g1: *const f64,
        g2: *const f64,
        g3: *const f64,
        g4: *const f64,
        g5: *const f64,
        g6: *const f64,
        g7: *const f64,
        g8: *const f64,
        s: *const f64,
        minus: *const i32,
        nb: i32,
        zb: *const f64,
        cg: *const f64,
        r0: f64,
        e0: f64,
    );
    fn nfor_set_cea(nr: i32, ne: i32, lr: *const f64, le: *const f64, t: *const f64, p: *const f64);
    fn nfor_mpi_step(
        rho: *const f64,
        mx: *const f64,
        my: *const f64,
        e: *const f64,
        nloc: i32,
        nx_l: i32,
        ny: i32,
        dx: f64,
        dy: f64,
        gamma: f64,
        cfl: f64,
        dt_cap: f64,
        bc_w: i32,
        bc_e: i32,
        bc_s: i32,
        bc_n: i32,
        bc_vals: *const f64,
        solid: *const i32,
        rho2: *mut f64,
        mx2: *mut f64,
        my2: *mut f64,
        e2: *mut f64,
        dt_out: *mut f64,
        ierr: *mut i32,
    );
    fn nfor_mpi_rank(rank: *mut i32);
    fn nfor_mpi_bcast4(buf: *mut f64, n4: i32, ierr: *mut i32);
    fn nfor_mpi_gather(
        rho: *const f64,
        mx: *const f64,
        my: *const f64,
        e: *const f64,
        nloc: i32,
        nx: i32,
        ny: i32,
        out: *mut f64,
        ierr: *mut i32,
    );
}

/// one mpi session per process: init once on first use, finalized at exit
/// through the process-lifetime leak (mpi forbids re-init after finalize, so
/// the session deliberately never drops).
pub struct MpiSession;

static INITED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

impl MpiSession {
    /// initialize mpi exactly once per process (idempotent).
    pub fn new() -> Result<Self, Error> {
        static ONCE: std::sync::Once = std::sync::Once::new();
        ONCE.call_once(|| {
            let mut ierr: i32 = 0;
            unsafe { nfor_mpi_init(&mut ierr) };
            INITED.store(ierr == 0, std::sync::atomic::Ordering::Release);
        });
        if !INITED.load(std::sync::atomic::Ordering::Acquire) {
            return Err(Error::InvalidArgs);
        }
        Ok(Self)
    }

    /// build this rank's decomposition; returns (rank, global x offset,
    /// local slab width).
    pub fn case_init(nx: usize, ny: usize) -> Result<(i32, i32, i32), Error> {
        let mut rank: i32 = 0;
        let mut i0: i32 = 0;
        let mut nx_local: i32 = 0;
        let mut ierr: i32 = 0;
        unsafe {
            nfor_mpi_case_init(
                nx as i32,
                ny as i32,
                &mut rank,
                &mut i0,
                &mut nx_local,
                &mut ierr,
            );
        }
        if ierr != 0 {
            return Err(Error::InvalidArgs);
        }
        Ok((rank, i0, nx_local))
    }
}

/// ship the tgas1 face-fit tables (once per process, before marching).
pub fn ship_tgas1_tables() {
    use crate::eqair::tgas1_tables;
    let tb = tgas1_tables();
    let flat2 = |a: &[(f64, f64)]| -> Vec<f64> { a.iter().flat_map(|(x, y)| [*x, *y]).collect() };
    let flat3 = |a: &[(f64, f64, f64)]| -> Vec<f64> {
        a.iter().flat_map(|(x, y, z)| [*x, *y, *z]).collect()
    };
    let flat2pad =
        |a: &[(f64, f64)]| -> Vec<f64> { a.iter().flat_map(|(x, y)| [*x, *y, 0.0]).collect() };
    let flat4 = |a: &[(f64, f64, f64, f64)]| -> Vec<f64> {
        a.iter().flat_map(|(w, x, y, z)| [*w, *x, *y, *z]).collect()
    };
    let g1 = flat2(&tb.g1);
    let g5 = flat2(&tb.g5);
    let g2 = flat2pad(&tb.g2);
    let g3 = flat3(&tb.g3);
    let g4 = flat3(&tb.g4);
    let g6 = flat2pad(&tb.g6);
    let g7 = flat3(&tb.g7);
    let g8 = flat3(&tb.g8);
    let s = flat4(&tb.s);
    let minus: Vec<i32> = tb.minus.iter().map(|&m| m as i32).collect();
    let zb: Vec<f64> = tb.z_breaks.iter().flat_map(|r| r.iter().copied()).collect();
    unsafe {
        nfor_set_tgas1(
            g1.as_ptr(),
            g2.as_ptr(),
            g3.as_ptr(),
            g4.as_ptr(),
            g5.as_ptr(),
            g6.as_ptr(),
            g7.as_ptr(),
            g8.as_ptr(),
            s.as_ptr(),
            minus.as_ptr(),
            12,
            zb.as_ptr(),
            tb.cold_gamm.as_ptr(),
            tb.rho0,
            tb.e0,
        );
    }
}

/// ship the cea table to the fortran side (once per process, before marching).
pub fn ship_cea_table() {
    // the rust table lives in eqair_cea as a private OnceLock; re-derive the
    // flat axes+values from the same include! source by calling its public
    // surface: p_t_at over every node would be roundabout. instead parse the
    // same tsv here (identical file, identical numbers).
    let raw = include_str!("../data/cea_eqair_table.tsv");
    let lines: Vec<&str> = raw.lines().collect();
    let first = lines[0].split_whitespace().collect::<Vec<_>>();
    let nr: usize = first[0].parse().unwrap_or(40);
    let ne: usize = first[1].parse().unwrap_or(60);
    let mut lr = Vec::with_capacity(nr);
    let mut le = Vec::with_capacity(ne);
    for s in lines[1].split_whitespace() {
        lr.push(s.parse::<f64>().unwrap_or(0.0));
    }
    for s in lines[2].split_whitespace() {
        le.push(s.parse::<f64>().unwrap_or(0.0));
    }
    let mut t = Vec::with_capacity(nr * ne);
    let mut p = Vec::with_capacity(nr * ne);
    // replicate the eqair_cea parser exactly: it starts at line index 4
    // (0-based), dropping the first data line, and pads any shortfall with
    // the (300.0, 101325.0) default. bit-identity with the rust table
    // requires the same quirk.
    for ln in lines.iter().skip(4).take(nr * ne) {
        let mut it = ln.split_whitespace();
        t.push(it.next().and_then(|s| s.parse().ok()).unwrap_or(300.0));
        p.push(it.next().and_then(|s| s.parse().ok()).unwrap_or(101325.0));
    }
    while t.len() < nr * ne {
        t.push(300.0);
        p.push(101325.0);
    }
    unsafe {
        nfor_set_cea(
            nr as i32,
            ne as i32,
            lr.as_ptr(),
            le.as_ptr(),
            t.as_ptr(),
            p.as_ptr(),
        );
    }
}

/// the mpi march driver for one euler_2d case: runs the full march on the
/// fortran side, gathering snapshots to rank 0 at the requested intervals.
/// only rank 0 returns output; other ranks return an empty log.
pub struct MpiMarchResult {
    /// rank 0 only: the final gathered state (full domain, row-major).
    pub state: Option<ConservedState2d>,
    /// rank 0 only: the run-log line fields.
    pub steps: usize,
    pub time: f64,
    pub wall: f64,
    pub dt_last: f64,
}

// the step surface mirrors advance2d_model plus the mpi gather contract.
#[allow(clippy::too_many_arguments)]
pub fn march_euler2d_mpi(
    state0: &ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    bc: &Boundaries2d,
    cfl: f64,
    final_time: f64,
    max_steps: usize,
    solid: &[bool],
) -> Result<MpiMarchResult, Error> {
    if !matches!(model, ThermoModel::EqAir) {
        // v1 gates the mpi march to the eqair capsule path; other closures
        // stay on the rust march until phase F2.
        return Err(Error::InvalidArgs);
    }
    // the process-global mpi session: init once, never finalize (the march
    // may be invoked repeatedly and mpi_init after finalize aborts).
    MpiSession::new()?;
    let (rank, i0, nx_local_i) = MpiSession::case_init(g.nx, g.ny).map_err(|e| {
        eprintln!("dbg: case_init failed nx={} ny={}", g.nx, g.ny);
        e
    })?;
    let nx_local = nx_local_i as usize;
    // ship the closure tables once per process: the march may be invoked
    // repeatedly (cli chunked mode) and re-shipping reallocates the fortran
    // cea allocatables every call.
    static SHIPPED: std::sync::Once = std::sync::Once::new();
    SHIPPED.call_once(|| {
        ship_cea_table();
        ship_tgas1_tables();
    });

    let n = g.nx * g.ny;
    // carve this rank's slab from the full-domain state (row-major j*nx+i)
    let mut rho_l = vec![0.0f64; nx_local * g.ny];
    let mut mx_l = vec![0.0f64; nx_local * g.ny];
    let mut my_l = vec![0.0f64; nx_local * g.ny];
    let mut e_l = vec![0.0f64; nx_local * g.ny];
    for j in 0..g.ny {
        for i in 0..nx_local {
            let k_src = j * g.nx + i0 as usize + i;
            let k_dst = j * nx_local + i;
            rho_l[k_dst] = state0.rho[k_src];
            mx_l[k_dst] = state0.mx[k_src];
            my_l[k_dst] = state0.my[k_src];
            e_l[k_dst] = state0.e[k_src];
        }
    }
    let solid_i32: Vec<i32> = {
        let mut s = vec![0i32; nx_local * g.ny];
        for j in 0..g.ny {
            for i in 0..nx_local {
                s[j * nx_local + i] = solid[j * g.nx + i0 as usize + i] as i32;
            }
        }
        s
    };

    // bc codes for fortran: 0=transmissive 1=inflow 2=outflow 3=slip 4=noslip.
    // each side ships [rho, u, v, p, ein_inflow]; ein is the closure-inverted
    // internal energy of the inflow state (ghost seeding for var 5) and is
    // only read when the side's code is inflow.
    let bc_code = |b: &crate::solver2d::Bc2d| -> (i32, [f64; 5]) {
        use crate::solver2d::Bc2d::*;
        match b {
            Transmissive | SupersonicOutflow => (0, [0.0; 5]),
            SupersonicInflow { rho, u, v, p } => {
                let ein = crate::eqair::energy_from_pressure(&[*rho], &[*p])
                    .ok()
                    .and_then(|v| v.first().copied())
                    .unwrap_or(0.0);
                (1, [*rho, *u, *v, *p, ein])
            }
            SlipWall => (3, [0.0; 5]),
            NoSlipWall => (4, [0.0; 5]),
            ProfileInflow { .. } => (1, [0.0; 5]),
        }
    };
    let (bw, bvw) = bc_code(&bc.west);
    let (be, bve) = bc_code(&bc.east);
    let (bs, bvs) = bc_code(&bc.south);
    let (bn, bvn) = bc_code(&bc.north);
    let _ = (bvw, bve, bvs, bvn); // v1: capsule bcs are w/s inflow, e/n outflow
    let bc_vals = [bvw, bve, bvs, bvn].concat();

    let dt_cap = final_time;
    let mut t: f64 = 0.0;
    let mut steps: usize = 0;
    let mut dt_last: f64 = 0.0;
    let wall_start = std::time::Instant::now();
    let mut ierr: i32 = 0;
    let mut rho2 = vec![0.0f64; rho_l.len()];
    let mut mx2 = vec![0.0f64; mx_l.len()];
    let mut my2 = vec![0.0f64; my_l.len()];
    let mut e2 = vec![0.0f64; e_l.len()];

    while t < final_time && steps < max_steps {
        let mut dt: f64 = 0.0;
        unsafe {
            nfor_mpi_step(
                rho_l.as_ptr(),
                mx_l.as_ptr(),
                my_l.as_ptr(),
                e_l.as_ptr(),
                rho_l.len() as i32,
                nx_local as i32,
                g.ny as i32,
                g.dx,
                g.dy,
                1.4,
                cfl,
                dt_cap,
                bw,
                be,
                bs,
                bn,
                bc_vals.as_ptr(),
                solid_i32.as_ptr(),
                rho2.as_mut_ptr(),
                mx2.as_mut_ptr(),
                my2.as_mut_ptr(),
                e2.as_mut_ptr(),
                &mut dt,
                &mut ierr,
            );
        }
        if ierr != 0 {
            eprintln!("dbg: nfor_mpi_step ierr={ierr} at step {steps}");
            return Err(Error::InvalidArgs);
        }
        std::mem::swap(&mut rho_l, &mut rho2);
        std::mem::swap(&mut mx_l, &mut mx2);
        std::mem::swap(&mut my_l, &mut my2);
        std::mem::swap(&mut e_l, &mut e2);
        t += dt;
        dt_last = dt;
        steps += 1;
    }
    let wall = wall_start.elapsed().as_secs_f64();

    // gather to rank 0 and rebuild the full-domain state there
    let mut gathered = vec![0.0f64; 4 * n];
    unsafe {
        nfor_mpi_gather(
            rho_l.as_ptr(),
            mx_l.as_ptr(),
            my_l.as_ptr(),
            e_l.as_ptr(),
            rho_l.len() as i32,
            g.nx as i32,
            g.ny as i32,
            gathered.as_mut_ptr(),
            &mut ierr,
        );
    }
    if ierr != 0 {
        return Err(Error::InvalidArgs);
    }
    let state = if rank == 0 {
        let mut st = ConservedState2d {
            rho: vec![0.0; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![0.0; n],
        };
        for k in 0..n {
            st.rho[k] = gathered[k * 4];
            st.mx[k] = gathered[k * 4 + 1];
            st.my[k] = gathered[k * 4 + 2];
            st.e[k] = gathered[k * 4 + 3];
        }
        Some(st)
    } else {
        None
    };
    Ok(MpiMarchResult {
        state,
        steps,
        time: t,
        wall,
        dt_last,
    })
}

/// whether this process was launched under mpirun with more than one rank
/// (env-only probe; never initializes mpi).
pub fn mpi_launched() -> bool {
    let size = std::env::var("OMPI_COMM_WORLD_SIZE")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1);
    size > 1 || std::env::var("NFOR_MPI").ok().as_deref() == Some("1")
}

/// rank 0 owns output: writes happen there only. initializes the session
/// on first use (safe under mpirun; plain launches return 0/true).
pub fn mpi_rank() -> Result<usize, Error> {
    if !mpi_launched() {
        return Ok(0);
    }
    MpiSession::new()?;
    let mut rank: i32 = 0;
    unsafe { nfor_mpi_rank(&mut rank) };
    Ok(rank.max(0) as usize)
}

/// release the session before process exit: mpirun treats a missing
/// finalize as abnormal termination and exits nonzero even when every
/// rank's work succeeded. collective, idempotent, and a no-op when the
/// session never initialized.
pub fn mpi_shutdown() {
    static DONE: std::sync::Once = std::sync::Once::new();
    DONE.call_once(|| {
        if mpi_launched() && INITED.load(std::sync::atomic::Ordering::Acquire) {
            let mut ierr: i32 = 0;
            unsafe { nfor_mpi_finalize(&mut ierr) };
        }
    });
}

/// chunked-body mode: rank 0 holds the projected full state; broadcast it
/// so every rank re-carves its slab from the same numbers.
pub fn mpi_bcast_state(st: &mut ConservedState2d) -> Result<(), Error> {
    MpiSession::new()?;
    let n = st.rho.len();
    let n4 = (4 * n) as i32;
    // one packed broadcast instead of four: four 32kib messages at the
    // ucx eager/rendezvous threshold corrupted neighbouring heap buffers
    // on this stack (openmpi 5.0.2).
    let mut buf: Vec<f64> = Vec::with_capacity(4 * n);
    let (r, mx, my, e) = (&st.rho, &st.mx, &st.my, &st.e);
    for k in 0..n {
        buf.push(r[k]);
        buf.push(mx[k]);
        buf.push(my[k]);
        buf.push(e[k]);
    }
    let mut ierr: i32 = 0;
    unsafe {
        nfor_mpi_bcast4(buf.as_mut_ptr(), n4, &mut ierr);
    }
    if ierr != 0 {
        return Err(Error::InvalidArgs);
    }
    for k in 0..n {
        st.rho[k] = buf[4 * k];
        st.mx[k] = buf[4 * k + 1];
        st.my[k] = buf[4 * k + 2];
        st.e[k] = buf[4 * k + 3];
    }
    Ok(())
}

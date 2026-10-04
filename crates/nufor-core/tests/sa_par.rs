//! the threaded sa march contract: identical physics at any thread count,
//! measured speedup on the plate geometry.

use nufor_core::{
    advance2d_sa_rk2, advance2d_sa_rk2_par, stretched_grid2d, Bc2d, BlasiusProfile, Boundaries2d,
    Clustering, ConservedState2d, InflowProfile, SaParams, StepConfig, TurbState,
};
use nufor_core::{prim_to_cons2d, wall_distance2d};

const GAMMA: f64 = 1.4;
const U_INF: f64 = 0.2;

/// the turbulent plate setup the validation uses, at a small size.
fn plate(
    nx: usize,
    ny: usize,
) -> (
    nufor_core::Grid2d,
    ConservedState2d,
    Boundaries2d,
    TurbState,
) {
    let g = stretched_grid2d(
        nx,
        0.0,
        1.2,
        ny,
        0.0,
        0.35,
        Clustering {
            first_cell: 5.18e-5,
            growth: 1.2,
        },
    )
    .unwrap();
    let n = nx * ny;
    let mu = 4.8e-7;
    let nu_tilde_inf = 3.0 * mu;
    let et = vec![1.0 / (GAMMA - 1.0) + 0.5 * U_INF * U_INF; n];
    let (mx, my, e) = prim_to_cons2d(&vec![1.0; n], &vec![U_INF; n], &vec![0.0; n], &et).unwrap();
    let st = ConservedState2d {
        rho: vec![1.0; n],
        mx,
        my,
        e,
    };
    let p = BlasiusProfile::new(U_INF, mu, -0.2)
        .unwrap()
        .anchored_at(0.0)
        .unwrap();
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: std::sync::Arc::new(InflowProfile::Blasius(p)),
            rho: 1.0,
            p: 1.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::NoSlipWall,
        north: Bc2d::SupersonicOutflow,
    };
    let turb = TurbState {
        nu_tilde: vec![nu_tilde_inf; n],
        d: wall_distance2d(&g, &bc, 0.35),
        params: SaParams {
            mu,
            pr: 0.72,
            pr_t: 0.9,
            nu_tilde_inf,
        },
    };
    (g, st, bc, turb)
}

#[test]
fn threaded_sa_march_is_bit_identical_to_serial() {
    // the same physical march at several thread counts must reproduce the
    // serial result exactly, cell for cell, in the mean flow and nu_tilde.
    let steps = 40;
    for &(nx, ny) in &[(48, 24), (96, 48)] {
        let ref_state = {
            let (g, st, bc, mut turb) = plate(nx, ny);
            let mut s = st;
            for _ in 0..steps {
                advance2d_sa_rk2(&mut s, &mut turb, &g, GAMMA, 0.3, true, &bc).unwrap();
            }
            (s, turb)
        };
        for &nt in &[2usize, 4, 6, 8] {
            let (g, st, bc, mut turb) = plate(nx, ny);
            let mut s = st;
            let cfg = StepConfig {
                gamma: GAMMA,
                cfl: 0.3,
                muscl: true,
                nthreads: nt,
            };
            for _ in 0..steps {
                advance2d_sa_rk2_par(&mut s, &mut turb, &g, &bc, cfg).unwrap();
            }
            assert_eq!(s.rho.len(), ref_state.0.rho.len());
            for k in 0..s.rho.len() {
                assert_eq!(
                    s.rho[k].to_bits(),
                    ref_state.0.rho[k].to_bits(),
                    "rho @ {k}"
                );
                assert_eq!(s.mx[k].to_bits(), ref_state.0.mx[k].to_bits(), "mx @ {k}");
                assert_eq!(s.my[k].to_bits(), ref_state.0.my[k].to_bits(), "my @ {k}");
                assert_eq!(s.e[k].to_bits(), ref_state.0.e[k].to_bits(), "e @ {k}");
                assert_eq!(
                    turb.nu_tilde[k].to_bits(),
                    ref_state.1.nu_tilde[k].to_bits(),
                    "nu_tilde @ {k}"
                );
            }
        }
    }
}

#[test]
#[ignore] // a timing run, not a correctness gate; invoke with -- --ignored
fn threaded_sa_speedup() {
    use std::time::Instant;
    let steps = 200;
    let (g, st, bc, mut turb) = plate(160, 60);
    let mut s = st.clone();
    let t0 = Instant::now();
    for _ in 0..steps {
        advance2d_sa_rk2(&mut s, &mut turb, &g, GAMMA, 0.3, true, &bc).unwrap();
    }
    let serial = t0.elapsed();
    for &nt in &[2usize, 4, 6, 8] {
        let (g2, st2, bc2, mut turb2) = plate(160, 60);
        let mut s2 = st2;
        let cfg = StepConfig {
            gamma: GAMMA,
            cfl: 0.3,
            muscl: true,
            nthreads: nt,
        };
        let t0 = Instant::now();
        for _ in 0..steps {
            advance2d_sa_rk2_par(&mut s2, &mut turb2, &g2, &bc2, cfg).unwrap();
        }
        let par = t0.elapsed();
        println!(
            "threads {nt}: serial {serial:?} par {par:?} speedup {:.2}x",
            serial.as_secs_f64() / par.as_secs_f64()
        );
    }
}

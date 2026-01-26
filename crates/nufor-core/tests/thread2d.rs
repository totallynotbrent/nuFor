//! threaded 2d step: the parallel result must be bit-identical to the serial one.

use nufor_core::{
    advance2d, advance2d_par, grid2d, prim_to_cons2d, Boundaries2d, ConservedState2d, Grid2d,
};

const GAMMA: f64 = 1.4;

fn sample_state(g: &Grid2d) -> ConservedState2d {
    // a smooth blob in pressure over a uniform base, safe and non-uniform.
    let n = g.nx * g.ny;
    let mut rho = vec![1.0; n];
    let mut u = vec![0.0; n];
    let v = vec![0.0; n];
    let mut p = vec![1.0; n];
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        for i in 0..g.nx {
            let x = g.centers_x[j * g.nx + i];
            let k = j * g.nx + i;
            p[k] = 1.0 + 0.4 * (-((x - 0.5).powi(2) + (y - 0.5).powi(2)) * 40.0).exp();
            u[k] = 0.1 * y;
            rho[k] = 1.0 + 0.1 * x;
        }
    }
    let et: Vec<f64> = p
        .iter()
        .zip(&u)
        .zip(&v)
        .map(|((pp, uu), vv)| *pp / (rho[0] * (GAMMA - 1.0)) + 0.5 * (uu * uu + vv * vv))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    ConservedState2d { rho, mx, my, e }
}

#[test]
fn parallel_matches_serial_exactly() {
    let g: Grid2d = grid2d(64, 64, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d::default();
    let mut serial = sample_state(&g);
    for _ in 0..5 {
        advance2d(&mut serial, &g, GAMMA, 0.5, true, &bc).unwrap();
    }
    for nthreads in [1usize, 2, 4, 8] {
        let mut par = sample_state(&g);
        for _ in 0..5 {
            advance2d_par(&mut par, &g, GAMMA, 0.5, true, &bc, nthreads).unwrap();
        }
        assert!(par.rho == serial.rho, "rho diverged at {nthreads} threads");
        assert!(par.mx == serial.mx, "mx diverged at {nthreads} threads");
        assert!(par.my == serial.my, "my diverged at {nthreads} threads");
        assert!(par.e == serial.e, "energy diverged at {nthreads} threads");
    }
}

//! 3d structured grid, state conversion, and ideal-gas eos.

use nufor_core::{
    cons_to_prim3d, eos_mach3d, eos_pressure3d, eos_sound_speed3d, grid3d, prim_to_cons3d,
    Bounds3d, Grid3d,
};

const GAMMA: f64 = 1.4;

#[test]
fn grid3d_has_the_right_geometry() {
    let g: Grid3d = grid3d(
        2,
        3,
        4,
        &Bounds3d {
            xmin: 0.0,
            xmax: 1.0,
            ymin: 0.0,
            ymax: 2.0,
            zmin: 0.0,
            zmax: 3.0,
        },
    )
    .unwrap();
    assert_eq!((g.nx, g.ny, g.nz), (2, 3, 4));
    assert!((g.dx - 0.5).abs() < 1e-12, "dx {}", g.dx);
    assert!((g.dy - 2.0 / 3.0).abs() < 1e-12);
    assert!((g.dz - 0.75).abs() < 1e-9);
    assert_eq!(g.centers_x.len(), 2 * 3 * 4);
    assert_eq!(g.faces_x.len(), 3);
    assert_eq!(g.faces_y.len(), 4);
    assert_eq!(g.faces_z.len(), 5);
    // the (0,0,0) cell is flat index 0: x=0.25, y=1/3, z=0.375.
    let k = 0;
    assert!((g.centers_x[k] - 0.25).abs() < 1e-12);
    assert!((g.centers_y[k] - 1.0 / 3.0).abs() < 1e-12);
    assert!((g.centers_z[k] - 0.375).abs() < 1e-9);
    // the far corner index (k=nz-1, j=ny-1, i=nx-1).
    let last = (g.nz - 1) * g.ny * g.nx + (g.ny - 1) * g.nx + (g.nx - 1);
    assert!((g.centers_x[last] - 0.75).abs() < 1e-12);
    assert!((g.centers_z[last] - 2.625).abs() < 1e-12);
}

#[test]
fn prim_cons_3d_round_trips_and_rejects_bad_density() {
    let n = 6;
    let rho = vec![1.0; n];
    let (u, v, w, et) = (vec![0.5; n], vec![-0.25; n], vec![0.1; n], vec![3.0; n]);
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &u, &v, &w, &et).unwrap();
    assert!(
        (mx[0] - 0.5).abs() < 1e-12 && (mz[0] - 0.1).abs() < 1e-12,
        "rho*u = mx"
    );
    let (u2, v2, w2, et2) = cons_to_prim3d(&rho, &mx, &my, &mz, &e).unwrap();
    for i in 0..n {
        assert!((u2[i] - u[i]).abs() < 1e-12, "u round trip");
        assert!((v2[i] - v[i]).abs() < 1e-12, "v round trip");
        assert!((w2[i] - w[i]).abs() < 1e-12, "w round trip");
        assert!((et2[i] - et[i]).abs() < 1e-12, "et round trip");
    }
    let bad = vec![0.0; n];
    let _ = bad;
    let result = prim_to_cons3d(&bad, &u, &v, &w, &et);
    assert!(result.is_err(), "zero density is rejected");
}

#[test]
fn eos3d_pressure_sound_and_mach_are_sane() {
    let n = 4;
    let rho = vec![1.0; n];
    let (u, v, w) = (vec![0.3; n], vec![0.0; n], vec![0.4; n]);
    let et = vec![2.5; n];
    let p = eos_pressure3d(GAMMA, &rho, &et, &u, &v, &w).unwrap();
    // p = (gamma-1) rho (et - 0.5(u^2+v^2+w^2)).
    let expect = 0.4 * (2.5 - 0.5 * (0.09 + 0.16));
    assert!((p[0] - expect).abs() < 1e-12, "p {}", p[0]);
    let a = eos_sound_speed3d(GAMMA, &p, &rho).unwrap();
    assert!((a[0] - (GAMMA * expect).sqrt()).abs() < 1e-9);
    let m = eos_mach3d(&u, &v, &w, &a).unwrap();
    let speed = (0.09f64 + 0.16).sqrt();
    assert!((m[0] - speed / a[0]).abs() < 1e-9);
}

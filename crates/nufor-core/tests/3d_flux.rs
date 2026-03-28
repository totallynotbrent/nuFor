//! 3d hllc flux sanity + the sod tube reproduced in three dimensions.

use nufor_core::{
    advance3d_rk2, grid3d, hllc_flux3, prim_to_cons3d, Bounds3d, ConservedState3d, FacePrim3,
    Grid3d,
};

const GAMMA: f64 = 1.4;

#[test]
fn hllc3d_uniform_flow_gives_the_exact_flux() {
    // uniform state (rho=1, u=0.5, v=-0.1, w=0.2, p=1) through an x-face.
    let l = FacePrim3 {
        rho: 1.0,
        u: 0.5,
        v: -0.1,
        w: 0.2,
        p: 1.0,
    };
    let f = hllc_flux3(GAMMA, l, l, 0);
    let e_t = 1.0 / (1.0 * (GAMMA - 1.0)) + 0.5 * (0.25 + 0.01 + 0.04);
    assert!((f.mass - 0.5).abs() < 1e-12, "mass flux {}", f.mass);
    assert!((f.mx - (0.25 + 1.0)).abs() < 1e-12, "mx flux {}", f.mx);
    // my = rho*u*v, mz = rho*u*w, e = u*(rho e_t + p).
    assert!((f.my - (0.5 * -0.1)).abs() < 1e-12);
    assert!((f.mz - (0.5 * 0.2)).abs() < 1e-12);
    assert!(
        (f.e - 0.5 * (e_t + 1.0)).abs() < 1e-9,
        "energy flux {}",
        f.e
    );
}

#[test]
fn sod_tube_reproduced_in_three_dimensions() {
    // a slice in x with uniform y and z: a pure 1d problem in a 3d box.
    let nx = 100;
    let (ny, nz) = (2, 2);
    let cell = nx * ny * nz;
    let g: Grid3d = grid3d(
        nx,
        ny,
        nz,
        &Bounds3d {
            xmin: 0.0,
            xmax: 1.0,
            ymin: 0.0,
            ymax: 0.5,
            zmin: 0.0,
            zmax: 0.5,
        },
    )
    .unwrap();
    let split = nx / 2;
    let mut rho = vec![0.0; cell];
    let mut pr = vec![0.0; cell];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let c = (k * ny + j) * nx + i;
                let left = i < split;
                rho[c] = if left { 1.0 } else { 0.125 };
                pr[c] = if left { 1.0 } else { 0.1 };
            }
        }
    }
    let (u, v, w) = (vec![0.0; cell], vec![0.0; cell], vec![0.0; cell]);
    let et: Vec<f64> = pr
        .iter()
        .zip(&rho)
        .map(|(p, r)| *p / (r * (GAMMA - 1.0)))
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &u, &v, &w, &et).unwrap();
    let mut st = ConservedState3d { rho, mx, my, mz, e };
    let mut t = 0.0;
    while t < 0.2 {
        let (dt, _) = advance3d_rk2(&mut st, &g, GAMMA, 0.5, true).unwrap();
        t += dt;
    }
    assert!(
        st.rho.iter().all(|&r| r.is_finite() && r > 0.0),
        "density stayed positive"
    );
    // the x-profile through the (j=0,k=0) slice should be the classic sod.
    let (mut rmax, mut rmin) = (0.0f64, f64::INFINITY);
    for i in 0..nx {
        let c = i; // (k=0, j=0) slice
        rmax = rmax.max(st.rho[c]);
        rmin = rmin.min(st.rho[c]);
    }
    assert!(
        (rmax - 1.0).abs() < 0.02,
        "left state preserved, got {rmax}"
    );
    assert!(
        (rmin - 0.125).abs() < 0.02,
        "right state preserved, got {rmin}"
    );
    // the star plateau between rarefaction tail and contact sits near rho=0.426
    // around x=0.6 at t=0.2 (the (k=0,j=0) slice is flat index i).
    let pb = st.rho[60];
    eprintln!("sod plateau at i=60: {pb:.4} (expect ~0.426)");
    assert!(
        (pb - 0.426).abs() < 0.05,
        "star plateau density, got {pb:.4}"
    );
    eprintln!("sod 3d: rmax {rmax:.4} rmin {rmin:.4}");
}

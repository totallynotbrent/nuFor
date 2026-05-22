//! unstructured finite-volume prototype: the sod tube on a quad mesh.

use nufor_core::{advance_ugrid, prim_to_cons2d, ConservedState2d, Ugrid};

const GAMMA: f64 = 1.4;

#[test]
fn sod_tube_reproduced_on_the_unstructured_path() {
    // a 100 x 4 cartesian quad "unstructured" mesh; the strip is uniform in y.
    let ug: Ugrid = Ugrid::cartesian_quads(100, 4);
    let cell = ug.cell_area.len();
    let mut rho = vec![0.0; cell];
    let mut pr = vec![1.0; cell];
    for k in 0..cell {
        let i = k % 100;
        rho[k] = if i < 50 { 1.0 } else { 0.125 };
        pr[k] = if i < 50 { 1.0 } else { 0.1 };
    }
    let (u, v) = (vec![0.0; cell], vec![0.0; cell]);
    let et: Vec<f64> = pr
        .iter()
        .zip(&rho)
        .map(|(p, r)| *p / (r * (GAMMA - 1.0)))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let mut t = 0.0;
    while t < 0.2 {
        let (dt, _) = advance_ugrid(&mut st, &ug, GAMMA, 0.3).unwrap();
        t += dt;
    }
    assert!(
        st.rho.iter().all(|&r| r.is_finite() && r > 0.0),
        "density stayed positive"
    );
    // sweep the x-profile along one y-row (cells 0..100 are row j=0).
    let (mut rmax, mut rmin) = (0.0f64, f64::INFINITY);
    for i in 0..100 {
        rmax = rmax.max(st.rho[i]);
        rmin = rmin.min(st.rho[i]);
    }
    assert!(
        (rmax - 1.0).abs() < 0.02,
        "left state preserved, got {rmax}"
    );
    assert!(
        (rmin - 0.125).abs() < 0.02,
        "right state preserved, got {rmin}"
    );
    // the star plateau near x=0.6 should read close to the classic 0.426.
    let pb = st.rho[58];
    eprintln!("ugrid sod plateau at i=58: {pb:.4} (expect ~0.426)");
    assert!(
        (pb - 0.426).abs() < 0.08,
        "star plateau density, got {pb:.4}"
    );
    eprintln!("ugrid sod: rmax {rmax:.4} rmin {rmin:.4}");
}

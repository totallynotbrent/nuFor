//! geometry + freestream contracts for the curvilinear body-fitted
//! grid. the mapping must be valid before any physics runs on it:
//! positive areas, face normals consistent, the inner boundary ON the
//! body, and a uniform eqair freestream held by the march.

use nufor_core::curvilinear::{advance2d_axi_curv, CurvGrid, Freestream};
use nufor_core::ConservedState2d;
use nufor_core::SphereCone;
use nufor_core::ThermoModel;

fn capsule() -> SphereCone {
    SphereCone {
        rn: 1.9,
        delta: 9.5663f64.to_radians(),
        rb: 2.365,
        xc: 3.9,
    }
}

#[test]
fn grid_geometry_is_valid() {
    let sc = capsule();
    let g = CurvGrid::over_sphere_cone(&sc, 96, 48, 6.0, 2.0);
    assert_eq!(g.cells.len(), 96 * 48);
    // every cell has positive area
    for c in &g.cells {
        assert!(c.area > 0.0, "non-positive cell area");
        assert!(c.r >= -1e-9, "cell below the axis");
    }
    // the inner boundary rows sit on the body meridian
    for i in 0..=g.nx {
        let (bx, br) = nufor_core::curvilinear::body_point(&sc, i as f64 / g.nx as f64);
        let (sx, sr) = (g.surface_x[i], g.surface_r[i]);
        assert!((sx - bx).abs() < 1e-9, "surface x mismatch at {i}");
        assert!((sr - br).abs() < 1e-9, "surface r mismatch at {i}");
    }
    // face normals are unit
    for f in &g.xfaces {
        assert!((f.nx * f.nx + f.nr * f.nr - 1.0).abs() < 1e-12);
        assert!(f.len > 0.0);
    }
    for f in &g.yfaces {
        assert!((f.nx * f.nx + f.nr * f.nr - 1.0).abs() < 1e-12);
        assert!(f.len > 0.0);
    }
}

#[test]
fn uniform_eqair_freestream_holds() {
    let sc = capsule();
    let (nx, ny) = (48, 24);
    let g = CurvGrid::over_sphere_cone(&sc, nx, ny, 6.0, 2.0);
    // a uniform M22 eqair freestream (the corrected ghost contract)
    let (rho, u, p) = (0.004f64, 6960.0f64, 287.0f64);
    let e_int = nufor_core::eqair_energy(&[rho], &[p]).unwrap()[0];
    let et = e_int + 0.5 * u * u;
    let n = nx * ny;
    let mut st = ConservedState2d {
        rho: vec![rho; n],
        mx: vec![rho * u; n],
        my: vec![0.0; n],
        e: vec![rho * et; n],
    };
    let model = ThermoModel::EqAir;
    for _ in 0..40 {
        let fs = Freestream { rho, u, p };
        let (_, resid) = advance2d_axi_curv(&mut st, &g, model, 0.2, fs).unwrap();
        assert!(resid.is_finite());
    }
    // a uniform axial velocity field is only a steady state far from
    // the body: near the nose the annular streamtube narrows and the
    // march correctly develops the compression. the freestream
    // contract applies to the outer rows, where the grid relaxes to
    // the far field.
    for j in (ny - 6)..ny {
        for i in (nx / 2)..nx {
            let k = j * nx + i;
            let dr = (st.rho[k] - rho).abs() / rho;
            assert!(dr < 0.01, "outer density drifted {dr:.3} at ({i},{j})");
        }
    }
}

#[test]
fn capsule_ogrid_marches() {
    let sc = capsule();
    let (nx, ny) = (96, 48);
    let g = CurvGrid::over_sphere_cone(&sc, nx, ny, 6.0, 2.0);
    let (rho, u, p) = (0.004f64, 6960.0f64, 287.0f64);
    let e_int = nufor_core::eqair_energy(&[rho], &[p]).unwrap()[0];
    let et = e_int + 0.5 * u * u;
    let n = nx * ny;
    let mut st = ConservedState2d {
        rho: vec![rho; n],
        mx: vec![rho * u; n],
        my: vec![0.0; n],
        e: vec![rho * et; n],
    };
    let fs = Freestream { rho, u, p };
    for step in 0..200 {
        {
            // pre-scan: cells whose internal energy went negative
            let mut bad = 0;
            for k in 0..n {
                let uu = st.mx[k] / st.rho[k];
                let vv = st.my[k] / st.rho[k];
                let ett = st.e[k] / st.rho[k];
                if ett - 0.5 * (uu * uu + vv * vv) <= 0.0 {
                    if bad < 3 {
                        eprintln!(
                            "neg e_int cell {k}: rho={:.3e} u={:.1} v={:.1} et={:.3e}",
                            st.rho[k], uu, vv, ett
                        );
                    }
                    bad += 1;
                }
            }
            if bad > 0 {
                eprintln!("step {step}: {bad} cells with e_int <= 0");
            }
        }
        match advance2d_axi_curv(&mut st, &g, ThermoModel::EqAir, 0.4, fs) {
            Ok((dt, resid)) => {
                assert!(resid.is_finite(), "resid NaN at {step}");
                if step % 20 == 0 {
                    eprintln!("step {step} t={:.5} resid={:.3e}", step as f64 * dt, resid);
                }
            }
            Err(e) => {
                let rr = st.rho.iter().cloned().fold(f64::INFINITY, f64::min);
                let ee = st.e.iter().cloned().fold(f64::INFINITY, f64::min);
                panic!("march err {e} at step {step}: min rho={rr:.5} min e={ee:.1}");
            }
        }
    }
}

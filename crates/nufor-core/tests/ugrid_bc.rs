//! gates for the boundary-group dispatch: a uniform flow on a tagged
//! mesh must hold exactly (wall reflection, inflow, outflow all
//! consistent with the interior), and the group bookkeeping must tag
//! only boundary faces.

use nufor_core::{advance_ugrid, prim_to_cons2d, ConservedState2d, Ugrid};

fn uniform_state(nc: usize, rho: f64, u: f64, v: f64, p: f64, gamma: f64) -> ConservedState2d {
    let et: Vec<f64> = vec![p / (rho * (gamma - 1.0)); nc];
    let (mx, my, e) = prim_to_cons2d(&vec![rho; nc], &vec![u; nc], &vec![v; nc], &et).unwrap();
    ConservedState2d {
        rho: vec![rho; nc],
        mx,
        my,
        e,
    }
}

/// a 4x4 quad mesh with every boundary face tagged wall, freestream set
/// to the interior: after one step the state must be bit-identical.
#[test]
fn uniform_flow_with_all_wall_groups_holds() {
    let mut ug = Ugrid::cartesian_quads(4, 4);
    // tag every boundary face as wall (group 1)
    let mut segs = std::collections::HashMap::new();
    for f in 0..ug.face_left.len() {
        if ug.face_right[f] < 0 {
            let (a, b) = ug.face_nodes[f];
            segs.insert((a, b), 1u32);
        }
    }
    ug.tag_boundary_faces(&segs);
    // quiescent gas: every flux row is zero at a wall, and the pressure
    // push cancels around each closed cell, so the state must hold
    // exactly. a nonzero wind would pile up against the walls, which
    // is real physics, not a consistency defect.
    let (rho, u, v, p, gamma) = (0.5, 0.0, 0.0, 2.0e5, 1.4);
    let mut st = uniform_state(16, rho, u, v, p, gamma);
    let before = st.clone();
    advance_ugrid(&mut st, &ug, gamma, 0.3, nufor_core::ThermoModel::Perfect { gamma }).unwrap();
    for c in 0..16 {
        assert!(
            (st.rho[c] - before.rho[c]).abs() < 1e-12,
            "rho drifted cell {c}: {} vs {}",
            st.rho[c],
            before.rho[c]
        );
        assert!(
            (st.mx[c] - before.mx[c]).abs() < 1e-9,
            "mx drifted cell {c}"
        );
        assert!(
            (st.my[c] - before.my[c]).abs() < 1e-9,
            "my drifted cell {c}"
        );
        assert!((st.e[c] - before.e[c]).abs() < 1e-6, "e drifted cell {c}");
    }
}

/// inflow and outflow groups set to the interior state must also hold a
/// uniform flow exactly, on a mesh whose walls are a mix of groups.
#[test]
fn uniform_flow_with_mixed_groups_holds() {
    let mut ug = Ugrid::cartesian_quads(6, 3);
    // west column = inflow (2), east column = outflow (3), top+bottom =
    // wall (1): tag by face position via the owning cell's center.
    let mut segs = std::collections::HashMap::new();
    for f in 0..ug.face_left.len() {
        if ug.face_right[f] >= 0 {
            continue;
        }
        let (a, b) = ug.face_nodes[f];
        // owning cell center tells us which wall this is
        let li = ug.face_left[f];
        let (cx, _) = ug.cell_center[li];
        let g = if cx < 0.2 {
            2
        } else if cx > 0.8 {
            3
        } else {
            1
        };
        segs.insert((a, b), g as u32);
    }
    ug.tag_boundary_faces(&segs);
    // freestream inflow matches the interior and the gas is quiescent
    // so the walls deflect nothing: exact hold.
    let (rho, u, v, p, gamma) = (1.2, 0.0, 0.0, 9.5e4, 1.4);
    ug.freestream = [rho, u, v, p];
    let mut st = uniform_state(18, rho, u, v, p, gamma);
    let before = st.clone();
    advance_ugrid(&mut st, &ug, gamma, 0.4, nufor_core::ThermoModel::Perfect { gamma }).unwrap();
    for c in 0..18 {
        assert!((st.rho[c] - before.rho[c]).abs() < 1e-12, "rho drifted");
        assert!((st.mx[c] - before.mx[c]).abs() < 1e-9, "mx drifted");
        assert!((st.my[c] - before.my[c]).abs() < 1e-9, "my drifted");
        assert!((st.e[c] - before.e[c]).abs() < 1e-6, "e drifted");
    }
}

/// a freestream inflow face against a different interior must push the
/// interior toward freestream, not the other way: cheap direction check.
#[test]
fn inflow_group_drives_interior_not_vice_versa() {
    let mut ug = Ugrid::cartesian_quads(4, 4);
    let mut segs = std::collections::HashMap::new();
    for f in 0..ug.face_left.len() {
        if ug.face_right[f] >= 0 {
            continue;
        }
        let (a, b) = ug.face_nodes[f];
        segs.insert((a, b), 2u32);
    }
    ug.tag_boundary_faces(&segs);
    ug.freestream = [2.0, 100.0, 0.0, 5.0e5];
    // interior starts as still air at a different pressure
    let mut st = uniform_state(16, 1.0, 0.0, 0.0, 1.0e5, 1.4);
    let (dt, _) = advance_ugrid(&mut st, &ug, 1.4, 0.2, nufor_core::ThermoModel::Perfect { gamma: 1.4 }).unwrap();
    assert!(dt > 0.0);
    // momentum must increase in +x somewhere (freestream blowing in)
    let mx_sum: f64 = st.mx.iter().sum();
    assert!(mx_sum > 0.0, "inflow did not push +x momentum in");
}

/// a uniform supersonic wind in an open duct (inflow left, outflow
/// right, slip walls top and bottom) must stay uniform: this is the
/// consistency property that makes body-fitted external flows work.
#[test]
fn uniform_wind_in_duct_holds() {
    let mut ug = Ugrid::cartesian_quads(8, 3);
    let mut segs = std::collections::HashMap::new();
    for f in 0..ug.face_left.len() {
        if ug.face_right[f] >= 0 {
            continue;
        }
        let (a, b) = ug.face_nodes[f];
        let li = ug.face_left[f];
        let (cx, _) = ug.cell_center[li];
        let g = if cx < 0.2 {
            2
        } else if cx > 0.8 {
            3
        } else {
            1
        };
        segs.insert((a, b), g as u32);
    }
    ug.tag_boundary_faces(&segs);
    let (rho, u, v, p, gamma) = (0.002, 2000.0, 0.0, 500.0, 1.4);
    ug.freestream = [rho, u, v, p];
    let nc = 24;
    let et_i = p / (rho * (gamma - 1.0));
    let et = vec![et_i + 0.5 * u * u; nc];
    let (mx, my, e) = prim_to_cons2d(&vec![rho; nc], &vec![u; nc], &vec![v; nc], &et).unwrap();
    let mut st = ConservedState2d {
        rho: vec![rho; nc],
        mx,
        my,
        e,
    };
    let before = st.clone();
    for _ in 0..200 {
        advance_ugrid(&mut st, &ug, gamma, 0.3, nufor_core::ThermoModel::Perfect { gamma }).unwrap();
    }
    let mut worst = 0.0f64;
    for c in 0..nc {
        worst = worst
            .max((st.rho[c] - before.rho[c]).abs())
            .max((st.mx[c] - before.mx[c]).abs())
            .max((st.e[c] - before.e[c]).abs());
    }
    println!("duct worst drift after 200 steps: {worst:.3e}");
    assert!(worst < 1e-7, "uniform wind drifted {worst}");
}

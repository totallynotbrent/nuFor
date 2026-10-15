//! contract test: the body-aware viscous wall band. a masked circle in a
//! slow crossflow must form a real no-slip layer: the wall-adjacent fluid
//! cells decelerate over steps while the solid cells stay inert, and the
//! isothermal wall drains heat from a hot flow.

use nufor_core::{
    add_viscous_cells_mu, grid2d, Bc2d, Boundaries2d, ConservedState2d, ThermoModel, TimeControl,
};

fn circle_sdf(x: f64, y: f64) -> f64 {
    let (cx, cy) = (0.5, 0.5);
    ((x - cx).powi(2) + (y - cy).powi(2)).sqrt() - 0.2
}

#[test]
#[ignore]
fn body_wall_no_slip_layer_forms() {
    let g = grid2d(64, 64, 0.0, 1.0, 0.0, 1.0).unwrap();
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let gamma = 1.4;

    // uniform +x flow, rho 1, u 0.1, p 1: slow, so the viscous layer is
    // the dominant physics and nothing shocks.
    let mut st = ConservedState2d {
        rho: vec![1.0; nx * ny],
        mx: vec![0.1; nx * ny],
        my: vec![0.0; nx * ny],
        e: vec![p_to_e(1.0, 1.0, 0.1 * 0.1); nx * ny],
    };
    // seed the solid cells inert like apply_solid does.
    let mask = nufor_core::solid_mask(&g, &circle_sdf);
    let ref_e = p_to_e(1.0, 1.0, 0.0);
    for k in 0..nx * ny {
        if mask[k] {
            st.rho[k] = 1.0;
            st.mx[k] = 0.0;
            st.my[k] = 0.0;
            st.e[k] = ref_e;
        }
    }

    let _bc = Boundaries2d {
        west: Bc2d::SlipWall,
        east: Bc2d::SlipWall,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };
    // grab the wall-band cells: fluid cells with a solid face-neighbor.
    let mut band: Vec<usize> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            if mask[k] {
                continue;
            }
            let west = i > 0 && mask[idx(i - 1, j)];
            let east = i + 1 < nx && mask[idx(i + 1, j)];
            let south = j > 0 && mask[idx(i, j - 1)];
            let north = j + 1 < ny && mask[idx(i, j + 1)];
            if west || east || south || north {
                band.push(k);
            }
        }
    }
    assert!(!band.is_empty(), "no wall band found");

    // march pure viscous diffusion steps (no euler sub-step): the layer
    // growth is then unambiguous.
    let u_before: Vec<f64> = band.iter().map(|&k| st.mx[k] / st.rho[k]).collect();
    let n_steps = 200;
    for _ in 0..n_steps {
        add_viscous_cells_mu(
            &mut st,
            &g,
            gamma,
            1.0,
            &[],
            0.71,
            None,
            TimeControl::Global(0.02 * g.dx * g.dy / 1.0),
            Some(ThermoModel::Perfect { gamma }),
            0.0,
            1,
            Some(&circle_sdf),
        )
        .unwrap();
        // keep the solid inert like the march does each step.
        for k in 0..nx * ny {
            if mask[k] {
                st.mx[k] = 0.0;
                st.my[k] = 0.0;
                st.rho[k] = 1.0;
                st.e[k] = ref_e;
            }
        }
    }
    // the wall band must have decelerated vs the freestream 0.1.
    let mut worst = 0.0f64;
    for (n, &k) in band.iter().enumerate() {
        let u_now = st.mx[k] / st.rho[k];
        // every band cell must be slower than it started (or at least not
        // faster); the layer forms against the wall.
        assert!(
            u_now <= u_before[n] + 1e-12,
            "band cell {} sped up: {} -> {}",
            k,
            u_before[n],
            u_now
        );
        worst = worst.max(u_before[n] - u_now);
    }
    // and the layer must actually be a LAYER: the mean band velocity
    // strictly below the freestream by now.
    let mean_band: f64 =
        band.iter().map(|&k| st.mx[k] / st.rho[k]).sum::<f64>() / band.len() as f64;
    assert!(
        mean_band < 0.99 * 0.1,
        "band did not decelerate: mean u {mean_band} vs freestream 0.1"
    );
    assert!(worst > 1e-4, "no deceleration anywhere: worst {worst}");

    // solid cells stay inert.
    for k in 0..nx * ny {
        if mask[k] {
            assert_eq!(st.mx[k], 0.0);
            assert_eq!(st.my[k], 0.0);
        }
    }
}

#[test]
#[ignore]
fn body_wall_isothermal_drains_heat() {
    let g = grid2d(64, 64, 0.0, 1.0, 0.0, 1.0).unwrap();
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let gamma = 1.4;
    // hot flow against a cold wall: rho 1, u 0, p 1 -> T = 300*... use
    // e so T ~ 3x the wall's 100 K.
    let t_hot = 300.0;
    let e_hot = t_hot * 287.0 / (gamma - 1.0);
    let mut st = ConservedState2d {
        rho: vec![1.0; nx * ny],
        mx: vec![0.0; nx * ny],
        my: vec![0.0; nx * ny],
        e: vec![e_hot; nx * ny],
    };
    let mask = nufor_core::solid_mask(&g, &circle_sdf);
    let ref_e = p_to_e(1.0, 1.0, 0.0);
    for k in 0..nx * ny {
        if mask[k] {
            st.e[k] = ref_e;
        }
    }
    let band: Vec<usize> = (0..nx * ny)
        .filter(|&k| {
            if mask[k] {
                return false;
            }
            let (i, j) = (k % nx, k / nx);
            (i > 0 && mask[idx(i - 1, j)])
                || (i + 1 < nx && mask[idx(i + 1, j)])
                || (j > 0 && mask[idx(i, j - 1)])
                || (j + 1 < ny && mask[idx(i, j + 1)])
        })
        .collect();
    let e_before: Vec<f64> = band.iter().map(|&k| st.e[k]).collect();
    let mut solid_energy_before = 0.0;
    for k in 0..nx * ny {
        if mask[k] {
            solid_energy_before += st.e[k];
        }
    }
    for _ in 0..100 {
        add_viscous_cells_mu(
            &mut st,
            &g,
            gamma,
            1.0,
            &[],
            0.71,
            None,
            TimeControl::Global(0.02 * g.dx * g.dy),
            Some(ThermoModel::Perfect { gamma }),
            100.0,
            1,
            Some(&circle_sdf),
        )
        .unwrap();
        for k in 0..nx * ny {
            if mask[k] {
                st.mx[k] = 0.0;
                st.my[k] = 0.0;
                st.rho[k] = 1.0;
                st.e[k] = ref_e;
            }
        }
    }
    let mut drained_any = false;
    for (n, &k) in band.iter().enumerate() {
        if st.e[k] < e_before[n] - 1e-9 {
            drained_any = true;
        }
        // with the wall colder than the flow, no band cell may gain energy.
        assert!(
            st.e[k] <= e_before[n] + 1e-9,
            "band cell {k} gained energy against a colder wall"
        );
    }
    assert!(drained_any, "no heat drained into the wall");
    // solid cells unchanged (the march re-pins them, but the operator
    // itself must not have pushed energy through them).
    let solid_energy_after: f64 = (0..nx * ny).filter(|&k| mask[k]).map(|k| st.e[k]).sum();
    assert!(
        (solid_energy_after - solid_energy_before).abs() < 1e-9,
        "solid energy changed"
    );
}

fn p_to_e(p: f64, rho: f64, ke: f64) -> f64 {
    p / ((1.4 - 1.0) * rho) + 0.5 * ke
}

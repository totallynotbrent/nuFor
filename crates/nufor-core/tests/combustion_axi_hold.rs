//! the axisymmetric combustion lane holds a moving state, and the
//! species transport keeps a uniform composition uniform in a duct.

use nufor_core::combustion_march::{advance2d_combustion_axi, CombustionState};
use nufor_core::combustion_thermo as ct;
use nufor_core::{Bc2d, Boundaries2d, ConservedState2d};

/// the burned chamber fill at cea's 70 bar equilibrium.
fn chamber_fill() -> [f64; ct::NS] {
    [
        0.01583304, 0.07194470, 0.77605355, 0.11605824, 0.00232724, 0.01719592, 0.00049223,
        0.00009298, 0.0,
    ]
}

/// a uniform moving state must survive the march. the combustion
/// closure's internal energy carries the nasa-7 formation offset, so a
/// positivity guard that tests it against zero brakes the flow on
/// every step.
#[test]
fn axi_combustion_holds_a_moving_state() {
    let (nx, ny) = (8, 8);
    let g = nufor_core::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
    let n = nx * ny;
    let fill = chamber_fill();
    let mut sp = CombustionState::with_fill(n, &fill);
    let rho = 3.6116f64;
    let p = 7.0e6f64;
    let u = 200.0f64;
    let e_int = nufor_core::combustion_march::h2o2_ic_energy_at(rho, p, &fill);
    let et = e_int + 0.5 * u * u;
    let mut state = ConservedState2d {
        rho: vec![rho; n],
        mx: vec![rho * u; n],
        my: vec![0.0; n],
        e: vec![rho * et; n],
    };
    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow { rho, u, v: 0.0, p },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };
    for _ in 0..40 {
        advance2d_combustion_axi(&mut state, &mut sp, &g, 0.4, false, &bc, 1, None).unwrap();
    }
    let (uu, _vv, ett) =
        nufor_core::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e).unwrap();
    let k = 0;
    let mut y = [0.0f64; ct::NS];
    for (s, dst) in y.iter_mut().enumerate() {
        *dst = sp.y[s][k];
    }
    let ei = ett[k] - 0.5 * uu[k] * uu[k];
    let p_cell = nufor_core::combustion_closure::pressure_at(state.rho[k], ei, &y);
    let t = nufor_core::combustion_closure::temperature_at(ei, &y);
    assert!(
        (uu[k] - u).abs() < 0.02 * u,
        "the flow must survive: u {} against {u}",
        uu[k]
    );
    assert!(
        (state.rho[k] - rho).abs() < 0.01 * rho,
        "the density must hold: {}",
        state.rho[k]
    );
    assert!(
        (p_cell - p).abs() < 0.02 * p,
        "the pressure must hold: {} Pa",
        p_cell
    );
    assert!(
        (t - 3676.53).abs() < 100.0,
        "the temperature must hold: {t} K"
    );
}

/// the axisymmetric species transport moves a radial gradient at the
/// rate the flow's radial flux does. the annular divergence already
/// holds the cell width, so the planar coefficient must not divide by
/// it a second time: the radial transport then runs at a fraction of
/// the flow's rate and the composition stratifies behind the flow.
#[test]
fn axi_species_transport_matches_the_radial_flux() {
    let (nx, ny) = (4, 8);
    let g = nufor_core::grid2d(nx, ny, 0.0, 0.04, 0.0, 0.08).unwrap();
    let n = nx * ny;
    let rho = 1.0f64;
    // a composition that falls off with the radius.
    let mut y = vec![vec![0.0f64; n]];
    for j in 0..ny {
        let r = (j as f64 + 0.5) * g.dy;
        for i in 0..nx {
            y[0][j * nx + i] = 0.2 * (1.0 - r / 0.08);
        }
    }
    let before = y[0].clone();
    // a radial mass flux, constant in the radius, so the divergence
    // term's radius weighting is what moves the composition.
    let fx: Vec<Vec<f64>> = (0..ny).map(|_| vec![0.0; nx + 1]).collect();
    let fy: Vec<Vec<f64>> = (0..nx).map(|_| vec![0.4; ny + 1]).collect();
    let dt = 1e-5;
    let dxs: Vec<f64> = (0..nx).map(|_| g.dx).collect();
    let dys: Vec<f64> = (0..ny).map(|_| g.dy).collect();
    nufor_core::species_transport::transport_species_generic(
        &mut y,
        &vec![rho; n],
        dt,
        &dxs,
        &dys,
        &fx,
        &fy,
        Some(&g.faces_y),
    )
    .unwrap();
    // the expected change: the annular divergence of the radial flux,
    // divided by the cell's radius.
    let mut worst = 0.0f64;
    for j in 1..ny - 1 {
        let rm = g.faces_y[j];
        let rp = g.faces_y[j + 1];
        let rc = 0.5 * (rm + rp);
        let dr = g.dy;
        let i = 1;
        let k = j * nx + i;
        // the upwind face composition: the flux is positive, so each
        // face carries the cell below it. face j sits between the
        // cells j-1 and j and takes the lower one; face j+1 sits
        // between j and j+1 and takes j.
        let y_below = before[k - nx];
        let y_here = before[k];
        let flux_lo = rm * y_below * 0.4;
        let flux_hi = rp * y_here * 0.4;
        let expected = -(flux_hi - flux_lo) / (rc * dr) * dt / rho;
        let got = y[0][k] - before[k];
        let rel = (got - expected).abs() / expected.abs().max(1e-30);
        worst = worst.max(rel);
    }
    assert!(
        worst < 1e-9,
        "the radial transport must match the annular divergence: {worst}"
    );
}

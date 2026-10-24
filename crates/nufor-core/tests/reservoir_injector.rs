//! the injector's boundary condition: a converging nozzle.
//!
//! an injector passes its own mass flow, and that flow follows the
//! nozzle relation between the tank state and the chamber's pressure.
//! the tests pin the shape of that relation and the shape of the exit
//! state, so a change to either shows up here.

use nufor_core::combustion_closure as cl;
use nufor_core::combustion_thermo as ct;
use nufor_core::{face_flow_area, injector_mass_flow, reservoir_real_gas_at, Bc2d, ThermoModel};

/// the tank state the tests use: 120 bar of the h2/o2 mixture.
const PT: f64 = 120.0e5;
const TT: f64 = 3676.53;

/// the burned chamber fill at cea's 70 bar equilibrium.
fn chamber_fill() -> Vec<f64> {
    vec![
        0.01583304, 0.07194470, 0.77605355, 0.11605824, 0.00232724, 0.01719592, 0.00049223,
        0.00009298, 0.0,
    ]
}

/// the tank's density at its own pressure and temperature.
fn tank_density() -> f64 {
    PT / (r_mix() * TT)
}

/// the mixture's gas constant at the tank's temperature.
fn r_mix() -> f64 {
    let mut y = [0.0f64; ct::NS];
    for (dst, src) in y.iter_mut().zip(chamber_fill()) {
        *dst = src;
    }
    cl::mix_gas_constant(&y)
}

/// the mixture's gamma at the tank's temperature.
fn gamma() -> f64 {
    let mut y = [0.0f64; ct::NS];
    for (dst, src) in y.iter_mut().zip(chamber_fill()) {
        *dst = src;
    }
    let cp = ct::mix_cp_mass(&y, TT);
    let rmix = r_mix();
    assert!(cp > rmix);
    cp / (cp - rmix)
}

/// a reservoir with an injector of the given area.
fn reservoir_bc(injector_area: Option<f64>) -> Bc2d {
    Bc2d::ReservoirInflowRealGas {
        pt: PT,
        tt: TT,
        // the tank's density at its own pressure and temperature.
        rho_t: tank_density(),
        fill: chamber_fill(),
        n_species: ct::NS,
        injector_area,
    }
}

/// the critical pressure ratio of the mixture.
fn critical_ratio() -> f64 {
    let g = gamma();
    (2.0 / (g + 1.0)).powf(g / (g - 1.0))
}

/// the choked mass flux of the tank state.
fn choked_flux() -> f64 {
    let g = gamma();
    tank_density()
        * (g * r_mix() * TT).sqrt()
        * (2.0 / (g + 1.0)).powf((g + 1.0) / (2.0 * (g - 1.0)))
}

/// the mass flow the injector passes at a chamber's pressure.
fn injector_flow(bc: &Bc2d, p_chamber: f64) -> f64 {
    let model = ThermoModel::Combustion9 {
        gamma_shaping: 1.195,
    };
    injector_mass_flow(bc, model, p_chamber).unwrap()
}

/// a choked injector holds its mass flow however far the chamber
/// drops below the critical pressure.
#[test]
fn a_choked_injector_holds_its_mass_flow() {
    let bc = reservoir_bc(Some(1.0e-3));
    let critical_p = critical_ratio() * PT;
    let want = choked_flux() * 1.0e-3;
    for p_chamber in [0.9 * critical_p, 0.7 * critical_p, 0.5 * critical_p] {
        let flow = injector_flow(&bc, p_chamber);
        assert!(
            (flow - want).abs() / want < 1e-6,
            "the injector's flow {flow} against the choked {want} at {p_chamber} pascal"
        );
    }
}

/// the unchoked flow rises with the chamber's pressure, because a
/// converging nozzle's mass flux grows as its exit pressure grows
/// until the exit reaches the sonic point.
#[test]
fn the_unchoked_flow_rises_with_the_chambers_pressure() {
    let bc = reservoir_bc(Some(1.0e-3));
    let near = injector_flow(&bc, 0.99 * PT);
    let mid = injector_flow(&bc, 0.9 * PT);
    let choked = injector_flow(&bc, 0.5 * critical_ratio() * PT);
    assert!(
        near < mid,
        "the unchoked feed must grow with the chamber: {near} against {mid}"
    );
    assert!(
        mid < choked,
        "and reach the choke below the critical pressure: {mid} against {choked}"
    );
}

/// a plain reservoir without an injector feeds the whole face.
#[test]
fn a_plain_reservoir_feeds_the_whole_face() {
    let bc = reservoir_bc(None);
    let model = ThermoModel::Combustion9 {
        gamma_shaping: 1.195,
    };
    let (rho, u, _) = reservoir_real_gas_at(&bc, model, 70.0e5, 3673.0).unwrap();
    assert!(rho > 0.0 && u > 0.0);
    // without an injector the ghost is the reservoir expanded to the
    // face's pressure, so the flow slows as the face's pressure rises.
    let (rho_hi, u_hi, _) = reservoir_real_gas_at(&bc, model, 0.99 * PT, 3673.0).unwrap();
    assert!(
        u_hi < u,
        "a plain reservoir must slow near the tank's pressure"
    );
    assert!(rho_hi > rho, "and it must be denser there");
}

/// the injector's area sets its mass flow: the area times the
/// choked flux of the tank state.
#[test]
fn the_injectors_area_sets_its_mass_flow() {
    let injector_area = 7.205373e-4;
    let bc = reservoir_bc(Some(injector_area));
    let model = ThermoModel::Combustion9 {
        gamma_shaping: 1.195,
    };
    let flow = injector_mass_flow(&bc, model, 70.0e5).unwrap();
    // the injector's mass flow is its area times the choked flux.
    assert!(
        (flow - 4.0078).abs() / 4.0078 < 0.01,
        "the injector's flow {flow} against the throat's 4.0078"
    );
}

/// a side without an injector has no mass flow of its own.
#[test]
fn a_side_without_an_injector_has_no_flow() {
    let bc = reservoir_bc(None);
    let model = ThermoModel::Combustion9 {
        gamma_shaping: 1.195,
    };
    assert!(injector_mass_flow(&bc, model, 70.0e5).is_none());
}

/// the face's flow area counts the fluid rows only.
///
/// a chamber's wall covers the rows above it, and those rows carry no
/// flow. counting them would spread an injector's mass flow over a
/// face the gas never sees, and the chamber would starve.
#[test]
fn the_face_area_skips_the_solid_rows() {
    let (nx, ny) = (4, 8);
    let g = nufor_core::grid2d(nx, ny, 0.0, 0.04, 0.0, 0.08).unwrap();
    // no solid rows: the whole face, the annulus from the axis out to
    // the outer radius.
    let whole = face_flow_area(&g, None, 0);
    let expected = std::f64::consts::PI * 0.08 * 0.08;
    assert!(
        (whole - expected).abs() / expected < 1e-12,
        "the whole face's area {whole} against {expected}"
    );
    // the top half is solid: the area drops to the inner annulus.
    let mut mask = vec![false; nx * ny];
    for j in ny / 2..ny {
        for i in 0..nx {
            mask[j * nx + i] = true;
        }
    }
    let inner = std::f64::consts::PI * 0.04 * 0.04;
    let half = face_flow_area(&g, Some(&mask), 0);
    assert!(
        (half - inner).abs() / inner < 1e-12,
        "the inner face's area {half} against {inner}"
    );
    // the solid rows must matter: without the mask the area is the
    // whole face again.
    assert!(half < whole, "the test needs the solid rows to matter");
}

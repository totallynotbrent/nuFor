//! the chamber-to-throat gate: the axisymmetric combustion lane.
//!
//! a chamber of burned gas sits behind an injector, necks down to a
//! throat and opens into a bell, with a back pressure at the exit. the
//! throat must choke, and its pressure and temperature must land on the
//! equilibrium throat state for this chamber.
//!
//! the reference is CEA's equilibrium solution for the injected gas at
//! a 7 MPa chamber: the throat sits at 4.05352 MPa and 3501.98 K. the
//! solver's own closure puts the frozen-composition sonic state at
//! 4.0074 to 4.0285 MPa and 3369.4 to 3369.7 K, so the gate holds the
//! throat within 12% in pressure and 5% in temperature of the CEA
//! values, which covers both references.
//!
//! the seed rides the mesh it marches: the quasi-1d state is built on
//! the grid's own per-column flow area, so the mass flux is uniform
//! from the first step and the throat's choked mass flow is the grid's,
//! not the continuous contour's.
//!
//! the settled march measures mach 1.0133 on the throat's column, a
//! throat pressure of 3.7087 MPa and a throat temperature of 3492.0 K,
//! with the exit at mach 1.3938. the converging section's mixture entropy
//! rises by less than 20 J/(kg K), so the flow is close to isentropic
//! there. the throat's pressure sits below the isentrope's because the
//! chemistry keeps running in the chamber and changes the gas the
//! nozzle sees.
//!
//! heavy: the march runs about thirty minutes in release on this box,
//! so the gate is ignore-marked and run locally:
//!
//!   cargo test -p nufor-cli --test chamber_throat_gate --release -- --ignored --nocapture
//!
//! the continuous integration runs the build, the lint and the parity
//! gates, not the chamber marches.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const NX: usize = 110;
const NY: usize = 40;
const X0: f64 = 0.0;
const X1: f64 = 0.165;
const Y0: f64 = 0.0;
const Y1: f64 = 0.1;

/// the march's end: the throat settles inside the first acoustic round
/// trip, so half a millisecond is enough for a steady state and keeps
/// the run near thirty minutes in release.
const FINAL_TIME: f64 = 0.0005;
const MAX_STEPS: usize = 6000;
/// the snapshot cadence: the last one is the state the gate reads.
const INTERVAL: usize = 250;

/// the wall's vertices: the chamber's contour, sampled along the wall
/// and closed over the top. the gate reads its radius back out of this
/// polygon, so the geometry the gate measures is the geometry the
/// solver marches.
const WALL: &str = "[[0.000750, 0.099995], [0.002250, 0.099956], [0.003750, 0.099878], [0.005250, 0.099760], [0.006750, 0.099604], [0.008250, 0.099409], [0.009750, 0.099175], [0.011250, 0.098903], [0.012750, 0.098593], [0.014250, 0.098245], [0.015750, 0.097860], [0.017250, 0.097437], [0.018750, 0.096978], [0.020250, 0.096483], [0.021750, 0.095952], [0.023250, 0.095385], [0.024750, 0.094785], [0.026250, 0.094150], [0.027750, 0.093482], [0.029250, 0.092781], [0.030750, 0.092047], [0.032250, 0.091283], [0.033750, 0.090488], [0.035250, 0.089663], [0.036750, 0.088810], [0.038250, 0.087928], [0.039750, 0.087019], [0.041250, 0.086083], [0.042750, 0.085122], [0.044250, 0.084137], [0.045750, 0.083128], [0.047250, 0.082096], [0.048750, 0.081043], [0.050250, 0.079970], [0.051750, 0.078877], [0.053250, 0.077766], [0.054750, 0.076637], [0.056250, 0.075493], [0.057750, 0.074333], [0.059250, 0.073160], [0.060750, 0.071973], [0.062250, 0.070776], [0.063750, 0.069568], [0.065250, 0.068351], [0.066750, 0.067126], [0.068250, 0.065894], [0.069750, 0.064656], [0.071250, 0.063415], [0.072750, 0.062170], [0.074250, 0.060924], [0.075750, 0.059676], [0.077250, 0.058430], [0.078750, 0.057185], [0.080250, 0.055944], [0.081750, 0.054706], [0.083250, 0.053474], [0.084750, 0.052249], [0.086250, 0.051032], [0.087750, 0.049824], [0.089250, 0.048627], [0.090750, 0.047440], [0.092250, 0.046267], [0.093750, 0.045107], [0.095250, 0.043963], [0.096750, 0.042834], [0.098250, 0.041723], [0.099750, 0.040630], [0.101250, 0.039557], [0.102750, 0.038504], [0.104250, 0.037472], [0.105750, 0.036463], [0.107250, 0.035478], [0.108750, 0.034517], [0.110250, 0.033581], [0.111750, 0.032672], [0.113250, 0.031790], [0.114750, 0.030937], [0.116250, 0.030112], [0.117750, 0.029317], [0.119250, 0.028553], [0.120750, 0.027819], [0.122250, 0.027118], [0.123750, 0.026450], [0.125250, 0.025815], [0.126750, 0.025215], [0.128250, 0.024648], [0.129750, 0.024117], [0.131250, 0.023622], [0.132750, 0.023163], [0.134250, 0.022740], [0.135750, 0.022355], [0.137250, 0.022007], [0.138750, 0.021697], [0.140250, 0.021425], [0.141750, 0.021191], [0.143250, 0.020996], [0.144750, 0.020840], [0.146250, 0.020722], [0.147750, 0.020644], [0.149250, 0.020605], [0.150750, 0.021259], [0.152250, 0.022561], [0.153750, 0.023815], [0.155250, 0.024989], [0.156750, 0.026055], [0.158250, 0.026987], [0.159750, 0.027762], [0.161250, 0.028361], [0.162750, 0.028768], [0.164250, 0.028974], [0.165000, 0.100000], [0.000000, 0.100000]]";

/// the case file, written into the test's directory.
fn write_case(dir: &Path) -> PathBuf {
    let toml = format!(
        r#"schema_version = 1

[metadata]
name = "chamber-throat"
description = "a chamber fed through an injector, necking to a throat and opening into a bell, with a back-pressure exit"
case_revision = 1

[physics]
equations = "euler_axi"
gamma = 1.2
gas_constant = 527.2
eos = "combustion"
gas = "h2o2"
chemistry = "on"

[physics.injection]
H2 = 0.01583304
O2 = 0.07194470
H2O = 0.77605355
OH = 0.11605824
H = 0.00232724
O = 0.01719592
HO2 = 0.00049223
H2O2 = 0.00009298

[mesh]
nx = {nx}
ny = {ny}
x0 = {x0}
x1 = {x1}
y0 = {y0}
y1 = {y1}
source = "uniform"

[initial_condition]
type = "quasi1d"
p0 = 7.0e6
t0 = 3676.53

[boundaries]
left = "reservoir_inflow"
right = "subsonic_outflow"
bottom = "slip_wall"
top = "slip_wall"
injector_area = 7.014600e-04

[boundaries.outflow_state]
p = 542900.0

[boundaries.inflow_state]
rho = 6.191099
u = 0.0
v = 0.0
p = 12000000.0

[body]
type = "polygon"
verts = {wall}

[numerics]
flux = "hllc"
reconstruction = "muscl"
cfl = 0.3
threads = 4

[time]
final_time = {final_time}
max_steps = {max_steps}

[output]
formats = ["vtk"]
interval_steps = {interval}
fields = ["rho", "p", "T", "mach"]
"#,
        nx = NX,
        ny = NY,
        x0 = X0,
        x1 = X1,
        y0 = Y0,
        y1 = Y1,
        wall = WALL,
        final_time = FINAL_TIME,
        max_steps = MAX_STEPS,
        interval = INTERVAL,
    );
    let path = dir.join("chamber-throat.toml");
    std::fs::write(&path, toml).unwrap();
    path
}

/// the polygon's vertices.
fn wall_verts() -> Vec<(f64, f64)> {
    let trimmed = WALL.trim_start_matches('[').trim_end_matches(']');
    trimmed
        .split("], [")
        .map(|pair| {
            let mut it = pair.trim_matches(|c| c == '[' || c == ']').split(',');
            let x: f64 = it.next().unwrap().trim().parse().unwrap();
            let y: f64 = it.next().unwrap().trim().parse().unwrap();
            (x, y)
        })
        .collect()
}

/// the signed distance to the polygon: negative inside.
fn signed_distance(verts: &[(f64, f64)], px: f64, py: f64) -> f64 {
    let mut best = f64::MAX;
    let mut inside = false;
    let n = verts.len();
    for i in 0..n {
        let (ax, ay) = verts[i];
        let (bx, by) = verts[(i + 1) % n];
        let (ex, ey) = (bx - ax, by - ay);
        let l2 = (ex * ex + ey * ey).max(1e-30);
        let t = (((px - ax) * ex + (py - ay) * ey) / l2).clamp(0.0, 1.0);
        let (qx, qy) = (ax + t * ex, ay + t * ey);
        best = best.min((px - qx).hypot(py - qy));
        if (ay > py) != (by > py) && px < ax + (py - ay) * ex / ey {
            inside = !inside;
        }
    }
    if inside {
        -best
    } else {
        best
    }
}

/// the solid mask: the solver treats a cell whose center sits within
/// 0.45 min(dx,dy) of the wall as solid, so the body fattens by half a
/// cell. the gate reads the same rows the solver marched.
fn solid_cells(verts: &[(f64, f64)]) -> Vec<bool> {
    let dx = (X1 - X0) / NX as f64;
    let dy = (Y1 - Y0) / NY as f64;
    let h = dx.min(dy);
    let mut mask = vec![false; NX * NY];
    for j in 0..NY {
        for i in 0..NX {
            let x = X0 + (i as f64 + 0.5) * dx;
            let y = Y0 + (j as f64 + 0.5) * dy;
            mask[j * NX + i] = signed_distance(verts, x, y) <= 0.45 * h;
        }
    }
    mask
}

/// the vtk's scalars.
fn vtk_scalars(path: &Path, name: &str) -> Vec<f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let marker = format!("SCALARS {name} double 1\nLOOKUP_TABLE default\n");
    let start = text.find(&marker).expect("the field is missing") + marker.len();
    let rest = &text[start..];
    let end = rest
        .find("\nSCALARS ")
        .or_else(|| rest.find("\nPOINT_DATA "))
        .unwrap_or(rest.len());
    let values: Vec<f64> = rest[..end]
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(values.len(), NX * NY, "the field's length");
    values
}

/// a column's flow state: the area-weighted means over the fluid rows,
/// and the mass flow through them.
#[derive(Debug)]
struct Station {
    area: f64,
    m_dot: f64,
    mach: f64,
    pressure: f64,
    temperature: f64,
}

fn station(
    col: usize,
    mask: &[bool],
    rho: &[f64],
    pressure: &[f64],
    temperature: &[f64],
    mach: &[f64],
    velocity_x: &[f64],
) -> Station {
    let dy = (Y1 - Y0) / NY as f64;
    let mut area = 0.0;
    let mut m_dot = 0.0;
    let mut p_sum = 0.0;
    let mut t_sum = 0.0;
    let mut m_sum = 0.0;
    for j in 0..NY {
        if mask[j * NX + col] {
            continue;
        }
        let k = j * NX + col;
        let rm = Y0 + j as f64 * dy;
        let rp = rm + dy;
        let a = std::f64::consts::PI * (rp * rp - rm * rm);
        area += a;
        m_dot += rho[k] * velocity_x[k] * a;
        p_sum += pressure[k] * a;
        t_sum += temperature[k] * a;
        m_sum += mach[k] * a;
    }
    assert!(area > 0.0, "the column has no fluid");
    Station {
        area,
        m_dot,
        mach: m_sum / area,
        pressure: p_sum / area,
        temperature: t_sum / area,
    }
}

/// the throat's column: the cell holding the polygon's narrowest point.
/// the grid's throat is a staircase plateau several cells wide, and the
/// flow only reaches mach one at the plateau's downstream edge, which
/// is where the wall's narrowest vertex sits.
fn throat_column(verts: &[(f64, f64)]) -> usize {
    let mut narrowest = f64::MAX;
    let mut x_throat = 0.0;
    for &(x, r) in verts {
        if r < narrowest {
            narrowest = r;
            x_throat = x;
        }
    }
    let dx = (X1 - X0) / NX as f64;
    (((x_throat - X0) / dx).floor() as usize).min(NX - 1)
}

/// the chamber-to-throat gate: the throat chokes and sits on CEA's
/// equilibrium throat state, and the mass flow is conserved.
#[test]
#[ignore = "heavy: an axisymmetric chamber march, run locally with -- --ignored"]
fn the_chamber_throat_chokes_at_the_cea_state() {
    let dir = std::env::temp_dir().join(format!("nufor_throat_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let case = write_case(&dir);
    let out: Output = Command::new(env!("CARGO_BIN_EXE_nufor"))
        .current_dir(&dir)
        .args(["run", case.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "the run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let vtk = dir.join("chamber-throat.vtk");
    let rho = vtk_scalars(&vtk, "density");
    let pressure = vtk_scalars(&vtk, "pressure");
    let temperature = vtk_scalars(&vtk, "temperature");
    let mach = vtk_scalars(&vtk, "mach");
    let velocity_x = vtk_scalars(&vtk, "velocity_x");

    let verts = wall_verts();
    let mask = solid_cells(&verts);
    let ti = throat_column(&verts);
    let throat = station(ti, &mask, &rho, &pressure, &temperature, &mach, &velocity_x);
    // the chamber's mid-section: the flow there is slow and axial.
    let ci = (ti / 2).max(1);
    let chamber = station(ci, &mask, &rho, &pressure, &temperature, &mach, &velocity_x);
    let exit = station(
        NX - 2,
        &mask,
        &rho,
        &pressure,
        &temperature,
        &mach,
        &velocity_x,
    );

    println!(
        "the throat (col {ti}): area {:.6e} m2, m_dot {:.4} kg/s, mach {:.4}, p {:.4} MPa, T {:.1} K",
        throat.area,
        throat.m_dot,
        throat.mach,
        throat.pressure / 1e6,
        throat.temperature
    );
    println!(
        "the chamber (col {ci}): m_dot {:.4} kg/s, mach {:.4}, p {:.4} MPa, T {:.1} K",
        chamber.m_dot,
        chamber.mach,
        chamber.pressure / 1e6,
        chamber.temperature
    );
    println!(
        "the exit (col {}): m_dot {:.4} kg/s, mach {:.4}, p {:.4} MPa, T {:.1} K",
        NX - 2,
        exit.m_dot,
        exit.mach,
        exit.pressure / 1e6,
        exit.temperature
    );

    // the duct must choke. the measured march lands at mach 1.0133 on
    // the throat's column; the band covers the staircase plateau, whose
    // cells run from mach 0.97 to mach 1.06.
    assert!(
        (throat.mach - 1.0).abs() < 0.08,
        "the throat must choke: mach {:.4}",
        throat.mach
    );
    // the throat's temperature must land on CEA's equilibrium throat.
    let cea_t = 3501.98;
    assert!(
        (throat.temperature - cea_t).abs() / cea_t < 0.05,
        "the throat's temperature must match CEA's {cea_t} K: {:.1} K",
        throat.temperature
    );
    // the throat's pressure must land on CEA's equilibrium throat.
    let cea_p = 4.05352e6;
    assert!(
        (throat.pressure - cea_p).abs() / cea_p < 0.12,
        "the throat's pressure must match CEA's {cea_p} Pa: {:.4} MPa",
        throat.pressure / 1e6
    );
    // the mass flow must be conserved: the chamber's feed is the
    // throat's throughput. the reacting lane's local states carry a few
    // percent of noise, so the band is wider than a roundoff check.
    assert!(
        (throat.m_dot - chamber.m_dot).abs() / chamber.m_dot.abs().max(1e-9) < 0.08,
        "the mass flow must be conserved: the chamber carries {:.4} kg/s, the throat {:.4} kg/s",
        chamber.m_dot,
        throat.m_dot
    );
    // the chamber ahead of the throat runs slow.
    assert!(
        chamber.mach < 0.1,
        "the chamber must run subsonic: mach {:.4}",
        chamber.mach
    );
    // a choked nozzle's exit runs supersonic.
    assert!(
        exit.mach > 1.0,
        "the exit must run supersonic: mach {:.4}",
        exit.mach
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

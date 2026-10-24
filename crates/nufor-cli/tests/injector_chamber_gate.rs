//! the injector feeds a chamber that holds its pressure.
//!
//! the gate drives the built binary on a real case: a chamber of
//! burned gas with an injector on its west face, a nozzle throat and a
//! bell, marched long enough for the chamber to fill. the chamber must
//! hold the pressure the seed gave it, the throat must stay near the
//! choked state, and no cell may collapse.
//!
//! this is the scenario that caught the ghost state. the ghost used to
//! carry the reservoir's own expansion, which runs colder than the
//! chamber's gas. the closure's internal energy carries the formation
//! offset, so the temperature difference was a megajoule-scale
//! enthalpy difference: the first column cooled, its pressure fell,
//! the ghost thinned, and the feed ran away. the ghost now takes the
//! chamber's own state.
//!
//! heavy: the march takes a few minutes, so the gate is ignore-marked
//! and run locally with `cargo test -p nufor-cli --test
//! injector_chamber_gate --release -- --ignored`. the continuous
//! integration runs the build, the lint and the parity gates, not the
//! chamber marches.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_nufor")
}

/// the case's chamber wall: a cylinder that necks down to a throat at
/// x = 0.15 and opens into a bell.
fn chamber_radius(x: f64) -> f64 {
    if x <= 0.15 {
        let s = x / 0.15;
        0.05 - 0.03 * (0.5 * (1.0 - (std::f64::consts::PI * s).cos()))
    } else {
        let s = ((x - 0.15) / 0.05).min(1.0);
        0.02 + 0.02 * (std::f64::consts::PI * s / 2.0).sin()
    }
}

/// the wall's vertices, from the chamber's radius.
fn wall_verts(n: usize) -> Vec<(f64, f64)> {
    let mut verts = Vec::with_capacity(n + 2);
    for i in 0..=n {
        let x = 0.20 * i as f64 / n as f64;
        verts.push((x, chamber_radius(x)));
    }
    verts.push((0.20, 0.10));
    verts.push((0.0, 0.10));
    verts
}

/// the case file, written into the test's directory.
fn write_case(dir: &Path) -> PathBuf {
    let verts = wall_verts(300);
    let verts = verts
        .iter()
        .map(|(x, y)| format!("[{x:.6}, {y:.6}]"))
        .collect::<Vec<String>>()
        .join(", ");
    let toml = format!(
        r#"schema_version = 1

[metadata]
name = "injector-chamber"
description = "a chamber fed through an injector, with a throat and a bell"
case_revision = 1

[physics]
equations = "euler_axi"
gamma = 1.2
gas_constant = 527.2
eos = "combustion"
gas = "h2o2"
chemistry = "off"

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
nx = 240
ny = 60
x0 = 0.0
x1 = 0.20
y0 = 0.0
y1 = 0.10
source = "uniform"

[initial_condition]
type = "quasi1d"
p0 = 7.0e6
t0 = 3676.53

[boundaries]
left = "reservoir_inflow"
right = "supersonic_outflow"
bottom = "slip_wall"
top = "slip_wall"
injector_area = 7.205373e-04

[boundaries.inflow_state]
rho = 6.191099
u = 0.0
v = 0.0
p = 12000000.0

[body]
type = "polygon"
verts = [{verts}]

[numerics]
flux = "hllc"
reconstruction = "muscl"
cfl = 0.05
threads = 1

[time]
final_time = 0.0002
max_steps = 1200

[output]
formats = ["vtk"]
interval_steps = 1200
fields = ["rho", "p", "T", "mach"]
"#
    );
    let path = dir.join("injector-chamber.toml");
    std::fs::write(&path, toml).unwrap();
    path
}

/// the vtk's scalars.
fn vtk_scalars(path: &Path, name: &str, nx: usize, ny: usize) -> Vec<f64> {
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
    assert_eq!(values.len(), nx * ny, "the field's length");
    values
}

/// the throat station's index: the narrowest radius in the grid.
fn throat_index(nx: usize) -> usize {
    let mut best = f64::MAX;
    let mut index = 0;
    for i in 0..nx {
        let x = (i as f64 + 0.5) * 0.20 / nx as f64;
        let a = chamber_radius(x) * chamber_radius(x);
        if a < best {
            best = a;
            index = i;
        }
    }
    index
}

/// the injector feeds a chamber that holds its pressure, and no cell
/// collapses while the march runs.
#[test]
#[ignore = "heavy: a chamber march, run locally with -- --ignored"]
fn an_injector_holds_a_chamber_at_its_pressure() {
    let dir = std::env::temp_dir().join(format!("nufor_injector_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let case = write_case(&dir);
    let out: Output = Command::new(bin())
        .current_dir(&dir)
        .args(["run", case.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "the run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let nx = 240;
    let ny = 60;
    let vtk = dir.join("injector-chamber.vtk");
    let pressure = vtk_scalars(&vtk, "pressure", nx, ny);
    let temperature = vtk_scalars(&vtk, "temperature", nx, ny);
    // no cell may collapse to the temperature floor.
    let coldest = temperature.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        coldest > 1000.0,
        "no cell may collapse to the floor: the coldest is {coldest} K"
    );
    // the chamber must hold its seed pressure. the chamber is the
    // volume upstream of the throat: the bell behind it has expanded,
    // and its pressure belongs to the nozzle, not to the chamber.
    let mut chamber = Vec::new();
    for i in 0..nx {
        let x = (i as f64 + 0.5) * 0.20 / nx as f64;
        if x > 0.15 {
            continue;
        }
        let r = chamber_radius(x);
        for j in 0..ny {
            let y = (j as f64 + 0.5) * 0.10 / ny as f64;
            if y < r {
                chamber.push(pressure[j * nx + i]);
            }
        }
    }
    let mean = chamber.iter().sum::<f64>() / chamber.len() as f64;
    assert!(
        (mean - 7.0e6).abs() / 7.0e6 < 0.10,
        "the chamber must hold its pressure: {mean} Pa"
    );
    // the throat must sit near the choked pressure.
    let ti = throat_index(nx);
    let r_throat = chamber_radius((ti as f64 + 0.5) * 0.20 / nx as f64);
    let mut throat = Vec::new();
    for j in 0..ny {
        let y = (j as f64 + 0.5) * 0.10 / ny as f64;
        if y < r_throat {
            throat.push(pressure[j * nx + ti]);
        }
    }
    let throat_mean = throat.iter().sum::<f64>() / throat.len() as f64;
    assert!(
        (throat_mean - 4.05352e6).abs() / 4.05352e6 < 0.1,
        "the throat must sit near the choked pressure: {throat_mean} Pa"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

//! a premix burns when the case says how hot it arrives.
//!
//! the gate drives the built binary on a real case: a constant-area
//! planar duct with the stoichiometric unburned mix entering at chamber
//! pressure. the mix must ignite, burn to a steady state, and reach the
//! exit burned.
//!
//! this is the scenario that caught the induction limit. an unburned
//! h2/o2 premix's own gas constant puts it near 800 K at 70 bar, and the
//! chain reaction is frozen there: the hydrogen's consumption rate
//! climbs by five orders of magnitude between 800 K and 1200 K. a case
//! that leaves the temperature unstated therefore sits frozen forever,
//! which looks like a working run and measures nothing.
//!
//! heavy: the hot march runs about forty minutes in release on this
//! box (the chemistry substep dominates once the gas is hot, 3.5 s per
//! step on a 120x20 grid), so the gate is ignore-marked and run locally
//! with `cargo test -p nufor-cli --test combustor_ignition_gate --release
//! -- --ignored`. the continuous integration runs the build, the lint
//! and the parity gates, not the combustor marches.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_nufor")
}

/// the case file, written into the test's directory. the inlet
/// temperature is the knob under test: at 800 K the mix must not burn,
/// at 1600 K it must.
fn write_case(dir: &Path, t_inlet: f64, name: &str) -> PathBuf {
    let toml = format!(
        r#"schema_version = 1

[metadata]
name = "{name}"
description = "a stoichiometric h2/o2 duct at chamber pressure"
case_revision = 1

[physics]
equations = "euler_2d"
gamma = 1.2
gas_constant = 692.3
eos = "combustion"
gas = "h2o2"
chemistry = "on"

[mesh]
nx = 120
ny = 20
x0 = 0.0
x1 = 0.6
y0 = 0.0
y1 = 0.1
source = "uniform"

[initial_condition]
type = "uniform"
rho = 12.64
u = 100.0
p = 7.0e6
t = {t_inlet:.1}

[boundaries]
left = "supersonic_inflow"
right = "subsonic_outflow"
bottom = "slip_wall"
top = "slip_wall"

[boundaries.inflow_state]
rho = 12.64
u = 100.0
v = 0.0
p = 7.0e6
t = {t_inlet:.1}

[boundaries.outflow_state]
rho = 12.64
u = 100.0
v = 0.0
p = 7.0e6

[numerics]
flux = "hllc"
reconstruction = "muscl"
cfl = 0.3
threads = 1

[time]
final_time = 0.001
max_steps = 700

[output]
formats = ["vtk"]
interval_steps = 700
fields = ["rho", "p", "T", "mach"]
"#
    );
    let path = dir.join(format!("{name}.toml"));
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

/// run a case and return the exit column's temperature and water
/// fraction.
fn run_case(dir: &Path, case: &Path, name: &str) -> (f64, f64) {
    let out: Output = Command::new(bin())
        .current_dir(dir)
        .args(["run", case.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "the run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let nx = 120;
    let ny = 20;
    let vtk = dir.join(format!("{name}.vtk"));
    let temperature = vtk_scalars(&vtk, "temperature", nx, ny);
    let water = vtk_scalars(&vtk, "y_H2O", nx, ny);
    // the exit column, at the duct's mid-height
    let j = ny / 2;
    let k = j * nx + nx - 1;
    (temperature[k], water[k])
}

/// a premix at the chain's floor stays cold: the mix must not burn
/// however long the march runs.
#[test]
#[ignore = "heavy: a combustor march, run locally with -- --ignored"]
fn a_cold_premix_never_burns() {
    let dir = std::env::temp_dir().join(format!("nufor_combustor_cold_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let case = write_case(&dir, 800.0, "combustor-cold");
    let (t_exit, water) = run_case(&dir, &case, "combustor-cold");
    assert!(
        t_exit < 1200.0,
        "the premix at 800 K must stay cold: the exit reads {t_exit} K"
    );
    assert!(
        water < 0.05,
        "the premix at 800 K must not make water: the exit reads {water}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

/// a premix above the chain's floor burns and reaches the exit burned.
#[test]
#[ignore = "heavy: a combustor march, run locally with -- --ignored"]
fn a_hot_premix_burns_to_the_exit() {
    let dir = std::env::temp_dir().join(format!("nufor_combustor_hot_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let case = write_case(&dir, 1600.0, "combustor-hot");
    let (t_exit, water) = run_case(&dir, &case, "combustor-hot");
    assert!(
        t_exit > 3000.0,
        "the premix at 1600 K must burn: the exit reads {t_exit} K"
    );
    assert!(
        water > 0.5,
        "the premix at 1600 K must make water: the exit reads {water}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

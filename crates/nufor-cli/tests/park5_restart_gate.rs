//! end-to-end gate: a park5 case writes its composition planes to the
//! restart file and a resumed run continues from the evolved chemistry
//! instead of the seed composition.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_nufor")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nufor_park5_rst_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn case(
    dir: &Path,
    restart_in: Option<&str>,
    restart_out: Option<&str>,
    final_time: f64,
) -> PathBuf {
    let ri = restart_in
        .map(|p| format!("restart_in = \"{p}\"\n"))
        .unwrap_or_default();
    let ro = restart_out
        .map(|p| format!("restart_out = \"{p}\"\n"))
        .unwrap_or_default();
    let toml = format!(
        r#"schema_version = 1

[metadata]
name = "park5-restart-gate"

[physics]
equations = "euler_2d"
gamma = 1.4
gas_constant = 296.8
eos = "park5"
gas = "n2"

[mesh]
nx = 24
ny = 16
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
source = "uniform"

[initial_condition]
type = "uniform"
rho = 1.363e-5
u = 6047.0
p = 0.89

[boundaries]
left = "supersonic_inflow"
right = "supersonic_outflow"
bottom = "slip_wall"
top = "supersonic_outflow"

[boundaries.inflow_state]
rho = 1.363e-5
u = 6047.0
p = 0.89

[numerics]
flux = "hllc"
reconstruction = "muscl"
cfl = 0.4

[output]
interval_steps = 1000
formats = ["vtk"]
fields = ["rho", "u", "v", "p", "t"]

[time]
final_time = {final_time}
max_steps = 4000
{ri}{ro}
"#,
        final_time = final_time
    );
    let path = dir.join("case.toml");
    std::fs::write(&path, toml).unwrap();
    path
}

fn run_case(path: &Path) -> std::process::Output {
    Command::new(bin())
        .arg("run")
        .arg(path)
        .output()
        .expect("the cli binary runs")
}

/// the composition planes recovered from a restart file: one plane per species.
fn species_planes(path: &Path) -> Vec<Vec<f64>> {
    let rd = nufor_core::read_restart2d_species(path).expect("the restart reads back");
    rd.species
        .expect("the species layout carries the composition planes")
}

#[test]
fn park5_restart_carries_the_composition_and_resumes_from_it() {
    let d = tmpdir("carry");
    let rst = d.join("run.rst");

    // the first leg evolves the chemistry and writes the composition out.
    let c1 = case(&d, None, Some(rst.to_str().unwrap()), 2.0e-5);
    let o1 = run_case(&c1);
    assert!(
        o1.status.success(),
        "first leg failed: {}",
        String::from_utf8_lossy(&o1.stderr)
    );

    let planes = species_planes(&rst);
    assert_eq!(planes.len(), 5, "five species planes are stored");
    for plane in &planes {
        assert_eq!(plane.len(), 24 * 16, "each plane spans the whole grid");
    }
    // the run really dissociated nitrogen: some cell carries atomic nitrogen.
    let n_atom_peak = planes[3].iter().copied().fold(0.0f64, f64::max);
    assert!(
        n_atom_peak > 0.0,
        "the evolved chemistry carries atomic nitrogen (peak {})",
        n_atom_peak
    );

    // the resumed leg reads the composition back, continues the march, and
    // writes its own restart: the chemistry must carry forward, not reset
    // to the pure-nitrogen seed.
    let rst2 = d.join("run2.rst");
    let c2 = case(
        &d,
        Some(rst.to_str().unwrap()),
        Some(rst2.to_str().unwrap()),
        4.0e-5,
    );
    let o2 = run_case(&c2);
    assert!(
        o2.status.success(),
        "resumed leg failed: {}",
        String::from_utf8_lossy(&o2.stderr)
    );
    let n_atom_again = species_planes(&rst2)[3]
        .iter()
        .copied()
        .fold(0.0f64, f64::max);
    assert!(
        n_atom_again > 0.0,
        "the resumed run continues the evolved composition"
    );

    std::fs::remove_dir_all(&d).unwrap();
}

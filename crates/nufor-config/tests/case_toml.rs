//! integration coverage for the case.toml schema against the real example file.

use nufor_config::{ConfigError, FluxScheme, InitialCondition, SCHEMA_VERSION};
use std::path::PathBuf;

fn example_case() -> PathBuf {
    // CARGO_MANIFEST_DIR is crates/nufor-config; the example lives under cases/.
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../cases/tutorial/case.toml")
}

#[test]
fn loads_example_case_file() {
    let config = nufor_config::load_case_config(example_case()).expect("example case must load");
    assert_eq!(config.schema_version, SCHEMA_VERSION);
    assert_eq!(config.metadata.name, "sod-shock-tube");
    assert_eq!(config.physics.gamma, 1.4);
    assert_eq!(config.mesh.nx, 100);
    assert_eq!(config.numerics.flux, FluxScheme::Hll);
    match config.initial_condition {
        InitialCondition::TwoState { left, right } => {
            assert_eq!(left.rho, 1.0);
            assert_eq!(right.rho, 0.125);
            assert_eq!(left.p, 1.0);
            assert_eq!(right.p, 0.1);
        }
        _ => panic!("example case must be a two-state shock tube"),
    }
}

#[test]
fn round_trip_preserves_config() {
    let config = nufor_config::load_case_config(example_case()).unwrap();
    let serialized = toml::to_string(&config).expect("config must serialize");
    let reparsed = nufor_config::parse_case_toml(&serialized).expect("round trip must parse");
    assert_eq!(reparsed, config);
}

#[test]
fn rejects_unsupported_schema_version() {
    let toml = std::fs::read_to_string(example_case()).unwrap();
    let bumped = toml.replace("schema_version = 1", "schema_version = 2");
    match nufor_config::parse_case_toml(&bumped) {
        Err(ConfigError::UnsupportedVersion { found, supported }) => {
            assert_eq!(found, 2);
            assert_eq!(supported, SCHEMA_VERSION);
        }
        other => panic!("expected UnsupportedVersion, got {other:?}"),
    }
}

#[test]
fn rejects_unknown_fields() {
    let toml = std::fs::read_to_string(example_case()).unwrap();
    let with_junk = format!("{toml}\nstray_key = \"oops\"\n");
    match nufor_config::parse_case_toml(&with_junk) {
        Err(ConfigError::Parse { message }) => {
            assert!(
                message.contains("stray_key") || message.contains("unknown field"),
                "parse error should mention the unknown key: {message}"
            );
        }
        other => panic!("expected Parse error, got {other:?}"),
    }
}

#[test]
fn missing_file_reports_io_error() {
    let missing = example_case().with_file_name("does-not-exist.toml");
    match nufor_config::load_case_config(missing) {
        Err(ConfigError::Io { .. }) => {}
        other => panic!("expected Io error, got {other:?}"),
    }
}

#[test]
fn rejects_wrong_type() {
    let toml = std::fs::read_to_string(example_case()).unwrap();
    let bad = toml.replace("flux = \"hll\"", "flux = true");
    assert!(matches!(
        nufor_config::parse_case_toml(&bad),
        Err(ConfigError::Parse { .. })
    ));
}

#[test]
fn invalid_case_reports_all_problems() {
    let toml = std::fs::read_to_string(example_case()).unwrap();
    let broken = toml
        .replace("gamma = 1.4", "gamma = 1.0")
        .replace("nx = 100", "nx = 1")
        .replace("cfl = 0.5", "cfl = 2.0");
    match nufor_config::parse_case_toml(&broken) {
        Err(ConfigError::Invalid { problems }) => {
            assert!(
                problems.len() >= 3,
                "expected multiple problems: {problems:?}"
            );
            assert!(problems.iter().any(|p| p.contains("gamma")));
            assert!(problems.iter().any(|p| p.contains("nx")));
            assert!(problems.iter().any(|p| p.contains("cfl")));
        }
        other => panic!("expected Invalid, got {other:?}"),
    }
}

#[test]
fn loads_and_validates_a_2d_blast_case() {
    let toml = r#"
schema_version = 1
[metadata]
name = "blast-2d"
[physics]
equations = "euler_2d"
gamma = 1.4
gas_constant = 287.0
[mesh]
nx = 80
ny = 80
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
[initial_condition]
type = "blast"
radius = 0.2
ambient = { rho = 1.0, u = 0.0, p = 1.0 }
fireball = { rho = 1.0, u = 0.0, p = 10.0 }
[boundaries]
left = "wall"
right = "wall"
[numerics]
flux = "hll"
reconstruction = "muscl"
cfl = 0.4
[time]
final_time = 0.15
max_steps = 10000
[output]
interval_steps = 100
formats = ["vtk"]
fields = ["rho"]
"#;
    let cfg = nufor_config::parse_case_toml(toml).expect("2d blast case must parse");
    assert_eq!(cfg.physics.equations, nufor_config::Equations::Euler2d);
    assert_eq!(cfg.mesh.ny, Some(80));
    match cfg.initial_condition {
        nufor_config::InitialCondition::Blast {
            radius, fireball, ..
        } => {
            assert_eq!(radius, 0.2);
            assert_eq!(fireball.p, 10.0);
        }
        _ => panic!("expected a blast initial condition"),
    }
    let problems = cfg.validate();
    assert!(
        problems.is_empty(),
        "blast case must validate: {problems:?}"
    );
}

#[test]
fn rejects_a_blast_ic_with_a_nonpositive_radius() {
    let toml = r#"
schema_version = 1
[metadata]
name = "reject-blast"
[physics]
equations = "euler_2d"
gamma = 1.4
gas_constant = 287.0
[mesh]
nx = 10
ny = 10
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
[initial_condition]
type = "blast"
radius = 0.0
ambient = { rho = 1.0, u = 0.0, p = 1.0 }
fireball = { rho = 1.0, u = 0.0, p = 10.0 }
[boundaries]
left = "wall"
right = "wall"
[numerics]
flux = "hll"
reconstruction = "first_order"
cfl = 0.5
[time]
final_time = 0.1
[output]
interval_steps = 100
formats = ["vtk"]
fields = ["rho"]
"#;
    let cfg = nufor_config::parse_case_toml(toml);
    match cfg {
        Err(nufor_config::ConfigError::Invalid { problems }) => assert!(
            problems.iter().any(|p| p.contains("radius")),
            "expected a radius complaint, got {problems:?}"
        ),
        other => panic!("expected a validation error, got {other:?}"),
    }
}

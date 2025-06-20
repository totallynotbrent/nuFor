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

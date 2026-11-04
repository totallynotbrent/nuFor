//! the zero-dimensional ignition case schema.
//!
//! a zerod case is a single reacting parcel in a rigid vessel. the
//! case names the mechanism and the parcel, and the march integrates
//! the parcel in time and writes the trace beside the case file.

use std::collections::HashMap;
use std::fmt;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

/// the schema version the crate understands.
pub const IGNITION_SCHEMA_VERSION: u32 = 1;

/// why an ignition case was rejected.
#[derive(Debug)]
pub enum IgnitionError {
    /// the file could not be read.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// the file is not valid TOML or does not match the schema.
    Parse { message: String },
    /// schema_version is present but not the supported version.
    UnsupportedVersion { found: u32, supported: u32 },
    /// the case parsed but broke a semantic rule.
    Invalid { problems: Vec<String> },
}

impl fmt::Display for IgnitionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IgnitionError::Io { path, source } => {
                write!(f, "cannot read case file {}: {source}", path.display())
            }
            IgnitionError::Parse { message } => write!(f, "invalid case toml: {message}"),
            IgnitionError::UnsupportedVersion { found, supported } => {
                write!(
                    f,
                    "unsupported ignition schema version {found}, this build supports {supported}"
                )
            }
            IgnitionError::Invalid { problems } => {
                write!(
                    f,
                    "ignition case validation failed: {}",
                    problems.join("; ")
                )
            }
        }
    }
}

impl std::error::Error for IgnitionError {}

/// a zero-dimensional ignition case.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IgnitionCase {
    /// must equal the supported schema version.
    pub schema_version: u32,
    /// the case name, used in the summary and the output files.
    pub name: String,
    /// the parcel's starting state and the march's settings.
    pub parcel: Parcel,
}

/// the reacting parcel and the march that drives it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parcel {
    /// the chemical mechanism: "gri30" today.
    pub mechanism: String,
    /// the fuel the parcel burns.
    pub fuel: String,
    /// the starting temperature, kelvin.
    pub initial_temperature_k: f64,
    /// the starting pressure, pascal.
    pub pressure_pa: f64,
    /// how long the march runs, seconds.
    pub final_time_s: f64,
    /// how finely the trace is sampled, seconds.
    pub sample_interval_s: f64,
    /// the largest step the march may take, seconds.
    pub max_step_s: f64,
    /// the temperature rise above the start that marks ignition.
    pub ignition_threshold_k: f64,
    /// the parcel's composition as mass fractions by species name.
    /// omitted means the stoichiometric unburned fuel and air.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub composition: Option<HashMap<String, f64>>,
}

/// every violated rule; empty means the case is valid.
pub fn validate(case: &IgnitionCase) -> Vec<String> {
    let mut problems = Vec::new();

    if case.name.trim().is_empty() {
        problems.push("name must not be empty".to_owned());
    }
    if case.parcel.mechanism != "gri30" {
        problems.push(format!(
            "parcel.mechanism must be \"gri30\" (got {:?})",
            case.parcel.mechanism
        ));
    }
    if case.parcel.initial_temperature_k <= 0.0 {
        problems.push("parcel.initial_temperature_k must be > 0".to_owned());
    }
    if case.parcel.pressure_pa <= 0.0 {
        problems.push("parcel.pressure_pa must be > 0".to_owned());
    }
    if case.parcel.final_time_s <= 0.0 {
        problems.push("parcel.final_time_s must be > 0".to_owned());
    }
    if case.parcel.sample_interval_s <= 0.0 {
        problems.push("parcel.sample_interval_s must be > 0".to_owned());
    }
    if case.parcel.max_step_s <= 0.0 {
        problems.push("parcel.max_step_s must be > 0".to_owned());
    }
    if case.parcel.sample_interval_s > case.parcel.final_time_s {
        problems.push("parcel.sample_interval_s must not exceed parcel.final_time_s".to_owned());
    }
    if case.parcel.ignition_threshold_k <= 0.0 {
        problems.push("parcel.ignition_threshold_k must be > 0".to_owned());
    }
    if let Some(composition) = &case.parcel.composition {
        let mut total = 0.0f64;
        for (name, value) in composition {
            if *value < 0.0 {
                problems.push(format!(
                    "composition[{name:?}] must not be negative (got {value})"
                ));
            }
            total += value;
        }
        if (total - 1.0).abs() > 1.0e-9 {
            problems.push(format!("composition masses must sum to one (got {total})"));
        }
    }

    problems
}

/// parse and validate an ignition case from a string.
pub fn parse_ignition_toml(input: &str) -> Result<IgnitionCase, IgnitionError> {
    let parsed: IgnitionCase = toml::from_str(input).map_err(|err| IgnitionError::Parse {
        message: err.to_string(),
    })?;
    if parsed.schema_version != IGNITION_SCHEMA_VERSION {
        return Err(IgnitionError::UnsupportedVersion {
            found: parsed.schema_version,
            supported: IGNITION_SCHEMA_VERSION,
        });
    }
    let problems = validate(&parsed);
    if problems.is_empty() {
        Ok(parsed)
    } else {
        Err(IgnitionError::Invalid { problems })
    }
}

/// read, parse, and validate an ignition case file.
pub fn load_ignition_case(path: impl AsRef<Path>) -> Result<IgnitionCase, IgnitionError> {
    let text = std::fs::read_to_string(path.as_ref()).map_err(|source| IgnitionError::Io {
        path: path.as_ref().to_path_buf(),
        source,
    })?;
    parse_ignition_toml(&text)
}

#[cfg(test)]
mod tests {
    //! the ignition case schema's own tests.

    use super::*;

    const GOOD: &str = r#"
schema_version = 1
name = "good case"

[parcel]
mechanism = "gri30"
fuel = "CH4"
initial_temperature_k = 2000.0
pressure_pa = 101325.0
final_time_s = 2.0e-3
sample_interval_s = 2.0e-6
max_step_s = 5.0e-7
ignition_threshold_k = 400.0
"#;

    #[test]
    fn a_good_case_parses() {
        let case = parse_ignition_toml(GOOD).expect("the case parses");
        assert_eq!(case.name, "good case");
        assert_eq!(case.parcel.mechanism, "gri30");
        assert_eq!(case.parcel.fuel, "CH4");
        assert!((case.parcel.initial_temperature_k - 2000.0).abs() < 1.0e-12);
        assert!(case.parcel.composition.is_none());
    }

    #[test]
    fn a_future_schema_version_is_refused() {
        let text = GOOD.replace("schema_version = 1", "schema_version = 99");
        let error = parse_ignition_toml(&text).expect_err("the version is refused");
        assert!(matches!(
            error,
            IgnitionError::UnsupportedVersion {
                found: 99,
                supported: IGNITION_SCHEMA_VERSION
            }
        ));
    }

    #[test]
    fn an_unknown_field_is_refused() {
        let text = GOOD.replace("name = \"good case\"", "name = \"good case\"\nmystery = 1");
        assert!(matches!(
            parse_ignition_toml(&text),
            Err(IgnitionError::Parse { .. })
        ));
    }

    #[test]
    fn a_negative_temperature_is_refused() {
        let text = GOOD.replace(
            "initial_temperature_k = 2000.0",
            "initial_temperature_k = -1.0",
        );
        let error = parse_ignition_toml(&text).expect_err("the temperature is refused");
        match error {
            IgnitionError::Invalid { problems } => {
                assert!(problems.iter().any(|p| p.contains("initial_temperature_k")));
            }
            other => panic!("unexpected error: {other}"),
        }
    }

    #[test]
    fn a_composition_that_does_not_sum_to_one_is_refused() {
        let text = format!("{GOOD}\n[parcel.composition]\nCH4 = 0.5\nO2 = 0.2\n");
        let error = parse_ignition_toml(&text).expect_err("the composition is refused");
        match error {
            IgnitionError::Invalid { problems } => {
                assert!(problems.iter().any(|p| p.contains("sum to one")));
            }
            other => panic!("unexpected error: {other}"),
        }
    }

    #[test]
    fn a_sample_interval_longer_than_the_horizon_is_refused() {
        let text = GOOD.replace("sample_interval_s = 2.0e-6", "sample_interval_s = 1.0e-2");
        let error = parse_ignition_toml(&text).expect_err("the interval is refused");
        match error {
            IgnitionError::Invalid { problems } => {
                assert!(problems.iter().any(|p| p.contains("must not exceed")));
            }
            other => panic!("unexpected error: {other}"),
        }
    }
}

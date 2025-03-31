//! Case-file configuration: the `case.toml` schema and its validation (PLAN
//! step 3, spec 40). The case file is the reproducibility root; it records
//! physics, mesh source, numerical methods, solver controls, output settings,
//! and metadata.
//!
//! Parsing goes through the `toml` crate with strict serde types
//! (`deny_unknown_fields`), so a case that does not match the schema fails
//! fast instead of silently running with a misspelled setting. Semantic rules
//! (ranges, termination criteria, physical positivity) run in
//! [`CaseConfig::validate`] after deserialization; see
//! `docs/formats/case-toml.md` for the human-readable schema reference.

pub mod schema;

mod validate;

pub use schema::{
    Boundaries, BoundaryKind, CaseConfig, Equations, FluxScheme, InitialCondition, Mesh, Metadata,
    Numerics, Output, OutputFormat, Physics, Reconstruction, ReferenceConditions, State,
    TimeControl, UnitSystem,
};

use std::fmt;
use std::io;
use std::path::{Path, PathBuf};

/// Schema version this crate understands. `case.toml` must declare exactly
/// this value in `schema_version`, so future incompatible layouts can be
/// detected instead of mis-parsed.
pub const SCHEMA_VERSION: u32 = 1;

/// Why a case file was rejected.
#[derive(Debug)]
pub enum ConfigError {
    /// The file could not be read.
    Io { path: PathBuf, source: io::Error },
    /// The file is not valid TOML or does not match the schema shape.
    Parse { message: String },
    /// `schema_version` is present but not [`SCHEMA_VERSION`].
    UnsupportedVersion { found: u32, supported: u32 },
    /// The file parsed but failed one or more semantic validation rules.
    Invalid { problems: Vec<String> },
}

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConfigError::Io { path, source } => {
                write!(f, "cannot read case file {}: {source}", path.display())
            }
            ConfigError::Parse { message } => write!(f, "invalid case.toml: {message}"),
            ConfigError::UnsupportedVersion { found, supported } => {
                write!(
                    f,
                    "unsupported case schema version {found}, this build supports {supported}"
                )
            }
            ConfigError::Invalid { problems } => {
                write!(f, "case validation failed: {}", problems.join("; "))
            }
        }
    }
}

impl std::error::Error for ConfigError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ConfigError::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Parse and validate a `case.toml` document from a string.
pub fn parse_case_toml(input: &str) -> Result<CaseConfig, ConfigError> {
    let parsed: CaseConfig = toml::from_str(input).map_err(|err| ConfigError::Parse {
        message: err.to_string(),
    })?;
    if parsed.schema_version != SCHEMA_VERSION {
        return Err(ConfigError::UnsupportedVersion {
            found: parsed.schema_version,
            supported: SCHEMA_VERSION,
        });
    }
    let problems = CaseConfig::validate(&parsed);
    if problems.is_empty() {
        Ok(parsed)
    } else {
        Err(ConfigError::Invalid { problems })
    }
}

/// Read, parse, and validate a `case.toml` file.
pub fn load_case_config(path: impl AsRef<Path>) -> Result<CaseConfig, ConfigError> {
    let text = std::fs::read_to_string(&path).map_err(|source| ConfigError::Io {
        path: path.as_ref().to_path_buf(),
        source,
    })?;
    parse_case_toml(&text)
}

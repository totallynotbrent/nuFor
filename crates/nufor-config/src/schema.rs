//! typed case.toml schema; strict structs reject unknown keys.

use serde::{Deserialize, Serialize};

/// top-level case file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CaseConfig {
    /// must equal the supported schema version.
    pub schema_version: u32,
    pub metadata: Metadata,
    pub physics: Physics,
    pub mesh: Mesh,
    pub initial_condition: InitialCondition,
    pub boundaries: Boundaries,
    pub numerics: Numerics,
    pub time: TimeControl,
    pub output: Output,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Metadata {
    /// short case name, used in logs and result files.
    pub name: String,
    /// free text: goal, expected result, history.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// manual revision counter for the case definition.
    #[serde(default = "default_case_revision")]
    pub case_revision: u32,
}

fn default_case_revision() -> u32 {
    1
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Physics {
    /// equation set to solve.
    pub equations: Equations,
    /// ratio of specific heats.
    pub gamma: f64,
    /// specific gas constant, J/(kg K).
    pub gas_constant: f64,
    /// unit convention for reported fields.
    #[serde(default)]
    pub unit_system: UnitSystem,
    /// nondimensionalization reference state.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference: Option<ReferenceConditions>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Equations {
    #[serde(rename = "euler_1d")]
    Euler1d,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum UnitSystem {
    #[default]
    Si,
    Nondim,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceConditions {
    pub rho: f64,
    pub p: f64,
    pub u: f64,
    /// length scale, m in SI.
    pub l: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mesh {
    /// number of cells in the 1D domain.
    pub nx: u32,
    /// domain left edge.
    pub x0: f64,
    /// domain right edge.
    pub x1: f64,
    /// how the mesh is produced; a file path once external meshes exist.
    #[serde(default = "default_mesh_source")]
    pub source: String,
}

fn default_mesh_source() -> String {
    "uniform".to_owned()
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum InitialCondition {
    Uniform {
        rho: f64,
        u: f64,
        p: f64,
    },
    /// left/right constant states, e.g. a shock tube.
    TwoState {
        left: State,
        right: State,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub rho: f64,
    /// x-velocity.
    pub u: f64,
    pub p: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Boundaries {
    pub left: BoundaryKind,
    pub right: BoundaryKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryKind {
    Wall,
    Inflow,
    Outflow,
    Periodic,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Numerics {
    /// riemann solver / approximate flux.
    pub flux: FluxScheme,
    /// spatial reconstruction order.
    pub reconstruction: Reconstruction,
    /// CFL number for the explicit time step.
    pub cfl: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FluxScheme {
    Hll,
    Hllc,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Reconstruction {
    FirstOrder,
    Muscl,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeControl {
    /// stop time; 0.0 with no other criterion is rejected.
    #[serde(default)]
    pub final_time: f64,
    /// hard cap on iterations.
    #[serde(default)]
    pub max_steps: u64,
    /// stop early when the residual drops below this.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub residual_target: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    /// write cadence in solver steps.
    pub interval_steps: u64,
    /// which writers run (CSV, VTK, HDF5).
    pub formats: Vec<OutputFormat>,
    /// fields written, e.g. rho, u, p.
    pub fields: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OutputFormat {
    Csv,
    Vtk,
    Hdf5,
}

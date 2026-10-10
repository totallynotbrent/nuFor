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
    /// immersed analytic body; absent means an empty domain.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub body: Option<BodySection>,
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
    /// molecular dynamic viscosity, required by the viscous/rans equation
    /// sets (euler ignores it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mu: Option<f64>,
    /// thermodynamic closure: perfect gas or equilibrium air.
    #[serde(default)]
    pub eos: Eos,
    /// molecular prandtl number for the heat flux (default 0.72, air).
    #[serde(
        default = "default_prandtl",
        skip_serializing_if = "is_default_prandtl"
    )]
    pub pr: f64,
    /// spalart-allmaras turbulence settings, required when equations is
    /// rans_2d_sa.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub turbulence: Option<Turbulence>,
    /// isothermal wall temperature (K) for the viscous wall boundary;
    /// omitted or non-positive means the adiabatic wall. the stagnation
    /// heating figures use a cold wall because the published entry
    /// correlations (fay-riddell) are posed at a fixed wall temperature.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_temperature: Option<f64>,
}

/// the sa turbulence block of a rans_2d_sa case.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Turbulence {
    /// freestream modified viscosity nu_tilde; the convention is 3*nu.
    pub nu_tilde_inf: f64,
    /// turbulent prandtl number for the eddy heat flux (default 0.9).
    #[serde(
        default = "default_prandtl_t",
        skip_serializing_if = "is_default_prandtl_t"
    )]
    pub pr_t: f64,
}

fn default_prandtl() -> f64 {
    0.72
}

fn default_prandtl_t() -> f64 {
    0.9
}

fn is_default_prandtl(v: &f64) -> bool {
    *v == 0.72
}

fn is_default_prandtl_t(v: &f64) -> bool {
    *v == 0.9
}

/// thermodynamic closure for the compressible equation sets.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Eos {
    /// constant gamma, constant R (the classical closure).
    #[default]
    #[serde(rename = "perfect")]
    Perfect,
    /// equilibrium air curve fits: p, T, a as functions of (rho, e).
    #[serde(rename = "eqair")]
    EqAir,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Equations {
    #[serde(rename = "euler_1d")]
    Euler1d,
    #[serde(rename = "euler_2d")]
    Euler2d,
    #[serde(rename = "euler_3d")]
    Euler3d,
    /// axisymmetric euler: the planar fluxes with the annular update
    /// (face radii and the p/r geometric source). the mesh's y axis is
    /// the radius; the south boundary is the symmetry axis.
    #[serde(rename = "euler_axi")]
    EulerAxi,
    #[serde(rename = "rans_2d_sa")]
    Rans2dSa,
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
    /// cells in y (required for euler_2d/3d).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ny: Option<u32>,
    /// y-domain lower edge.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub y0: Option<f64>,
    /// y-domain upper edge.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub y1: Option<f64>,
    /// cells in z (required for euler_3d).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nz: Option<u32>,
    /// z-domain lower edge.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub z0: Option<f64>,
    /// z-domain upper edge.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub z1: Option<f64>,
    /// how the mesh is produced: "uniform" from nx/x0/x1, "file" loading a
    /// list of cell centers from `path`, or "clustered" stretching the y (and
    /// optionally x) axis toward a wall with `first_cell`/`growth`.
    #[serde(default = "default_mesh_source")]
    pub source: String,
    /// cell-center mesh file read when source = "file".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    /// mesh-file dialect when source = "file": "vtk" (rectilinear/structured
    /// points) or "coordinates" (per-axis face lists). auto-detected when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<String>,
    /// first-cell height on a clustered axis (source = "clustered").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub first_cell: Option<f64>,
    /// geometric growth ratio on a clustered axis (source = "clustered").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub growth: Option<f64>,
    /// where source = "clustered" refines: "wall" (bottom), "top", or
    /// "channel" (both sides); defaults to "wall".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster: Option<String>,
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
    /// an over-pressured sphere in ambient surroundings (the classic blast);
    /// natural for the 2D and 3D equations.
    Blast {
        /// radius of the fireball.
        radius: f64,
        /// ambient state outside the fireball.
        ambient: State,
        /// high-pressure fireball state.
        fireball: State,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub rho: f64,
    /// x-velocity.
    pub u: f64,
    /// y-velocity, when the context has one (1d states ignore it).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub v: f64,
    pub p: f64,
}

fn is_zero(v: &f64) -> bool {
    *v == 0.0
}

/// an immersed analytic body. the meridian is evaluated per cell; the
/// mask (and cut cells where the solver supports them) come from the
/// same signed distance. `inside` convention: negative distance is
/// inside the solid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum BodySection {
    /// a sphere-cone body of revolution: a spherical nose cap of nose
    /// radius `rn` tangent to a conical flank of half-angle
    /// `delta_deg` (degrees), flaring to base radius `rb`. `xc` is the
    /// axial station of the sphere center, so the nose tip sits at
    /// `xc - rn`.
    SphereCone {
        /// nose radius R_N.
        rn: f64,
        /// cone half-angle in degrees.
        delta_deg: f64,
        /// base radius R_b.
        rb: f64,
        /// axial station of the sphere center.
        xc: f64,
    },
    /// a circle (a cylinder in 2d planar, a sphere in axisymmetric).
    Sphere { r: f64, cx: f64, cy: f64 },
    /// a closed polygon from its vertices, ordered around the boundary.
    Polygon {
        /// vertex list, [[x, y], ...] in mesh units.
        verts: Vec<(f64, f64)>,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Boundaries {
    pub left: BoundaryKind,
    pub right: BoundaryKind,
    /// extra 2d sides (top/bottom aliases) honored by the 2d solvers; the
    /// 1d runner ignores them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top: Option<BoundaryKind>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bottom: Option<BoundaryKind>,
    /// the inflow profile when a side is `profile_inflow`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inflow_profile: Option<InflowProfileSpec>,
    /// the fixed state a `supersonic_inflow` side carries.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inflow_state: Option<State>,
}

/// how a profile-inflow side prescribes its boundary layer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case", tag = "type")]
pub enum InflowProfileSpec {
    /// the blasius laminar layer; `leading_edge` is the virtual x station
    /// where the layer starts, upstream of or at the inflow face.
    Blasius {
        /// the virtual leading-edge x station (same units as mesh.x0).
        leading_edge: f64,
    },
    /// a tabulated profile: y, u, v columns read from `path`.
    Table {
        /// the profile file: one `y u v` row per line, y strictly increasing.
        path: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryKind {
    Wall,
    Inflow,
    Outflow,
    Periodic,
    /// inflow carrying a boundary-layer profile (see inflow_profile).
    ProfileInflow,
    /// a fixed supersonic inflow state (see inflow_state).
    SupersonicInflow,
    /// an inviscid slip wall: the ghost mirrors the normal velocity.
    /// the axisymmetric south boundary uses this as the symmetry axis.
    SlipWall,
    /// supersonic outflow: the ghost copies the interior state.
    SupersonicOutflow,
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
    /// advance each cell at its own stability bound instead of one global
    /// increment: a steady-state acceleration, not time-accurate.
    #[serde(default)]
    pub local_time_stepping: bool,
    /// worker threads for the threaded marches; 1 (the default) keeps the
    /// serial kernels.
    #[serde(default)]
    pub threads: usize,
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

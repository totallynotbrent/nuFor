//! conservative 1d euler time-marching: ghost-cell boundaries, finite-volume update, residual, logging.
//!
//! the heavy numerics (cfl, state conversion, hll flux) call into the fortran kernels; this module
//! owns the control flow: assembling face states, the conservative update, the residual, and the log.

use crate::cfl::cfl_dt;
use crate::flux::hll_flux;
use crate::state::cons_to_prim;
use crate::Error;

type FaceStates = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

/// boundary condition applied at one end of the 1D domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// open end; the ghost copies the interior state so waves leave freely.
    Transmissive,
    /// solid wall; the ghost mirrors the interior state with the normal velocity negated.
    Reflective,
}

/// why a solver run stopped, reported so a caller can diagnose without reading the log.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminationReason {
    /// the residual dropped to or below the requested tolerance.
    Converged,
    /// the accumulated simulated time reached t_end.
    TimeEnd,
    /// the step budget ran out before any other stop rule fired.
    MaxSteps,
    /// a non-physical state (non-finite density or non-positive pressure) appeared mid-run.
    BlowUp,
}

/// the three conservative variables over the cells.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservedState {
    /// density rho per cell.
    pub rho: Vec<f64>,
    /// momentum rho*u per cell.
    pub m: Vec<f64>,
    /// total energy rho*e_t per cell.
    pub e: Vec<f64>,
}

/// one recorded row of a solve: running time, time step, and residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EulerLog {
    /// step index, starting at one.
    pub step: usize,
    /// accumulated simulation time after this step.
    pub time: f64,
    /// the time step taken.
    pub dt: f64,
    /// the largest |dU| across the cells that this step produced.
    pub residual: f64,
}

/// the outcome of a solver run.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerResult {
    /// copy of the caller's state after the run.
    pub state: ConservedState,
    /// one row per step taken, oldest first.
    pub log: Vec<EulerLog>,
    /// the residual of the final step.
    pub residual: f64,
    /// the accumulated simulation time.
    pub time: f64,
    /// how many steps were taken.
    pub steps: usize,
    /// true when the residual reached the requested tolerance.
    pub converged: bool,
    /// why the run stopped.
    pub reason: TerminationReason,
}

/// the result of a physical-validity scan over a state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalCheck {
    /// true when every cell has finite, non-negative density and positive pressure.
    pub ok: bool,
    /// the smallest density seen across the cells.
    pub min_rho: f64,
    /// the smallest pressure seen where the density was usable.
    pub min_p: f64,
    /// the first cell index that failed the scan, if any.
    pub bad_cell: Option<usize>,
}

/// the run configuration: the gas, step control, boundaries, and stopping rules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EulerConfig {
    /// ratio of specific heats.
    pub gamma: f64,
    /// cfl number for the explicit step.
    pub cfl: f64,
    /// uniform cell width.
    pub dx: f64,
    /// boundary condition at the left end.
    pub left: Boundary,
    /// boundary condition at the right end.
    pub right: Boundary,
    /// maximum number of steps.
    pub max_steps: usize,
    /// stop when the simulated time reaches this.
    pub t_end: f64,
    /// stop when the residual drops to or below this.
    pub tol: f64,
}

/// ghost-cell conserved state for a boundary given the interior primitive state.
fn ghost_cons(boundary: Boundary, r: f64, u: f64, et: f64) -> (f64, f64, f64) {
    let up = match boundary {
        Boundary::Transmissive => u,
        Boundary::Reflective => -u,
    };
    (r, r * up, r * et)
}

/// assemble conserved left/right states at every face (n+1 of them) from interior cells and ghosts.
fn build_faces(
    state: &ConservedState,
    u: &[f64],
    et: &[f64],
    left: Boundary,
    right: Boundary,
) -> FaceStates {
    let n = state.rho.len();
    let gl = ghost_cons(left, state.rho[0], u[0], et[0]);
    let gr = ghost_cons(right, state.rho[n - 1], u[n - 1], et[n - 1]);
    let (mut lr, mut lm, mut le) = (
        Vec::with_capacity(n + 1),
        Vec::with_capacity(n + 1),
        Vec::with_capacity(n + 1),
    );
    let (mut rr, mut rm, mut re) = (
        Vec::with_capacity(n + 1),
        Vec::with_capacity(n + 1),
        Vec::with_capacity(n + 1),
    );
    for j in 0..=n {
        if j == 0 {
            lr.push(gl.0);
            lm.push(gl.1);
            le.push(gl.2);
        } else {
            lr.push(state.rho[j - 1]);
            lm.push(state.m[j - 1]);
            le.push(state.e[j - 1]);
        }
        if j == n {
            rr.push(gr.0);
            rm.push(gr.1);
            re.push(gr.2);
        } else {
            rr.push(state.rho[j]);
            rm.push(state.m[j]);
            re.push(state.e[j]);
        }
    }
    (lr, lm, le, rr, rm, re)
}

/// one explicit conservation step: HLL fluxes at every face, CFL time step,
/// conservative update, and the max|dU| residual. Advances the state in place.
pub fn advance(
    state: &mut ConservedState,
    gamma: f64,
    cfl: f64,
    dx: f64,
    left: Boundary,
    right: Boundary,
) -> Result<(f64, f64), Error> {
    let n = state.rho.len();
    if n < 2 || state.rho.len() != state.m.len() || state.rho.len() != state.e.len() {
        return Err(Error::InvalidArgs);
    }
    let step = cfl_dt(gamma, cfl, dx, &state.rho, &state.m, &state.e)?;
    let dt = step.dt;
    let (u, et) = cons_to_prim(&state.rho, &state.m, &state.e)?;
    let (lr, lm, le, rr, rm, re) = build_faces(state, &u, &et, left, right);
    let flux = hll_flux(gamma, &lr, &lm, &le, &rr, &rm, &re)?;
    let dtdx = dt / dx;
    let mut resid: f64 = 0.0;
    for i in 0..n {
        let drho = dtdx * (flux.rho[i + 1] - flux.rho[i]);
        let dm = dtdx * (flux.m[i + 1] - flux.m[i]);
        let de = dtdx * (flux.e[i + 1] - flux.e[i]);
        state.rho[i] -= drho;
        state.m[i] -= dm;
        state.e[i] -= de;
        resid = resid.max(drho.abs()).max(dm.abs()).max(de.abs());
    }
    Ok((dt, resid))
}

/// scan a conserved state for non-finite or non-physical cells.
///
/// a cell is valid when its density is finite and positive and its derived
/// pressure is finite and positive; the scan reports the running minima and the
/// first bad cell so a caller can either continue or stop cleanly.
pub fn check_physical(state: &ConservedState, gamma: f64) -> PhysicalCheck {
    let mut out = PhysicalCheck {
        ok: true,
        min_rho: f64::INFINITY,
        min_p: f64::INFINITY,
        bad_cell: None,
    };
    for i in 0..state.rho.len() {
        let r = state.rho[i];
        let (m, e) = (state.m[i], state.e[i]);
        out.min_rho = out.min_rho.min(r);
        let r_ok = r.is_finite() && r > 0.0;
        if r_ok {
            let u = m / r;
            let et = e / r;
            let p = (gamma - 1.0) * r * (et - 0.5 * u * u);
            out.min_p = out.min_p.min(p);
            if !(p.is_finite() && p > 0.0) {
                out.ok = false;
                out.bad_cell = out.bad_cell.or(Some(i));
            }
        } else {
            out.ok = false;
            out.bad_cell = out.bad_cell.or(Some(i));
        }
        if !(m.is_finite() && e.is_finite()) {
            out.ok = false;
            out.bad_cell = out.bad_cell.or(Some(i));
        }
    }
    out
}

/// run an explicit solve until the residual drops to tol, t_end elapses, or max_steps is used up.
pub fn euler_solve(state: &mut ConservedState, cfg: &EulerConfig) -> Result<EulerResult, Error> {
    let mut log: Vec<EulerLog> = Vec::new();
    let mut time = 0.0;
    let mut residual = f64::INFINITY;
    let mut reason = TerminationReason::MaxSteps;
    while log.len() < cfg.max_steps && time < cfg.t_end && residual > cfg.tol {
        let (dt, r) = advance(state, cfg.gamma, cfg.cfl, cfg.dx, cfg.left, cfg.right)?;
        let chk = check_physical(state, cfg.gamma);
        if !chk.ok {
            reason = TerminationReason::BlowUp;
            residual = r;
            break;
        }
        residual = r;
        time += dt;
        log.push(EulerLog {
            step: log.len() + 1,
            time,
            dt,
            residual: r,
        });
        if residual <= cfg.tol {
            reason = TerminationReason::Converged;
        } else if time >= cfg.t_end {
            reason = TerminationReason::TimeEnd;
        }
    }
    let steps = log.len();
    Ok(EulerResult {
        state: state.clone(),
        log,
        residual,
        time,
        steps,
        converged: reason == TerminationReason::Converged,
        reason,
    })
}

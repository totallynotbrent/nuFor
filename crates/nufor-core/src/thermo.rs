//! thermodynamic closure dispatch: one call site per solver stage that
//! needs pressure and sound speed, perfect gas or equilibrium air.
//!
//! the eqair branch lands with the curve fits (Srinivasan-Tannehill
//! fits, coefficients pending from the research pass); until then
//! every closure call runs the perfect-gas path. the dispatch keeps
//! the solver code model-agnostic: the march passes a `Thermo` and
//! never branches on the model itself.

use crate::eos2d::{eos_pressure2d, eos_sound_speed2d};
use crate::Error;

/// which thermodynamic closure a run uses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermoModel {
    /// constant gamma, constant R.
    Perfect { gamma: f64 },
    /// equilibrium air curve fits evaluated from (rho, e).
    EqAir,
}

/// pressure for a whole field under the chosen closure.
pub fn pressure(
    model: ThermoModel,
    rho: &[f64],
    e: &[f64],
    u: &[f64],
    v: &[f64],
) -> Result<Vec<f64>, Error> {
    match model {
        ThermoModel::Perfect { gamma } => {
            // the argument is internal energy per mass; velocity is
            // unused here but kept so both closures share one signature.
            let _ = (u, v);
            let n = rho.len();
            if n == 0 || e.len() != n || gamma <= 1.0 {
                return Err(Error::InvalidArgs);
            }
            let mut p = vec![0.0; n];
            for i in 0..n {
                let r = rho[i];
                if !r.is_finite() || r <= 0.0 {
                    return Err(Error::InvalidArgs);
                }
                p[i] = (gamma - 1.0) * r * e[i];
            }
            Ok(p)
        }
        ThermoModel::EqAir => crate::eqair_cea::pressure(rho, e, u, v),
    }
}

/// gas temperature for a whole field under the chosen closure. for
/// perfect gas this is p/(rho*R). for equilibrium air the fits carry
/// p directly but no T, so T is recovered from the perfect-gas form
/// with the effective gas constant the fits imply: R_eff = p/(rho*T)
/// requires T, which equilibrium air defines through the fit's own
/// partition — the practical closure is T = p/(rho*R0) with R0 the
/// cold-air value (the fits' e0 = R0*T0 anchor), which is exact in
/// the cold limit and the standard engineering approximation up the
/// Z ladder (dissociation changes the mean molecular weight; the
/// difference shows up as a correction the vibrational model would
/// carry, and for transport-property temperatures this accuracy
/// suffices).
pub fn temperature(model: ThermoModel, rho: &[f64], p: &[f64]) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    if n == 0 || p.len() != n {
        return Err(Error::InvalidArgs);
    }
    match model {
        ThermoModel::Perfect { .. } => {
            let mut out = vec![0.0; n];
            for i in 0..n {
                let r = rho[i];
                if !r.is_finite() || r <= 0.0 {
                    return Err(Error::InvalidArgs);
                }
                out[i] = p[i] / (r * 287.05);
            }
            Ok(out)
        }
        ThermoModel::EqAir => crate::eqair_cea::temperature(rho, p),
    }
}

pub fn sound_speed(
    model: ThermoModel,
    rho: &[f64],
    et: &[f64],
    u: &[f64],
    v: &[f64],
) -> Result<Vec<f64>, Error> {
    match model {
        ThermoModel::Perfect { gamma } => {
            let p = eos_pressure2d(gamma, rho, et, u, v)?;
            // the march historically clamped pressure at its wave-speed
            // use; transient cells can carry a tiny negative there, so
            // floor the field the same way instead of failing the step.
            let p_floored: Vec<f64> = p.iter().map(|x| x.max(1e-12)).collect();
            eos_sound_speed2d(gamma, &p_floored, rho)
        }
        ThermoModel::EqAir => {
            let n = rho.len();
            if n == 0 || et.len() != n || u.len() != n || v.len() != n {
                return Err(Error::InvalidArgs);
            }
            let mut e = vec![0.0; n];
            for i in 0..n {
                e[i] = et[i] - 0.5 * (u[i] * u[i] + v[i] * v[i]);
            }
            // the fused table read returns the cheap gamma_eff speed; the
            // legacy FD-based `eqair::sound_speed_from_energy` costs 3 fd
            // lookups per cell and misreads equilibrium relaxation as
            // acoustic stiffness, so route the cfl cap through the table.
            Ok(crate::eqair_cea::sound_speed(rho, &e, u, v)?)
        }
    }
}

/// the closure from the case config.
pub fn model_from_config(cfg: &nufor_config::CaseConfig) -> ThermoModel {
    use nufor_config::Eos;
    match cfg.physics.eos {
        Eos::Perfect => ThermoModel::Perfect {
            gamma: cfg.physics.gamma,
        },
        Eos::EqAir => ThermoModel::EqAir,
    }
}

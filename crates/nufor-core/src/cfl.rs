//! explicit time-step control from the cfl condition.

use crate::error::from_code;
use crate::ffi;
use crate::util::checked_len;
use crate::Error;

/// largest characteristic speed and the cfl step it yields.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CflStep {
    /// largest characteristic speed |u| + a over the cells.
    pub max_speed: f64,
    /// global explicit time step cfl*dx/s_max.
    pub dt: f64,
}

/// stable explicit time step for a uniform grid of width dx from the CFL condition.
pub fn cfl_dt(
    gamma: f64,
    cfl: f64,
    dx: f64,
    rho: &[f64],
    m: &[f64],
    e: &[f64],
) -> Result<CflStep, Error> {
    let n = checked_len(rho, m, e)?;
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    if !cfl.is_finite() || !(0.0 < cfl && cfl <= 1.0) {
        return Err(Error::InvalidArgs);
    }
    if !dx.is_finite() || dx <= 0.0 {
        return Err(Error::InvalidArgs);
    }
    let mut step = CflStep {
        max_speed: 0.0,
        dt: 0.0,
    };
    let mut err = 0i32;
    // SAFETY: equal-length non-empty slices sized by checked_len; the kernel touches [0,n) of each.
    unsafe {
        ffi::nfor_cfl_dt(
            gamma,
            cfl,
            dx,
            n as i32,
            rho.as_ptr(),
            m.as_ptr(),
            e.as_ptr(),
            &mut step.max_speed,
            &mut step.dt,
            &mut err,
        );
    }
    from_code(err)?;
    Ok(step)
}

//! ideal-gas equation of state helpers backed by the fortran kernels.

use crate::error::from_code;
use crate::ffi;
use crate::util::{checked_len, checked_len2};
use crate::Error;

/// ideal-gas pressure p = (gamma-1)*rho*e_int, computed in the fortran kernel.
pub fn eos_pressure(gamma: f64, rho: &[f64], et: &[f64], u: &[f64]) -> Result<Vec<f64>, Error> {
    let n = checked_len(rho, et, u)?;
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut p = vec![0.0; n];
    let mut err = 0i32;
    // SAFETY: equal-length non-empty slices; the kernel indexes [0,n) of each.
    unsafe {
        ffi::nfor_eos_pressure(
            gamma,
            n as i32,
            rho.as_ptr(),
            et.as_ptr(),
            u.as_ptr(),
            p.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok(p)
}

/// ideal-gas sound speed a = sqrt(gamma*p/rho); needs positive density and pressure.
pub fn eos_sound_speed(gamma: f64, rho: &[f64], p: &[f64]) -> Result<Vec<f64>, Error> {
    let n = checked_len2(rho, p)?;
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut a = vec![0.0; n];
    let mut err = 0i32;
    // SAFETY: as in eos_pressure; buffers below are exactly n entries.
    unsafe {
        ffi::nfor_eos_sound_speed(
            gamma,
            n as i32,
            rho.as_ptr(),
            p.as_ptr(),
            a.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok(a)
}

/// mach number M = |u|/a; the kernel rejects a non-positive sound speed.
pub fn eos_mach(u: &[f64], a: &[f64]) -> Result<Vec<f64>, Error> {
    let n = checked_len2(u, a)?;
    let mut mach = vec![0.0; n];
    let mut err = 0i32;
    // SAFETY: as in eos_pressure; buffers below are exactly n entries.
    unsafe {
        ffi::nfor_eos_mach(
            n as i32,
            u.as_ptr(),
            a.as_ptr(),
            mach.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok(mach)
}

/// ideal-gas temperature from p = rho*R*T; r, density, and pressure must be positive.
pub fn eos_temperature(r: f64, rho: &[f64], p: &[f64]) -> Result<Vec<f64>, Error> {
    let n = checked_len2(rho, p)?;
    if !r.is_finite() || r <= 0.0 {
        return Err(Error::InvalidArgs);
    }
    let mut t = vec![0.0; n];
    let mut err = 0i32;
    // SAFETY: as in eos_pressure; buffers below are exactly n entries.
    unsafe {
        ffi::nfor_eos_temperature(
            r,
            n as i32,
            rho.as_ptr(),
            p.as_ptr(),
            t.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok(t)
}

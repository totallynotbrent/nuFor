//! the hll approximate riemann solver with davis wave-speed estimates.

use crate::error::from_code;
use crate::ffi;
use crate::util::checked_len6;
use crate::Error;

/// HLL flux components at n faces: mass, momentum, and total-energy flux.
#[derive(Debug, Clone, PartialEq)]
pub struct HllFlux {
    /// mass flux rho*u per face.
    pub rho: Vec<f64>,
    /// momentum flux rho*u^2 + p per face.
    pub m: Vec<f64>,
    /// total-energy flux u*(E + p) per face.
    pub e: Vec<f64>,
}

/// HLL flux for n 1D euler faces using davis wave-speed estimates, no Roe average.
pub fn hll_flux(
    gamma: f64,
    rho_l: &[f64],
    m_l: &[f64],
    e_l: &[f64],
    rho_r: &[f64],
    m_r: &[f64],
    e_r: &[f64],
) -> Result<HllFlux, Error> {
    let n = checked_len6(rho_l, m_l, e_l, rho_r, m_r, e_r)?;
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut flux = HllFlux {
        rho: vec![0.0; n],
        m: vec![0.0; n],
        e: vec![0.0; n],
    };
    let mut err = 0i32;
    // SAFETY: all six slices are equal-length with n entries; the kernel stops at the first bad state.
    unsafe {
        ffi::nfor_hll_flux(
            gamma,
            n as i32,
            rho_l.as_ptr(),
            m_l.as_ptr(),
            e_l.as_ptr(),
            rho_r.as_ptr(),
            m_r.as_ptr(),
            e_r.as_ptr(),
            flux.rho.as_mut_ptr(),
            flux.m.as_mut_ptr(),
            flux.e.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok(flux)
}

//! conversions between conserved and primitive state, handled in the fortran kernel.

use crate::error::from_code;
use crate::ffi;
use crate::util::checked_len;
use crate::Error;

/// conserved state from primitives: m = rho*u, E = rho*e_t; density must be positive.
pub fn prim_to_cons(rho: &[f64], u: &[f64], et: &[f64]) -> Result<(Vec<f64>, Vec<f64>), Error> {
    let n = checked_len(rho, u, et)?;
    let mut m = vec![0.0; n];
    let mut e = vec![0.0; n];
    let mut err = 0i32;
    // SAFETY: equal-length non-empty slices; the kernel indexes [0,n) of each.
    unsafe {
        ffi::nfor_prim_to_cons(
            n as i32,
            rho.as_ptr(),
            u.as_ptr(),
            et.as_ptr(),
            m.as_mut_ptr(),
            e.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok((m, e))
}

/// primitives from conservative state, the inverse of prim_to_cons; density must be positive.
pub fn cons_to_prim(rho: &[f64], m: &[f64], e: &[f64]) -> Result<(Vec<f64>, Vec<f64>), Error> {
    let n = checked_len(rho, m, e)?;
    let mut u = vec![0.0; n];
    let mut et = vec![0.0; n];
    let mut err = 0i32;
    // SAFETY: as in prim_to_cons; buffers below are exactly n entries.
    unsafe {
        ffi::nfor_cons_to_prim(
            n as i32,
            rho.as_ptr(),
            m.as_ptr(),
            e.as_ptr(),
            u.as_mut_ptr(),
            et.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok((u, et))
}

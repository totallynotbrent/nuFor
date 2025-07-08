//! uniform 1d control-volume geometry built in the fortran kernel.

use crate::error::from_code;
use crate::ffi;
use crate::Error;
use std::os::raw::c_int;

/// uniform 1D control-volume geometry built in the fortran kernel.
#[derive(Debug, Clone, PartialEq)]
pub struct Grid1d {
    /// cell-center coordinates, one per control volume.
    pub centers: Vec<f64>,
    /// interface coordinates, one more than cells; the domain closes exactly.
    pub faces: Vec<f64>,
    /// uniform cell width (x_max - x_min)/n.
    pub dx: f64,
}

/// builds a uniform 1D grid of n cells over [xmin, xmax]; n >= 2 and xmax > xmin.
pub fn grid1d(n: usize, xmin: f64, xmax: f64) -> Result<Grid1d, Error> {
    if !(n >= 2 && n <= c_int::MAX as usize)
        || !(xmin.is_finite() && xmax.is_finite())
        || xmax <= xmin
    {
        return Err(Error::InvalidArgs);
    }
    let mut centers = vec![0.0; n];
    let mut faces = vec![0.0; n + 1];
    let mut dx: f64 = 0.0;
    let mut err: c_int = 0;
    // SAFETY: both buffers are exactly as long as the kernel fills them.
    unsafe {
        ffi::nfor_grid1d_init(
            n as c_int,
            xmin,
            xmax,
            centers.as_mut_ptr(),
            faces.as_mut_ptr(),
            &mut dx,
            &mut err,
        );
    }
    from_code(err)?;
    Ok(Grid1d { centers, faces, dx })
}

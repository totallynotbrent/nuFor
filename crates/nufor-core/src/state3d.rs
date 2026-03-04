//! 3d conserved state and its conversion to primitive variables.

use crate::Error;

/// four flat cell arrays returned together from a conserved conversion.
type Quad = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);
/// the five conserved variables over the cells of a 3d grid, flat row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservedState3d {
    pub rho: Vec<f64>,
    pub mx: Vec<f64>,
    pub my: Vec<f64>,
    pub mz: Vec<f64>,
    pub e: Vec<f64>,
}

/// primitive -> conserved (mx, my, mz, e), rejecting non-positive density.
pub fn prim_to_cons3d(
    rho: &[f64],
    u: &[f64],
    v: &[f64],
    w: &[f64],
    et: &[f64],
) -> Result<Quad, Error> {
    let n = rho.len();
    if n == 0 || u.len() != n || v.len() != n || w.len() != n || et.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut mx = vec![0.0; n];
    let mut my = vec![0.0; n];
    let mut mz = vec![0.0; n];
    let mut e = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        mx[i] = r * u[i];
        my[i] = r * v[i];
        mz[i] = r * w[i];
        e[i] = r * et[i];
    }
    Ok((mx, my, mz, e))
}

/// conserved -> primitive (u, v, w, total specific energy).
pub fn cons_to_prim3d(
    rho: &[f64],
    mx: &[f64],
    my: &[f64],
    mz: &[f64],
    e: &[f64],
) -> Result<Quad, Error> {
    let n = rho.len();
    if n == 0 || mx.len() != n || my.len() != n || mz.len() != n || e.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut u = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut w = vec![0.0; n];
    let mut et = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        u[i] = mx[i] / r;
        v[i] = my[i] / r;
        w[i] = mz[i] / r;
        et[i] = e[i] / r;
    }
    Ok((u, v, w, et))
}

//! 2d conserved state and its conversions to primitive variables.

use crate::Error;

/// three flat cell arrays returned together from a state conversion.
type Triple = (Vec<f64>, Vec<f64>, Vec<f64>);
/// the four conserved variables over the cells of a 2d grid, flat row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservedState2d {
    /// density rho per cell.
    pub rho: Vec<f64>,
    /// x-momentum rho*u per cell.
    pub mx: Vec<f64>,
    /// y-momentum rho*v per cell.
    pub my: Vec<f64>,
    /// total energy rho*e_t per cell.
    pub e: Vec<f64>,
}

/// primitive -> conserved, rejecting non-positive or non-finite density.
pub fn prim_to_cons2d(rho: &[f64], u: &[f64], v: &[f64], et: &[f64]) -> Result<Triple, Error> {
    let n = rho.len();
    if n == 0 || u.len() != n || v.len() != n || et.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut mx = vec![0.0; n];
    let mut my = vec![0.0; n];
    let mut e = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        mx[i] = r * u[i];
        my[i] = r * v[i];
        e[i] = r * et[i];
    }
    Ok((mx, my, e))
}

/// conserved -> primitive (u, v, total specific energy), rejecting non-positive density.
pub fn cons_to_prim2d(rho: &[f64], mx: &[f64], my: &[f64], e: &[f64]) -> Result<Triple, Error> {
    let n = rho.len();
    if n == 0 || mx.len() != n || my.len() != n || e.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut u = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut et = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        u[i] = mx[i] / r;
        v[i] = my[i] / r;
        et[i] = e[i] / r;
    }
    Ok((u, v, et))
}

/// a physical-validity report over a 2d conserved state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalCheck2d {
    /// true when every cell has finite, positive density and finite conserved values.
    pub ok: bool,
    /// the smallest density seen.
    pub min_rho: f64,
    /// the first bad cell index, if any.
    pub bad_cell: Option<usize>,
}

/// scan a 2d conserved state for non-finite or non-positive density.
///
/// pressure is derived elsewhere (see eos); this catches density and
/// continuity failures, which are the ones that corrupt the momentum/energy split.
pub fn check_physical2d(state: &ConservedState2d) -> PhysicalCheck2d {
    let mut chk = PhysicalCheck2d {
        ok: true,
        min_rho: f64::INFINITY,
        bad_cell: None,
    };
    for i in 0..state.rho.len() {
        let r = state.rho[i];
        chk.min_rho = chk.min_rho.min(r);
        let good = r.is_finite()
            && r > 0.0
            && state.mx[i].is_finite()
            && state.my[i].is_finite()
            && state.e[i].is_finite();
        if !good {
            chk.ok = false;
            chk.bad_cell = chk.bad_cell.or(Some(i));
        }
    }
    chk
}

//! ideal-gas relations for the 3d euler equations.

use crate::Error;

/// pressure from density, total specific energy, and velocity.
pub fn eos_pressure3d(
    gamma: f64,
    rho: &[f64],
    et: &[f64],
    u: &[f64],
    v: &[f64],
    w: &[f64],
) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    if n == 0 || et.len() != n || u.len() != n || v.len() != n || w.len() != n || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut p = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        p[i] = (gamma - 1.0) * r * (et[i] - 0.5 * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
    }
    Ok(p)
}

/// sound speed a = sqrt(gamma p / rho).
pub fn eos_sound_speed3d(gamma: f64, p: &[f64], rho: &[f64]) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    if n == 0 || p.len() != n || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut a = vec![0.0; n];
    for i in 0..n {
        let (pr, r) = (p[i], rho[i]);
        if !r.is_finite() || r <= 0.0 || !pr.is_finite() || pr <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        a[i] = (gamma * pr / r).sqrt();
    }
    Ok(a)
}

/// mach number |v|/a from velocity and sound speed.
pub fn eos_mach3d(u: &[f64], v: &[f64], w: &[f64], a: &[f64]) -> Result<Vec<f64>, Error> {
    let n = u.len();
    if n == 0 || v.len() != n || w.len() != n || a.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut mach = vec![0.0; n];
    for i in 0..n {
        if a[i] <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        mach[i] = (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]).sqrt() / a[i];
    }
    Ok(mach)
}

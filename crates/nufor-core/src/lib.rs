//! rust boundary over the fortran numerical library: marshals arguments, checks the boundary, and maps error codes.

use std::ffi::c_char;
use std::os::raw::c_int;

/// error codes from the fortran kernels; must match nuforkernels.f90.
mod codes {
    use crate::c_int;

    pub const OK: c_int = 0;
    pub const E_ARGS: c_int = 1;
    pub const E_DATA: c_int = 2;
}

/// structured boundary error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// invalid argument (count below one or arrays not equal length).
    InvalidArgs,
    /// the kernel reported a numerical failure on valid inputs.
    KernelFailure,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::InvalidArgs => write!(f, "invalid arguments across the FFI boundary"),
            Error::KernelFailure => write!(f, "fortran kernel reported a numerical failure"),
        }
    }
}

impl std::error::Error for Error {}

fn from_code(code: c_int) -> Result<(), Error> {
    match code {
        codes::OK => Ok(()),
        codes::E_ARGS => Err(Error::InvalidArgs),
        codes::E_DATA => Err(Error::KernelFailure),
        _ => Err(Error::KernelFailure),
    }
}

// raw C ABI into the statically linked fortran library (bind(C) names).
mod ffi {
    use super::*;

    #[link(name = "nuforkernels", kind = "static")]
    extern "C" {
        pub fn nfor_version(ver: *mut c_char, ver_len: c_int, err: *mut c_int);
        pub fn nfor_saxpy(count: c_int, alpha: f64, x: *const f64, y: *mut f64, err: *mut c_int);
        pub fn nfor_grid1d_init(
            n: c_int,
            xmin: f64,
            xmax: f64,
            centers: *mut f64,
            faces: *mut f64,
            dx: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_prim_to_cons(
            n: c_int,
            rho: *const f64,
            u: *const f64,
            et: *const f64,
            m: *mut f64,
            e: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_cons_to_prim(
            n: c_int,
            rho: *const f64,
            m: *const f64,
            e: *const f64,
            u: *mut f64,
            et: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_eos_pressure(
            gamma: f64,
            n: c_int,
            rho: *const f64,
            et: *const f64,
            u: *const f64,
            p: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_eos_sound_speed(
            gamma: f64,
            n: c_int,
            rho: *const f64,
            p: *const f64,
            a: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_eos_mach(
            n: c_int,
            u: *const f64,
            a: *const f64,
            mach: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_eos_temperature(
            r: f64,
            n: c_int,
            rho: *const f64,
            p: *const f64,
            t: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_hll_flux(
            gamma: f64,
            n: c_int,
            rho_l: *const f64,
            m_l: *const f64,
            e_l: *const f64,
            rho_r: *const f64,
            m_r: *const f64,
            e_r: *const f64,
            f_rho: *mut f64,
            f_m: *mut f64,
            f_e: *mut f64,
            err: *mut c_int,
        );
        pub fn nfor_cfl_dt(
            gamma: f64,
            cfl: f64,
            dx: f64,
            n: c_int,
            rho: *const f64,
            m: *const f64,
            e: *const f64,
            s_max: *mut f64,
            dt: *mut f64,
            err: *mut c_int,
        );
    }
}

/// two slices must agree on a positive length.
fn checked_len2(a: &[f64], b: &[f64]) -> Result<usize, Error> {
    if a.is_empty() || a.len() != b.len() || a.len() > c_int::MAX as usize {
        return Err(Error::InvalidArgs);
    }
    Ok(a.len())
}

/// three slices must agree on a positive length.
fn checked_len(a: &[f64], b: &[f64], c: &[f64]) -> Result<usize, Error> {
    if a.len() != c.len() {
        return Err(Error::InvalidArgs);
    }
    checked_len2(a, b)
}

/// six state slices must agree on a positive length.
fn checked_len6(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
    e: &[f64],
    f: &[f64],
) -> Result<usize, Error> {
    if a.len() != d.len() || a.len() != e.len() || a.len() != f.len() {
        return Err(Error::InvalidArgs);
    }
    checked_len(a, b, c)
}

/// version string reported by the fortran kernel library.
pub fn version() -> String {
    let mut buf = [0 as c_char; 64];
    let mut err: c_int = 0;
    // SAFETY: buf is a writable 64-byte buffer; the kernel writes at most its length.
    unsafe {
        ffi::nfor_version(buf.as_mut_ptr(), buf.len() as c_int, &mut err);
    }
    debug_assert_eq!(err, codes::OK, "nfor_version failed");
    let bytes: Vec<u8> = buf
        .iter()
        .map(|&c| c as u8)
        .take_while(|&b| b != 0)
        .collect();
    String::from_utf8_lossy(&bytes).into_owned()
}

/// y = alpha*x + y elementwise (BLAS-like saxpy) over equal-length slices.
pub fn saxpy(alpha: f64, x: &[f64], y: &mut [f64]) -> Result<(), Error> {
    if x.len() != y.len() {
        return Err(Error::InvalidArgs);
    }
    if x.len() > c_int::MAX as usize {
        return Err(Error::InvalidArgs);
    }
    if x.is_empty() {
        return Ok(());
    }
    let mut err: c_int = 0;
    // SAFETY: both slices are non-empty, equal-length, contiguous f64 arrays.
    unsafe {
        ffi::nfor_saxpy(
            x.len() as c_int,
            alpha,
            x.as_ptr(),
            y.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)
}

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
    if n < 2 || n > c_int::MAX as usize || !(xmin.is_finite() && xmax.is_finite()) || xmax <= xmin {
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

/// conserved state from primitives: m = rho*u, E = rho*e_t; density must be positive.
pub fn prim_to_cons(rho: &[f64], u: &[f64], et: &[f64]) -> Result<(Vec<f64>, Vec<f64>), Error> {
    let n = checked_len(rho, u, et)?;
    let mut m = vec![0.0; n];
    let mut e = vec![0.0; n];
    let mut err: c_int = 0;
    // SAFETY: equal-length non-empty slices; the kernel indexes [0,n) of each.
    unsafe {
        ffi::nfor_prim_to_cons(
            n as c_int,
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
    let mut err: c_int = 0;
    // SAFETY: as in prim_to_cons; buffers below are exactly n entries.
    unsafe {
        ffi::nfor_cons_to_prim(
            n as c_int,
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

/// ideal-gas pressure p = (gamma-1)*rho*e_int, computed in the fortran kernel.
pub fn eos_pressure(gamma: f64, rho: &[f64], et: &[f64], u: &[f64]) -> Result<Vec<f64>, Error> {
    let n = checked_len(rho, et, u)?;
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut p = vec![0.0; n];
    let mut err: c_int = 0;
    // SAFETY: equal-length non-empty slices; the kernel indexes [0,n) of each.
    unsafe {
        ffi::nfor_eos_pressure(
            gamma,
            n as c_int,
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
    let mut err: c_int = 0;
    // SAFETY: as in `eos_pressure`; buffers below are exactly `n` entries.
    unsafe {
        ffi::nfor_eos_sound_speed(
            gamma,
            n as c_int,
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
    let mut err: c_int = 0;
    // SAFETY: as in `eos_pressure`; buffers below are exactly `n` entries.
    unsafe {
        ffi::nfor_eos_mach(
            n as c_int,
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
    let mut err: c_int = 0;
    // SAFETY: as in `eos_pressure`; buffers below are exactly `n` entries.
    unsafe {
        ffi::nfor_eos_temperature(
            r,
            n as c_int,
            rho.as_ptr(),
            p.as_ptr(),
            t.as_mut_ptr(),
            &mut err,
        );
    }
    from_code(err)?;
    Ok(t)
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
    let mut err: c_int = 0;
    // SAFETY: all six slices are equal-length with n entries; the kernel stops at the first bad state.
    unsafe {
        ffi::nfor_hll_flux(
            gamma,
            n as c_int,
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
    if !cfl.is_finite() || cfl <= 0.0 || cfl > 1.0 {
        return Err(Error::InvalidArgs);
    }
    if !dx.is_finite() || dx <= 0.0 {
        return Err(Error::InvalidArgs);
    }
    let mut step = CflStep {
        max_speed: 0.0,
        dt: 0.0,
    };
    let mut err: c_int = 0;
    // SAFETY: equal-length non-empty slices sized by `checked_len`; the kernel
    // touches indices [0, n) of each and stops at the first invalid state.
    unsafe {
        ffi::nfor_cfl_dt(
            gamma,
            cfl,
            dx,
            n as c_int,
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

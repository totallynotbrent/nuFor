//! Rust half of the mixed-language boundary (PLAN step 2).
//!
//! All numerics live in the Fortran kernel library; this crate only marshals
//! arguments, checks the boundary, and maps error codes. ABI rules follow the
//! FFI policy in `docs/research/ffi-boundary.md` (spec 39): bind(C), contiguous
//! arrays with explicit lengths, no Fortran derived types, structured errors.

use std::ffi::c_char;
use std::os::raw::c_int;

/// Error codes produced by the Fortran kernels (must match `nuforkernels.f90`).
mod codes {
    use crate::c_int;

    pub const OK: c_int = 0;
    pub const E_ARGS: c_int = 1;
    pub const E_DATA: c_int = 2;
}

/// Structured boundary error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// Invalid argument (count below one, arrays not equal length).
    InvalidArgs,
    /// The kernel reported a numerical failure for valid inputs.
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

// Raw C ABI into the statically linked Fortran library (bind(C) names).
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
    }
}

/// Two slices must agree on a positive length representable in `c_int`.
fn checked_len2(a: &[f64], b: &[f64]) -> Result<usize, Error> {
    if a.is_empty() || a.len() != b.len() || a.len() > c_int::MAX as usize {
        return Err(Error::InvalidArgs);
    }
    Ok(a.len())
}

/// Three slices must agree on a positive length representable in `c_int`.
fn checked_len(a: &[f64], b: &[f64], c: &[f64]) -> Result<usize, Error> {
    if a.len() != c.len() {
        return Err(Error::InvalidArgs);
    }
    checked_len2(a, b)
}

/// Six state slices must agree on a positive length representable in `c_int`.
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

/// Version string reported by the Fortran kernel library.
pub fn version() -> String {
    let mut buf = [0 as c_char; 64];
    let mut err: c_int = 0;
    // SAFETY: `buf` is a writable 64-byte buffer and its length is passed
    // explicitly; the kernel writes at most that many bytes plus a NUL terminator.
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

/// `y = alpha * x + y` elementwise (BLAS-like saxp) over equal-length slices.
///
/// Contiguous slices map directly to the C array convention; the Fortran kernel
/// owns the arithmetic and reports a structured error code on bad input.
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
    // SAFETY: both slices are non-empty, equal length, and each points to a
    // valid contiguous f64 array; the kernel only touches indices [0, len).
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

/// Uniform 1D control-volume geometry (spec 18), built in the Fortran kernel.
#[derive(Debug, Clone, PartialEq)]
pub struct Grid1d {
    /// Cell-center coordinates, one per control volume.
    pub centers: Vec<f64>,
    /// Interface coordinates, one more than cells; the domain closes exactly.
    pub faces: Vec<f64>,
    /// Uniform cell width `(x_max - x_min) / n`.
    pub dx: f64,
}

/// Builds a uniform 1D grid of `n` cells over `[xmin, xmax]` in the Fortran
/// kernel layer. `n` must be at least 2 and `xmax` must exceed `xmin`.
pub fn grid1d(n: usize, xmin: f64, xmax: f64) -> Result<Grid1d, Error> {
    if n < 2 || n > c_int::MAX as usize || !(xmin.is_finite() && xmax.is_finite()) || xmax <= xmin {
        return Err(Error::InvalidArgs);
    }
    let mut centers = vec![0.0; n];
    let mut faces = vec![0.0; n + 1];
    let mut dx: f64 = 0.0;
    let mut err: c_int = 0;
    // SAFETY: both buffers have the exact lengths the kernel fills
    // (`centers` n entries, `faces` n+1); it writes nothing beyond them.
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

/// Conservative state from primitives: `m = rho * u`, `E = rho * e_t`.
///
/// `e_t` is the total specific energy; pressure recovery is the EOS step's
/// job. Density must be positive or the kernel reports a numerical failure.
pub fn prim_to_cons(rho: &[f64], u: &[f64], et: &[f64]) -> Result<(Vec<f64>, Vec<f64>), Error> {
    let n = checked_len(rho, u, et)?;
    let mut m = vec![0.0; n];
    let mut e = vec![0.0; n];
    let mut err: c_int = 0;
    // SAFETY: equal-length non-empty slices sized by `checked_len`; the kernel
    // touches indices [0, n) of each and stops at the first bad density.
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

/// Primitives from conservative state, the inverse of [`prim_to_cons`]:
/// `u = m / rho`, `e_t = E / rho`. Density must be positive.
pub fn cons_to_prim(rho: &[f64], m: &[f64], e: &[f64]) -> Result<(Vec<f64>, Vec<f64>), Error> {
    let n = checked_len(rho, m, e)?;
    let mut u = vec![0.0; n];
    let mut et = vec![0.0; n];
    let mut err: c_int = 0;
    // SAFETY: as in `prim_to_cons`; buffers below are exactly `n` entries.
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

/// Ideal-gas pressure from primitives (spec 25): `p = (gamma - 1) * rho * e_int`
/// with `e_int = e_t - u^2/2`, computed in the Fortran kernel.
///
/// `gamma` must exceed 1. Non-positive density or internal energy is an
/// inadmissible Euler state, so the kernel reports a numerical failure.
pub fn eos_pressure(gamma: f64, rho: &[f64], et: &[f64], u: &[f64]) -> Result<Vec<f64>, Error> {
    let n = checked_len(rho, et, u)?;
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(Error::InvalidArgs);
    }
    let mut p = vec![0.0; n];
    let mut err: c_int = 0;
    // SAFETY: equal-length non-empty slices sized by `checked_len`; the kernel
    // touches indices [0, n) of each and rejects the first invalid state.
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

/// Ideal-gas sound speed (spec 25): `a = sqrt(gamma * p / rho)`.
///
/// `gamma` must exceed 1 and the state must have positive density and
/// pressure so the speed of sound stays real and positive.
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

/// Mach number from flow velocity and sound speed (spec 25): `M = |u| / a`.
///
/// The kernel rejects a non-positive sound speed, which would make the Mach
/// number meaningless.
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

/// Ideal-gas temperature from the thermal law `p = rho * R * T` (spec 25).
///
/// `r` is the specific gas constant from the case config in J/(kg K) and must
/// be positive; density and pressure must be positive too.
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

/// HLL numerical flux (spec 52, 90) of the 1D Euler equations for `n` face
/// Riemann problems taken at once.
///
/// Left states `rho_l, m_l, e_l` and right states `rho_r, m_r, e_r` are
/// conserved variables (`m = rho * u`, `e` total energy per volume); the
/// kernel recovers primitives, pressure, and sound speed through the
/// ideal-gas EOS and returns the three HLL flux components. Wave speeds are
/// the Davis estimates, so no Roe average is needed (see
/// `docs/numerics/flux-hll.md`). `gamma` must exceed 1 and every state must
/// be admissible (positive density and internal energy), or the kernel
/// reports a structured failure.
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
    // SAFETY: equal-length non-empty slices sized by `checked_len6`; every
    // buffer has exactly `n` entries and the kernel stops at the first bad
    // state before writing to any output slice.
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

/// HLL flux components at `n` faces: mass, momentum, and total-energy flux.
#[derive(Debug, Clone, PartialEq)]
pub struct HllFlux {
    /// Mass flux `rho * u` per face.
    pub rho: Vec<f64>,
    /// Momentum flux `rho * u^2 + p` per face.
    pub m: Vec<f64>,
    /// Total-energy flux `u * (E + p)` per face.
    pub e: Vec<f64>,
}

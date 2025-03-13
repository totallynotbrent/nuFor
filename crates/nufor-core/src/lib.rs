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
    }
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

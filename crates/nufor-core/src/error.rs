//! error type surfaced across the rust boundary, mapping kernel codes to a rust Error.

use std::os::raw::c_int;

/// error codes from the fortran kernels; must match nuforkernels.f90.
pub(crate) mod codes {
    use super::c_int;

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

/// maps a kernel return code onto the structured error.
pub(crate) fn from_code(code: c_int) -> Result<(), Error> {
    match code {
        codes::OK => Ok(()),
        codes::E_ARGS => Err(Error::InvalidArgs),
        codes::E_DATA => Err(Error::KernelFailure),
        _ => Err(Error::KernelFailure),
    }
}

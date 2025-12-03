//! rust boundary over the fortran numerical library: marshals arguments, checks the boundary, and maps error codes.
//!
//! the crate is split into focused modules: error, ffi, shared slice checks, and one
//! module per numerical area (grid, state conversion, eos, hll flux, cfl, and the solver).

use std::ffi::c_char;
use std::os::raw::c_int;

mod cfl;
mod eos;
mod eos2d;
mod error;
mod exact;
mod ffi;
mod flux;
mod grid;
mod grid2d;
mod h5;
mod output;
mod restart;
mod solver;
mod state;
mod state2d;
mod util;

pub use cfl::{cfl_dt, CflStep};
pub use eos::{eos_mach, eos_pressure, eos_sound_speed, eos_temperature};
pub use eos2d::{eos_mach2d, eos_pressure2d, eos_sound_speed2d};
pub use error::Error;
pub use exact::{riemann, ExactSolution, PrimState};
pub use flux::{hll_flux, HllFlux};
pub use grid::{grid1d, Grid1d};
pub use grid2d::{grid2d, Grid2d};
pub use h5::{read_h5, write_h5};
pub use output::{write_csv, write_vtk, OutputState};
pub use restart::{read_restart, write_restart, RestartData};
pub use solver::{
    advance, check_physical, euler_solve, Boundary, ConservedState, EulerConfig, EulerLog,
    EulerResult, PhysicalCheck, TerminationReason,
};
pub use state::{cons_to_prim, prim_to_cons};
pub use state2d::{
    check_physical2d, cons_to_prim2d, prim_to_cons2d, ConservedState2d, PhysicalCheck2d,
};

use error::{codes, from_code};

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
    if x.len() != y.len() || x.len() > c_int::MAX as usize {
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

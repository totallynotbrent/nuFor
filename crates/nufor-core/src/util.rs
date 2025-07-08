//! shared slice-length checks that keep the arrays across the boundary consistent.

use crate::error::Error;
use std::os::raw::c_int;

/// two slices must agree on a positive length.
pub(crate) fn checked_len2(a: &[f64], b: &[f64]) -> Result<usize, Error> {
    if a.is_empty() || a.len() != b.len() || a.len() > c_int::MAX as usize {
        return Err(Error::InvalidArgs);
    }
    Ok(a.len())
}

/// three slices must agree on a positive length.
pub(crate) fn checked_len(a: &[f64], b: &[f64], c: &[f64]) -> Result<usize, Error> {
    if a.len() != c.len() {
        return Err(Error::InvalidArgs);
    }
    checked_len2(a, b)
}

/// six state slices must agree on a positive length.
pub(crate) fn checked_len6(
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

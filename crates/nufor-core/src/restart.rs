//! restart: save and reload the full run state as a compact binary file.
//!
//! the format is a magic header carrying an explicit schema version, the block
//! shape, and the conserved arrays as little-endian f64, so a restart recovers
//! the solution bit-for-bit. the reader is strict: it validates the schema
//! version, the exact file size, and the header values, so a corrupt or
//! future-schema file fails loudly instead of being silently misread. the
//! legacy pre-versioned layout is still accepted and reported as version 0.

use std::path::Path;

use crate::solver::ConservedState;
use crate::Error;

/// the current restart schema version.
pub const RESTART_VERSION: u8 = 1;
/// magic that prefixes the current, versioned layout.
const MAGIC: &[u8; 4] = b"NUFR";
/// magic of the legacy, pre-versioned layout, still readable as version 0.
const LEGACY_MAGIC: &[u8; 6] = b"NUFR1\0";
/// magic of the 2d layout (rho, mx, my, e over an nx*ny grid).
const MAGIC2D: &[u8; 4] = b"NUF2";
/// an upper bound on n so a corrupt size field cannot drive a huge allocation.
const MAX_CELLS: u64 = 100_000_000;
/// reserved header bytes for n, gamma, time, and step.
const HDR: usize = 32;
/// header bytes for the 2d layout: n, nx, ny, gamma, time, step.
const HDR2D: usize = 48;

/// the full state recovered from a restart file.
#[derive(Debug, Clone, PartialEq)]
pub struct RestartData {
    /// the conserved state.
    pub state: ConservedState,
    /// ratio of specific heats.
    pub gamma: f64,
    /// the simulated time the run reached.
    pub time: f64,
    /// the number of steps taken.
    pub step: u64,
    /// the restart schema version the file was written with.
    pub version: u8,
}

/// the full 2d state recovered from a restart file.
#[derive(Debug, Clone, PartialEq)]
pub struct RestartData2d {
    /// the conserved 2d state (rho, mx, my, e, flat row-major over nx*ny).
    pub state: crate::state2d::ConservedState2d,
    /// grid columns.
    pub nx: usize,
    /// grid rows.
    pub ny: usize,
    /// ratio of specific heats.
    pub gamma: f64,
    /// the simulated time the run reached.
    pub time: f64,
    /// the number of steps taken.
    pub step: u64,
}

/// writes a full 2d state to a restart file in the NUF2 layout.
pub fn write_restart2d(
    path: &Path,
    state: &crate::state2d::ConservedState2d,
    nx: usize,
    ny: usize,
    gamma: f64,
    time: f64,
    step: u64,
) -> Result<(), Error> {
    let n = state.rho.len();
    if n < 2
        || n != state.mx.len()
        || n != state.my.len()
        || n != state.e.len()
        || n != nx * ny
        || nx < 2
        || ny < 2
    {
        return Err(Error::InvalidArgs);
    }
    if (n as u64) > MAX_CELLS {
        return Err(Error::InvalidArgs);
    }
    let mut out = Vec::with_capacity(MAGIC2D.len() + HDR2D + 4 * n * 8);
    out.extend_from_slice(MAGIC2D);
    out.extend_from_slice(&(n as u64).to_le_bytes());
    out.extend_from_slice(&(nx as u64).to_le_bytes());
    out.extend_from_slice(&(ny as u64).to_le_bytes());
    out.extend_from_slice(&gamma.to_le_bytes());
    out.extend_from_slice(&time.to_le_bytes());
    out.extend_from_slice(&step.to_le_bytes());
    for a in [&state.rho, &state.mx, &state.my, &state.e] {
        for v in a {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    std::fs::write(path, &out).map_err(|_| Error::KernelFailure)
}

/// reads a 2d restart file back, validating magic, schema shape, and size.
pub fn read_restart2d(path: &Path) -> Result<RestartData2d, Error> {
    let bytes = std::fs::read(path).map_err(|_| Error::KernelFailure)?;
    if bytes.len() < MAGIC2D.len() + HDR2D || &bytes[..4] != MAGIC2D {
        return Err(Error::InvalidArgs);
    }
    let rd = &bytes[MAGIC2D.len()..MAGIC2D.len() + HDR2D];
    let n = u64::from_le_bytes(rd[0..8].try_into().unwrap()) as usize;
    let nx = u64::from_le_bytes(rd[8..16].try_into().unwrap()) as usize;
    let ny = u64::from_le_bytes(rd[16..24].try_into().unwrap()) as usize;
    let gamma = f64::from_le_bytes(rd[24..32].try_into().unwrap());
    let time = f64::from_le_bytes(rd[32..40].try_into().unwrap());
    if !(2..=MAX_CELLS as usize).contains(&n)
        || nx < 2
        || ny < 2
        || n != nx * ny
        || !gamma.is_finite()
        || !time.is_finite()
        || gamma <= 1.0
    {
        return Err(Error::InvalidArgs);
    }
    if bytes.len() != MAGIC2D.len() + HDR2D + 4 * n * 8 {
        return Err(Error::InvalidArgs);
    }
    let mut cur = MAGIC2D.len() + HDR2D;
    let mut take = |count: usize| -> Vec<f64> {
        let mut v = Vec::with_capacity(count);
        for _ in 0..count {
            let slice: [u8; 8] = bytes[cur..cur + 8].try_into().unwrap();
            cur += 8;
            v.push(f64::from_le_bytes(slice));
        }
        v
    };
    let rho = take(n);
    let mx = take(n);
    let my = take(n);
    let e = take(n);
    Ok(RestartData2d {
        state: crate::state2d::ConservedState2d { rho, mx, my, e },
        nx,
        ny,
        gamma,
        time,
        step: 0,
    })
}

/// writes a full state to a restart file in the current versioned format.
pub fn write_restart(
    path: &Path,
    state: &ConservedState,
    gamma: f64,
    time: f64,
    step: u64,
) -> Result<(), Error> {
    let n = state.rho.len();
    if n < 2 || state.rho.len() != state.m.len() || state.rho.len() != state.e.len() {
        return Err(Error::InvalidArgs);
    }
    if (n as u64) > MAX_CELLS {
        return Err(Error::InvalidArgs);
    }
    let mut out = Vec::with_capacity(MAGIC.len() + 1 + HDR + 3 * n * 8);
    out.extend_from_slice(MAGIC);
    out.push(RESTART_VERSION);
    out.extend_from_slice(&(n as u64).to_le_bytes());
    out.extend_from_slice(&gamma.to_le_bytes());
    out.extend_from_slice(&time.to_le_bytes());
    out.extend_from_slice(&step.to_le_bytes());
    for a in [&state.rho, &state.m, &state.e] {
        for v in a {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    std::fs::write(path, &out).map_err(|_| Error::KernelFailure)
}

/// reads a restart file back into a full state, validating schema and shape.
pub fn read_restart(path: &Path) -> Result<RestartData, Error> {
    let bytes = std::fs::read(path).map_err(|_| Error::KernelFailure)?;

    // detect the layout: versioned magic plus an explicit version byte, or the
    // legacy magic; anything else is not a restart file.
    let (version, base) = if bytes.len() >= LEGACY_MAGIC.len() && &bytes[..6] == LEGACY_MAGIC {
        (0u8, 6usize)
    } else if bytes.len() >= 5 && &bytes[..4] == MAGIC {
        (bytes[4], 5usize)
    } else {
        return Err(Error::InvalidArgs);
    };

    // a file from a newer schema version than this build cannot be trusted.
    if version > RESTART_VERSION {
        return Err(Error::IncompatibleVersion);
    }

    let want = base.checked_add(HDR).ok_or(Error::InvalidArgs)?;
    if bytes.len() < want {
        return Err(Error::InvalidArgs);
    }
    let rd = &bytes[base..base + HDR];
    let n = u64::from_le_bytes(rd[0..8].try_into().unwrap());
    if !(2..=MAX_CELLS).contains(&n) {
        return Err(Error::InvalidArgs);
    }
    let gamma = f64::from_le_bytes(rd[8..16].try_into().unwrap());
    let time = f64::from_le_bytes(rd[16..24].try_into().unwrap());
    let step = u64::from_le_bytes(rd[24..32].try_into().unwrap());
    let n = n as usize;

    // strict: the file must be exactly the expected size, no truncation or
    // trailing bytes that would hide a mismatched layout.
    let expected = base + HDR + 3 * n * 8;
    if bytes.len() != expected {
        return Err(Error::InvalidArgs);
    }
    // a well-formed file still must carry meaningful physics.
    if !(gamma.is_finite() && time.is_finite() && gamma > 1.0) {
        return Err(Error::InvalidArgs);
    }

    let mut cur = base + HDR;
    let mut take = |count: usize| -> Vec<f64> {
        let mut v = Vec::with_capacity(count);
        for _ in 0..count {
            let slice: [u8; 8] = bytes[cur..cur + 8].try_into().unwrap();
            cur += 8;
            v.push(f64::from_le_bytes(slice));
        }
        v
    };
    let rho = take(n);
    let m = take(n);
    let e = take(n);
    Ok(RestartData {
        state: ConservedState { rho, m, e },
        gamma,
        time,
        step,
        version,
    })
}

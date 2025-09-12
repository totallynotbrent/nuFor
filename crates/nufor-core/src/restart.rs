//! restart: save and reload the full run state as a compact binary file.
//!
//! the format is a tiny header plus the conserved arrays as little-endian f64,
//! so a restart recovers the solution bit-for-bit; the header validates the
//! block sizes on load so a corrupt file fails cleanly.

use std::path::Path;

use crate::solver::ConservedState;
use crate::Error;

const MAGIC: &[u8; 6] = b"NUFR1\0";
const MAX_CELLS: u64 = 100_000_000;

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
}

/// writes a full state to a restart file.
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
    let mut out = Vec::with_capacity(MAGIC.len() + 32 + 3 * n * 8);
    out.extend_from_slice(MAGIC);
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

/// reads a restart file back into a full state.
pub fn read_restart(path: &Path) -> Result<RestartData, Error> {
    let bytes = std::fs::read(path).map_err(|_| Error::KernelFailure)?;
    if bytes.len() < MAGIC.len() + 32 {
        return Err(Error::InvalidArgs);
    }
    if &bytes[..MAGIC.len()] != MAGIC {
        return Err(Error::InvalidArgs);
    }
    let rd = &bytes[MAGIC.len()..];
    let n = u64::from_le_bytes(rd[0..8].try_into().unwrap());
    if !(2..=MAX_CELLS).contains(&n) {
        return Err(Error::InvalidArgs);
    }
    let gamma = f64::from_le_bytes(rd[8..16].try_into().unwrap());
    let time = f64::from_le_bytes(rd[16..24].try_into().unwrap());
    let step = u64::from_le_bytes(rd[24..32].try_into().unwrap());
    let n = n as usize;
    if bytes.len() < MAGIC.len() + 32 + 3 * n * 8 {
        return Err(Error::InvalidArgs);
    }
    let mut cur = MAGIC.len() + 32;
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
    })
}

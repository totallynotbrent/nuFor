//! plain-text structured output for the 1d state: a csv table and a vtk file.
//!
//! the csv is a human/script-friendly table of the cells; the vtk file is a
//! legacy-format structured-points block (ascii) that paraview and gnuplot
//! can read for plotting snapshots.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::eos::eos_pressure;
use crate::Error;

/// the conserved state plus cell centers and the uniform cell width.
pub struct OutputState<'a> {
    /// cell-center coordinates.
    pub centers: &'a [f64],
    /// density per cell.
    pub rho: &'a [f64],
    /// momentum per cell.
    pub m: &'a [f64],
    /// total energy per cell.
    pub e: &'a [f64],
    /// ratio of specific heats.
    pub gamma: f64,
}

fn checked(state: &OutputState) -> Result<usize, Error> {
    let n = state.centers.len();
    if n < 2 || state.rho.len() != n || state.m.len() != n || state.e.len() != n {
        return Err(Error::InvalidArgs);
    }
    Ok(n)
}

/// writes a csv table with one row per cell: x, rho, m, e, u, p.
pub fn write_csv(path: &Path, state: &OutputState) -> Result<(), Error> {
    let n = checked(state)?;
    let u: Vec<f64> = state
        .m
        .iter()
        .zip(state.rho)
        .map(|(&m, &r)| m / r)
        .collect();
    let et: Vec<f64> = state
        .e
        .iter()
        .zip(state.rho)
        .map(|(&e, &r)| e / r)
        .collect();
    let p = eos_pressure(state.gamma, state.rho, &et, &u)?;
    let mut out = File::create(path).map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "x,rho,m,e,u,p").map_err(|_| Error::InvalidArgs)?;
    for i in 0..n {
        writeln!(
            out,
            "{},{},{},{},{},{}",
            state.centers[i], state.rho[i], state.m[i], state.e[i], u[i], p[i]
        )
        .map_err(|_| Error::InvalidArgs)?;
    }
    Ok(())
}

/// writes a legacy vtk structured-points file with scalar rho, m, and e arrays.
pub fn write_vtk(path: &Path, state: &OutputState) -> Result<(), Error> {
    let n = checked(state)?;
    let dx = if state.centers.len() > 1 {
        state.centers[1] - state.centers[0]
    } else {
        1.0
    };
    let origin = state.centers[0] - 0.5 * dx;
    let mut out = File::create(path).map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "# vtk DataFile Version 3.0").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "nufor 1d snapshot").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "ASCII").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "DATASET STRUCTURED_POINTS").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "DIMENSIONS {n} 1 1").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "ORIGIN {origin} 0 0").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "SPACING {dx} 0 0").map_err(|_| Error::InvalidArgs)?;
    writeln!(out, "POINT_DATA {n}").map_err(|_| Error::InvalidArgs)?;
    for (name, arr) in [
        ("rho", state.rho),
        ("momentum", state.m),
        ("energy", state.e),
    ] {
        writeln!(out, "SCALARS {name} double 1").map_err(|_| Error::InvalidArgs)?;
        writeln!(out, "LOOKUP_TABLE default").map_err(|_| Error::InvalidArgs)?;
        for v in arr {
            writeln!(out, "{v}").map_err(|_| Error::InvalidArgs)?;
        }
    }
    Ok(())
}

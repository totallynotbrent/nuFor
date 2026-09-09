//! plain-text structured output for the 1d state: a csv table and a vtk file.
//!
//! the csv is a human/script-friendly table of the cells; the vtk file is a
//! legacy-format structured-points block (ascii) that paraview and gnuplot
//! can read for plotting snapshots.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::eos::eos_pressure;
use crate::grid2d::Grid2d;
use crate::grid3d::Grid3d;
use crate::state2d::ConservedState2d;
use crate::state3d::ConservedState3d;
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

/// a derived pressure field. p = (gamma-1)(e - kinetic), with kinetic from the
/// momentum components, so we need no primitive conversion up front.
fn pressure2d(rho: &[f64], mx: &[f64], my: &[f64], e: &[f64], gamma: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let kin = 0.5 * (mx[k] * mx[k] + my[k] * my[k]) / rho[k];
            (gamma - 1.0) * (e[k] - kin)
        })
        .collect()
}

fn pressure3d(
    rho: &[f64],
    mx: &[f64],
    my: &[f64],
    mz: &[f64],
    e: &[f64],
    gamma: f64,
    n: usize,
) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let kin = 0.5 * (mx[k] * mx[k] + my[k] * my[k] + mz[k] * mz[k]) / rho[k];
            (gamma - 1.0) * (e[k] - kin)
        })
        .collect()
}

/// write a 2d structured-points vtk block (ascii) from a resolved 2d state.
pub fn write_vtk2d(
    path: &Path,
    g: &Grid2d,
    st: &ConservedState2d,
    gamma: f64,
) -> Result<(), Error> {
    let n = g.nx * g.ny;
    let p = pressure2d(&st.rho, &st.mx, &st.my, &st.e, gamma, n);
    let mut s = String::with_capacity(n * 48);
    s.push_str("# vtk DataFile Version 3.0\nnuFor 2d snapshot\nASCII\nDATASET STRUCTURED_POINTS\n");
    s.push_str(&format!(
        "DIMENSIONS {} {} 1\nORIGIN {} {} 0\nSPACING {} {} 1\n",
        g.nx, g.ny, g.xmin, g.ymin, g.dx, g.dy
    ));
    s.push_str(&format!("POINT_DATA {}\n", n));
    append_scalar(&mut s, "density", &st.rho);
    append_scalar(&mut s, "pressure", &p);
    append_scalar(&mut s, "energy", &st.e);
    write_all(path, &s)
}

/// write a 3d structured-points vtk block (ascii) from a resolved 3d state.
pub fn write_vtk3d(
    path: &Path,
    g: &Grid3d,
    st: &ConservedState3d,
    gamma: f64,
) -> Result<(), Error> {
    let n = g.nx * g.ny * g.nz;
    let p = pressure3d(&st.rho, &st.mx, &st.my, &st.mz, &st.e, gamma, n);
    let mut s = String::with_capacity(n * 48);
    s.push_str("# vtk DataFile Version 3.0\nnuFor 3d snapshot\nASCII\nDATASET STRUCTURED_POINTS\n");
    s.push_str(&format!(
        "DIMENSIONS {} {} {}\nORIGIN {} {} {}\nSPACING {} {} {}\n",
        g.nx, g.ny, g.nz, g.xmin, g.ymin, g.zmin, g.dx, g.dy, g.dz
    ));
    s.push_str(&format!("POINT_DATA {}\n", n));
    append_scalar(&mut s, "density", &st.rho);
    append_scalar(&mut s, "pressure", &p);
    append_scalar(&mut s, "energy", &st.e);
    write_all(path, &s)
}

fn append_scalar(s: &mut String, name: &str, vals: &[f64]) {
    s.push_str(&format!("SCALARS {name} double 1\nLOOKUP_TABLE default\n"));
    let mut k = 0;
    for v in vals {
        s.push_str(&format!("{v:.6e} "));
        k += 1;
        if k % 6 == 0 {
            s.push('\n');
        }
    }
    s.push('\n');
}

fn write_all(path: &Path, s: &str) -> Result<(), Error> {
    let mut f = match File::create(path) {
        Ok(f) => f,
        Err(_) => return Err(Error::InvalidArgs),
    };
    if f.write_all(s.as_bytes()).is_err() {
        return Err(Error::InvalidArgs);
    }
    Ok(())
}

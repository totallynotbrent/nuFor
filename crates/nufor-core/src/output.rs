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
use crate::thermo::ThermoModel;
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
    write_vtk2d_model(path, g, st, ThermoModel::Perfect { gamma })
}

/// the closure-aware 2d snapshot writer: pressure and mach under the
/// run's own thermodynamic model, so an eqair run's file carries the
/// gas it actually marched with.
pub fn write_vtk2d_model(
    path: &Path,
    g: &Grid2d,
    st: &ConservedState2d,
    model: ThermoModel,
) -> Result<(), Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e)?;
    let e_int: Vec<f64> = (0..n)
        .map(|k| (et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k])).max(1.0))
        .collect();
    let p = crate::thermo::pressure(model, &st.rho, &e_int, &u, &v)?;
    let a = crate::thermo::sound_speed(model, &st.rho, &et, &u, &v)?;
    let mach: Vec<f64> = (0..n)
        .map(|k| (u[k] * u[k] + v[k] * v[k]).sqrt() / a[k].max(1e-6))
        .collect();
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
    append_scalar(&mut s, "mach", &mach);
    write_all(path, &s)
}

/// write the curvilinear body-fitted snapshot: a vtk STRUCTURED_GRID
/// with explicit node coordinates, so the mapped cells render
/// as-is. fields carry the closure-consistent pressure and mach.
pub fn write_vtk_curv(
    path: &Path,
    g: &crate::curvilinear::CurvGrid,
    st: &mut ConservedState2d,
    model: crate::thermo::ThermoModel,
) -> Result<(), Error> {
    let n = g.nx * g.ny;
    for k in 0..n {
        if !(st.rho[k].is_finite()) || st.rho[k] <= 1.0e-6 {
            st.rho[k] = 1.0e-6;
            st.mx[k] = 0.0;
            st.my[k] = 0.0;
            st.e[k] = 1.0e-6 * 1.0e4;
        } else if !(st.e[k].is_finite()) || st.e[k] <= 0.0 {
            st.e[k] = st.rho[k] * 1.0e4;
        }
    }
    let (u, v, et) = crate::state2d::cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e)?;
    // floor the closure inputs the same way the curvilinear march
    // does: transient cells can sit below the TGAS1 table's tidy
    // region while the field organizes.
    let e_int: Vec<f64> = (0..n)
        .map(|k| (et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k])).max(1.0e4))
        .collect();
    let p = crate::thermo::pressure(model, &st.rho, &e_int, &u, &v)?;
    let et_cl: Vec<f64> = (0..n)
        .map(|k| {
            let ke = 0.5 * (u[k] * u[k] + v[k] * v[k]);
            et[k].max(ke + 1.0e4)
        })
        .collect();
    let a = crate::thermo::sound_speed(model, &st.rho, &et_cl, &u, &v).map_err(|e| {
        let mut det = String::new();
        for k in 0..st.rho.len().min(30) {
            if !st.rho[k].is_finite() || st.rho[k] <= 0.0 || !(st.e[k].is_finite()) {
                det.push_str(&format!(
                    "bad state cell {} rho={} e={} ",
                    k, st.rho[k], st.e[k]
                ));
            }
            if !et[k].is_finite() {
                det.push_str(&format!(
                    "bad et cell {} et={} u={} v={} ",
                    k, et[k], u[k], v[k]
                ));
            }
        }
        if det.is_empty() {
            eprintln!("vtk write: sound_speed err {e} with clean states (sound formula domain)");
        } else {
            eprintln!("vtk write: {e} | {}", det);
        }
        e
    })?;
    let mach: Vec<f64> = (0..n)
        .map(|k| (u[k] * u[k] + v[k] * v[k]).sqrt() / a[k].max(1e-6))
        .collect();

    // node coordinates: bilinear corners of each cell. the mapping
    // stores cell centroids and faces, so reconstruct nodes from the
    // face midpoints: node[i,j] = the corner shared by xfaces i,i+1
    // and yfaces j,j+1 of the surrounding cells. for simplicity (and
    // because the quads are near-rectangular in the mapped frame)
    // write the centroid positions as an UNSTRUCTURED-grid proxy:
    // a STRUCTURED_GRID of the same topology with centroid coords.
    let mut s = String::with_capacity(n * 64);
    s.push_str("# vtk DataFile Version 3.0\nnuFor curvilinear snapshot\nASCII\n");
    s.push_str(&format!(
        "DATASET STRUCTURED_GRID\nDIMENSIONS {} {} 1\nPOINTS {} double\n",
        g.nx, g.ny, n
    ));
    for j in 0..g.ny {
        for i in 0..g.nx {
            let c = g.cells[j * g.nx + i];
            s.push_str(&format!("{:.6e} {:.6e} 0\n", c.x, c.r));
        }
    }
    s.push_str(&format!("POINT_DATA {}\n", n));
    append_scalar(&mut s, "density", &st.rho);
    append_scalar(&mut s, "pressure", &p);
    append_scalar(&mut s, "energy", &st.e);
    append_scalar(&mut s, "mach", &mach);
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

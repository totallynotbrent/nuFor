//! rectilinear (possibly non-uniform, axis-aligned) structured grids.
//!
//! a rectilinear grid is specified by per-axis coordinate arrays holding the
//! `n+1` face positions (or, equivalently, the `n` cell widths). this is the
//! natural shape of a VTK `RECTILINEAR_GRID`/`STRUCTURED_POINTS` file and of
//! the plain "coordinates" files the case format reads. it generalizes the
//! uniform `grid2d`/`grid3d` constructors by replacing the single `dx`/`dy`/
//! `dz` with per-axis spacing.

use crate::grid2d::Grid2d;
use crate::grid3d::{Bounds3d, Grid3d};
use crate::Error;

/// build a 2d rectilinear grid from per-axis face coordinates.
///
/// `fx` has `nx+1` entries (x faces, left to right) and `fy` has `ny+1`
/// (y faces, bottom to top). both must be strictly increasing and finite.
pub fn rectilinear_grid2d(fx: &[f64], fy: &[f64]) -> Result<Grid2d, Error> {
    if fx.len() < 2 || fy.len() < 2 {
        return Err(Error::InvalidArgs);
    }
    if !strictly_increasing(fx) || !strictly_increasing(fy) {
        return Err(Error::InvalidArgs);
    }
    let (nx, ny) = (fx.len() - 1, fy.len() - 1);
    let (xmin, xmax) = (fx[0], fx[nx]);
    let (ymin, ymax) = (fy[0], fy[ny]);
    let mut centers_x = Vec::with_capacity(nx * ny);
    let mut centers_y = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            centers_x.push(0.5 * (fx[i] + fx[i + 1]));
            centers_y.push(0.5 * (fy[j] + fy[j + 1]));
        }
    }
    let dx = widths(fx);
    let dy = widths(fy);
    Ok(Grid2d {
        nx,
        ny,
        xmin,
        xmax,
        ymin,
        ymax,
        dx: minimum(&dx),
        dy: minimum(&dy),
        dxs: dx,
        dys: dy,
        centers_x,
        centers_y,
        faces_x: fx.to_vec(),
        faces_y: fy.to_vec(),
    })
}

/// build a 3d rectilinear grid from per-axis face coordinates.
pub fn rectilinear_grid3d(fx: &[f64], fy: &[f64], fz: &[f64]) -> Result<Grid3d, Error> {
    if fx.len() < 2 || fy.len() < 2 || fz.len() < 2 {
        return Err(Error::InvalidArgs);
    }
    if !strictly_increasing(fx) || !strictly_increasing(fy) || !strictly_increasing(fz) {
        return Err(Error::InvalidArgs);
    }
    let (nx, ny, nz) = (fx.len() - 1, fy.len() - 1, fz.len() - 1);
    let b = Bounds3d {
        xmin: fx[0],
        xmax: fx[nx],
        ymin: fy[0],
        ymax: fy[ny],
        zmin: fz[0],
        zmax: fz[nz],
    };
    let mut centers_x = Vec::with_capacity(nx * ny * nz);
    let mut centers_y = Vec::with_capacity(nx * ny * nz);
    let mut centers_z = Vec::with_capacity(nx * ny * nz);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                centers_x.push(0.5 * (fx[i] + fx[i + 1]));
                centers_y.push(0.5 * (fy[j] + fy[j + 1]));
                centers_z.push(0.5 * (fz[k] + fz[k + 1]));
            }
        }
    }
    let (dx, dy, dz) = (widths(fx), widths(fy), widths(fz));
    Ok(Grid3d {
        nx,
        ny,
        nz,
        xmin: b.xmin,
        xmax: b.xmax,
        ymin: b.ymin,
        ymax: b.ymax,
        zmin: b.zmin,
        zmax: b.zmax,
        dx: minimum(&dx),
        dy: minimum(&dy),
        dz: minimum(&dz),
        dxs: dx,
        dys: dy,
        dzs: dz,
        centers_x,
        centers_y,
        centers_z,
        faces_x: fx.to_vec(),
        faces_y: fy.to_vec(),
        faces_z: fz.to_vec(),
    })
}

/// the per-cell width along one axis, `n = cells`, inputs are `n+1` faces.
fn widths(faces: &[f64]) -> Vec<f64> {
    faces.windows(2).map(|w| w[1] - w[0]).collect()
}

fn minimum(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::INFINITY, f64::min)
}

fn strictly_increasing(v: &[f64]) -> bool {
    v.iter().all(|x| x.is_finite()) && v.windows(2).all(|w| w[1] > w[0])
}

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

/// face coordinates for one axis with geometric clustering: cells grow from
/// `h0` by ratio `r` while they fit, and the remainder fills with equal
/// width cells at least as wide as the last geometric one, so the spacing
/// never jumps downward. `r <= 1` recovers a uniform grid.
fn clustered_faces(a0: f64, a1: f64, n: usize, h0: f64, r: f64) -> Result<Vec<f64>, Error> {
    if n < 1 || h0 <= 0.0 || r <= 0.0 {
        return Err(Error::InvalidArgs);
    }
    let span = a1 - a0;
    if span <= 0.0 {
        return Err(Error::InvalidArgs);
    }
    if n == 1 || (r - 1.0).abs() < 1e-12 {
        let w = span / n as f64;
        return Ok((0..=n).map(|k| a0 + k as f64 * w).collect());
    }
    // cumulative width of m geometrically growing cells.
    let cum = |m: usize| h0 * (r.powi(m as i32) - 1.0) / (r - 1.0);
    let width = |k: usize| h0 * r.powi(k as i32);
    // the largest geometric prefix that fits and leaves tail cells no
    // narrower than the last geometric cell (monotone spacing).
    let mut m = 0usize;
    for cand in 1..n {
        if cum(cand) > span {
            break;
        }
        let tail = (span - cum(cand)) / (n - cand) as f64;
        if tail + 1e-12 >= width(cand - 1) {
            m = cand;
        }
    }
    if m == 0 {
        let w = span / n as f64;
        return Ok((0..=n).map(|k| a0 + k as f64 * w).collect());
    }
    let mut faces = Vec::with_capacity(n + 1);
    faces.push(a0);
    for k in 0..m {
        faces.push(faces[k] + width(k));
    }
    let tail = (span - cum(m)) / (n - m) as f64;
    for k in 0..n - m {
        faces.push(faces[m + k] + tail);
    }
    faces.truncate(n + 1);
    faces[n] = a1;
    Ok(faces)
}

/// how one axis is clustered toward a wall: the first-cell height and the
/// geometric growth ratio of successive cells.
#[derive(Debug, Clone, Copy)]
pub struct Clustering {
    /// first-cell height at the refined side.
    pub first_cell: f64,
    /// growth ratio of each next cell.
    pub growth: f64,
}

/// a wall-normal clustered 2d grid for boundary-layer work: uniform in x,
/// y clustered toward the bottom wall (the plate).
pub fn stretched_grid2d(
    nx: usize,
    x0: f64,
    x1: f64,
    ny: usize,
    y0: f64,
    y1: f64,
    c: Clustering,
) -> Result<Grid2d, Error> {
    let fx: Vec<f64> = (0..=nx)
        .map(|i| x0 + (x1 - x0) * i as f64 / nx as f64)
        .collect();
    let fy = clustered_faces(y0, y1, ny, c.first_cell, c.growth)?;
    rectilinear_grid2d(&fx, &fy)
}

/// a channel grid clustered toward both walls: the half-height is clustered
/// bottom-up and mirrored, so the first cell at each wall is `first_cell`
/// with growth `growth`. `ny` must be even so the channel halves split cleanly.
pub fn channel_grid2d(
    nx: usize,
    x0: f64,
    x1: f64,
    ny: usize,
    y0: f64,
    y1: f64,
    c: Clustering,
) -> Result<Grid2d, Error> {
    if ny < 2 || ny % 2 != 0 {
        return Err(Error::InvalidArgs);
    }
    let fx: Vec<f64> = (0..=nx)
        .map(|i| x0 + (x1 - x0) * i as f64 / nx as f64)
        .collect();
    let half = clustered_faces(y0, 0.5 * (y0 + y1), ny / 2, c.first_cell, c.growth)?;
    let mut fy = half.clone();
    // mirror the half-cluster about the channel centerline.
    let ymid = 0.5 * (y0 + y1);
    for &y in half.iter().rev().skip(1) {
        fy.push(2.0 * ymid - y);
    }
    fy.truncate(ny + 1);
    fy[ny] = y1;
    rectilinear_grid2d(&fx, &fy)
}

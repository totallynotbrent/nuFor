//! structured 2d cartesian grid: cell centers, cell widths, and axis-aligned faces.

use crate::Error;

/// a uniform structured grid over [xmin,xmax] x [ymin,ymax] with nx x ny cells.
///
/// arrays are row-major (j*ny? no: j*nx+i), stored flat with nx*ny entries;
/// faces are the nx+1 vertical lines and ny+1 horizontal lines that close the
/// domain.
#[derive(Debug, Clone)]
pub struct Grid2d {
    /// number of cells in the x direction.
    pub nx: usize,
    /// number of cells in the y direction.
    pub ny: usize,
    /// lower x bound of the domain.
    pub xmin: f64,
    /// upper x bound of the domain.
    pub xmax: f64,
    /// lower y bound of the domain.
    pub ymin: f64,
    /// upper y bound of the domain.
    pub ymax: f64,
    /// uniform cell width in x (minimum over the axis when non-uniform).
    pub dx: f64,
    /// uniform cell width in y (minimum over the axis when non-uniform).
    pub dy: f64,
    /// per-cell x widths (nx entries); equals [dx; nx] on a uniform grid.
    pub dxs: Vec<f64>,
    /// per-cell y widths (ny entries); equals [dy; ny] on a uniform grid.
    pub dys: Vec<f64>,
    /// x coordinate of each cell center, flat row-major (j*nx + i).
    pub centers_x: Vec<f64>,
    /// y coordinate of each cell center, flat row-major (j*nx + i).
    pub centers_y: Vec<f64>,
    /// x coordinates of the nx+1 vertical faces.
    pub faces_x: Vec<f64>,
    /// y coordinates of the ny+1 horizontal faces.
    pub faces_y: Vec<f64>,
}

/// build a structured cartesian grid over the given domain.
pub fn grid2d(
    nx: usize,
    ny: usize,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
) -> Result<Grid2d, Error> {
    if nx < 1 || ny < 1 {
        return Err(Error::InvalidArgs);
    }
    let finite = [xmin, xmax, ymin, ymax].iter().all(|v| v.is_finite());
    if !finite || xmax <= xmin || ymax <= ymin {
        return Err(Error::InvalidArgs);
    }
    let dx = (xmax - xmin) / nx as f64;
    let dy = (ymax - ymin) / ny as f64;
    let mut centers_x = Vec::with_capacity(nx * ny);
    let mut centers_y = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            centers_x.push(xmin + (i as f64 + 0.5) * dx);
            centers_y.push(ymin + (j as f64 + 0.5) * dy);
        }
    }
    let faces_x = (0..=nx).map(|i| xmin + i as f64 * dx).collect();
    let faces_y = (0..=ny).map(|j| ymin + j as f64 * dy).collect();
    Ok(Grid2d {
        nx,
        ny,
        xmin,
        xmax,
        ymin,
        ymax,
        dx,
        dy,
        dxs: vec![dx; nx],
        dys: vec![dy; ny],
        centers_x,
        centers_y,
        faces_x,
        faces_y,
    })
}

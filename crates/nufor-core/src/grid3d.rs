//! structured 3d cartesian grid: cell centers, cell widths, and axis faces.

use crate::Error;

/// a uniform structured grid over [xmin,xmax]x[ymin,ymax]x[zmin,zmax].
///
/// arrays are flat with nx*ny*nz entries, index (k*ny + j)*nx + i so x varies
/// fastest; faces are the nx+1, ny+1, nz+1 planes that close the domain.
#[derive(Debug, Clone)]
pub struct Grid3d {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub zmin: f64,
    pub zmax: f64,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub centers_x: Vec<f64>,
    pub centers_y: Vec<f64>,
    pub centers_z: Vec<f64>,
    pub faces_x: Vec<f64>,
    pub faces_y: Vec<f64>,
    pub faces_z: Vec<f64>,
}

/// the six domain bounds of a 3d box, grouped so the grid constructor stays tidy.
#[derive(Debug, Clone, Copy)]
pub struct Bounds3d {
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub zmin: f64,
    pub zmax: f64,
}

/// build a structured cartesian grid over the given domain.
pub fn grid3d(nx: usize, ny: usize, nz: usize, b: &Bounds3d) -> Result<Grid3d, Error> {
    let (xmin, xmax, ymin, ymax, zmin, zmax) = (b.xmin, b.xmax, b.ymin, b.ymax, b.zmin, b.zmax);
    if nx < 1 || ny < 1 || nz < 1 {
        return Err(Error::InvalidArgs);
    }
    let finite = [xmin, xmax, ymin, ymax, zmin, zmax]
        .iter()
        .all(|v| v.is_finite());
    if !finite || xmax <= xmin || ymax <= ymin || zmax <= zmin {
        return Err(Error::InvalidArgs);
    }
    let (dx, dy, dz) = (
        (xmax - xmin) / nx as f64,
        (ymax - ymin) / ny as f64,
        (zmax - zmin) / nz as f64,
    );
    let cells = nx * ny * nz;
    let (mut cx, mut cy, mut cz) = (
        Vec::with_capacity(cells),
        Vec::with_capacity(cells),
        Vec::with_capacity(cells),
    );
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                cx.push(xmin + (i as f64 + 0.5) * dx);
                cy.push(ymin + (j as f64 + 0.5) * dy);
                cz.push(zmin + (k as f64 + 0.5) * dz);
            }
        }
    }
    Ok(Grid3d {
        nx,
        ny,
        nz,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
        dx,
        dy,
        dz,
        centers_x: cx,
        centers_y: cy,
        centers_z: cz,
        faces_x: (0..=nx).map(|i| xmin + i as f64 * dx).collect(),
        faces_y: (0..=ny).map(|j| ymin + j as f64 * dy).collect(),
        faces_z: (0..=nz).map(|k| zmin + k as f64 * dz).collect(),
    })
}

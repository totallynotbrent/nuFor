//! wall distance for the structured 2d grid: the perpendicular distance from
//! each cell center to the nearest solid-wall boundary plane, the standard
//! cheap field for near-wall turbulence models on cartesian meshes.

use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;

/// distance from each cell center to the nearest wall among the solid sides.
///
/// for the cartesian grid a wall side is a plane, so the perpendicular
/// distance is exact; a domain with no solid sides has no wall, and the
/// caller supplies the freestream distance instead (the sa model needs a
/// finite d; infinity would zero the destruction term through (nu_tilde/d)^2
/// only in the limit, so the caller passes something like the domain size).
pub fn wall_distance2d(g: &Grid2d, bc: &Boundaries2d, far: f64) -> Vec<f64> {
    let n = g.nx * g.ny;
    let is_wall = |b: crate::solver2d::Bc2d| {
        matches!(
            b,
            crate::solver2d::Bc2d::SlipWall | crate::solver2d::Bc2d::NoSlipWall
        )
    };
    let mut out = vec![far; n];
    // distances to each wall plane; only solid sides participate.
    let west = is_wall(bc.west).then_some(g.xmin);
    let east = is_wall(bc.east).then_some(g.xmax);
    let south = is_wall(bc.south).then_some(g.ymin);
    let north = is_wall(bc.north).then_some(g.ymax);
    if west.is_none() && east.is_none() && south.is_none() && north.is_none() {
        return out;
    }
    for j in 0..g.ny {
        for i in 0..g.nx {
            let c = j * g.nx + i;
            let (x, y) = (g.centers_x[c], g.centers_y[c]);
            let mut best = far;
            if let Some(xw) = west {
                best = best.min(x - xw);
            }
            if let Some(xe) = east {
                best = best.min(xe - x);
            }
            if let Some(ys) = south {
                best = best.min(y - ys);
            }
            if let Some(yn) = north {
                best = best.min(yn - y);
            }
            out[c] = best;
        }
    }
    out
}

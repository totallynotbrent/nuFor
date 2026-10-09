//! cut-cell geometry fractions from a signed distance field: the
//! first half of the cut-cell upgrade (days 2-3 of the sprint).
//!
//! a cut cell straddling a solid boundary carries a fluid volume
//! fraction alpha_vol and open face fractions alpha_face[4]; the
//! conservative update on the partial cell scales each face flux by
//! its open fraction and divides by the volume fraction. the fractions
//! come from sampling the SDF: the mean-value theorem route, where the
//! fraction of a segment lying in the fluid equals the fraction of its
//! sample points with positive SDF (the midpoint rule over N samples
//! per face). accuracy is one order below the sample count supports,
//! which for a first-order cut treatment means 2 samples per face and
//! 2x2 per cell are enough; the defaults here use 4 for margin.

/// the fluid fraction of a 1d segment [0, len] sampled at n midpoints,
/// SDF positive inside the fluid.
fn segment_fraction(sdf: &dyn Fn(f64) -> f64, len: f64, n: usize) -> f64 {
    if len <= 0.0 {
        return 0.0;
    }
    let mut open = 0usize;
    for k in 0..n {
        let t = (k as f64 + 0.5) * len / n as f64;
        if sdf(t) > 0.0 {
            open += 1;
        }
    }
    open as f64 / n as f64
}

/// the fluid fraction of one cell, by the product of segment
/// quadratures along each axis (a 2x2 midpoint tensor product).
fn cell_fraction(sdf: &dyn Fn(f64, f64) -> f64, x0: f64, x1: f64, y0: f64, y1: f64) -> f64 {
    let (nx, ny) = (2usize, 2usize);
    let mut open = 0usize;
    for j in 0..ny {
        let y = y0 + (j as f64 + 0.5) * (y1 - y0) / ny as f64;
        for i in 0..nx {
            let x = x0 + (i as f64 + 0.5) * (x1 - x0) / nx as f64;
            if sdf(x, y) > 0.0 {
                open += 1;
            }
        }
    }
    open as f64 / (nx * ny) as f64
}

/// the cut-cell geometry of one cell: fluid volume fraction and open
/// fractions of the west, east, south, north faces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutFractions {
    pub vol: f64,
    pub west: f64,
    pub east: f64,
    pub south: f64,
    pub north: f64,
}

impl CutFractions {
    /// solid cell (all zero) or full cell (all one) as the degenerate
    /// cases; the fraction layout matches the staircase mask exactly
    /// when vol is 0 or 1.
    pub fn full() -> Self {
        CutFractions {
            vol: 1.0,
            west: 1.0,
            east: 1.0,
            south: 1.0,
            north: 1.0,
        }
    }

    pub fn solid() -> Self {
        CutFractions {
            vol: 0.0,
            west: 0.0,
            east: 0.0,
            south: 0.0,
            north: 0.0,
        }
    }
}

/// compute the fractions for the cell with corners (x0,y0)-(x1,y1).
pub fn cut_fractions(
    sdf: &dyn Fn(f64, f64) -> f64,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
) -> CutFractions {
    let west = segment_fraction(&|t| sdf(x0, y0 + t), y1 - y0, 4);
    let east = segment_fraction(&|t| sdf(x1, y0 + t), y1 - y0, 4);
    let south = segment_fraction(&|t| sdf(x0 + t, y0), x1 - x0, 4);
    let north = segment_fraction(&|t| sdf(x0 + t, y1), x1 - x0, 4);
    let vol = cell_fraction(sdf, x0, x1, y0, y1);
    // a cell with zero fluid volume passes nothing through its faces:
    // a vol-0 sliver with open apertures would divide the update by
    // ~0. the geometry and the update must agree that it is solid.
    if vol <= 0.0 {
        return CutFractions::solid();
    }
    CutFractions {
        vol,
        west,
        east,
        south,
        north,
    }
}

/// the corner-interpolation fractions: sample the SDF at the cell
/// corners, place edge crossings by linear interpolation, and read the
/// volume and face fractions off the interface polygon. exact for a
/// planar interface, O(h^2) interface position for smooth walls, the
/// recipe the method survey recommends over point sampling.
#[allow(clippy::too_many_arguments)]
pub fn cut_fractions_corner(
    sdf: &dyn Fn(f64, f64) -> f64,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
) -> CutFractions {
    let bl = sdf(x0, y0);
    let br = sdf(x1, y0);
    let tr = sdf(x1, y1);
    let tl = sdf(x0, y1);
    if bl <= 0.0 && br <= 0.0 && tr <= 0.0 && tl <= 0.0 {
        return CutFractions::solid();
    }
    if bl > 0.0 && br > 0.0 && tr > 0.0 && tl > 0.0 {
        return CutFractions::full();
    }
    // the fluid fraction of an edge with endpoint samples (a at the
    // start, b at the end): both same sign is all-or-nothing; on a
    // mixed edge the crossing sits at t = a/(a-b) from the start, so
    // the fluid span is t when the start is fluid and 1-t when it is
    // solid.
    let edge_frac = |a: f64, b: f64| -> f64 {
        if a * b > 0.0 {
            if a > 0.0 {
                1.0
            } else {
                0.0
            }
        } else if (a - b).abs() < 1e-30 {
            0.5
        } else {
            let t = (a / (a - b)).clamp(0.0, 1.0);
            if a > 0.0 {
                t
            } else {
                1.0 - t
            }
        }
    };
    let south = edge_frac(bl, br);
    let north = edge_frac(tl, tr);
    let west = edge_frac(bl, tl);
    let east = edge_frac(br, tr);

    // the fluid polygon: corners with phi > 0 plus the edge crossings,
    // in order; shoelace gives the area fraction.
    let (dx, dy) = (x1 - x0, y1 - y0);
    let corners = [(x0, y0, bl), (x1, y0, br), (x1, y1, tr), (x0, y1, tl)];
    let mut poly: Vec<(f64, f64)> = Vec::with_capacity(6);
    for w in 0..4 {
        let (cx, cy, cv) = corners[w];
        let (nx, ny, nv) = corners[(w + 1) % 4];
        if cv > 0.0 {
            poly.push((cx, cy));
        }
        if cv * nv < 0.0 {
            let t = cv / (cv - nv);
            poly.push((cx + t * (nx - cx), cy + t * (ny - cy)));
        }
    }
    let mut area2 = 0.0f64;
    for w in 0..poly.len() {
        let (ax, ay) = poly[w];
        let (bx, by) = poly[(w + 1) % poly.len()];
        area2 += ax * by - bx * ay;
    }
    let vol = (0.5 * area2 / (dx * dy)).abs();
    if vol <= 1e-12 {
        return CutFractions::solid();
    }
    if west <= 0.0 && east <= 0.0 && south <= 0.0 && north <= 0.0 {
        return CutFractions::solid();
    }
    CutFractions {
        vol,
        west,
        east,
        south,
        north,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_plane_cuts_exact() {
        // a body occupying x > 0.5: the SDF is 0.5 - x.
        let sdf = |x: f64, _y: f64| 0.5 - x;
        // a cell spanning [0.4, 0.6] x [0, 1]: half fluid.
        let f = cut_fractions(&sdf, 0.4, 0.6, 0.0, 1.0);
        assert!((f.vol - 0.5).abs() < 1e-12, "vol {vol}", vol = f.vol);
        assert!((f.north - 0.5).abs() < 1e-12);
        assert!((f.south - 0.5).abs() < 1e-12);
        assert_eq!(f.west, 1.0);
        assert_eq!(f.east, 0.0);
    }

    #[test]
    fn corner_builder_half_plane_exact() {
        let sdf = |x: f64, _y: f64| 0.5 - x;
        let f = cut_fractions_corner(&sdf, 0.4, 0.6, 0.0, 1.0);
        assert!((f.vol - 0.5).abs() < 1e-12, "vol {vol}", vol = f.vol);
        assert!((f.north - 0.5).abs() < 1e-12);
        assert!((f.south - 0.5).abs() < 1e-12);
        assert_eq!(f.west, 1.0);
        assert_eq!(f.east, 0.0);
    }

    #[test]
    fn corner_builder_diagonal_exact() {
        // a 45-degree interface crossing the left edge at y = 0.2 and
        // the top edge at x = 0.8 (the exact diagonal through corners
        // is degenerate): the fluid corner triangle has area
        // 0.5 * 0.8 * 0.8 = 0.32, and the open apertures are the left
        // and top spans.
        let sdf = |x: f64, y: f64| y - x - 0.2;
        let f = cut_fractions_corner(&sdf, 0.0, 1.0, 0.0, 1.0);
        assert!((f.vol - 0.32).abs() < 1e-12, "vol {vol}", vol = f.vol);
        assert_eq!(f.south, 0.0);
        assert_eq!(f.east, 0.0);
        assert!((f.west - 0.8).abs() < 1e-12);
        assert!((f.north - 0.8).abs() < 1e-12);
    }

    #[test]
    fn corner_builder_thin_sliver_is_consistent() {
        // a sliver so thin the interface passes between the face
        // midpoints and the far corners: the volume never reaches 0
        // with open apertures (the midpoint builder's failure mode).
        // a plane at distance 0.001 inside a unit cell from the west
        // face: vol ~ 0.001, west aperture 1.
        let sdf = |x: f64, _y: f64| 0.001 - x;
        let f = cut_fractions_corner(&sdf, 0.0, 1.0, 0.0, 1.0);
        assert!(f.vol > 0.0, "a thin sliver keeps a positive volume");
        assert_eq!(f.west, 1.0);
        assert_eq!(f.east, 0.0);
    }

    #[test]
    fn circle_fractions_sum_sensibly() {
        // a circle radius 1 at the origin, fluid outside.
        let sdf = |x: f64, y: f64| (x * x + y * y).sqrt() - 1.0;
        // the cell [0.6, 0.8]^2 straddles the circle's first-quadrant
        // arc (corners r = 0.85 inside, r = 1.13 outside). the volume
        // fraction must lie between the min and max face fractions.
        let f = cut_fractions(&sdf, 0.6, 0.8, 0.6, 0.8);
        let lo = f.west.min(f.east).min(f.south).min(f.north);
        let hi = f.west.max(f.east).max(f.south).max(f.north);
        assert!(f.vol >= lo - 1e-12 && f.vol <= hi + 1e-12, "{f:?}");
        assert!(f.vol > 0.0 && f.vol < 1.0);
    }
}

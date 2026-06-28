//! line probes of a 2d scalar field through bilinear interpolation.
//!
//! a probe samples a cell-centered scalar field along a straight segment and
//! returns the value at evenly spaced points, so a chart can show how a
//! quantity varies across a flow rather than only as a colormap.

use crate::{Error, Grid2d};

/// sample `data` (row-major, j*nx+i) along the segment (x0,y0)->(x1,y1).
///
/// returns `n` pairs of (distance along the segment, interpolated value),
/// evenly spaced from the start to the end. points outside the domain are
/// clamped to the boundary cell so the probe never reads past the grid.
pub fn probe_line(
    data: &[f64],
    g: &Grid2d,
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    n: usize,
) -> Result<Vec<(f64, f64)>, Error> {
    if data.len() != g.nx * g.ny {
        return Err(Error::InvalidArgs);
    }
    if n < 2 {
        return Err(Error::InvalidArgs);
    }
    let (len_x, len_y) = (x1 - x0, y1 - y0);
    let seg = (len_x * len_x + len_y * len_y).sqrt();
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let t = k as f64 / (n - 1) as f64;
        let (x, y) = (x0 + t * len_x, y0 + t * len_y);
        let v = sample(data, g, x, y);
        out.push((t * seg, v));
    }
    Ok(out)
}

/// bilinearly interpolate `data` at the point (x, y), clamping to the grid.
fn sample(data: &[f64], g: &Grid2d, x: f64, y: f64) -> f64 {
    // fractional cell index through the center lattice; clamp at the boundary.
    let fx = ((x - g.xmin) / g.dx - 0.5).clamp(0.0, (g.nx - 1) as f64);
    let fy = ((y - g.ymin) / g.dy - 0.5).clamp(0.0, (g.ny - 1) as f64);
    let i0 = fx.floor() as usize;
    let j0 = fy.floor() as usize;
    let i1 = (i0 + 1).min(g.nx - 1);
    let j1 = (j0 + 1).min(g.ny - 1);
    let (tx, ty) = (fx - i0 as f64, fy - j0 as f64);
    let (f00, f10) = (data[j0 * g.nx + i0], data[j0 * g.nx + i1]);
    let (f01, f11) = (data[j1 * g.nx + i0], data[j1 * g.nx + i1]);
    let top = f00 + tx * (f10 - f00);
    let bot = f01 + tx * (f11 - f01);
    top + ty * (bot - top)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid2d;

    #[test]
    fn linear_field_is_interpolated_exactly() {
        // f = x + 2y on the unit cell-center lattice.
        let g = grid2d(8, 8, 0.0, 1.0, 0.0, 1.0).unwrap();
        let data: Vec<f64> = g
            .centers_x
            .iter()
            .zip(&g.centers_y)
            .map(|(&x, &y)| x + 2.0 * y)
            .collect();
        let probe = probe_line(&data, &g, 0.2, 0.25, 0.8, 0.75, 11).unwrap();
        assert_eq!(probe.len(), 11);
        for &(s, v) in &probe {
            // segment vector is (0.6, 0.5); t = distance / |segment|.
            let t = s / 0.61f64.sqrt();
            let (x, y) = (0.2 + 0.6 * t, 0.25 + 0.5 * t);
            let want = x + 2.0 * y;
            assert!((v - want).abs() < 1e-9, "s={s} v={v} want={want}");
        }
    }

    #[test]
    fn horizontal_row_through_the_midplane() {
        let g = grid2d(16, 16, 0.0, 1.0, 0.0, 1.0).unwrap();
        let data: Vec<f64> = (0..g.nx * g.ny).map(|k| k as f64).collect();
        let probe = probe_line(&data, &g, 0.05, 0.5, 0.95, 0.5, 5).unwrap();
        // symmetric field, so the centre point should sit near the mean.
        let mid = probe[probe.len() / 2].1;
        let lo = probe.first().unwrap().1;
        let hi = probe.last().unwrap().1;
        assert!(mid > lo && mid < hi);
        // constant y: values rise with x because the field increases with k.
        assert!(probe[1].1 > probe[0].1);
    }

    #[test]
    fn out_of_domain_point_clamps_to_the_edge() {
        let g = grid2d(4, 4, 0.0, 1.0, 0.0, 1.0).unwrap();
        let data = vec![5.0; 16];
        let probe = probe_line(&data, &g, -0.5, -0.5, 1.5, 1.5, 3).unwrap();
        for &(_, v) in &probe {
            assert_eq!(v, 5.0);
        }
    }
}

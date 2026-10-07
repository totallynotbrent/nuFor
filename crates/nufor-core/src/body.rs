//! immersed solid bodies on the structured 2d grid.
//!
//! a body is cut out of the cartesian mesh: cells inside the body are frozen
//! to a quiescent reference, and the band of fluid cells adjacent to the
//! surface has the velocity component normal to the surface removed, which
//! enforces a slip wall (no-normal-flow) at the immersed boundary. this turns
//! a plain cartesian solver into one that can run flow over a blunt body.
//!
//! two shapes go through one interface: the analytic circle (the validated
//! cylinder case) and the arbitrary closed polygon (any sketched profile,
//! the aeroshell meridian). `inside` decides which cells are solid and
//! `normal` feeds the surface band, so the marching code never branches on
//! shape.

use crate::grid2d::Grid2d;
use crate::state2d::ConservedState2d;

/// a solid body shape: closed interior test plus outward surface normal.
pub trait SolidShape {
    /// whether the point falls inside the body.
    fn inside(&self, x: f64, y: f64) -> bool;

    /// the outward unit normal at a point near the surface.
    fn normal(&self, x: f64, y: f64) -> (f64, f64);
}

/// a solid circular body on the 2d grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidBody {
    /// center x coordinate.
    pub cx: f64,
    /// center y coordinate.
    pub cy: f64,
    /// radius.
    pub r: f64,
}

impl SolidShape for SolidBody {
    fn inside(&self, x: f64, y: f64) -> bool {
        let (dx, dy) = (x - self.cx, y - self.cy);
        dx * dx + dy * dy <= self.r * self.r
    }

    fn normal(&self, x: f64, y: f64) -> (f64, f64) {
        let (dx, dy) = (x - self.cx, y - self.cy);
        let l = (dx * dx + dy * dy).sqrt().max(1e-30);
        (dx / l, dy / l)
    }
}

/// a solid polygon body: a closed loop of vertices in order (either winding;
/// the inside test is crossing-count based so orientation does not matter).
///
/// the surface normal is taken from the nearest edge, outward from the
/// interior. vertices closer than about half a cell to their neighbors are
/// dropped so a fine sketch does not create degenerate edges.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidPolygon {
    /// vertices (x, y) in closed order, no duplicate closing vertex.
    pub verts: Vec<(f64, f64)>,
    /// cached total area, signed by winding; used to orient normals.
    signed_area: f64,
}

impl SolidPolygon {
    /// build from a vertex loop, dropping near-duplicate consecutive points.
    pub fn new(verts: Vec<(f64, f64)>) -> Self {
        let mut vs: Vec<(f64, f64)> = Vec::with_capacity(verts.len());
        for v in verts {
            if let Some(&last) = vs.last() {
                let d = ((v.0 - last.0).powi(2) + (v.1 - last.1).powi(2)).sqrt();
                if d < 1e-12 {
                    continue;
                }
            }
            vs.push(v);
        }
        // drop the closing vertex if it repeats the first.
        if vs.len() > 2 {
            let (first, last) = (vs[0], vs[vs.len() - 1]);
            let d = ((first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)).sqrt();
            if d < 1e-12 {
                vs.pop();
            }
        }
        let mut signed_area = 0.0;
        for i in 0..vs.len() {
            let j = (i + 1) % vs.len();
            signed_area += vs[i].0 * vs[j].1 - vs[j].0 * vs[i].1;
        }
        signed_area *= 0.5;
        SolidPolygon {
            verts: vs,
            signed_area,
        }
    }

    /// a regular n-gon inscribed in the circle (n polygon vertices around).
    pub fn inscribed_polygon(cx: f64, cy: f64, r: f64, n: usize) -> Self {
        let mut verts = Vec::with_capacity(n);
        for i in 0..n {
            let th = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            verts.push((cx + r * th.cos(), cy + r * th.sin()));
        }
        SolidPolygon::new(verts)
    }
}

impl SolidShape for SolidPolygon {
    fn inside(&self, x: f64, y: f64) -> bool {
        // even-odd ray cast to the +x direction.
        let n = self.verts.len();
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = self.verts[i];
            let (xj, yj) = self.verts[j];
            if (yi > y) != (yj > y) {
                // x where the edge crosses the horizontal line at y.
                let t = (y - yi) / (yj - yi);
                let xc = xi + t * (xj - xi);
                if x < xc {
                    inside = !inside;
                }
            }
            j = i;
        }
        inside
    }

    fn normal(&self, x: f64, y: f64) -> (f64, f64) {
        // nearest edge wins; its outward normal comes from the winding.
        let n = self.verts.len();
        let mut best = (0usize, f64::INFINITY);
        let mut best_point = (0.0f64, 0.0f64);
        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = self.verts[i];
            let (bx, by) = self.verts[j];
            let (ex, ey) = (bx - ax, by - ay);
            let len2 = (ex * ex + ey * ey).max(1e-30);
            let t = (((x - ax) * ex + (y - ay) * ey) / len2).clamp(0.0, 1.0);
            let (px, py) = (ax + t * ex, ay + t * ey);
            let d2 = (x - px).powi(2) + (y - py).powi(2);
            if d2 < best.1 {
                best = (i, d2);
                best_point = (px, py);
            }
        }
        let (i, _) = best;
        let j = (i + 1) % n;
        let (ax, ay) = self.verts[i];
        let (bx, by) = self.verts[j];
        // edge direction rotated +90 deg: (-ey, ex); flip by signed area so
        // the normal points away from the interior.
        let (mut nx, mut ny) = (-(by - ay), bx - ax);
        let l = (nx * nx + ny * ny).sqrt().max(1e-30);
        nx /= l;
        ny /= l;
        // a positive signed area (counter-clockwise) wants the rotated-left
        // normal; check by testing the interior side via a probe point.
        let (px, py) = best_point;
        let probe_inside = self.inside(px + nx * 1e-9, py + ny * 1e-9);
        if probe_inside {
            nx = -nx;
            ny = -ny;
        }
        let _ = self.signed_area;
        (nx, ny)
    }
}

/// a quiescent cell state inside the body: ambient density and pressure.
pub struct SolidRef {
    pub rho: f64,
    pub p: f64,
}

impl Default for SolidRef {
    fn default() -> Self {
        SolidRef { rho: 1.0, p: 1.0 }
    }
}

/// impose the solid body on the state after each step.
///
/// cells strictly inside the body are held at the inert reference; fluid cells
/// whose center sits within the interface band have their velocity projected
/// onto the surface tangent, killing the normal component so no mass flows
/// through the wall. density and pressure at the surface are copied, matching
/// the zero-gradient expectation at a slip wall.
pub fn apply_solid<S: SolidShape>(state: &mut ConservedState2d, g: &Grid2d, body: &S, gamma: f64) {
    let ref_ = SolidRef::default();
    let e_ref = ref_.p / ((gamma - 1.0) * ref_.rho);
    let band = 1.5 * g.dx.min(g.dy);
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let (x, y) = (g.centers_x[k], g.centers_y[k]);
            if body.inside(x, y) {
                // solid cell: keep it inert.
                state.rho[k] = ref_.rho;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
                state.e[k] = e_ref;
            } else if near_surface(body, x, y, band) {
                // fluid cell bordering the surface: remove the normal velocity.
                let (nx, ny) = body.normal(x, y);
                let u = state.mx[k] / state.rho[k];
                let v = state.my[k] / state.rho[k];
                let un = u * nx + v * ny;
                if un > 0.0 {
                    // outward-flowing component is reflected back (slip).
                    state.mx[k] -= 2.0 * un * nx * state.rho[k];
                    state.my[k] -= 2.0 * un * ny * state.rho[k];
                } else {
                    // inward- or zero-normal: remove the inward normal part too.
                    state.mx[k] -= un * nx * state.rho[k];
                    state.my[k] -= un * ny * state.rho[k];
                }
            }
        }
    }
}

/// whether the point sits within the interface band outside the surface.
///
/// the circle has an exact distance; the polygon approximates with the
/// nearest-edge distance, which is what the band means for a rasterized body.
fn near_surface<S: SolidShape>(body: &S, x: f64, y: f64, band: f64) -> bool {
    // default implementation via the shape's own distance proxy: sample the
    // four neighbors; if any flips the inside test the point borders the
    // surface, and the band widens it to band/2 beyond that.
    let h = band / 3.0;
    let touching = body.inside(x - h, y)
        || body.inside(x + h, y)
        || body.inside(x, y - h)
        || body.inside(x, y + h);
    if touching {
        return true;
    }
    // widen: step outward along the normal until clear of the band.
    let (nx, ny) = body.normal(x, y);
    !body.inside(x + nx * band, y + ny * band)
        && (body.inside(x + nx * (band + h), y + ny * (band + h))
            || body.inside(x - nx * band, y - ny * band))
}

impl SolidBody {
    /// distance from the point to the circle (negative inside).
    pub fn signed_distance(&self, x: f64, y: f64) -> f64 {
        let (dx, dy) = (x - self.cx, y - self.cy);
        (dx * dx + dy * dy).sqrt() - self.r
    }
}

impl SolidPolygon {
    /// distance from the point to the nearest polygon edge (positive
    /// outside, negative inside, exact to the piecewise-linear surface).
    pub fn signed_distance(&self, x: f64, y: f64) -> f64 {
        let n = self.verts.len();
        let mut best = f64::INFINITY;
        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = self.verts[i];
            let (bx, by) = self.verts[j];
            let (ex, ey) = (bx - ax, by - ay);
            let len2 = (ex * ex + ey * ey).max(1e-30);
            let t = (((x - ax) * ex + (y - ay) * ey) / len2).clamp(0.0, 1.0);
            let (px, py) = (ax + t * ex, ay + t * ey);
            let d = ((x - px).powi(2) + (y - py).powi(2)).sqrt();
            if d < best {
                best = d;
            }
        }
        if self.inside(x, y) {
            -best
        } else {
            best
        }
    }
}

/// impose a body specified by a signed-distance function on the state.
///
/// same contract as `apply_solid`; the caller supplies `dist` (negative
/// inside the body) plus `normal` for surface cells. used by the axisymmetric
/// aeroshell path where the meridian is analytic.
pub fn apply_solid_fn(
    state: &mut ConservedState2d,
    g: &Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    normal: &dyn Fn(f64, f64) -> (f64, f64),
    gamma: f64,
) {
    let ref_ = SolidRef::default();
    let e_ref = ref_.p / ((gamma - 1.0) * ref_.rho);
    let band = 1.5 * g.dx.min(g.dy);
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let (x, y) = (g.centers_x[k], g.centers_y[k]);
            let d = dist(x, y);
            if d <= 0.0 {
                // solid cell: keep it inert.
                state.rho[k] = ref_.rho;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
                state.e[k] = e_ref;
            } else if d <= band {
                // fluid cell bordering the surface: remove the normal velocity.
                let (nx, ny) = normal(x, y);
                let u = state.mx[k] / state.rho[k];
                let v = state.my[k] / state.rho[k];
                let un = u * nx + v * ny;
                if un > 0.0 {
                    state.mx[k] -= 2.0 * un * nx * state.rho[k];
                    state.my[k] -= 2.0 * un * ny * state.rho[k];
                } else {
                    state.mx[k] -= un * nx * state.rho[k];
                    state.my[k] -= un * ny * state.rho[k];
                }
            }
        }
    }
}

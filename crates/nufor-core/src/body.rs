//! immersed solid bodies on the structured 2d grid.
//!
//! a circular body is cut out of the cartesian mesh: cells inside the body are
//! frozen to a quiescent reference, and the band of fluid cells adjacent to the
//! surface has the velocity component normal to the surface removed, which
//! enforces a slip wall (no-normal-flow) at the immersed boundary. this turns
//! a plain cartesian solver into one that can run flow over a blunt body.

use crate::grid2d::Grid2d;
use crate::state2d::ConservedState2d;

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

impl SolidBody {
    /// whether the point falls inside the body.
    pub fn inside(&self, x: f64, y: f64) -> bool {
        let (dx, dy) = (x - self.cx, y - self.cy);
        dx * dx + dy * dy <= self.r * self.r
    }

    /// the outward unit normal at a point outside the body.
    pub fn normal(&self, x: f64, y: f64) -> (f64, f64) {
        let (dx, dy) = (x - self.cx, y - self.cy);
        let l = (dx * dx + dy * dy).sqrt().max(1e-30);
        (dx / l, dy / l)
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
pub fn apply_solid(state: &mut ConservedState2d, g: &Grid2d, body: &SolidBody, gamma: f64) {
    let ref_ = SolidRef::default();
    let e_ref = ref_.p / ((gamma - 1.0) * ref_.rho);
    let band = 1.5 * g.dx.min(g.dy);
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let (x, y) = (g.centers_x[k], g.centers_y[k]);
            let (dx, dy) = (x - body.cx, y - body.cy);
            let d = (dx * dx + dy * dy).sqrt();
            if d <= body.r {
                // solid cell: keep it inert.
                state.rho[k] = ref_.rho;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
                state.e[k] = e_ref;
            } else if d - body.r <= band {
                // fluid cell bordering the surface: remove the normal velocity.
                let (nx, ny) = (dx / d, dy / d);
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

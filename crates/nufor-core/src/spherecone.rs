//! the sphere-cone meridian: an analytic body-of-revolution profile
//! (spherical nose cap tangent to a conical flank) used by the body mask
//! and the axisymmetric runs.
//!
//! the shape is the capsule forebody from brent's trade study, an
//! apollo-style sphere-cone: a spherical nose cap of radius R_N whose arc
//! grows downstream until it is tangent to a conical flank of half-angle
//! delta, which flares out to the base radius R_b at the base plane. the
//! meridian is the upper half; the body of revolution is the full shape.
//! the nose sits at `xc - rn`, the flank is tangent at `x_tangent`, and
//! the base plane closes the shape at `x_base`.
//!
//! sketch (x downstream, y the radius):
//!
//!   y ^            base
//!     |           .-'
//!     |    flank /
//!     |  .-----'   tangent point
//!     | .'  spherical cap (R_N)
//!     |'
//!     +-----------------> x

/// the sphere-cone parameters, nose forward: x measured downstream from
/// the nose, y is the radius from the symmetry axis.
#[derive(Debug, Clone, Copy)]
pub struct SphereCone {
    /// nose radius R_N.
    pub rn: f64,
    /// cone half-angle delta (radians).
    pub delta: f64,
    /// base radius R_b.
    pub rb: f64,
    /// axial station of the sphere center (the nose sits at xc - rn).
    pub xc: f64,
}

impl SphereCone {
    /// the reference capsule geometry: R_N = 1.90 m, base radius 2.365 m at
    /// x_base = 4.5 m, tangent at x = 1.584 m. returns the shape with the
    /// nose at x = 0.
    pub fn brent_shell() -> Self {
        // delta from the tangent-geometry closure: the sphere arc reaches
        // the tangent point (xt, y_t) with slope tan(delta), then the flank
        // runs to (xb, rb). with t = tan(delta):
        //   rb - y_t = (xb - xt) * t,  y_t = rn / sqrt(1 + t^2)
        // => rb - rn / sqrt(1+t^2) = (xb - xt) * t, solved by bisection.
        let (rn, rb, xt, xb) = (1.90, 2.365, 1.584, 4.5);
        let solve = |t: f64| rb - rn / (1.0 + t * t).sqrt() - (xb - xt) * t;
        let (mut lo, mut hi) = (1e-4f64, 2.0f64);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if solve(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t = 0.5 * (lo + hi);
        let delta = t.atan();
        let xc = xt - rn * t / (1.0 + t * t).sqrt();
        SphereCone { rn, delta, rb, xc }
    }

    /// the axial station where the cap meets the flank (the tangent point).
    pub fn x_tangent(&self) -> f64 {
        self.xc + self.rn * self.delta.sin()
    }

    /// the radius at the tangent point.
    pub fn y_tangent(&self) -> f64 {
        self.rn * self.delta.cos()
    }

    /// the axial station of the base plane.
    pub fn x_base(&self) -> f64 {
        self.x_tangent() + (self.rb - self.y_tangent()) / self.delta.tan()
    }

    /// the meridian radius at axial station x (the upper surface).
    pub fn radius_at(&self, x: f64) -> f64 {
        let xt = self.x_tangent();
        if x <= xt {
            // spherical cap: circle centered at (xc, 0).
            let arg = (self.rn * self.rn - (x - self.xc) * (x - self.xc)).max(0.0);
            arg.sqrt()
        } else {
            // conical flank: radius grows downstream at tan(delta).
            self.y_tangent() + (x - xt) * self.delta.tan()
        }
    }

    /// the meridian as points: cap arc then flank, nose to base.
    pub fn meridian_points(&self, x_base: f64, n_cap: usize) -> Vec<(f64, f64)> {
        let x_nose = self.xc - self.rn;
        let xt = self.x_tangent();
        let mut pts = Vec::with_capacity(n_cap + 8);
        for i in 0..=n_cap {
            let x = x_nose + (xt - x_nose) * (i as f64) / (n_cap as f64);
            pts.push((x, self.radius_at(x)));
        }
        pts.push((x_base, self.radius_at(x_base)));
        pts
    }

    /// the closed body-of-revolution meridian polygon: outer surface from
    /// nose to base corner, the base edge across, the mirrored surface
    /// back to the nose.
    pub fn body_polygon(&self, x_base: f64, n_cap: usize) -> Vec<(f64, f64)> {
        let upper = self.meridian_points(x_base, n_cap);
        let mut pts = upper.clone();
        // base corner down to the axis.
        pts.push((x_base, 0.0));
        // mirror the upper surface back along the axis to the nose.
        for &(x, y) in upper.iter().rev() {
            pts.push((x, -y));
        }
        // close at the nose on the axis.
        pts.push((self.xc - self.rn, 0.0));
        pts
    }
}

/// the sphere-cone as a signed-distance body: negative inside, positive
/// outside, with the outward normal from the analytic surface.
pub struct SphereConeSdf {
    pub sc: SphereCone,
}

impl SphereConeSdf {
    /// signed distance to the body of revolution, negative inside.
    ///
    /// the interior is: nose <= x <= base, |y| <= radius(x). distances are
    /// exact on the cap (circle) and the flank (line) where they govern.
    pub fn dist(&self, x: f64, y: f64) -> f64 {
        let sc = &self.sc;
        let xt = sc.x_tangent();
        let xb = sc.x_base();
        let nose = sc.xc - sc.rn;

        // distance to the cap circle (governs for x < xt).
        let d_cap = ((x - sc.xc).powi(2) + y * y).sqrt() - sc.rn;
        // distance to the flank line: y = y_t + (x - xt) * tan(d); the
        // outward side is above the line. normalized line distance,
        // positive above (outside).
        let (td, y_t) = (sc.delta.tan(), sc.y_tangent());
        let a = td;
        let c = y_t - a * xt;
        let d_flank = -(a * x - y + c) / (a * a + 1.0).sqrt();
        // distance to the base plane (outward is +x past xb).
        let d_base = x - xb;

        let inside = x > nose && x < xb && y.abs() < sc.radius_at(x);
        let d = if x < xt {
            // the spherical cap governs.
            d_cap
        } else if x <= xb {
            // the flank governs the whole band from the tangent station
            // to the base; its sign handles both sides of the line.
            d_flank
        } else if y.abs() > sc.rb {
            // past the base and outside the base radius: the corner.
            ((x - xb).powi(2) + (y.abs() - sc.rb).powi(2)).sqrt()
        } else {
            // past the base within the radius: the base plane.
            d_base
        };
        if inside {
            // interior: negative distance to the nearest face; the band
            // never queries this branch, only the sign matters plus a
            // sane magnitude near the surface.
            let mut m = d.abs();
            if x < xt {
                m = m.min(d_cap.abs());
            } else {
                m = m.min(d_flank.abs());
            }
            m = m.min(d_base.abs());
            -m
        } else {
            d
        }
    }

    /// the outward unit normal on the surface nearest the point.
    pub fn normal(&self, x: f64, y: f64) -> (f64, f64) {
        let sc = &self.sc;
        let xt = sc.x_tangent();
        let xb = sc.x_base();
        if x > xb {
            // base face: outward is +x.
            (1.0, 0.0)
        } else if x < xt {
            // cap: radial from the sphere center (outward).
            let (dx, dy) = (x - sc.xc, y);
            let l = (dx * dx + dy * dy).sqrt().max(1e-30);
            (dx / l, dy / l)
        } else {
            // flank: the outward normal of the upper surface points away
            // from the body interior: upstream and up, (-sin d, cos d).
            let d = sc.delta;
            let sgn = if y >= 0.0 { 1.0 } else { -1.0 };
            (-d.sin(), sgn * d.cos())
        }
    }
}

//! the blasius flat-plate similarity solution: the laminar boundary layer
//! the profile-inflow boundary condition prescribes and the first
//! quantitative validation target (cf = 2 f''(0)/sqrt(re_x)).
//!
//! the similarity ode f''' + 0.5 f f'' = 0 with f(0) = f'(0) = 0 and
//! f'(inf) = 1 is integrated once at construction with fixed-step rk4 over
//! eta in [0, 12] (h = 2e-3, f'(12) = 1 to roundoff with the literature
//! shooting value f''(0) = 0.3320573362151963), and the stored rows are
//! interpolated monotonically. a physical (x, y) evaluation rescales eta
//! with the local similarity variable y*sqrt(u_inf/(nu*x)).

use crate::Error;

/// the integrated similarity table: rows of (eta, f, f', f'').
#[derive(Debug, Clone, PartialEq)]
pub struct BlasiusTable {
    etas: Vec<f64>,
    fs: Vec<f64>,
    fps: Vec<f64>,
    fpps: Vec<f64>,
}

/// the shoot value every textbook quotes; rk4 with h=2e-3 lands f'(12) = 1
/// to roundoff from it, and the table reproduces the canonical f' rows.
const FPP0: f64 = 0.3320573362151963;
const ETA_MAX: f64 = 12.0;
const H: f64 = 2e-3;

impl BlasiusTable {
    /// integrate the similarity ode once.
    pub fn integrate() -> Self {
        let n = (ETA_MAX / H) as usize;
        let mut etas = Vec::with_capacity(n + 1);
        let mut fs = Vec::with_capacity(n + 1);
        let mut fps = Vec::with_capacity(n + 1);
        let mut fpps = Vec::with_capacity(n + 1);
        let (mut f, mut fp, mut fpp) = (0.0f64, 0.0f64, FPP0);
        let mut eta = 0.0f64;
        for _ in 0..=n {
            etas.push(eta);
            fs.push(f);
            fps.push(fp);
            fpps.push(fpp);
            // one rk4 step of (f, f', f'') with f''' = -0.5 f f''.
            let k1 = (fp, fpp, -0.5 * f * fpp);
            let (f2, fp2, fpp2) = (
                f + 0.5 * H * k1.0,
                fp + 0.5 * H * k1.1,
                fpp + 0.5 * H * k1.2,
            );
            let k2 = (fp2, fpp2, -0.5 * f2 * fpp2);
            let (f3, fp3, fpp3) = (
                f + 0.5 * H * k2.0,
                fp + 0.5 * H * k2.1,
                fpp + 0.5 * H * k2.2,
            );
            let k3 = (fp3, fpp3, -0.5 * f3 * fpp3);
            let (f4, fp4, fpp4) = (f + H * k3.0, fp + H * k3.1, fpp + H * k3.2);
            let k4 = (fp4, fpp4, -0.5 * f4 * fpp4);
            f += (H / 6.0) * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
            fp += (H / 6.0) * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
            fpp += (H / 6.0) * (k1.2 + 2.0 * k2.2 + 2.0 * k3.2 + k4.2);
            eta += H;
        }
        BlasiusTable {
            etas,
            fs,
            fps,
            fpps,
        }
    }

    /// f'(eta) by monotone interpolation: u/u_inf for the layer.
    pub fn fp(&self, eta: f64) -> f64 {
        let i = self.index(eta);
        let (e0, e1) = (self.etas[i], self.etas[i + 1]);
        let w = if e1 > e0 { (eta - e0) / (e1 - e0) } else { 0.0 };
        let v = self.fps[i] + w * (self.fps[i + 1] - self.fps[i]);
        v.clamp(0.0, 1.0)
    }

    /// f(eta) by monotone interpolation.
    pub fn f(&self, eta: f64) -> f64 {
        let i = self.index(eta);
        let (e0, e1) = (self.etas[i], self.etas[i + 1]);
        let w = if e1 > e0 { (eta - e0) / (e1 - e0) } else { 0.0 };
        self.fs[i] + w * (self.fs[i + 1] - self.fs[i])
    }

    /// f''(eta) by monotone interpolation (the wall shear at eta=0).
    pub fn fpp(&self, eta: f64) -> f64 {
        let i = self.index(eta);
        let (e0, e1) = (self.etas[i], self.etas[i + 1]);
        let w = if e1 > e0 { (eta - e0) / (e1 - e0) } else { 0.0 };
        self.fpps[i] + w * (self.fpps[i + 1] - self.fpps[i])
    }

    /// the bracketing index for eta, clamped to the table range.
    fn index(&self, eta: f64) -> usize {
        if eta <= 0.0 {
            return 0;
        }
        if eta >= ETA_MAX {
            return self.etas.len() - 2;
        }
        // uniform step, so the bracket is a direct scale.
        ((eta / H) as usize).min(self.etas.len() - 2)
    }
}

/// the wall shear coefficient constant: cf = 2 f''(0)/sqrt(re_x).
pub const CF_CONST: f64 = 2.0 * FPP0;

/// the blasius layer over a plate: freestream speed u_inf, kinematic
/// viscosity nu, a virtual leading edge at x_le (the station where the
/// layer starts, which may sit upstream of the domain), and the x station
/// where an inflow plane evaluates the profile.
#[derive(Debug, Clone, PartialEq)]
pub struct BlasiusProfile {
    table: BlasiusTable,
    /// freestream speed.
    pub u_inf: f64,
    /// kinematic viscosity.
    pub nu: f64,
    /// the virtual leading-edge x station.
    pub x_le: f64,
    /// the x station an inflow plane anchors the profile at.
    pub x_in: f64,
}

impl BlasiusProfile {
    /// build the layer; nu and u_inf must be positive.
    pub fn new(u_inf: f64, nu: f64, x_le: f64) -> Result<Self, Error> {
        if u_inf <= 0.0 || nu <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        Ok(BlasiusProfile {
            table: BlasiusTable::integrate(),
            u_inf,
            nu,
            x_le,
            x_in: x_le,
        })
    }

    /// anchor the inflow-plane evaluation at x (the distance from the
    /// virtual leading edge to the inflow face).
    pub fn anchored_at(mut self, x_in: f64) -> Self {
        self.x_in = x_in;
        self
    }

    /// the local similarity variable eta = y sqrt(u_inf/(nu (x - x_le))).
    fn eta(&self, x: f64, y: f64) -> f64 {
        let dx = (x - self.x_le).max(1e-12);
        y * (self.u_inf / (self.nu * dx)).sqrt()
    }

    /// the streamwise velocity u(x, y).
    pub fn u(&self, x: f64, y: f64) -> f64 {
        self.u_inf * self.table.fp(self.eta(x, y))
    }

    /// the wall-normal velocity v(x, y) from the similarity form.
    pub fn v(&self, x: f64, y: f64) -> f64 {
        let eta = self.eta(x, y);
        let dx = (x - self.x_le).max(1e-12);
        // v/u_inf = 0.5 (eta f' - f)/sqrt(re_x), re_x = u_inf dx / nu.
        let re_x = (self.u_inf * dx / self.nu).sqrt();
        0.5 * self.u_inf * (eta * self.table.fp(eta) - self.table.f(eta)) / re_x
    }

    /// the analytic wall shear tau_w(x) = mu u_inf f''(0) sqrt(u_inf/(nu dx)).
    pub fn tau_w(&self, x: f64) -> f64 {
        let dx = (x - self.x_le).max(1e-12);
        self.nu * self.u_inf * self.u_inf * FPP0 * (self.u_inf / (self.nu * dx)).sqrt()
    }

    /// cf(x) = 2 tau_w/(rho u_inf^2) with the local freestream density.
    pub fn cf(&self, x: f64, rho_inf: f64) -> f64 {
        2.0 * self.tau_w(x) / (rho_inf * self.u_inf * self.u_inf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// canonical f' rows from the textbook blasius table.
    const CANON: &[(f64, f64)] = &[
        (0.0, 0.0),
        (0.5, 0.165968),
        (1.0, 0.329780),
        (2.0, 0.629771),
        (3.0, 0.846050),
        (4.0, 0.955518),
        (5.0, 0.991551),
        (6.0, 0.998977),
    ];

    #[test]
    fn table_reproduces_the_canonical_rows() {
        let t = BlasiusTable::integrate();
        for &(eta, want) in CANON {
            let got = t.fp(eta);
            assert!(
                (got - want).abs() < 1e-4,
                "f'({eta}) = {got:.6}, want {want:.6}"
            );
        }
        // the far field closes on 1.
        assert!((t.fp(12.0) - 1.0).abs() < 1e-9);
        // monotone in eta.
        let mut prev = -1.0;
        for k in (0..t.etas.len()).step_by(50) {
            assert!(t.fps[k] >= prev);
            prev = t.fps[k];
        }
    }

    #[test]
    fn wall_values_match_the_literature() {
        let t = BlasiusTable::integrate();
        assert!((t.fpp(0.0) - 0.3320573362151963).abs() < 1e-9);
        assert!((t.f(5.0) - 3.283274).abs() < 2e-3);
        // cf constant: 2 f''(0) = 0.6641.
        assert!((CF_CONST - 0.664115).abs() < 1e-5);
    }

    #[test]
    fn profile_evaluates_a_physical_layer() {
        let p = BlasiusProfile::new(0.2, 1e-4, 0.0).unwrap();
        // at the wall u = 0; far above the layer u = u_inf.
        assert!(p.u(1.0, 0.0).abs() < 1e-12);
        assert!((p.u(1.0, 0.3) - 0.2).abs() < 1e-6);
        // the wall shear decays like 1/sqrt(x): cf(4x) = cf(x)/2.
        let c1 = p.cf(0.5, 1.0);
        let c4 = p.cf(2.0, 1.0);
        assert!((c4 - 0.5 * c1).abs() < 1e-12);
        // v is small and positive (the displacement effect).
        let v = p.v(1.0, 0.05);
        assert!(v > 0.0 && v < 0.01 * p.u_inf * 10.0);
    }

    #[test]
    fn rejects_nonphysical_args() {
        assert!(BlasiusProfile::new(0.0, 1e-4, 0.0).is_err());
        assert!(BlasiusProfile::new(0.2, 0.0, 0.0).is_err());
        assert!(BlasiusProfile::new(-0.2, 1e-4, 0.0).is_err());
    }
}

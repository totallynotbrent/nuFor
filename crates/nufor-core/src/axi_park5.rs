//! the park5 mixture closure for one axisymmetric flux sweep.
//!
//! the closure kernel is owned by the fortran core and speaks in
//! batches. the face loop used to send the two sides of every
//! interface through it separately, so one ffi crossing and roughly
//! nine buffer allocations rode along on each face. this module lifts
//! the call to the sweep level: the left and right states of every
//! face in the sweep are interleaved into one batch, the fortran
//! kernel answers once, and the per-face values are picked back out.
//!
//! the kernel is per-cell exact — each cell's pressure and sound
//! speed derive from its own density, energy and composition with no
//! neighbour term — so raising the batch size cannot change a single
//! cell's bits. the fallback that walks the rust reference stays for
//! the error path only, as it does everywhere else the shipping lane
//! meets the fortran core.

use crate::hllc2d::FacePrim;
use crate::species_fortran::NS;
use crate::Error;

/// the mixture closure over a whole sweep's faces: pressure and sound
/// speed per cell, laid out left/right per face — cell `2*f` is the
/// left side of face `f`, cell `2*f+1` its right side.
pub(crate) struct FaceClosure {
    p: Vec<f64>,
    a: Vec<f64>,
}

impl FaceClosure {
    /// the face prim pair for face `f` under the park5 closure.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pair(
        &self,
        f: usize,
        rl: f64,
        ul: f64,
        vl: f64,
        xl: f64,
        rr: f64,
        ur: f64,
        vr: f64,
        xr: f64,
    ) -> (FacePrim, FacePrim) {
        (
            FacePrim {
                rho: rl,
                u: ul,
                v: vl,
                p: self.p[2 * f],
                a: self.a[2 * f],
                e: xl,
            },
            FacePrim {
                rho: rr,
                u: ur,
                v: vr,
                p: self.p[2 * f + 1],
                a: self.a[2 * f + 1],
                e: xr,
            },
        )
    }
}

/// the batched park5 closure for one sweep. `rho_*`, `e_*` and the per
/// species rows are indexed by face; the two sides of face `f` are the
/// left and right reconstructed states arriving with it. the composition
/// rows may be absent altogether (a park5 model without a species field
/// pads with zeros, exactly as the per-face call did).
pub(crate) fn sweep_closure(
    rho_l: &[f64],
    rho_r: &[f64],
    e_l: &[f64],
    e_r: &[f64],
    y_l: &[Vec<f64>],
    y_r: &[Vec<f64>],
) -> Result<FaceClosure, Error> {
    let nf = rho_l.len();
    if rho_r.len() != nf || e_l.len() != nf || e_r.len() != nf {
        return Err(Error::InvalidArgs);
    }
    let m = 2 * nf;
    let mut rho = vec![0.0f64; m];
    let mut e = vec![0.0f64; m];
    let mut planes: Vec<Vec<f64>> = vec![vec![0.0f64; m]; NS];
    for f in 0..nf {
        rho[2 * f] = rho_l[f];
        rho[2 * f + 1] = rho_r[f];
        e[2 * f] = e_l[f];
        e[2 * f + 1] = e_r[f];
        for (s, plane) in planes.iter_mut().enumerate() {
            if let Some(row) = y_l.get(s) {
                plane[2 * f] = row[f];
            }
            if let Some(row) = y_r.get(s) {
                plane[2 * f + 1] = row[f];
            }
        }
    }
    match crate::species_fortran::closure(&rho, &e, &planes) {
        Ok(f) => Ok(FaceClosure { p: f.p, a: f.a }),
        // the rust reference stays the yardstick on the error path only.
        Err(_) => {
            let mut p = vec![0.0f64; m];
            let mut a = vec![0.0f64; m];
            for k in 0..m {
                let mut y = [0.0f64; NS];
                for (s, dst) in y.iter_mut().enumerate() {
                    *dst = planes[s][k];
                }
                p[k] = crate::species_closure::pressure_at(rho[k], e[k], &y);
                a[k] = crate::species_closure::sound_speed_at(rho[k], e[k], &y);
            }
            Ok(FaceClosure { p, a })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the batched sweep closure must agree with the rust reference for
    /// every cell, and the batch must not couple neighbouring cells.
    #[test]
    fn sweep_closure_matches_the_rust_reference() {
        let nf = 6;
        let mut rho_l = Vec::new();
        let mut rho_r = Vec::new();
        let mut e_l = Vec::new();
        let mut e_r = Vec::new();
        let mut y_l = vec![Vec::new(); NS];
        let mut y_r = vec![Vec::new(); NS];
        for f in 0..nf {
            let rho = 0.5 + 0.37 * f as f64;
            rho_l.push(rho);
            rho_r.push(rho * 1.11);
            e_l.push(900.0 + 130.0 * f as f64);
            e_r.push(950.0 + 90.0 * f as f64);
            for s in 0..NS {
                let frac = 0.2 + 0.1 * (s + f) as f64;
                y_l[s].push(frac);
                y_r[s].push(frac * 0.93);
            }
        }
        let closure = sweep_closure(&rho_l, &rho_r, &e_l, &e_r, &y_l, &y_r).expect("closure");
        for f in 0..nf {
            for (side, rho, e, ys) in [
                (0usize, rho_l[f], e_l[f], &y_l),
                (1usize, rho_r[f], e_r[f], &y_r),
            ] {
                let mut y = [0.0f64; NS];
                for (s, dst) in y.iter_mut().enumerate() {
                    *dst = ys[s][f];
                }
                let k = 2 * f + side;
                let p_ref = crate::species_closure::pressure_at(rho, e, &y);
                let a_ref = crate::species_closure::sound_speed_at(rho, e, &y);
                let tol = 1.0e-12 * p_ref.abs().max(1.0);
                assert!(
                    (closure.p[k] - p_ref).abs() <= tol,
                    "pressure differs at face {f} side {side}: {} vs {p_ref}",
                    closure.p[k]
                );
                let atol = 1.0e-12 * a_ref.abs().max(1.0);
                assert!(
                    (closure.a[k] - a_ref).abs() <= atol,
                    "sound speed differs at face {f} side {side}: {} vs {a_ref}",
                    closure.a[k]
                );
            }
        }
    }

    /// a mismatched input length is an argument error, not a silent read.
    #[test]
    fn sweep_closure_rejects_mismatched_lengths() {
        let err = sweep_closure(&[1.0, 2.0], &[1.0], &[1.0, 2.0], &[1.0, 2.0], &[], &[]);
        assert!(matches!(err, Err(Error::InvalidArgs)));
    }
}

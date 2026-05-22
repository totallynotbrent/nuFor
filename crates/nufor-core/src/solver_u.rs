//! cell-centered unstructured finite volume: rotate the hllc flux by each face
//! normal and accumulate the divergence face by face.

use crate::eos2d::eos_pressure2d;
use crate::error::Error;
use crate::hllc2d::{hllc_flux, FacePrim};
use crate::state2d::{cons_to_prim2d, ConservedState2d};
use crate::ugrid::Ugrid;

/// the hllc flux of a pair of primitive states rotated onto a face normal.
fn flux_on_normal(gamma: f64, left: FacePrim, right: FacePrim, nx: f64, ny: f64) -> [f64; 4] {
    let un = |u: f64, v: f64| u * nx + v * ny;
    let ut = |u: f64, v: f64| -u * ny + v * nx;
    let l = FacePrim {
        rho: left.rho,
        u: un(left.u, left.v),
        v: ut(left.u, left.v),
        p: left.p,
    };
    let r = FacePrim {
        rho: right.rho,
        u: un(right.u, right.v),
        v: ut(right.u, right.v),
        p: right.p,
    };
    let f = hllc_flux(gamma, l, r, 0);
    // rotate the face-frame flux back: normal component is [1], tangent [2].
    [f.mass, f.mx * nx - f.my * ny, f.mx * ny + f.my * nx, f.e]
}

/// one conservative step for a cell-centered unstructured state.
pub fn advance_ugrid(
    state: &mut ConservedState2d,
    ug: &Ugrid,
    gamma: f64,
    cfl: f64,
) -> Result<(f64, f64), Error> {
    let nc = ug.cell_area.len();
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    // a robust per-cell length scale for anisotropic cells: the inscribed width
    // area / (half perimeter), so a tall thin cell sizes by its narrow side.
    let mut h = f64::INFINITY;
    let mut perim = vec![0.0; nc];
    for f in 0..ug.face_left.len() {
        let len = (ug.face_dx[f] * ug.face_dx[f] + ug.face_dy[f] * ug.face_dy[f]).sqrt();
        perim[ug.face_left[f]] += len;
        if let Some(rc) = (ug.face_right[f] >= 0).then(|| ug.face_right[f] as usize) {
            perim[rc] += len;
        }
    }
    for (area, p) in ug.cell_area.iter().zip(&perim) {
        h = h.min(area / (0.5 * p).max(1e-12));
    }
    let mut smax = 0.0f64;
    for c in 0..nc {
        let a = (gamma * p[c] / state.rho[c].max(1e-12)).sqrt();
        smax = smax.max(u[c].abs().max(v[c].abs()) + a);
    }
    let dt = cfl * h / smax.max(1e-12);

    // per-cell divergence accumulators (the negative of the outflow).
    let (mut dr, mut dmx, mut dmy, mut de) =
        (vec![0.0; nc], vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]);
    for f in 0..ug.face_left.len() {
        let nxl = ug.face_dx[f];
        let nyl = ug.face_dy[f];
        let len = (nxl * nxl + nyl * nyl).sqrt();
        let (nx, ny) = (nxl / len, nyl / len);
        let l = ug.face_left[f];
        let r = ug.face_right[f];
        let li = l;
        let ri = if r >= 0 { Some(r as usize) } else { None };
        // transmissive boundary: the far state mirrors the owning cell.
        let (rr, rux, rv, rp) = if let Some(rc) = ri {
            (state.rho[rc], u[rc], v[rc], p[rc])
        } else {
            (state.rho[li], u[li], v[li], p[li])
        };
        let flux = flux_on_normal(
            gamma,
            FacePrim {
                rho: state.rho[li],
                u: u[li],
                v: v[li],
                p: p[li],
            },
            FacePrim {
                rho: rr,
                u: rux,
                v: rv,
                p: rp,
            },
            nx,
            ny,
        );
        // the flux leaves the left cell (outward normal) and enters the right.
        dr[li] -= len * flux[0];
        dmx[li] -= len * flux[1];
        dmy[li] -= len * flux[2];
        de[li] -= len * flux[3];
        if let Some(rc) = ri {
            dr[rc] += len * flux[0];
            dmx[rc] += len * flux[1];
            dmy[rc] += len * flux[2];
            de[rc] += len * flux[3];
        }
    }
    let mut resid = 0.0f64;
    for c in 0..nc {
        let s = dt / ug.cell_area[c];
        state.rho[c] += s * dr[c];
        state.mx[c] += s * dmx[c];
        state.my[c] += s * dmy[c];
        state.e[c] += s * de[c];
        resid = resid
            .max(s * dr[c].abs())
            .max(s * dmx[c].abs())
            .max(s * de[c].abs());
    }
    Ok((dt, resid))
}

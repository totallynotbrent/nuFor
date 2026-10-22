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
        a: left.a,
        e: left.e,
    };
    let r = FacePrim {
        rho: right.rho,
        u: un(right.u, right.v),
        v: ut(right.u, right.v),
        p: right.p,
        a: right.a,
        e: right.e,
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
    let (mut dr_acc, mut dmx_acc, mut dmy_acc, mut de_acc) =
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

        // a slip wall takes the exact flux [0, p n, 0]: no riemann
        // solve, no mirror state, so a uniform field holds bit-exactly
        // and the wall pressure is the cell pressure by construction.
        if ri.is_none() && ug.face_group[f] == 1 {
            dr_acc[li] -= 0.0;
            dmx_acc[li] -= p[li] * nxl;
            dmy_acc[li] -= p[li] * nyl;
            de_acc[li] -= 0.0;
            continue;
        }

        // boundary groups: 2 freestream inflow, 3 outflow, 0
        // transmissive mirror.
        let (rr, rux, rv, rp) = if let Some(rc) = ri {
            (state.rho[rc], u[rc], v[rc], p[rc])
        } else {
            match ug.face_group[f] {
                2 => (
                    ug.freestream[0],
                    ug.freestream[1],
                    ug.freestream[2],
                    ug.freestream[3],
                ),
                3 => (state.rho[li], u[li], v[li], p[li]),
                _ => (state.rho[li], u[li], v[li], p[li]),
            }
        };
        let flux = flux_on_normal(
            gamma,
            FacePrim::perfect(state.rho[li], u[li], v[li], p[li], gamma),
            FacePrim::perfect(rr, rux, rv, rp, gamma),
            nx,
            ny,
        );
        // the flux leaves the left cell (outward normal) and enters the right.
        dr_acc[li] -= len * flux[0];
        dmx_acc[li] -= len * flux[1];
        dmy_acc[li] -= len * flux[2];
        de_acc[li] -= len * flux[3];
        if let Some(rc) = ri {
            dr_acc[rc] += len * flux[0];
            dmx_acc[rc] += len * flux[1];
            dmy_acc[rc] += len * flux[2];
            de_acc[rc] += len * flux[3];
        }
    }
    // positivity floors tied to the case freestream (set once by the
    // caller from the IC, never re-derived from the marching state: a
    // depressed cell would ratchet the floor down each step). cells
    // are clamped to 1e-3 of the freestream density and internal
    // energy, far below any physical state but well above vacuum.
    let rho_floor = (ug.freestream[0] * 1e-3).max(1e-12);
    let e_floor = (ug.freestream[3] / (gamma - 1.0) * 1e-3).max(1e-12);
    let mut resid = 0.0f64;
    for c in 0..nc {
        let s = dt / ug.cell_area[c];
        state.rho[c] += s * dr_acc[c];
        state.mx[c] += s * dmx_acc[c];
        state.my[c] += s * dmy_acc[c];
        state.e[c] += s * de_acc[c];
        resid = resid
            .max(s * dr_acc[c].abs())
            .max(s * dmx_acc[c].abs())
            .max(s * de_acc[c].abs());
        // positivity floors: a first-order cell-average update can push
        // a cell slightly negative in the shock layer at extreme mach.
        // the floors come from the case freestream, set once by the
        // caller, so a sod tube at p ~ 0.1 is clamped just as little
        // as a reentry flow at p ~ 1e5.
        if state.rho[c] < rho_floor {
            state.rho[c] = rho_floor;
        }
        let ke = 0.5 * (state.mx[c] * state.mx[c] + state.my[c] * state.my[c]) / state.rho[c];
        if state.e[c] - ke < e_floor {
            state.e[c] = ke + e_floor;
        }
    }
    Ok((dt, resid))
}

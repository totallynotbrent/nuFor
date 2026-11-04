//! laminar viscous terms on the unstructured (dual-control-volume) mesh.
//!
//! mirrors the structured `viscous2d` kernels on the Ugrid face loop:
//! green-gauss cell gradients of u, v, and T drive a Newtonian stress +
//! conduction flux contracted with each face's outward area vector, added
//! with the same sign convention as the inviscid term in `advance_ugrid`
//! (momentum leaving the left cell subtracts; entering the right adds).
//! wall-group faces use one-sided no-slip gradients toward the wall state.

use crate::state2d::{cons_to_prim2d, specific_kinetic_energy};
use crate::ugrid::Ugrid;
use crate::viscous2d::sutherland_mu;
use crate::ConservedState2d;

/// neighbor-mean face hint of a cell-centered field (copy on boundary).
fn face_mean(ug: &Ugrid, phi: &[f64]) -> Vec<f64> {
    (0..ug.face_left.len())
        .map(|f| {
            let l = ug.face_left[f];
            match (ug.face_right[f] >= 0).then(|| ug.face_right[f] as usize) {
                Some(r) => 0.5 * (phi[l] + phi[r]),
                None => phi[l],
            }
        })
        .collect()
}

/// green-gauss cell gradient of a primitive field, via face-mean hints.
#[allow(clippy::needless_range_loop)]
fn gg_grad(ug: &Ugrid, phi_face: &[f64]) -> Vec<(f64, f64)> {
    let nc = ug.cell_area.len();
    let mut gx = vec![0.0f64; nc];
    let mut gy = vec![0.0f64; nc];
    for f in 0..ug.face_left.len() {
        let l = ug.face_left[f];
        let (nx, ny) = (ug.face_dx[f], ug.face_dy[f]);
        let pf = phi_face[f];
        gx[l] += pf * nx;
        gy[l] += pf * ny;
        if let Some(r) = (ug.face_right[f] >= 0).then(|| ug.face_right[f] as usize) {
            gx[r] -= pf * nx;
            gy[r] -= pf * ny;
        }
    }
    for c in 0..nc {
        let a = ug.cell_area[c].max(1e-15);
        gx[c] /= a;
        gy[c] /= a;
    }
    gx.into_iter().zip(gy).collect()
}

/// viscous tendencies (mass, x-momentum, y-momentum, energy) per cell for
/// one step of size dt, NOT yet divided by cell area. interior faces use the
/// neighbor-averaged green-gauss gradient; faces in `wall_group` use a
/// one-sided gradient toward the no-slip wall state over the face length.
#[allow(clippy::too_many_arguments)]
pub fn viscous_ugrid_tendency(
    state: &ConservedState2d,
    ug: &Ugrid,
    dt: f64,
    mu0: f64,
    t0: f64,
    s_ref: f64,
    pr: f64,
    cp: f64,
    wall_group: u32,
    wall_temp: f64,
    wall_u: f64,
    wall_v: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let nc = ug.cell_area.len();
    let zero = || {
        (
            vec![0.0f64; nc],
            vec![0.0f64; nc],
            vec![0.0f64; nc],
            vec![0.0f64; nc],
        )
    };
    let Ok((u, v, _et)) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e) else {
        return zero();
    };
    // temperature proportional to specific internal energy (constant cp):
    // adequate for the gradient since cp factors out of the conduction term.
    let t: Vec<f64> = (0..nc)
        .map(|c| {
            let ke = specific_kinetic_energy(state.rho[c], state.mx[c], state.my[c]);
            (state.e[c] / state.rho[c] - ke).max(1.0) / cp
        })
        .collect();

    let gu = gg_grad(ug, &face_mean(ug, &u));
    let gv = gg_grad(ug, &face_mean(ug, &v));
    let gt = gg_grad(ug, &face_mean(ug, &t));

    let nf = ug.face_left.len();
    let mut acc = zero();
    for f in 0..nf {
        let li = ug.face_left[f];
        let ri = (ug.face_right[f] >= 0).then(|| ug.face_right[f] as usize);
        let len = (ug.face_dx[f] * ug.face_dx[f] + ug.face_dy[f] * ug.face_dy[f]).sqrt();
        let nx = ug.face_dx[f] / len; // outward normal of the left cell
        let ny = ug.face_dy[f] / len;
        let is_wall = ug.face_group[f] == wall_group;

        // face gradient: averaged cell gradients inside, one-sided to the
        // wall state on wall faces.
        let (gux, guy, gvx, gvy, gtx, gty) = match ri {
            Some(r) => (
                0.5 * (gu[li].0 + gu[r].0),
                0.5 * (gu[li].1 + gu[r].1),
                0.5 * (gv[li].0 + gv[r].0),
                0.5 * (gv[li].1 + gv[r].1),
                0.5 * (gt[li].0 + gt[r].0),
                0.5 * (gt[li].1 + gt[r].1),
            ),
            None if is_wall => (
                (wall_u - u[li]) / len,
                0.0,
                (wall_v - v[li]) / len,
                0.0,
                (wall_temp - t[li]) / len,
                0.0,
            ),
            // a boundary face not in the wall group is adiabatic + slip:
            // zero flux crosses it, which conserves energy globally.
            None => continue,
        };

        let tf = match ri {
            Some(r) => 0.5 * (t[li] + t[r]),
            None => 0.5 * (t[li] + wall_temp),
        };
        let mu = sutherland_mu(mu0, tf.max(1.0), t0, s_ref);
        let k = mu * cp / pr; // conductivity

        let div = gux + gvy;
        let txx = 2.0 * gux - 2.0 / 3.0 * div;
        let tyy = 2.0 * gvy - 2.0 / 3.0 * div;
        let txy = guy + gvx;

        // viscous flux through the face (force/unit length on the left cell):
        // momentum = tau . n, energy = (tau.u).n + k gradT . n. the left cell
        // LOSES this, the right GAINS it (minus sign on the left accumulator).
        let fm_x = mu * (txx * nx + txy * ny);
        let fm_y = mu * (txy * nx + tyy * ny);
        // fourier heat flux is q = -k gradT, so the energy LEAVING the left
        // cell through the face is (tau.u).n + q.n = (tau.u).n - k gradT.n.
        let fe = mu * (txx * nx + txy * ny) * u[li] + mu * (txy * nx + tyy * ny) * v[li]
            - k * (gtx * nx + gty * ny);

        acc.1[li] -= len * fm_x;
        acc.2[li] -= len * fm_y;
        acc.3[li] -= len * fe;
        if let Some(r) = ri {
            acc.1[r] += len * fm_x;
            acc.2[r] += len * fm_y;
            acc.3[r] += len * fe;
        }
    }
    // scale by dt now so callers apply area division directly.
    let _ = dt;
    acc
}

//! hllc flux for the 2d euler equations.
//!
//! computes the three-wave hllc flux across a face given the left/right
//! primitive states and the face axis (0 = vertical, 1 = horizontal); the grid
//! is axis-aligned so the rotation into the face frame is a coordinate swap.

/// the primitive state on one side of a face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FacePrim {
    /// density.
    pub rho: f64,
    /// x-velocity.
    pub u: f64,
    /// y-velocity.
    pub v: f64,
    /// pressure.
    pub p: f64,
}

/// an inviscid physical flux vector (mass, mx, my, energy) across a face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flux4 {
    pub mass: f64,
    pub mx: f64,
    pub my: f64,
    pub e: f64,
}

/// hllc flux across a face. `axis` is 0 for a vertical face (normal +x) and 1
/// for a horizontal face (normal +y).
pub fn hllc_flux(gamma: f64, left: FacePrim, right: FacePrim, axis: usize) -> Flux4 {
    // clamp to a small positive floor so a reconstruction overshoot of rho or
    // p near a shock cannot poison the sound-speed sqrt with NaN.
    let rl = left.rho.max(1e-12);
    let ul = left.u;
    let vl = left.v;
    let pl = left.p.max(1e-12);
    let rr = right.rho.max(1e-12);
    let ur = right.u;
    let vr = right.v;
    let pr = right.p.max(1e-12);
    let (un_l, ut_l) = if axis == 0 { (ul, vl) } else { (vl, ul) };
    let (un_r, ut_r) = if axis == 0 { (ur, vr) } else { (vr, ur) };
    let al = (gamma * pl / rl).sqrt();
    let ar = (gamma * pr / rr).sqrt();

    // pvrs estimate of the intermediate pressure, used only to sharpen the wave
    // speeds so a rarefaction head is not under-ranked by a low-density state.
    let rho_bar = 0.5 * (rl + rr);
    let a_bar = 0.5 * (al + ar);
    let p_star = (0.5 * (pl + pr) - 0.5 * (un_r - un_l) * rho_bar * a_bar).max(0.0);
    let ql = if p_star <= pl {
        1.0
    } else {
        (1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / pl - 1.0)).sqrt()
    };
    let qr = if p_star <= pr {
        1.0
    } else {
        (1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / pr - 1.0)).sqrt()
    };
    let sl = (un_l - al * ql).min(un_r - ar * qr);
    let sr = (un_l + al * ql).max(un_r + ar * qr);
    let denom = rl * (sl - un_l) - rr * (sr - un_r);
    let sstar = (pr - pl + rl * un_l * (sl - un_l) - rr * un_r * (sr - un_r)) / denom;

    // physical flux in the face frame: [mass, normal-mom, tangent-mom, energy];
    // energy is (rho e_t + p) so it carries the full total-enthalpy flux.
    let phys = |r: f64, un: f64, ut: f64, p: f64, e_tot: f64| {
        let en = e_tot + p;
        [r * un, r * un * un + p, r * un * ut, un * en]
    };
    let et_l = pl / (rl * (gamma - 1.0)) + 0.5 * (un_l * un_l + ut_l * ut_l);
    let et_r = pr / (rr * (gamma - 1.0)) + 0.5 * (un_r * un_r + ut_r * ut_r);
    // total energy (not specific), which the flux turns into enthalpy.
    let e_l = rl * et_l;
    let e_r = rr * et_r;
    let f_l = phys(rl, un_l, ut_l, pl, e_l);
    let f_r = phys(rr, un_r, ut_r, pr, e_r);

    // star conserved state (toro) for a side, and the star flux via rankine--
    // hugoniot: f* = f + s (u* - u).
    let ustar = |r: f64, un: f64, ut: f64, p: f64, e: f64, s_side: f64| -> [f64; 4] {
        let fac = r * (s_side - un) / (s_side - sstar);
        [
            fac,
            fac * sstar,
            fac * ut,
            fac * (e / r + (sstar - un) * (sstar + p / (r * (s_side - un)))),
        ]
    };
    let u_l = [rl, rl * un_l, rl * ut_l, e_l];
    let u_r = [rr, rr * un_r, rr * ut_r, e_r];
    let us_l = ustar(rl, un_l, ut_l, pl, e_l, sl);
    let us_r = ustar(rr, un_r, ut_r, pr, e_r, sr);
    let fstar_l = [
        f_l[0] + sl * (us_l[0] - u_l[0]),
        f_l[1] + sl * (us_l[1] - u_l[1]),
        f_l[2] + sl * (us_l[2] - u_l[2]),
        f_l[3] + sl * (us_l[3] - u_l[3]),
    ];
    let fstar_r = [
        f_r[0] + sr * (us_r[0] - u_r[0]),
        f_r[1] + sr * (us_r[1] - u_r[1]),
        f_r[2] + sr * (us_r[2] - u_r[2]),
        f_r[3] + sr * (us_r[3] - u_r[3]),
    ];

    // pick the flux from the location of the contact wave relative to the face.
    let f = match () {
        _ if sl > 0.0 => f_l,
        _ if sr < 0.0 => f_r,
        _ if sstar > 0.0 => fstar_l,
        _ => fstar_r,
    };

    // rotate the face-frame flux back to global (mass, mx, my, energy).
    if axis == 0 {
        Flux4 {
            mass: f[0],
            mx: f[1],
            my: f[2],
            e: f[3],
        }
    } else {
        Flux4 {
            mass: f[0],
            mx: f[2],
            my: f[1],
            e: f[3],
        }
    }
}

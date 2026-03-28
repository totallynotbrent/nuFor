//! hllc flux for the 3d euler equations across an axis-aligned face.
//!
//! the 3d flux has five components (mass, three momenta, energy); the third
//! velocity is a passive tangent along the face, carried through unchanged.

/// the primitive state on one side of a 3d face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FacePrim3 {
    pub rho: f64,
    pub u: f64,
    pub v: f64,
    pub w: f64,
    pub p: f64,
}

/// an inviscid physical flux vector (mass, mx, my, mz, energy) across a face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flux5 {
    pub mass: f64,
    pub mx: f64,
    pub my: f64,
    pub mz: f64,
    pub e: f64,
}

/// hllc flux across a face. axis: 0 normal +x, 1 normal +y, 2 normal +z.
pub fn hllc_flux3(gamma: f64, left: FacePrim3, right: FacePrim3, axis: usize) -> Flux5 {
    let rl = left.rho.max(1e-12);
    let rr = right.rho.max(1e-12);
    let pl = left.p.max(1e-12);
    let pr = right.p.max(1e-12);
    // rotate into the face frame: (un, t1, t2) about the chosen normal axis.
    let (n, t1, t2) = match axis {
        0 => (left.u, left.v, left.w),
        1 => (left.v, left.u, left.w),
        _ => (left.w, left.u, left.v),
    };
    let (nr, t1r, t2r) = match axis {
        0 => (right.u, right.v, right.w),
        1 => (right.v, right.u, right.w),
        _ => (right.w, right.u, right.v),
    };
    let (un_l, ut1_l, ut2_l) = (n, t1, t2);
    let (un_r, ut1_r, ut2_r) = (nr, t1r, t2r);
    let al = (gamma * pl / rl).sqrt();
    let ar = (gamma * pr / rr).sqrt();
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

    let phys = |r: f64, un: f64, t1: f64, t2: f64, p: f64, e_tot: f64| -> [f64; 5] {
        let en = e_tot + p;
        [r * un, r * un * un + p, r * un * t1, r * un * t2, un * en]
    };
    let et_l = pl / (rl * (gamma - 1.0)) + 0.5 * (un_l * un_l + ut1_l * ut1_l + ut2_l * ut2_l);
    let et_r = pr / (rr * (gamma - 1.0)) + 0.5 * (un_r * un_r + ut1_r * ut1_r + ut2_r * ut2_r);
    let (e_l, e_r) = (rl * et_l, rr * et_r);
    let f_l = phys(rl, un_l, ut1_l, ut2_l, pl, e_l);
    let f_r = phys(rr, un_r, ut1_r, ut2_r, pr, e_r);

    let ustar = |r: f64, un: f64, t1: f64, t2: f64, p: f64, e: f64, s_side: f64| -> [f64; 5] {
        let fac = r * (s_side - un) / (s_side - sstar);
        [
            fac,
            fac * sstar,
            fac * t1,
            fac * t2,
            fac * (e / r + (sstar - un) * (sstar + p / (r * (s_side - un)))),
        ]
    };
    let u_l = [rl, rl * un_l, rl * ut1_l, rl * ut2_l, e_l];
    let u_r = [rr, rr * un_r, rr * ut1_r, rr * ut2_r, e_r];
    let us_l = ustar(rl, un_l, ut1_l, ut2_l, pl, e_l, sl);
    let us_r = ustar(rr, un_r, ut1_r, ut2_r, pr, e_r, sr);
    let fstar_l = [
        f_l[0] + sl * (us_l[0] - u_l[0]),
        f_l[1] + sl * (us_l[1] - u_l[1]),
        f_l[2] + sl * (us_l[2] - u_l[2]),
        f_l[3] + sl * (us_l[3] - u_l[3]),
        f_l[4] + sl * (us_l[4] - u_l[4]),
    ];
    let fstar_r = [
        f_r[0] + sr * (us_r[0] - u_r[0]),
        f_r[1] + sr * (us_r[1] - u_r[1]),
        f_r[2] + sr * (us_r[2] - u_r[2]),
        f_r[3] + sr * (us_r[3] - u_r[3]),
        f_r[4] + sr * (us_r[4] - u_r[4]),
    ];
    let f = match () {
        _ if sl > 0.0 => f_l,
        _ if sr < 0.0 => f_r,
        _ if sstar > 0.0 => fstar_l,
        _ => fstar_r,
    };
    let (mx, my, mz) = match axis {
        0 => (f[1], f[2], f[3]),
        1 => (f[2], f[1], f[3]),
        _ => (f[2], f[3], f[1]),
    };
    Flux5 {
        mass: f[0],
        mx,
        my,
        mz,
        e: f[4],
    }
}

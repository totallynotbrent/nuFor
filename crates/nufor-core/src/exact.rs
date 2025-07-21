//! exact solution of the 1D ideal-gas Riemann problem (Toro, ch. 4).
//!
//! used to verify the numerical scheme on shock-tube cases like sod and lax:
//! given left and right primitive states it solves for the star-region pressure
//! by bisection, then samples the exact density/velocity/pressure anywhere.

/// a primitive ideal-gas state (density, velocity, pressure).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrimState {
    /// density.
    pub rho: f64,
    /// velocity.
    pub u: f64,
    /// pressure.
    pub p: f64,
}

/// solution of a 1D riemann problem at a point in space-time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExactSolution {
    /// star-region (contact) pressure.
    pub p_star: f64,
    /// contact velocity.
    pub u_star: f64,
    /// exact state at the sampled (x, t).
    pub state: PrimState,
}

/// the velocity-jump function for one side used to solve the star pressure.
fn velocity_jump(p: f64, rho: f64, a: f64, p0: f64, g: f64) -> f64 {
    let gm1 = g - 1.0;
    let gp1 = g + 1.0;
    if p > p0 {
        // shock: the normal-shock relation for the velocity jump.
        (p - p0) * (2.0 / (gp1 * rho * (p + gm1 / gp1 * p0) + 1e-300)).sqrt()
    } else {
        // rarefaction: isentropic fan relations.
        2.0 * a / gm1 * ((p / p0).powf(gm1 / (2.0 * g)) - 1.0)
    }
}

/// exact solution of a riemann problem: solves the star pressure, then samples.
pub fn riemann(left: PrimState, right: PrimState, gamma: f64, x: f64, t: f64) -> ExactSolution {
    let a_l = (gamma * left.p / left.rho).sqrt();
    let a_r = (gamma * right.p / right.rho).sqrt();

    // f(p) = fl(p) + fr(p) + (u_r - u_l); the root is the star pressure.
    let f = |p: f64| {
        velocity_jump(p, left.rho, a_l, left.p, gamma)
            + velocity_jump(p, right.rho, a_r, right.p, gamma)
            + (right.u - left.u)
    };

    // bisection on a bracket that always contains the star pressure.
    let (mut lo, mut hi) = (1e-12, 1e6);
    let flo = f(lo);
    let fhi0 = f(hi);
    let mut guard = 0;
    while fhi0 * flo > 0.0 && guard < 300 {
        hi *= 4.0;
        guard += 1;
    }
    let mut p_star = 0.5 * (lo + hi);
    for _ in 0..400 {
        let fm = f(p_star);
        if fm == 0.0 {
            break;
        }
        if flo * fm < 0.0 {
            hi = p_star;
        } else {
            lo = p_star;
        }
        p_star = 0.5 * (lo + hi);
    }

    // contact velocity follows directly from the two sides' jumps.
    let fl = velocity_jump(p_star, left.rho, a_l, left.p, gamma);
    let fr = velocity_jump(p_star, right.rho, a_r, right.p, gamma);
    let u_star = 0.5 * (left.u + right.u) + 0.5 * (fr - fl);

    let state = sample(left, right, gamma, p_star, u_star, x, t);
    ExactSolution {
        p_star,
        u_star,
        state,
    }
}

/// star-region density for one side: isentropic for a rarefaction, the shock ratio otherwise.
fn star_density(p_star: f64, p0: f64, rho0: f64, g: f64) -> f64 {
    let gm1 = g - 1.0;
    let gp1 = g + 1.0;
    if p_star > p0 {
        rho0 * ((p_star / p0 + gm1 / gp1) / (gm1 / gp1 * p_star / p0 + 1.0))
    } else {
        rho0 * (p_star / p0).powf(1.0 / g)
    }
}

/// sample the exact solution at (x, t) through the wave pattern.
fn sample(
    left: PrimState,
    right: PrimState,
    g: f64,
    p_star: f64,
    u_star: f64,
    x: f64,
    t: f64,
) -> PrimState {
    let gm1 = g - 1.0;
    let gp1 = g + 1.0;
    let a_l = (g * left.p / left.rho).sqrt();
    let a_r = (g * right.p / right.rho).sqrt();
    let xi = x / t;
    let rho_l_star = star_density(p_star, left.p, left.rho, g);
    let rho_r_star = star_density(p_star, right.p, right.rho, g);

    if xi < u_star {
        // left of the contact: a shock or a rarefaction fan leads in from the left state.
        if p_star > left.p {
            let s = left.u - a_l * ((gp1 / (2.0 * g)) * p_star / left.p + gm1 / (2.0 * g)).sqrt();
            if xi < s {
                left
            } else {
                PrimState {
                    rho: rho_l_star,
                    u: u_star,
                    p: p_star,
                }
            }
        } else {
            let a_l_star = a_l * (p_star / left.p).powf(gm1 / (2.0 * g));
            let head = left.u - a_l;
            let tail = u_star - a_l_star;
            if xi < head {
                left
            } else if xi > tail {
                PrimState {
                    rho: rho_l_star,
                    u: u_star,
                    p: p_star,
                }
            } else {
                // inside the fan the isentrope gives the state.
                let u = 2.0 / gp1 * (xi + a_l + gm1 / 2.0 * left.u);
                let a = a_l + gm1 / 2.0 * (left.u - u);
                let p = left.p * (a / a_l).powf(2.0 * g / gm1);
                let rho = left.rho * (p / left.p).powf(1.0 / g);
                PrimState { rho, u, p }
            }
        }
    } else {
        // right of the contact: a shock or a rarefaction fan leads in from the right state.
        if p_star > right.p {
            let s = right.u + a_r * ((gp1 / (2.0 * g)) * p_star / right.p + gm1 / (2.0 * g)).sqrt();
            if xi > s {
                right
            } else {
                PrimState {
                    rho: rho_r_star,
                    u: u_star,
                    p: p_star,
                }
            }
        } else {
            let a_r_star = a_r * (p_star / right.p).powf(gm1 / (2.0 * g));
            let head = right.u + a_r;
            let tail = u_star + a_r_star;
            if xi > head {
                right
            } else if xi < tail {
                PrimState {
                    rho: rho_r_star,
                    u: u_star,
                    p: p_star,
                }
            } else {
                let u = 2.0 / gp1 * (xi - a_r + gm1 / 2.0 * right.u);
                let a = a_r - gm1 / 2.0 * (u - right.u);
                let p = right.p * (a / a_r).powf(2.0 * g / gm1);
                let rho = right.rho * (p / right.p).powf(1.0 / g);
                PrimState { rho, u, p }
            }
        }
    }
}

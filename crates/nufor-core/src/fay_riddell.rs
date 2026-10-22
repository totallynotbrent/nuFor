//! fay-riddell stagnation-point heating: the closed form used as the
//! stagnation-heating validation anchor. perfect-gas implementation with sutherland
//! transport and lewis number one (no dissociation enthalpy term);
//! the equilibrium-air extension rides the species march when it lands.

/// sutherland viscosity for n2 (s = 111 k, mu_ref 1.7168e-5 pa s at 273.15 k).
pub fn mu_n2(t: f64) -> f64 {
    const S: f64 = 111.0;
    const MU_REF: f64 = 1.7168e-5;
    const T_REF: f64 = 273.15;
    MU_REF * (t / T_REF).powf(1.5) * (T_REF + S) / (t + S)
}

/// the stagnation state behind a normal shock, perfect gas:
/// (p_0, t_0, h_0) at the nose from the freestream mach and statics.
pub fn stagnation_state(
    mach: f64,
    p_inf: f64,
    t_inf: f64,
    gas_const: f64,
    gamma: f64,
) -> (f64, f64, f64) {
    let m2 = mach * mach;
    // total pressure and temperature via the rayleigh/pitot chain:
    // normal-shock total pressure ratio (nasa sp-113 form), then
    // isentropic total state of the post-shock flow.
    let p2_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (m2 - 1.0);
    let m2_post = (1.0 + (gamma - 1.0) / 2.0 * m2) / (gamma * m2 - (gamma - 1.0) / 2.0);
    let m_post = m2_post.sqrt();
    let p02_p2 = (1.0 + (gamma - 1.0) / 2.0 * m_post * m_post).powf(gamma / (gamma - 1.0));
    let p0 = p_inf * p2_p1 * p02_p2;
    let t0 = t_inf * (1.0 + (gamma - 1.0) / 2.0 * m2);
    let cp = gamma * gas_const / (gamma - 1.0);
    let h0 = cp * t0;
    (p0, t0, h0)
}

/// the stagnation-point edge state: post-normal-shock statics used by
/// the boundary-layer edge (p_e = p_0 for the stagnation point, t_e =
/// the post-shock static temperature).
pub fn post_shock_statics(
    mach: f64,
    p_inf: f64,
    t_inf: f64,
    gas_const: f64,
    gamma: f64,
) -> (f64, f64, f64) {
    let m2 = mach * mach;
    let p2_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (m2 - 1.0);
    let rho2_rho1 = (gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0);
    let p_e = p_inf * p2_p1;
    let t_e = t_inf * p2_p1 / rho2_rho1;
    let rho_e = p_e / (gas_const * t_e);
    (p_e, t_e, rho_e)
}

/// the edge-and-wall inputs the closed form needs, all si.
pub struct FayRiddellEdge {
    pub p_inf: f64,
    pub p_e: f64,
    pub rho_e: f64,
    pub t_e: f64,
    pub h0e: f64,
    pub t_w: f64,
    pub gas_const: f64,
    pub gamma: f64,
    pub pr: f64,
    pub r_n: f64,
}

/// fay-riddell stagnation-point cold-wall heat flux, w/m^2:
/// q = 0.763 pr^-0.6 (rho_e mu_e beta)^0.5 (h_0e - h_w), with the
/// newtonian stagnation velocity gradient
/// beta = sqrt(2 (p_e - p_inf) / rho_e) / r_n.
/// lewis number is one, so the dissociation-enthalpy bracket is absent.
pub fn fay_riddell_q(s: &FayRiddellEdge) -> f64 {
    let mu_e = mu_n2(s.t_e);
    let beta = (2.0 * (s.p_e - s.p_inf) / s.rho_e).sqrt() / s.r_n;
    let cp = s.gamma * s.gas_const / (s.gamma - 1.0);
    let h_w = cp * s.t_w;
    0.763 * s.pr.powf(-0.6) * (s.rho_e * mu_e * beta).sqrt() * (s.h0e - h_w)
}

/// planar (2d cylinder) stagnation-line cold-wall heat flux, w/m^2:
/// q = 0.57 pr^-0.6 sqrt(rho_e mu_e beta) (h0e - h_w), with the same
/// newtonian gradient beta = sqrt(2 (p_e - p_inf)/rho_e)/r_n as the
/// sphere. sources: beckwith & cohen eq. (16) p. 409 (ntrs 19660082088);
/// tauber review eq. (36a) citing cohen nasa tr r-118. the ratio to
/// the sphere form is 0.57/0.763 = 0.747; reshotko's cold-wall
/// 1/sqrt(2) = 0.707 alternative is the systematic lower bound (see
/// plans/fay-riddell-2d-research.md section 4 for the disagreement).
pub fn fay_riddell_q_2d(s: &FayRiddellEdge) -> f64 {
    let mu_e = mu_n2(s.t_e);
    let beta = (2.0 * (s.p_e - s.p_inf) / s.rho_e).sqrt() / s.r_n;
    let cp = s.gamma * s.gas_const / (s.gamma - 1.0);
    let h_w = cp * s.t_w;
    0.57 * s.pr.powf(-0.6) * (s.rho_e * mu_e * beta).sqrt() * (s.h0e - h_w)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the rankine-hugoniot tables at m=20, gamma=1.4: p2/p1 = 465.5,
    /// rho2/rho1 = 5.926, t2/t1 = 78.55. the stagnation total state
    /// follows the nasa chain; both are textbook-locked.
    #[test]
    fn normal_shock_matches_tables() {
        let (p_e, t_e, rho_e) = post_shock_statics(20.0, 1.0, 220.0, 296.8, 1.4);
        let rho_inf = 1.0 / (296.8 * 220.0);
        let p_ratio = p_e / 1.0;
        let rho_ratio = rho_e / rho_inf;
        let t_ratio = t_e / 220.0;
        assert!((p_ratio - 465.5).abs() / 465.5 < 5e-3, "p ratio {p_ratio}");
        assert!(
            (rho_ratio - 5.926).abs() / 5.926 < 5e-3,
            "rho ratio {rho_ratio}"
        );
        assert!((t_ratio - 78.55).abs() / 78.55 < 5e-3, "t ratio {t_ratio}");
    }

    /// scaling sanity: q grows with the square root of nose radius
    /// inverse, density, and edge enthalpy; units land in w/m^2.
    #[test]
    fn scales_and_units() {
        let p_inf = 0.89;
        let t_inf = 220.0;
        let r = 296.8;
        let (p_e, t_e, rho_e) = post_shock_statics(20.0, p_inf, t_inf, r, 1.4);
        let (_, _, h0e) = stagnation_state(20.0, p_inf, t_inf, r, 1.4);
        let mk = |r_n: f64| FayRiddellEdge {
            p_inf,
            p_e,
            rho_e,
            t_e,
            h0e,
            t_w: 1000.0,
            gas_const: r,
            gamma: 1.4,
            pr: 0.71,
            r_n,
        };
        let q1 = fay_riddell_q(&mk(1.0));
        let q2 = fay_riddell_q(&mk(4.0));
        // 4x radius -> 2x less heating
        let ratio = q1 / q2;
        assert!((ratio - 2.0).abs() < 1e-9, "r_n scaling {ratio}");
        // magnitude: this classic case sits near 1e5 w/m^2 cold-wall
        assert!(q1 > 1e4 && q1 < 1e6, "q1 {q1} w/m^2 out of envelope");
    }

    /// the planar (cylinder) gate from plans/fay-riddell-2d-research.md
    /// section 6(a): m 20 perfect n2, r_n 0.01 m, tw 300 k -> q_2d =
    /// 89.6 w/cm^2 on the 0.57 form and 84.8 on the 1/sqrt(2) bound,
    /// with the ratio to the axisymmetric form inside [0.70, 0.75].
    #[test]
    fn planar_gate_matches_research_doc() {
        let m = 20.0f64;
        let (p_e, t_e, rho_e) = post_shock_statics(m, 0.89, 220.0, 296.8, 1.4);
        let cp = 1.4 * 296.8 / (1.4 - 1.0);
        // the stagnation enthalpy is the freestream total: cp t0 with
        // t0 = t_inf (1 + (gamma-1)/2 m^2), not the post-shock state.
        let t0 = 220.0 * (1.0 + 0.5 * (1.4 - 1.0) * m * m);
        let h0e = cp * t0;
        let s = FayRiddellEdge {
            p_inf: 0.89,
            p_e,
            rho_e,
            t_e,
            h0e,
            r_n: 0.01,
            t_w: 300.0,
            pr: 0.71,
            gamma: 1.4,
            gas_const: 296.8,
        };
        let q_axi = fay_riddell_q(&s);
        let q_2d = fay_riddell_q_2d(&s);
        let ratio = q_2d / q_axi;
        assert!(
            (ratio - 0.747).abs() < 5e-3,
            "2d/axi ratio {ratio}, expected 0.57/0.763 = 0.747"
        );
        assert!(
            ratio > 0.70 && ratio < 0.75,
            "ratio {ratio} outside the published spread [0.70, 0.75]"
        );
        // the doc's worked example is w/cm^2; our case has r_n 0.01 m and
        // p_inf 0.89 pa like the gate case, so q in w/cm^2 lands near
        // the doc's 89.6 within the perfect-gas edge-state differences.
        let q_2d_wcm2 = q_2d / 1e4;
        assert!(
            q_2d_wcm2 > 80.0 && q_2d_wcm2 < 100.0,
            "q_2d {q_2d_wcm2} w/cm^2 outside the doc's worked-example band"
        );
    }
}

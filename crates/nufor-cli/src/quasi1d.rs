//! the isentropic area-mach solution of a chamber-throat contour.
//!
//! a choked chamber started from a uniform fill rings for a long time
//! before the standing acoustic wave settles. seeding the case with the
//! quasi-1d steady state lets the march be checked for HOLDING that
//! state, which is the physics the case is about.

/// one station of the quasi-1d solution.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct Station {
    /// the station's x coordinate.
    pub x: f64,
    /// the local radius of the flow tube.
    pub r: f64,
    /// the local mach number.
    pub mach: f64,
    /// the local pressure.
    pub p: f64,
    /// the local temperature.
    pub t: f64,
    /// the local density.
    pub rho: f64,
    /// the local axial velocity.
    pub u: f64,
}

/// the mach number at a given area ratio on the requested branch.
pub fn mach_at_area_ratio(ar: f64, gamma: f64, supersonic: bool) -> f64 {
    let ar = ar.max(1.0);
    let gm = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    let (mut lo, mut hi) = if supersonic { (1.0, 40.0) } else { (1e-6, 1.0) };
    for _ in 0..100 {
        let m = 0.5 * (lo + hi);
        let f = (1.0 / m) * ((2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * m * m)).powf(gm);
        // the area ratio falls with mach on the subsonic branch and
        // rises with mach on the supersonic one.
        if supersonic {
            if f > ar {
                hi = m;
            } else {
                lo = m;
            }
        } else if f > ar {
            lo = m;
        } else {
            hi = m;
        }
    }
    0.5 * (lo + hi)
}

/// the isentropic state at a given mach number.
pub fn station_state(
    mach: f64,
    gamma: f64,
    p0: f64,
    t0: f64,
    gas_const: f64,
) -> (f64, f64, f64, f64) {
    let a = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    let p = p0 * a.powf(-gamma / (gamma - 1.0));
    let t = t0 / a;
    let rho = p / (gas_const * t);
    let u = mach * (gamma * gas_const * t).sqrt();
    (rho, u, p, t)
}

/// the quasi-1d solution along a contour, choked at the smallest area.
///
/// the contour runs from the chamber to the throat; the smallest area
/// carries mach one, and the upstream stations take the subsonic
/// branch of the area-mach relation.
pub fn quasi1d_solve(
    contour: &[(f64, f64)],
    gamma: f64,
    p0: f64,
    t0: f64,
    gas_const: f64,
) -> Vec<Station> {
    if contour.len() < 2 {
        return Vec::new();
    }
    let areas: Vec<f64> = contour
        .iter()
        .map(|(_, r)| std::f64::consts::PI * r * r)
        .collect();
    let a_min = areas.iter().cloned().fold(f64::INFINITY, f64::min);
    // the throat station: the smallest area, the last one on a
    // converging contour.
    let mut i_throat = 0usize;
    for (k, a) in areas.iter().enumerate() {
        if *a < areas[i_throat] {
            i_throat = k;
        }
    }
    let mut out = Vec::with_capacity(contour.len());
    for (k, (x, r)) in contour.iter().enumerate() {
        // upstream of the throat the flow is subsonic, downstream
        // supersonic; at the throat itself the flow is sonic.
        let supersonic = k > i_throat;
        let mach = if k == i_throat {
            1.0
        } else {
            mach_at_area_ratio(areas[k] / a_min, gamma, supersonic)
        };
        let (rho, u, p, t) = station_state(mach, gamma, p0, t0, gas_const);
        out.push(Station {
            x: *x,
            r: *r,
            mach,
            p,
            t,
            rho,
            u,
        });
    }
    out
}

/// the contour's radius at a given x, linearly interpolated.
/// the slope dr/dx of the flow-tube wall at `x`, from the contour's
/// two bracketing vertices. the axisymmetric seed needs it: the
/// streamlines of a converging duct converge inward, so the radial
/// velocity is u * (r/R) * dR/dx and a seed with v = 0 violates
/// continuity by the factor R'.
pub fn contour_slope(contour: &[(f64, f64)], x: f64) -> f64 {
    if contour.len() < 2 {
        return 0.0;
    }
    let last = contour[contour.len() - 1];
    if x <= contour[0].0 {
        let (x0, r0) = contour[0];
        let (x1, r1) = contour[1];
        return if (x1 - x0).abs() > 1e-15 {
            (r1 - r0) / (x1 - x0)
        } else {
            0.0
        };
    }
    if x >= last.0 {
        let (x0, r0) = contour[contour.len() - 2];
        let (x1, r1) = last;
        return if (x1 - x0).abs() > 1e-15 {
            (r1 - r0) / (x1 - x0)
        } else {
            0.0
        };
    }
    for k in 1..contour.len() {
        let (x0, r0) = contour[k - 1];
        let (x1, r1) = contour[k];
        if x <= x1 {
            return if (x1 - x0).abs() > 1e-15 {
                (r1 - r0) / (x1 - x0)
            } else {
                0.0
            };
        }
    }
    0.0
}

pub fn contour_radius(contour: &[(f64, f64)], x: f64) -> f64 {
    if contour.is_empty() {
        return 0.0;
    }
    if x <= contour[0].0 {
        return contour[0].1;
    }
    let last = contour[contour.len() - 1];
    if x >= last.0 {
        return last.1;
    }
    for k in 1..contour.len() {
        let (x0, r0) = contour[k - 1];
        let (x1, r1) = contour[k];
        if x <= x1 {
            let f = (x - x0) / (x1 - x0);
            return r0 + f * (r1 - r0);
        }
    }
    last.1
}

/// the contour the case's polygon describes: the radius at each x.
///
/// the polygon's boundary is walked edge by edge and the gas side is
/// the lower branch, so a chamber polygon that closes over the top
/// still yields the chamber's own wall contour.
pub fn contour_from_polygon(verts: &[(f64, f64)], samples: usize) -> Vec<(f64, f64)> {
    let chain = contour_chain(verts);
    if chain.len() < 2 {
        return Vec::new();
    }
    let xmin = chain[0].0;
    let xmax = chain[chain.len() - 1].0;
    let mut out = Vec::with_capacity(samples);
    for k in 0..samples {
        let x = xmin + (xmax - xmin) * (k as f64) / ((samples - 1).max(1) as f64);
        let r = contour_radius(&chain, x);
        out.push((x, r));
    }
    out
}

/// the wall contour of the polygon: its longest chain of vertices with
/// a rising x. a chamber polygon closes with a short chain over the
/// top, so the long rising chain is the wall the gas runs against.
pub fn contour_chain(verts: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let n = verts.len();
    if n < 2 {
        return Vec::new();
    }
    // the wall contour of the polygon: the longest run of vertices
    // with a rising x whose radius varies smoothly. a polygon that
    // closes over the top of the chamber ends its wall chain with a
    // jump to the outer boundary, and that jump must not enter the
    // chain: the closing vertex sits at the box corner, and seeding
    // the duct with its radius blows the exit area up by a factor of
    // twenty.
    let mut best: Vec<(f64, f64)> = Vec::new();
    for s in 0..n {
        let mut chain = vec![verts[s]];
        for k in 1..n {
            let cur = verts[(s + k) % n];
            let prev = chain[chain.len() - 1];
            if cur.0 <= prev.0 {
                break;
            }
            // a radius jump of more than a factor of two between
            // neighbours is the closing edge, not the wall.
            let ratio = cur.1 / prev.1.max(1e-12);
            if ratio > 2.0 || ratio < 0.5 {
                break;
            }
            chain.push(cur);
        }
        if chain.len() > best.len() {
            best = chain;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the wall contour of a chamber polygon must not pick up the
    /// closing vertex at the box corner: that vertex sits at the outer
    /// radius, and seeding the duct with it blows the exit area up by
    /// a factor of twenty.
    #[test]
    fn contour_chain_stops_at_the_closing_edge() {
        let mut poly = vec![(0.00075, 0.05)];
        for i in 1..100 {
            let x = 0.00075 + i as f64 * 0.0015;
            poly.push((x, 0.05 - 0.03 * (x / 0.15)));
        }
        poly.push((0.15, 0.1));
        poly.push((0.0, 0.1));
        let chain = contour_chain(&poly);
        let last = chain[chain.len() - 1];
        assert!(
            last.1 < 0.03,
            "the chain must end on the wall, not the box corner: {last:?}"
        );
        assert!(
            chain.len() > 90,
            "the chain must cover the wall: {}",
            chain.len()
        );
    }

    /// the slope of the wall contour: the axisymmetric seed's radial
    /// velocity comes from it.
    #[test]
    fn contour_slope_follows_the_wall() {
        let contour = vec![(0.0, 0.10), (0.1, 0.05), (0.2, 0.02)];
        let s = contour_slope(&contour, 0.05);
        assert!((s + 0.5).abs() < 1e-9, "slope {s}");
        let s = contour_slope(&contour, 0.15);
        assert!((s + 0.3).abs() < 1e-9, "slope {s}");
    }

    /// the area-mach relation must return mach one at an area ratio of
    /// one on both branches.
    #[test]
    fn area_ratio_one_is_sonic_on_both_branches() {
        let g = 1.2;
        assert!((mach_at_area_ratio(1.0, g, false) - 1.0).abs() < 1e-6);
        assert!((mach_at_area_ratio(1.0, g, true) - 1.0).abs() < 1e-6);
    }

    /// the area-mach relation must be monotone: a bigger area carries
    /// a smaller mach on the subsonic branch and a bigger one on the
    /// supersonic branch.
    #[test]
    fn area_mach_is_monotone_on_each_branch() {
        let g = 1.2;
        let mut prev = mach_at_area_ratio(1.0, g, false);
        for ar in [1.5, 2.0, 4.0, 8.0, 16.0] {
            let m = mach_at_area_ratio(ar, g, false);
            assert!(m < prev, "subsonic mach must fall as the area grows");
            prev = m;
        }
        let mut prev = mach_at_area_ratio(1.0, g, true);
        for ar in [1.5, 2.0, 4.0, 8.0, 16.0] {
            let m = mach_at_area_ratio(ar, g, true);
            assert!(m > prev, "supersonic mach must rise as the area grows");
            prev = m;
        }
    }

    /// the isentropic relations must recover the stagnation state as
    /// the mach goes to zero.
    #[test]
    fn station_state_recovers_the_stagnation_state() {
        let (rho, u, p, t) = station_state(1e-4, 1.2, 7.0e6, 3676.53, 527.2);
        assert!((p - 7.0e6).abs() / 7.0e6 < 1e-6);
        assert!((t - 3676.53).abs() / 3676.53 < 1e-6);
        assert!(u < 1.0, "the stagnation velocity is zero");
        assert!((rho - 7.0e6 / (527.2 * 3676.53)).abs() < 1e-3);
    }

    /// a converging contour must choke at its smallest area and carry
    /// a subsonic branch upstream.
    #[test]
    fn converging_contour_chokes_at_its_minimum_area() {
        let contour: Vec<(f64, f64)> = (0..20)
            .map(|k| {
                let x = k as f64 * 0.01;
                (x, 0.1 - 0.004 * k as f64)
            })
            .collect();
        let st = quasi1d_solve(&contour, 1.2, 7.0e6, 3676.53, 527.2);
        assert_eq!(st.len(), contour.len());
        let i_throat = st
            .iter()
            .enumerate()
            .min_by(|a, b| {
                a.1.r
                    .partial_cmp(&b.1.r)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(k, _)| k)
            .unwrap_or(0);
        assert!(
            (st[i_throat].mach - 1.0).abs() < 1e-9,
            "the smallest area carries mach one"
        );
        assert!(
            st[0].mach < 1.0 && st[0].mach > 0.0,
            "the chamber is subsonic, got {}",
            st[0].mach
        );
        // the pressure must fall from the chamber toward the throat.
        assert!(st[0].p > st[i_throat].p);
    }

    /// the polygon contour must recover the gas-side radius of a
    /// chamber polygon that closes over the top.
    #[test]
    fn polygon_contour_recovers_the_gas_side_radius() {
        // the gas sits under the sloped wall; the polygon closes with a
        // flat top edge and a flat floor.
        let verts = vec![(0.0, 0.1), (1.0, 0.2), (1.0, 0.1), (0.0, 0.1)];
        let c = contour_from_polygon(&verts, 11);
        assert_eq!(c.len(), 11);
        assert!((c[0].1 - 0.1).abs() < 1e-9, "the inlet radius");
        assert!((c[10].1 - 0.2).abs() < 1e-9, "the exit radius");
        let mid = contour_radius(&c, 0.5);
        assert!((mid - 0.15).abs() < 1e-6, "the mid radius, got {mid}");
    }
}

//! the cell closure and the legacy face closure have to agree.
//!
//! cells read the cea table through `eqair_cea::p_a_at`, while the
//! legacy face fits read `eqair::eqair_pressure_at`. a divergence
//! between the two shows up as a pressure jump at every interface
//! they meet, so the two closures are held to a shared tolerance over
//! the region both cover.
//!
//! the measured behaviour over the campaign envelope (specific internal
//! energy from 1.4e5 to 3e7 J/kg, density from 1e-4 to 1 kg/m3) is a
//! median of about 5 % and a ninetieth percentile of about 10 %. above
//! 3e7 J/kg the legacy fit walks away from the cea table and the median
//! gap rises past 20 %, which the gate records rather than hides.
//!
//! this lives in the crate because the legacy closure is private to
//! it, so an example or an integration test cannot reach it.

#[cfg(test)]
mod tests {
    use crate::eqair;
    use crate::eqair_cea;

    /// the closure the cells read.
    fn cell_pressure(rho: f64, e: f64) -> f64 {
        eqair_cea::p_a_at(rho, e).0
    }

    /// the closure the legacy face fits read.
    fn face_pressure(rho: f64, e: f64) -> f64 {
        eqair::eqair_pressure_at(rho, e)
    }

    /// relative gap between the two closures at one state.
    fn gap(rho: f64, e: f64) -> Option<f64> {
        let p_cell = cell_pressure(rho, e);
        let p_face = face_pressure(rho, e);
        if !(p_cell.is_finite() && p_face.is_finite()) || p_cell <= 0.0 || p_face <= 0.0 {
            return None;
        }
        Some((p_cell - p_face).abs() / p_cell.max(p_face))
    }

    /// sample the campaign envelope on a logarithmic lattice.
    fn campaign_samples(e_lo: f64, e_hi: f64) -> Vec<f64> {
        let mut gaps = Vec::new();
        for i in 0..30 {
            let rho = 1e-4 * 10f64.powf(4.0 * (i as f64) / 29.0);
            for j in 0..24 {
                let e = e_lo * (e_hi / e_lo).powf((j as f64) / 23.0);
                if let Some(g) = gap(rho, e) {
                    gaps.push(g);
                }
            }
        }
        gaps
    }

    fn percentile(sorted: &[f64], q: f64) -> f64 {
        let idx = (((sorted.len() as f64) * q) as usize).min(sorted.len() - 1);
        sorted[idx]
    }

    #[test]
    fn cell_and_face_closures_agree_over_the_campaign_envelope() {
        let mut gaps = campaign_samples(1.43e5, 3.0e7);
        assert!(gaps.len() > 500, "sampled only {} states", gaps.len());
        gaps.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median = percentile(&gaps, 0.5);
        let p90 = percentile(&gaps, 0.9);
        let max = gaps[gaps.len() - 1];

        assert!(median <= 0.08, "median gap rose to {:.2}%", median * 100.0);
        assert!(
            p90 <= 0.15,
            "ninetieth-percentile gap rose to {:.2}%",
            p90 * 100.0
        );
        assert!(max <= 0.30, "worst gap rose to {:.2}%", max * 100.0);
    }

    #[test]
    fn the_hot_corner_divergence_stays_inside_its_documented_band() {
        // above 3e7 J/kg the legacy fit is known to walk away from the
        // cea table. this gate keeps that departure bounded so a change
        // in either closure shows up here instead of in a campaign run.
        let mut gaps = campaign_samples(3.0e7, 7.0e7);
        assert!(!gaps.is_empty(), "no hot-corner states sampled");
        gaps.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median = percentile(&gaps, 0.5);
        assert!(
            median <= 0.35,
            "hot-corner median gap rose to {:.2}%",
            median * 100.0
        );
    }
}

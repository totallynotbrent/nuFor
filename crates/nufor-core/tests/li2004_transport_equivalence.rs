//! the li-2004 lane's transport numbers must not move.
//!
//! the shared implementation in `transport_general.rs` replaced the
//! formulas that lived inside `combustion_transport.rs`. the two are
//! numerically different expressions of the same algebra, so this
//! module compares them over a grid covering the temperatures and
//! pressures the hydrogen cases use, and demands agreement to the
//! last bits of the doubles.
//!
//! the reference implementation below is copied from the module as it
//! stood before the refactor. it is intentionally repetitive: it is
//! the yardstick, not new code.

use nufor_core::combustion_thermo::{M, NS};
use nufor_core::combustion_transport::{
    binary_diffusivity, effective_diffusivity, mixture_conductivity, mixture_viscosity,
    species_conductivity, species_viscosity,
};

/// the pre-refactor collision-integral fits.
fn omega_mu_ref(tstar: f64) -> f64 {
    let t = tstar.max(0.3);
    1.16145 / t.powf(0.14874)
        + 0.52487 * (-0.7732 * t).exp()
        + 2.16178 * (-2.43787 * t).exp()
        + 6.7353e-2 / (t * t).max(1e-12) * (-3.0 * t).exp().min(1.0)
}

/// the pre-refactor diffusion collision-integral fit.
fn omega_d_ref(tstar: f64) -> f64 {
    let t = tstar.max(0.3);
    1.06036 / t.powf(0.15610)
        + 0.19300 * (-0.47635 * t).exp()
        + 1.03587 * (-1.52996 * t).exp()
        + 1.76474 * (-3.89411 * t).exp()
}

/// the pre-refactor viscosity.
fn species_viscosity_ref(i: usize, t: f64) -> f64 {
    let eps: [f64; NS] = [38.0, 107.4, 572.4, 80.0, 145.0, 80.0, 107.4, 107.4, 97.53];
    let sig: [f64; NS] = [2.92, 3.458, 2.605, 2.75, 2.05, 2.75, 3.458, 3.458, 3.621];
    let tstar = t / eps[i];
    let m_g: f64 = M[i] * 1.0e3;
    2.669e-6 * (m_g * t).sqrt() / (sig[i] * sig[i]) / omega_mu_ref(tstar)
}

/// the pre-refactor conductivity.
fn species_conductivity_ref(i: usize, t: f64, cp_i: f64) -> f64 {
    let mu = species_viscosity_ref(i, t);
    mu * (cp_i + 1.25 * 8.314462 / M[i])
}

/// the pre-refactor binary diffusivity.
fn binary_diffusivity_ref(i: usize, j: usize, t: f64, p: f64) -> f64 {
    let eps: [f64; NS] = [38.0, 107.4, 572.4, 80.0, 145.0, 80.0, 107.4, 107.4, 97.53];
    let sig: [f64; NS] = [2.92, 3.458, 2.605, 2.75, 2.05, 2.75, 3.458, 3.458, 3.621];
    let eps_ij: f64 = (eps[i] * eps[j]).sqrt();
    let sig_ij: f64 = 0.5 * (sig[i] + sig[j]);
    let tstar = t / eps_ij;
    let m_g: [f64; 2] = [M[i] * 1.0e3, M[j] * 1.0e3];
    let mred: f64 = 1.0 / m_g[0] + 1.0 / m_g[1];
    let p_atm = (p / 101_325.0).max(1e-6);
    0.0018583 * t.powf(1.5) * mred.sqrt() / (p_atm * sig_ij * sig_ij * omega_d_ref(tstar)) * 1.0e-4
}

/// the pre-refactor wilke viscosity.
fn mixture_viscosity_ref(y: &[f64; NS], t: f64) -> f64 {
    let inv_m: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if inv_m <= 0.0 {
        return species_viscosity_ref(nufor_core::combustion_thermo::IN2, t);
    }
    let mut x = [0.0f64; NS];
    for i in 0..NS {
        x[i] = (y[i] / M[i]) / inv_m;
    }
    let mut mu_i = [0.0f64; NS];
    for (i, entry) in mu_i.iter_mut().enumerate() {
        *entry = species_viscosity_ref(i, t);
    }
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * mu_i[i];
        let mut phi_sum = 0.0;
        for j in 0..NS {
            if x[j] < 1e-12 {
                continue;
            }
            let mr: f64 = M[i] / M[j];
            let phi = (1.0 + (mu_i[i] / mu_i[j]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        return mu_i[nufor_core::combustion_thermo::IN2];
    }
    num / den
}

/// the pre-refactor wilke conductivity.
fn mixture_conductivity_ref(y: &[f64; NS], t: f64, cps: &[f64; NS]) -> f64 {
    let inv_m: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if inv_m <= 0.0 {
        return 0.0;
    }
    let mut x = [0.0f64; NS];
    for i in 0..NS {
        x[i] = (y[i] / M[i]) / inv_m;
    }
    let mut k_i = [0.0f64; NS];
    for i in 0..NS {
        k_i[i] = species_conductivity_ref(i, t, cps[i]);
    }
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * k_i[i];
        let mut phi_sum = 0.0;
        for j in 0..NS {
            if x[j] < 1e-12 {
                continue;
            }
            let mr: f64 = M[i] / M[j];
            let phi = (1.0 + (k_i[i] / k_i[j]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        return 0.0;
    }
    num / den
}

/// the pre-refactor effective diffusivity.
fn effective_diffusivity_ref(y: &[f64; NS], t: f64, p: f64) -> [f64; NS] {
    let mut x = [0.0f64; NS];
    let mut inv_m = 0.0;
    for i in 0..NS {
        inv_m += y[i] / M[i];
    }
    if inv_m > 1e-30 {
        for i in 0..NS {
            x[i] = (y[i] / M[i]) / inv_m;
        }
    }
    let mut d = [0.0f64; NS];
    for (i, entry) in d.iter_mut().enumerate() {
        let mut inv = 0.0f64;
        for (j, &xj) in x.iter().enumerate() {
            if j != i && xj > 1e-30 {
                inv += xj / binary_diffusivity_ref(i, j, t, p);
            }
        }
        *entry = if inv > 1e-30 { 1.0 / inv } else { 0.0 };
    }
    d
}

/// the compositions the grid walks: stoichiometric mixtures, rich
/// premixes, dilute quenches, and inert-dominated blends.
fn compositions() -> Vec<[f64; NS]> {
    use nufor_core::combustion_thermo::{IH2, IH2O, IN2, IO2};
    let mut out = Vec::new();
    let stoich = [(IH2, 2.0f64), (IO2, 1.0), (IN2, 3.76)];
    let mut y = [0.0f64; NS];
    for &(i, v) in stoich.iter() {
        y[i] = v;
    }
    let sum: f64 = y.iter().sum();
    for v in y.iter_mut() {
        *v /= sum;
    }
    out.push(y);

    let mut y = [0.0f64; NS];
    y[IH2] = 0.5;
    y[IO2] = 0.2;
    y[IN2] = 0.2;
    y[IH2O] = 0.1;
    out.push(y);

    let mut y = [0.0f64; NS];
    y[IH2] = 0.05;
    y[IO2] = 0.95;
    out.push(y);

    let mut y = [0.0f64; NS];
    y[IH2O] = 0.7;
    y[IN2] = 0.3;
    out.push(y);

    out
}

/// every pure-species viscosity agrees bit for bit.
#[test]
fn pure_viscosity_is_bit_identical() {
    for i in 0..NS {
        for &t in &[200.0f64, 300.0, 900.0, 1800.0, 3500.0] {
            let new = species_viscosity(i, t);
            let old = species_viscosity_ref(i, t);
            assert_eq!(
                new.to_bits(),
                old.to_bits(),
                "species {i} at {t} k: {new:e} vs {old:e}"
            );
        }
    }
}

/// every pure-species conductivity agrees bit for bit.
#[test]
fn pure_conductivity_is_bit_identical() {
    for i in 0..NS {
        for &t in &[200.0f64, 300.0, 900.0, 1800.0, 3500.0] {
            let cp = 1000.0 + 100.0 * i as f64;
            let new = species_conductivity(i, t, cp);
            let old = species_conductivity_ref(i, t, cp);
            assert_eq!(
                new.to_bits(),
                old.to_bits(),
                "species {i} at {t} k: {new:e} vs {old:e}"
            );
        }
    }
}

/// every binary diffusivity agrees bit for bit.
#[test]
fn binary_diffusivity_is_bit_identical() {
    for i in 0..NS {
        for j in 0..NS {
            for &t in &[250.0f64, 800.0, 2500.0] {
                for &p in &[101_325.0f64, 506_625.0, 2_026_500.0] {
                    let new = binary_diffusivity(i, j, t, p);
                    let old = binary_diffusivity_ref(i, j, t, p);
                    assert_eq!(
                        new.to_bits(),
                        old.to_bits(),
                        "pair {i}-{j} at {t} k, {p} pa: {new:e} vs {old:e}"
                    );
                }
            }
        }
    }
}

/// every mixture viscosity agrees bit for bit.
#[test]
fn mixture_viscosity_is_bit_identical() {
    for y in compositions() {
        for &t in &[300.0f64, 1200.0, 3000.0] {
            let new = mixture_viscosity(&y, t);
            let old = mixture_viscosity_ref(&y, t);
            assert_eq!(
                new.to_bits(),
                old.to_bits(),
                "mix at {t} k: {new:e} vs {old:e}"
            );
        }
    }
}

/// every mixture conductivity agrees bit for bit.
#[test]
fn mixture_conductivity_is_bit_identical() {
    for y in compositions() {
        let cps: [f64; NS] = std::array::from_fn(|i| 1000.0 + 100.0 * i as f64);
        for &t in &[300.0f64, 1200.0, 3000.0] {
            let new = mixture_conductivity(&y, t, &cps);
            let old = mixture_conductivity_ref(&y, t, &cps);
            assert_eq!(
                new.to_bits(),
                old.to_bits(),
                "mix at {t} k: {new:e} vs {old:e}"
            );
        }
    }
}

/// every effective diffusivity agrees bit for bit.
#[test]
fn effective_diffusivity_is_bit_identical() {
    for y in compositions() {
        for &t in &[300.0f64, 1500.0] {
            for &p in &[101_325.0f64, 1_013_250.0] {
                let new = effective_diffusivity(&y, t, p);
                let old = effective_diffusivity_ref(&y, t, p);
                for i in 0..NS {
                    assert_eq!(
                        new[i].to_bits(),
                        old[i].to_bits(),
                        "species {i} at {t} k: {:e} vs {:e}",
                        new[i],
                        old[i]
                    );
                }
            }
        }
    }
}

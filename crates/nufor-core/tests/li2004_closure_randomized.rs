//! the generic slice closure must agree with the legacy nine-species
//! closure bit for bit on any composition, not only the sampled ones.
//!
//! the chamber migrates its closure calls onto the scheme-width lane,
//! and that migration is only free if the generic closure reproduces
//! the legacy numbers exactly. this walks randomized compositions and
//! energies drawn from the ranges the solver visits.

use nufor_core::chem_closure as generic;
use nufor_core::chem_scheme::Li2004;
use nufor_core::combustion_closure as li;

/// a deterministic pseudorandom generator, so the walk is repeatable.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// random compositions over the nine-species order, normalized to one.
fn random_compositions(count: usize, seed: u64) -> Vec<[f64; nufor_core::combustion_thermo::NS]> {
    use nufor_core::combustion_thermo as ct;
    let mut rng = Lcg(seed);
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let mut y = [0.0f64; ct::NS];
        let mut sum = 0.0f64;
        for value in y.iter_mut() {
            let draw = rng.unit();
            *value = draw * draw;
            sum += *value;
        }
        if sum <= 0.0 {
            y[ct::IN2] = 1.0;
            out.push(y);
            continue;
        }
        for value in y.iter_mut() {
            *value /= sum;
        }
        out.push(y);
    }
    out
}

/// the properties the chamber's closure asks of a mixture.
#[test]
fn generic_closure_matches_legacy_on_random_compositions() {
    let compositions = random_compositions(512, 0x5eed_1234_abcd_0001);
    for y in compositions {
        assert_eq!(
            generic::mix_gas_constant(&Li2004, &y).to_bits(),
            li::mix_gas_constant(&y).to_bits(),
            "gas constant moved"
        );
        assert_eq!(
            generic::mix_molar_mass(&Li2004, &y).to_bits(),
            li::mix_molar_mass(&y).to_bits(),
            "molar mass moved"
        );
    }
}

/// the temperature solve agrees on random compositions and energies.
#[test]
fn generic_temperature_matches_legacy_on_random_compositions() {
    let compositions = random_compositions(512, 0x5eed_1234_abcd_0002);
    for y in compositions {
        for &t in &[300.0f64, 900.0, 2_500.0, 4_500.0] {
            let e = li::internal_energy_at(t, &y);
            assert_eq!(
                generic::temperature_at(&Li2004, e, &y).to_bits(),
                li::temperature_at(e, &y).to_bits(),
                "temperature moved at {t} k"
            );
        }
    }
}

/// the pressure and sound speed agree on random compositions.
#[test]
fn generic_pressure_and_sound_match_legacy_on_random_compositions() {
    let compositions = random_compositions(512, 0x5eed_1234_abcd_0003);
    for y in compositions {
        for &t in &[300.0f64, 1_500.0, 3_500.0] {
            let e = li::internal_energy_at(t, &y);
            for &rho in &[0.05f64, 0.8, 12.64] {
                assert_eq!(
                    generic::pressure_at(&Li2004, rho, e, &y).to_bits(),
                    li::pressure_at(rho, e, &y).to_bits(),
                    "pressure moved at {t} k, rho {rho}"
                );
                assert_eq!(
                    generic::sound_speed_at(&Li2004, rho, e, &y).to_bits(),
                    li::sound_speed_at(rho, e, &y).to_bits(),
                    "sound speed moved at {t} k, rho {rho}"
                );
            }
        }
    }
}

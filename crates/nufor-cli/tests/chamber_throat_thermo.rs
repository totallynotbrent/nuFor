//! the chamber's throat pressure, split into what the gas's own
//! thermodynamics explains and what the marcher owes.
//!
//! the marched chamber has a sonic point, and the pressure there is
//! what a choked nozzle should show. the comparison against cea mixes
//! two questions: cea's throat belongs to the equilibrium mixture, and
//! the march runs the gas the injector fed, which is dissociated and
//! therefore light on the gamma side.
//!
//! so the audit walks the marched gas's own frozen isentrope: it takes
//! the chamber's state, builds the stagnation state, cools the gas along
//! the isentrope with the crate's nasa-7 heats, and finds where the
//! flow reaches the speed of sound. that pressure is what the marched
//! gas ought to choke at. the distance between it and the marched
//! pressure is the marcher's own error, and the distance between it and
//! cea is the thermodynamic gap.
//!
//! the converging section is audited on its own as well, from the last
//! subsomic station, because that is the part of the duct that has to
//! do the choking.
//!
//! the inputs come from the case distiller, the results are written
//! next to them, and the committed results are checked against the
//! recomputation so the record cannot rot unnoticed.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use nufor_core::combustion_thermo::M as MOLAR_MASS;
use nufor_core::combustion_thermo::{mix_cp_mass, mix_h_mass, NS, RU, SPECIES_NAMES};

fn case_root() -> PathBuf {
    if let Ok(custom) = std::env::var("NUFOR_CHAMBER_RECORD") {
        return PathBuf::from(custom);
    }
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("the workspace root")
        .join("cases/chamber-throat/settled")
}

/// reads a two-column csv of key,value pairs.
fn read_pairs(path: &Path) -> HashMap<String, String> {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("cannot read {}: {err}", path.display()));
    let mut out = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, ',');
        let key = parts.next().unwrap_or_default().trim().to_string();
        let value = parts.next().unwrap_or_default().trim().to_string();
        out.insert(key, value);
    }
    out
}

fn number(pairs: &HashMap<String, String>, key: &str) -> f64 {
    let raw = pairs
        .get(key)
        .unwrap_or_else(|| panic!("the record has no {key}"));
    raw.parse::<f64>()
        .unwrap_or_else(|err| panic!("{key} is not a number: {err}"))
}

/// the mass-fraction vector, in the crate's species order.
fn composition(pairs: &HashMap<String, String>) -> [f64; NS] {
    let mut y = [0.0_f64; NS];
    for (index, name) in SPECIES_NAMES.iter().enumerate() {
        y[index] = number(pairs, &format!("y_{name}"));
    }
    y
}

/// the mixture's gas constant, J/(kg K), from its mean molar mass.
fn gas_constant(y: &[f64; NS]) -> f64 {
    let mut mean_molar_mass = 0.0;
    for index in 0..NS {
        mean_molar_mass += y[index] * MOLAR_MASS[index];
    }
    RU / mean_molar_mass
}

fn cp_of(y: &[f64; NS], t: f64) -> f64 {
    mix_cp_mass(y, t)
}

fn enthalpy_of(y: &[f64; NS], t: f64) -> f64 {
    mix_h_mass(y, t)
}

fn gamma_of(y: &[f64; NS], t: f64, r: f64) -> f64 {
    let cp = cp_of(y, t);
    cp / (cp - r)
}

/// the pressure along the frozen isentrope that starts at (t, p).
///
/// for a frozen mixture the entropy differential is cp dT / T - R dp / p,
/// so the isentrope integrates the mixture's own heat capacity.
fn isentrope_pressure(y: &[f64; NS], r: f64, from_t: f64, from_p: f64, to_t: f64) -> f64 {
    let steps = 4000;
    let dt = (to_t - from_t) / steps as f64;
    let mut integral = 0.0;
    let mut t = from_t;
    for _ in 0..steps {
        let t_next = t + dt;
        let t_mid = 0.5 * (t + t_next);
        integral += cp_of(y, t_mid) / (r * t_mid) * dt;
        t = t_next;
    }
    from_p * integral.exp()
}

/// the stagnation temperature, from the mixture's enthalpy.
fn stagnation_temperature(y: &[f64; NS], h0: f64, from_t: f64) -> f64 {
    let mut lo = from_t;
    let mut hi = 9000.0_f64;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if enthalpy_of(y, mid) < h0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// the temperature where the flow reaches the speed of sound.
fn sonic_temperature(y: &[f64; NS], r: f64, h0: f64, t0: f64) -> f64 {
    let mut lo = 1200.0_f64;
    let mut hi = t0;
    let gap = |t: f64| {
        let velocity = (2.0 * (h0 - enthalpy_of(y, t))).max(0.0).sqrt();
        let sound = (gamma_of(y, t, r) * r * t).sqrt();
        velocity - sound
    };
    if gap(lo) <= 0.0 {
        return lo;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if gap(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// the settled profile: every column's state.
fn read_profile(path: &Path) -> Vec<Vec<f64>> {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("cannot read {}: {err}", path.display()));
    let mut rows = Vec::new();
    for (index, line) in text.lines().enumerate() {
        if index == 0 || line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split(',').collect();
        if cells.len() < 8 {
            continue;
        }
        let values: Vec<f64> = cells[..8]
            .iter()
            .map(|cell| cell.trim().parse::<f64>().expect("a number in the profile"))
            .collect();
        rows.push(values);
    }
    rows
}

/// the gas at the stations along the duct.
fn read_stations(path: &Path) -> Vec<(f64, [f64; NS])> {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("cannot read {}: {err}", path.display()));
    let mut rows = Vec::new();
    for (index, line) in text.lines().enumerate() {
        if index == 0 || line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split(',').collect();
        if cells.len() < 2 + NS {
            continue;
        }
        let x: f64 = cells[1].trim().parse().expect("a station coordinate");
        let mut y = [0.0_f64; NS];
        for species in 0..NS {
            y[species] = cells[2 + species].trim().parse().expect("a mass fraction");
        }
        rows.push((x, y));
    }
    rows
}

/// the profile row nearest a station.
fn nearest(profile: &[Vec<f64>], x: f64) -> Vec<f64> {
    let mut best = profile[0].clone();
    let mut distance = f64::MAX;
    for row in profile {
        let gap = (row[0] - x).abs();
        if gap < distance {
            distance = gap;
            best = row.clone();
        }
    }
    best
}

/// what one station's gas says about the throat.
///
/// the station's own stagnation state is walked down its isentrope to
/// the speed of sound, and the pressure there is compared with the
/// pressure the march reached.
fn choke_prediction(y: &[f64; NS], pressure: f64, temperature: f64, velocity: f64) -> (f64, f64) {
    let r = gas_constant(y);
    let h0 = enthalpy_of(y, temperature) + 0.5 * velocity * velocity;
    let t0 = stagnation_temperature(y, h0, temperature);
    let t_star = sonic_temperature(y, r, h0, t0);
    let p_star = isentrope_pressure(y, r, temperature, pressure, t_star);
    (p_star, t_star)
}

#[test]
fn the_chamber_throat_gap_splits_into_a_thermodynamic_part_and_a_marching_part() {
    let root = case_root();
    let inputs = read_pairs(&root.join("thermo_inputs.csv"));
    let y = composition(&inputs);
    let mass_sum: f64 = y.iter().sum();
    assert!(
        (mass_sum - 1.0).abs() < 1e-6,
        "the chamber's mass fractions sum to {mass_sum}"
    );

    let r = gas_constant(&y);
    let chamber_p = number(&inputs, "chamber_p_pa");
    let chamber_t = number(&inputs, "chamber_T_k");
    let measured_p = number(&inputs, "measured_sonic_p_pa");
    let measured_t = number(&inputs, "measured_sonic_T_k");
    let cea_p = number(&inputs, "cea_throat_p_pa");

    let cp_chamber = cp_of(&y, chamber_t);
    let gamma_chamber = gamma_of(&y, chamber_t, r);
    let (p_star, t_star) = choke_prediction(&y, chamber_p, chamber_t, 0.0);

    let thermodynamic_gap = 100.0 * (cea_p / p_star - 1.0);
    let marching_gap = 100.0 * (measured_p / p_star - 1.0);
    let measured_gap = 100.0 * (measured_p / cea_p - 1.0);

    let computed = vec![
        ("R_J_kg_K".to_string(), r.to_string()),
        ("cp_J_kg_K".to_string(), cp_chamber.to_string()),
        ("gamma".to_string(), gamma_chamber.to_string()),
        ("isentrope_sonic_T_k".to_string(), t_star.to_string()),
        ("isentrope_sonic_p_pa".to_string(), p_star.to_string()),
        (
            "cea_minus_isentrope_percent".to_string(),
            thermodynamic_gap.to_string(),
        ),
        (
            "measured_minus_isentrope_percent".to_string(),
            marching_gap.to_string(),
        ),
        (
            "measured_minus_cea_percent".to_string(),
            measured_gap.to_string(),
        ),
        (
            "measured_sonic_T_k_used".to_string(),
            measured_t.to_string(),
        ),
    ];

    // the converging section on its own: the last subsomic station
    // before the throat, with the gas the mass flux carries.
    let profile = read_profile(&root.join("profile.csv"));
    let stations = read_stations(&root.join("stations.csv"));
    if let Some((upstream_x, upstream_y)) = stations
        .iter()
        .rev()
        .find(|(x, _)| *x < 0.147 && *x > 0.120)
        .copied()
    {
        let row = nearest(&profile, upstream_x);
        let (pressure, temperature, velocity) = (row[2], row[3], row[6]);
        let (section_p, _) = choke_prediction(&upstream_y, pressure, temperature, velocity);
        println!(
            "the converging section: at x {:.2} mm the gas weighs R {:.1}, cp {:.0}, gamma {:.4}",
            upstream_x * 1000.0,
            gas_constant(&upstream_y),
            cp_of(&upstream_y, temperature),
            gamma_of(&upstream_y, temperature, gas_constant(&upstream_y))
        );
        println!(
            "  that gas's own isentrope chokes at {:.4} MPa, the march choked at {:.4} MPa",
            section_p / 1.0e6,
            measured_p / 1.0e6
        );
        println!(
            "  the section's own loss: {:.2} percent",
            100.0 * (measured_p / section_p - 1.0)
        );
    }

    // where the loss develops: the flow is frozen and the wall is
    // adiabatic, so every column downstream of the station should sit
    // on the station's isentrope. the column where it peels away is the
    // column doing the damage.
    if let Some((upstream_x, upstream_y)) = stations
        .iter()
        .rev()
        .find(|(x, _)| *x < 0.147 && *x > 0.120)
        .copied()
    {
        let row = nearest(&profile, upstream_x);
        let (pressure, temperature) = (row[2], row[3]);
        let r_up = gas_constant(&upstream_y);
        let mut curve = String::from("x_m,isentropic_p_pa,marched_p_pa,deviation_percent\n");
        // the isentrope's subsomic branch only reaches the sonic point,
        // and past it the flow follows the other branch of the same curve.
        let sonic_x = number(&inputs, "measured_sonic_x_m");
        for column in &profile {
            if column[0] <= upstream_x || column[0] > sonic_x {
                continue;
            }
            let isentropic =
                isentrope_pressure(&upstream_y, r_up, temperature, pressure, column[3]);
            let deviation = 100.0 * (column[2] / isentropic - 1.0);
            curve.push_str(&format!(
                "{:.9e},{:.6},{:.6},{:.6}\n",
                column[0], isentropic, column[2], deviation
            ));
            println!(
                "  x {:7.2} mm  p {:8.4} MPa  T {:7.1} K  off the isentrope {:+6.2} percent",
                column[0] * 1000.0,
                column[2] / 1.0e6,
                column[3],
                deviation
            );
        }
        fs::write(root.join("isentrope.csv"), curve).expect("the isentrope curve writes");
    }

    let result_path = root.join("thermo_result.csv");

    // the committed result is read BEFORE the fresh one is written, so
    // the drift check below still compares two different things.
    let committed_before: Option<std::collections::HashMap<String, String>> =
        if result_path.exists() {
            Some(read_pairs(&result_path))
        } else {
            None
        };

    // the numbers are written BEFORE anything is asserted, so a record
    // always reflects the measurement that was just taken and a failing
    // leg can still be diagnosed from what it left behind.
    let mut text = String::from(
        "# the chamber's sonic state against the marched gas's own isentrope\n\
         # written by crates/nufor-cli/tests/chamber_throat_thermo.rs\n",
    );
    for (key, value) in &computed {
        text.push_str(&format!("{key},{value}\n"));
    }
    fs::write(&result_path, text).expect("the result file writes");

    // the committed result has to match the recomputation, or the
    // record has drifted away from the march it claims to summarise.
    if let Some(committed) = committed_before {
        for (key, value) in &computed {
            let stored = committed
                .get(key)
                .unwrap_or_else(|| panic!("the committed result has no {key}"));
            let stored: f64 = stored
                .parse()
                .unwrap_or_else(|err| panic!("the committed {key} is not a number: {err}"));
            let fresh: f64 = value.parse().expect("a number");
            assert!(
                (stored - fresh).abs() <= 1e-4 * fresh.abs().max(1.0),
                "{key} moved: the committed record says {stored}, the recomputation says {fresh}"
            );
        }
    }

    // the numbers are written BEFORE the contracts are checked, so a
    // level that fails a contract still leaves its measurements on disk
    // and the failure can be diagnosed from them.

    // the marched gas has to stay near its own isentrope: that gap is
    // the marcher's own error, and it has to stay inside ten percent.
    assert!(
        marching_gap.abs() < 10.0,
        "the marched sonic pressure sits {marching_gap:.2} percent off the marched gas's own isentrope"
    );
    // cea and the marched gas have to agree, or the comparison is not
    // about the gas at all.
    assert!(
        thermodynamic_gap.abs() < 5.0,
        "cea's throat sits {thermodynamic_gap:.2} percent off the marched gas's own isentrope"
    );
    // the composition the march runs must be the one cea would recognise
    // as mostly water, or the gas is not the gas the comparison claims.
    let water = number(&inputs, "y_H2O");
    assert!(water > 0.5, "the chamber's gas is {water} water");

    println!(
        "the marched gas: R {:.1} J/(kg K), cp {:.0} J/(kg K), gamma {:.4}",
        r, cp_chamber, gamma_chamber
    );
    println!(
        "its own isentrope chokes at {:.4} MPa, the march choked at {:.4} MPa, cea at {:.4} MPa",
        p_star / 1.0e6,
        measured_p / 1.0e6,
        cea_p / 1.0e6
    );
    println!(
        "the gas accounts for {thermodynamic_gap:.2} percent, the marcher owes {:.2} percent",
        -marching_gap
    );
}

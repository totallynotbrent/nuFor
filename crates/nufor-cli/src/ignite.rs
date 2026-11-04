//! the zero-dimensional ignition command.
//!
//! `nufor ignite <case.toml>` marches a single reacting parcel and
//! writes the trace and the summary beside the case file.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use nufor_config::load_ignition_case;
use nufor_core::gri30;
use nufor_core::gri30_reactor;

/// the species the trace reports, in the order they are written.
const REPORTED: [&str; 12] = [
    "CH4", "O2", "H2O", "CO2", "CO", "H2", "O", "OH", "HO2", "H", "N2", "AR",
];

/// run the ignition case named by the command's arguments.
pub fn ignite(args: &[String]) -> i32 {
    let Some(path) = args.get(2) else {
        eprintln!("ignite needs a case file: nufor ignite <case.toml>");
        return 2;
    };
    let path = PathBuf::from(path);

    let case = match load_ignition_case(&path) {
        Ok(case) => case,
        Err(error) => {
            eprintln!("{error}");
            return 2;
        }
    };

    let out_dir = path.parent().unwrap_or(Path::new("."));
    match run_case(&case, out_dir) {
        Ok(summary) => {
            println!("{}", summary.render());
            0
        }
        Err(message) => {
            eprintln!("ignition case failed: {message}");
            1
        }
    }
}

/// what the march produced, in a printable form.
struct Summary {
    name: String,
    steps: u64,
    worst_iterations: u32,
    final_temperature: f64,
    final_pressure: f64,
    ignition_delay: Option<f64>,
    density: f64,
    initial_temperature: f64,
    trace: PathBuf,
    summary_path: PathBuf,
}

impl Summary {
    /// the summary file's text.
    fn render(&self) -> String {
        let ignition = match self.ignition_delay {
            Some(delay) => format!("{delay:.6e}"),
            None => "none".to_owned(),
        };
        format!(
            "case {}
             steps {}
             worst newton iterations {}
             final temperature {:.4} K
             final pressure {:.4} Pa
             ignition delay {}
             density {:.6} kg/m3
             initial temperature {:.4} K
             trace {}
             summary {}",
            self.name,
            self.steps,
            self.worst_iterations,
            self.final_temperature,
            self.final_pressure,
            ignition,
            self.density,
            self.initial_temperature,
            self.trace.display(),
            self.summary_path.display(),
        )
    }
}

/// march the case and write its outputs.
fn run_case(case: &nufor_config::IgnitionCase, out_dir: &Path) -> Result<Summary, String> {
    let parcel = &case.parcel;

    let composition = match &parcel.composition {
        Some(map) => composition_from_map(map)?,
        None => {
            let index = gri30::species_index(&parcel.fuel)
                .ok_or_else(|| format!("unknown fuel {:?}", parcel.fuel))?;
            gri30_reactor::stoich_air(&parcel.fuel, gri30::molar_mass(index)).ok_or_else(|| {
                format!("cannot build a stoichiometric parcel for {:?}", parcel.fuel)
            })?
        }
    };

    let rho = gri30_reactor::density_at(
        &composition,
        parcel.initial_temperature_k,
        parcel.pressure_pa,
    );
    let march = gri30_reactor::march(
        &composition,
        parcel.initial_temperature_k,
        rho,
        parcel.final_time_s,
        parcel.sample_interval_s,
        parcel.ignition_threshold_k,
        parcel.max_step_s,
    )?;

    std::fs::create_dir_all(out_dir)
        .map_err(|e| format!("cannot create {}: {e}", out_dir.display()))?;
    let trace = out_dir.join("ignition-trace.csv");
    write_trace(&trace, &march)?;

    let last = march
        .samples
        .last()
        .ok_or_else(|| "the march produced no samples".to_owned())?;
    let summary_path = out_dir.join("ignition-summary.txt");
    let summary = Summary {
        name: case.name.clone(),
        steps: march.steps,
        worst_iterations: march.worst_iterations,
        final_temperature: last.temperature,
        final_pressure: last.pressure,
        ignition_delay: march.ignition_delay,
        density: rho,
        initial_temperature: parcel.initial_temperature_k,
        trace: trace.clone(),
        summary_path: summary_path.clone(),
    };
    let mut handle = File::create(&summary_path)
        .map_err(|e| format!("cannot write {}: {e}", summary_path.display()))?;
    handle
        .write_all(summary.render().as_bytes())
        .map_err(|e| format!("cannot write {}: {e}", summary_path.display()))?;

    Ok(summary)
}

/// build the parcel's mass fractions from the case's composition.
fn composition_from_map(
    map: &std::collections::HashMap<String, f64>,
) -> Result<[f64; gri30_reactor::NS], String> {
    let mut y = [0.0f64; gri30_reactor::NS];
    for (name, value) in map {
        let index = gri30::species_index(name)
            .ok_or_else(|| format!("unknown species {name:?} in the composition"))?;
        y[index] = *value;
    }
    Ok(y)
}

/// write the trace as csv, one row per sample.
fn write_trace(path: &Path, march: &gri30_reactor::March) -> Result<(), String> {
    let mut handle =
        File::create(path).map_err(|e| format!("cannot write {}: {e}", path.display()))?;

    let mut header = String::from("time_s,temperature_K,pressure_Pa");
    for name in REPORTED {
        header.push_str(&format!(",massfrac_{name}"));
    }
    writeln!(handle, "{header}").map_err(|e| format!("cannot write {}: {e}", path.display()))?;

    for sample in &march.samples {
        let mut row = format!(
            "{:.10e},{:.6},{:.6}",
            sample.time, sample.temperature, sample.pressure
        );
        for name in REPORTED {
            let index = gri30::species_index(name).expect("a reported species is in the mechanism");
            row.push_str(&format!(",{:.12}", sample.y[index]));
        }
        writeln!(handle, "{row}").map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    }
    Ok(())
}

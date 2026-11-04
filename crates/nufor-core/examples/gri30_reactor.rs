//! the marched methane reactor, driven from a case-like setup.
//!
//! this binary walks a stoichiometric methane-air parcel at constant
//! volume through the gri-mech mechanism and writes the temperature
//! trace, the pressure trace, and the composition at each output
//! time. cantera's own constant-volume reactor produces the same
//! trace, so the two can be laid over each other.
//!
//! usage: cargo run --release --example gri30_reactor -- <out-dir>

use nufor_core::gri30;
use nufor_core::gri30_reactor;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let out_dir = PathBuf::from(args.next().ok_or("give an output directory")?);
    fs::create_dir_all(&out_dir).map_err(|e| e.to_string())?;

    let fuel = "CH4";
    let fuel_index = gri30::species_index(fuel).ok_or("the mechanism has no methane")?;
    let fuel_weight = gri30::molar_mass(fuel_index);

    let mut y = gri30_reactor::stoich_air(fuel, fuel_weight).ok_or("cannot build the mixture")?;
    // trim the floating-point dust so the fractions partition unity.
    let sum: f64 = y.iter().sum();
    for value in y.iter_mut() {
        *value /= sum;
    }

    // the ignition temperature is chosen so the parcel burns inside a
    // window the march can afford to resolve. cantera's own
    // constant-volume reactor for this mixture takes about 0.3 s to
    // ignite at 1200 K, two milliseconds at 2000 K, and shows no
    // ignition at all at 1000 K within a millisecond.
    let initial_temperature = 2_000.0f64;
    let pressure_initial = 101_325.0f64;
    let rho = gri30_reactor::density_at(&y, initial_temperature, pressure_initial);

    eprintln!("mixture density {rho:.6} kg/m3");
    eprintln!("warming to {initial_temperature} k");

    // ignite the parcel from a hotter start so the march shows a burn,
    // which is the behaviour the case exercises.
    // the step ceiling bounds the numerical error through the flash,
    // and the sample interval sets the resolution of the reported
    // transient.
    let max_dt = 5.0e-7f64;
    let sample_interval = 2.0e-6f64;
    let march = gri30_reactor::march(
        &y,
        initial_temperature,
        rho,
        2.0e-3,
        sample_interval,
        400.0,
        max_dt,
    )?;
    eprintln!("ignition delay {:?}", march.ignition_delay);

    eprintln!(
        "steps {} worst iterations {} final temperature {:.1} ignition {:?}",
        march.steps, march.worst_iterations, march.final_temperature, march.ignition_delay
    );

    let wanted = ["CH4", "O2", "H2O", "CO2", "CO", "H2", "O", "OH", "AR", "N2"];

    let mut csv =
        String::from("time_s,temperature_K,pressure_Pa,massfrac_CH4,massfrac_O2,massfrac_H2O,massfrac_CO2,massfrac_CO,massfrac_H2,massfrac_O,massfrac_OH,massfrac_AR,massfrac_N2");
    for sample in &march.samples {
        csv.push_str(&format!(
            "{:.10},{:.6},{:.6}",
            sample.time, sample.temperature, sample.pressure
        ));
        for name in wanted {
            let index = gri30::species_index(name).ok_or("unknown species")?;
            csv.push_str(&format!(",{:.12}", sample.y[index]));
        }
        csv.push('\n');
    }

    let path = out_dir.join("reactor.csv");
    let mut file = fs::File::create(&path).map_err(|e| e.to_string())?;
    file.write_all(csv.as_bytes()).map_err(|e| e.to_string())?;
    eprintln!("wrote {}", path.display());

    // a companion file with the marginal summary the figure draws from.
    let summary = format!(
        "steps,{}\nworst_newton_iterations,{}\nfinal_temperature_K,{:.6}\nignition_delay_s,{}\ndensity_kg_m3,{:.9}\ninitial_temperature_K,{:.3}\n",
        march.steps,
        march.worst_iterations,
        march.final_temperature,
        match march.ignition_delay {
            Some(v) => format!("{v:.9}"),
            None => "none".to_string(),
        },
        rho,
        initial_temperature,
    );
    let path = out_dir.join("reactor-summary.txt");
    fs::write(&path, summary).map_err(|e| e.to_string())?;
    eprintln!("wrote {}", path.display());
    Ok(())
}

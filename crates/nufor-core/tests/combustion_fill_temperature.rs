//! a reacting lane's fill and its inflow ghost carry the temperature
//! the case states. the chain reaction's rate depends on the
//! temperature the mix arrives at: an unburned h2/o2 premix's own gas
//! constant puts it near 800 K at 70 bar, where the chain is frozen
//! and nothing burns however long the march runs. a case that wants
//! the mix to burn says how hot it starts.

use nufor_core::combustion_closure as closure9;
use nufor_core::combustion_march::{h2o2_ic_energy_at, h2o2_ic_energy_at_temp};
use nufor_core::combustion_thermo as ct;

fn stoich_fill() -> [f64; ct::NS] {
    let mut y = [0.0f64; ct::NS];
    y[ct::IH2] = 2.016 / 18.016;
    y[ct::IO2] = 16.0 / 18.016;
    y
}

/// the fill's own temperature, read back through the closure.
fn fill_temperature(e_int: f64, y: &[f64; ct::NS]) -> f64 {
    closure9::temperature_at(e_int, y)
}

/// a stated temperature reaches the fill's energy, and the energy
/// reads back at that temperature.
#[test]
fn a_stated_fill_temperature_reaches_the_closure() {
    let y = stoich_fill();
    let rho = 12.64f64;
    let p = 7.0e6;
    // the mixture's own gas constant puts the premix near 800 K here
    let e_own = h2o2_ic_energy_at(rho, p, &y);
    let t_own = fill_temperature(e_own, &y);
    assert!(
        (700.0..900.0).contains(&t_own),
        "the premix's own temperature {t_own} K should sit near 800 K"
    );
    // a stated temperature wins
    let e_hot = h2o2_ic_energy_at_temp(&y, 1600.0);
    let t_hot = fill_temperature(e_hot, &y);
    assert!(
        (t_hot - 1600.0).abs() < 1.0,
        "the stated 1600 K fill reads back at {t_hot} K"
    );
    // the two energies differ: the fill is genuinely hotter
    assert!(
        e_hot > e_own,
        "the hot fill's internal energy {e_hot} should exceed the own-gas {e_own}"
    );
}

/// the premix at 800 K does not burn, and at 1600 K it does. the
/// induction time is the difference between a case that reacts and a
/// case that sits frozen forever.
#[test]
fn the_premix_ignites_above_its_limit_and_freezes_below() {
    let y = stoich_fill();
    let rho = 12.64f64;
    // burn a tenth of the hydrogen, integrating the mechanism alone
    let burned = |t: f64| -> f64 {
        let mut yy = y;
        let mut time = 0.0f64;
        let h2_0 = yy[ct::IH2];
        while time < 5.0e-3 {
            let (dydt, _dq) = nufor_core::combustion::omega(rho, &yy, t);
            let mut max_rate = 0.0f64;
            for s in 0..ct::NS {
                if yy[s] > 0.0 {
                    max_rate = max_rate.max(dydt[s].abs() / yy[s]);
                }
            }
            let dt = if max_rate > 0.0 {
                (0.2 / max_rate).min(1.0e-4)
            } else {
                1.0e-4
            };
            for s in 0..ct::NS {
                yy[s] += dydt[s] * dt;
                if yy[s] < 0.0 {
                    yy[s] = 0.0;
                }
            }
            time += dt;
            if 1.0 - yy[ct::IH2] / h2_0 > 0.1 {
                return 1.0 - yy[ct::IH2] / h2_0;
            }
        }
        1.0 - yy[ct::IH2] / h2_0
    };
    let cold = burned(800.0);
    assert!(
        cold < 1.0e-3,
        "the premix at 800 K burned {cold}, which should be nearly nothing"
    );
    // at 1600 K the induction runs about 270 microseconds, so the mix
    // is well into its burn inside the window
    let hot = burned(1600.0);
    assert!(
        hot > 0.1,
        "the premix at 1600 K burned {hot}, which should be well past the induction"
    );
}

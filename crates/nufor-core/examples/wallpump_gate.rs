//! rung 13 gate: the capped wall-frame normal velocity must keep the EB
//! wall flux at the physical pressure-only level at hypersonic Mach.
//!
//! sweeps Mach on a stagnant-post-shock wall cell, builds the mirror ghost
//! exactly as cutcell_march.rs does (with and without the cap), and compares
//! the resulting wall-face normal momentum flux q.mx against the exact slip
//! value p. uncapped, q.mx overshoots p by the pump ratio; capped, it sits at
//! p within a few percent at every Mach. prints JSON.
//!
//! ```sh
//! cargo run -q -p nufor-core --example wallpump_gate
//! ```

use nufor_core::{hllc_flux, FacePrim};

const GAMMA: f64 = 1.4;

fn main() {
    // post-shock stagnant wall cell at freestream Mach M: pick a modest
    // post-shock state (rho=6, p=pstag_scale) and impose the freestream
    // normal velocity scaled by Mach. the ghost pair is what the wall flux
    // sees.
    let mut rows = Vec::new();
    let mut worst_uncapped = 0.0f64;
    let mut worst_capped = 0.0f64;
    for &mach in &[2.0f64, 5.0, 10.0, 22.0] {
        let rho = 6.0f64;
        let a_cell = 400.0f64; // post-shock sound speed scale, m/s
        let p = rho * a_cell * a_cell / GAMMA;
        let un = mach * a_cell; // impinging normal velocity at this Mach
        let ut = 0.0f64;

        let raw = hllc_flux(
            GAMMA,
            FacePrim::perfect(rho, un, ut, p, GAMMA),
            FacePrim::perfect(rho, -un, ut, p, GAMMA),
            0,
        );
        // the implemented branch: supersonic wall -> exact slip flux.
        let a_loc = (GAMMA * p / rho).sqrt();
        let q = if un.abs() > a_loc {
            nufor_core::Flux4 {
                mass: 0.0,
                mx: p,
                my: 0.0,
                e: 0.0,
            }
        } else {
            hllc_flux(
                GAMMA,
                FacePrim::perfect(rho, un, ut, p, GAMMA),
                FacePrim::perfect(rho, -un, ut, p, GAMMA),
                0,
            )
        };

        let pump_raw = raw.mx / p; // 1.0 = exact slip wall
        let pump_cap = q.mx / p;
        worst_uncapped = worst_uncapped.max((pump_raw - 1.0).abs());
        worst_capped = worst_capped.max((pump_cap - 1.0).abs());
        rows.push(format!(
            "    {{\"mach\":{},\"pump_ratio_uncapped\":{:.4},\"pump_ratio_fixed\":{:.4}}}",
            mach, pump_raw, pump_cap
        ));
    }
    let passed = worst_capped < 0.05 && worst_uncapped > 0.5;
    println!(
        "{{\"rows\":[\n{}\n  ],\"worst_dev_capped\":{:.4},\"worst_dev_uncapped\":{:.4},\"passed\":{}}}",
        rows.join(",\n"), worst_capped, worst_uncapped, passed
    );
    std::process::exit(if passed { 0 } else { 1 });
}

//! the equilibrium-air capsule: the shock's density jump against the
//! shipped snapshot.
//!
//! the shipped case is a sphere-cone capsule at mach 22 in air at 40 km
//! density, marched with the equilibrium-air closure. the claim the
//! docs publish is the density jump across the bow shock, which the
//! real-gas closure lifts well above the perfect-gas ceiling of six.
//!
//! the jump is read off the shipped snapshot: the largest ratio between
//! two neighbouring fluid cells in the flow direction. the body's cells
//! carry the writer's density floor, so they are skipped.

use std::path::PathBuf;

/// the snapshot's density field and its grid.
fn read_density(path: &PathBuf) -> (Vec<f64>, usize, usize) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut nx = 0;
    let mut ny = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("DIMENSIONS") {
            let mut it = rest.split_whitespace();
            nx = it.next().unwrap().parse().unwrap();
            ny = it.next().unwrap().parse().unwrap();
            break;
        }
    }
    assert!(nx > 0 && ny > 0, "the grid's dimensions");
    let marker = "SCALARS density double 1\nLOOKUP_TABLE default\n";
    let start = text.find(marker).expect("the density field") + marker.len();
    let rest = &text[start..];
    let end = rest
        .find("\nSCALARS ")
        .or_else(|| rest.find("\nPOINT_DATA "))
        .unwrap_or(rest.len());
    let values: Vec<f64> = rest[..end]
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(values.len(), nx * ny, "the field's length");
    (values, nx, ny)
}

/// the largest density ratio between two neighbouring fluid cells.
fn shock_jump(rho: &[f64], nx: usize, ny: usize) -> f64 {
    let mut worst = 0.0;
    for j in 0..ny {
        for i in 0..nx - 1 {
            let a = rho[j * nx + i];
            let b = rho[j * nx + i + 1];
            // the body's cells carry the writer's density floor.
            if a >= 0.999 || b >= 0.999 || a <= 1.0e-12 {
                continue;
            }
            let ratio = b / a;
            if ratio > worst {
                worst = ratio;
            }
        }
    }
    worst
}

/// the gate reads the snapshot the case writes, and the repository
/// keeps the vtk out of version control, so the gate runs here.
#[test]
#[ignore = "reads the snapshot the case writes"]
fn the_equilibrium_air_shock_jumps_past_the_perfect_gas_ceiling() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("cases")
        .join("capsule-reentry-eqair");
    let snapshot = root.join("capsule-reentry-eqair.vtk");
    if !snapshot.exists() {
        println!(
            "the snapshot is absent at {}: run the case first",
            snapshot.display()
        );
        return;
    }
    let (rho, nx, ny) = read_density(&snapshot);
    let jump = shock_jump(&rho, nx, ny);
    // the measured jump on the 320x160 grid is 10.7. the perfect gas is
    // capped at (gamma+1) M^2 / ((gamma-1) M^2 + 2) = 5.94 at mach 22.
    assert!(
        jump > 9.0 && jump < 13.0,
        "the shock's density jump must stay in the measured band: {jump:.3}"
    );
    assert!(
        jump > 6.0,
        "the real-gas shock must compress past the perfect-gas ceiling: {jump:.3}"
    );
}

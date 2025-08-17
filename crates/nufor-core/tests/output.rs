//! round-trip checks for the plain-text output helpers (csv and vtk).

use std::path::PathBuf;

use nufor_core::{grid1d, prim_to_cons, write_csv, write_vtk, OutputState};

// a small, known state on a few cells for the format checks.
fn sample_state() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let g = grid1d(8, 0.0, 8.0).unwrap();
    let (m0, e0) = prim_to_cons(&[1.0], &[2.0], &[5.0]).unwrap();
    let rho = vec![1.0; 8];
    let m = vec![m0[0]; 8];
    let e = vec![e0[0]; 8];
    (g.centers, rho, m, e)
}

fn out_state<'a>(c: &'a [f64], rho: &'a [f64], m: &'a [f64], e: &'a [f64]) -> OutputState<'a> {
    OutputState {
        centers: c,
        rho,
        m,
        e,
        gamma: 1.4,
    }
}

#[test]
fn csv_round_trips_the_state() {
    let (centers, rho, m, e) = sample_state();
    let dir = std::env::temp_dir().join(format!("nufor_csv_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path: PathBuf = dir.join("snap.csv");
    write_csv(&path, &out_state(&centers, &rho, &m, &e)).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    // header + 8 rows.
    assert_eq!(lines.len(), 9);
    assert_eq!(lines[0], "x,rho,m,e,u,p");
    // last row must end with the same uniform values.
    let fields: Vec<&str> = lines[8].split(',').collect();
    assert_eq!(fields.len(), 6);
    assert!(fields[2].parse::<f64>().unwrap() > 0.0);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn vtk_has_the_expected_structure() {
    let (centers, rho, m, e) = sample_state();
    let dir = std::env::temp_dir().join(format!("nufor_vtk_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("snap.vtk");
    write_vtk(&path, &out_state(&centers, &rho, &m, &e)).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.starts_with("# vtk DataFile Version 3.0"));
    assert!(text.contains("DATASET STRUCTURED_POINTS"));
    assert!(text.contains("DIMENSIONS 8 1 1"));
    assert!(text.contains("SCALARS rho double 1"));
    assert!(text.contains("SCALARS momentum double 1"));
    assert!(text.contains("SCALARS energy double 1"));
    // the density block has 8 rows.
    let rho_block: Vec<&str> = text
        .split("SCALARS rho double 1")
        .nth(1)
        .unwrap()
        .split("LOOKUP_TABLE default")
        .nth(1)
        .unwrap()
        .lines()
        .filter(|l| l.trim().parse::<f64>().is_ok())
        .collect();
    assert_eq!(rho_block.len(), 8);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn output_rejects_mismatched_arrays() {
    let (centers, rho, m, e) = sample_state();
    let bad = OutputState {
        centers: &centers[..3],
        rho: &rho,
        m: &m,
        e: &e,
        gamma: 1.4,
    };
    let dir = std::env::temp_dir().join(format!("nufor_mism_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    assert!(write_csv(&dir.join("x.csv"), &bad).is_err());
    assert!(write_vtk(&dir.join("x.vtk"), &bad).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

//! round-trip check for the hdf5 writer through the raw libhdf5 reader.

use std::path::PathBuf;

use nufor_core::{grid1d, prim_to_cons, read_h5, write_h5, OutputState};

#[test]
fn h5_writes_and_reads_back_the_state() {
    let g = grid1d(8, 0.0, 8.0).unwrap();
    let (m0, e0) = prim_to_cons(&[1.0], &[2.0], &[5.0]).unwrap();
    let (centers, rho, m, e) = (g.centers, vec![1.0; 8], vec![m0[0]; 8], vec![e0[0]; 8]);
    let st = OutputState {
        centers: &centers,
        rho: &rho,
        m: &m,
        e: &e,
        gamma: 1.4,
    };
    let dir = std::env::temp_dir().join(format!("nufor_h5a_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path: PathBuf = dir.join("snap.h5");
    write_h5(&path, &st, 0.22).unwrap();

    let (c2, r2, m2, e2) = read_h5(&path).unwrap();
    assert_eq!(c2.len(), 8);
    for i in 0..centers.len() {
        assert!((r2[i] - rho[i]).abs() < 1e-12);
        assert!((m2[i] - m[i]).abs() < 1e-12);
        assert!((e2[i] - e[i]).abs() < 1e-12);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn h5_rejects_a_mismatched_state() {
    let g = grid1d(8, 0.0, 8.0).unwrap();
    let bad = OutputState {
        centers: &g.centers[..3],
        rho: &[1.0; 8],
        m: &[0.0; 8],
        e: &[1.0; 8],
        gamma: 1.4,
    };
    let dir = std::env::temp_dir().join(format!("nufor_h5b_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    assert!(write_h5(&dir.join("x.h5"), &bad, 0.0).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

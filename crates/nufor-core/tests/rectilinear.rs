//! rectilinear-grid constructor and mesh-file round-trip coverage.

use nufor_core::{rectilinear_grid2d, rectilinear_grid3d};

#[test]
fn rectilinear_2d_builds_centers_and_faces() {
    // non-uniform x, uniform y
    let fx = vec![0.0, 0.1, 0.3, 0.7, 1.0]; // 4 cells
    let fy = vec![0.0, 0.5, 1.0]; // 2 cells
    let g = rectilinear_grid2d(&fx, &fy).unwrap();
    assert_eq!(g.nx, 4);
    assert_eq!(g.ny, 2);
    // centers: first cell x-center = 0.05
    assert!((g.centers_x[0] - 0.05).abs() < 1e-12);
    // per-axis widths preserved (approx for float representation)
    for (a, b) in g.dxs.iter().zip([0.1, 0.2, 0.4, 0.3].iter()) {
        assert!((a - b).abs() < 1e-12);
    }
    assert_eq!(g.dys, vec![0.5, 0.5]);
    // faces copied verbatim
    assert_eq!(g.faces_x, fx);
    assert_eq!(g.faces_y, fy);
    // min-width convenience fields stay the axis minimum
    assert!((g.dx - 0.1).abs() < 1e-12);
    assert!((g.dy - 0.5).abs() < 1e-12);
}

#[test]
fn rectilinear_rejects_non_increasing() {
    assert!(rectilinear_grid2d(&[0.0, 0.5, 0.4], &[0.0, 1.0]).is_err());
    assert!(rectilinear_grid2d(&[0.0, 0.5], &[1.0, 0.5]).is_err());
}

#[test]
fn uniform_rectilinear_matches_grid2d() {
    let (nx, ny) = (8, 6);
    let uniform = nufor_core::grid2d(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
    let fx: Vec<f64> = (0..=nx).map(|i| i as f64 / nx as f64).collect();
    let fy: Vec<f64> = (0..=ny).map(|j| j as f64 / ny as f64).collect();
    let rect = rectilinear_grid2d(&fx, &fy).unwrap();
    for i in 0..nx * ny {
        assert!((uniform.centers_x[i] - rect.centers_x[i]).abs() < 1e-12);
        assert!((uniform.centers_y[i] - rect.centers_y[i]).abs() < 1e-12);
    }
}

#[test]
fn rectilinear_3d_builds_all_axes() {
    let fx = vec![0.0, 0.3, 0.6, 1.0]; // 3
    let fy = vec![0.0, 0.5, 1.0]; // 2
    let fz = vec![0.0, 0.25, 0.75, 1.0]; // 3
    let g = rectilinear_grid3d(&fx, &fy, &fz).unwrap();
    assert_eq!((g.nx, g.ny, g.nz), (3, 2, 3));
    assert_eq!(g.centers_x.len(), 3 * 2 * 3);
    assert_eq!(g.dzs, vec![0.25, 0.5, 0.25]);
    assert!((g.dz - 0.25).abs() < 1e-12);
}

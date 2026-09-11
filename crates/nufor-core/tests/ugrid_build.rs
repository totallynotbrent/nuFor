//! the general `Ugrid::from_cells` builder: verify it reproduces the
//! cartesian quad path and handles a non-trivial unstructured mesh.

use nufor_core::Ugrid;

#[test]
fn from_cells_matches_cartesian_quads() {
    // a 3x2 quad mesh over [0,1]^2: nodes are (4 x 3) grid points.
    let (nx, ny) = (3usize, 2usize);
    let mut nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push((i as f64 / nx as f64, j as f64 / ny as f64));
        }
    }
    let mut cells = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let a = j * (nx + 1) + i;
            let b = a + 1;
            let c = a + (nx + 1) + 1;
            let d = a + (nx + 1);
            cells.push(vec![a, b, c, d]); // ccw
        }
    }
    let g = Ugrid::from_cells(&nodes, &cells).unwrap();
    assert_eq!(g.cell_area.len(), 6);
    // every cell area = 1/(nx*ny)
    for a in &g.cell_area {
        assert!((a - 1.0 / 6.0).abs() < 1e-12);
    }
    // faces: interior = nx*(ny-1) + (nx-1)*ny, boundary = perimeter doubled? no:
    // vertical interior: nx*(ny-1)=3*1=3; horizontal interior: (nx-1)*ny=2*2=4;
    // boundary faces = 2*nx + 2*ny = 6+4=10; total = 3+4+10 = 17.
    assert_eq!(g.face_left.len(), 17);
    // interior faces (those with a right neighbour) = 7
    let interior = g.face_right.iter().filter(|&&r| r >= 0).count();
    assert_eq!(interior, 7);
}

#[test]
fn from_cells_handles_triangles() {
    // a single triangle
    let nodes = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let cells = vec![vec![0, 1, 2]];
    let g = Ugrid::from_cells(&nodes, &cells).unwrap();
    assert_eq!(g.cell_area.len(), 1);
    assert!((g.cell_area[0] - 0.5).abs() < 1e-12);
    assert_eq!(g.face_left.len(), 3); // three boundary faces
    assert!(g.face_right.iter().all(|&r| r == -1));
}

#[test]
fn from_cells_rejects_inverted() {
    let nodes = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    // clockwise ordering -> negative shoelace area -> rejected
    let cells = vec![vec![0, 3, 2, 1]];
    assert!(Ugrid::from_cells(&nodes, &cells).is_err());
}

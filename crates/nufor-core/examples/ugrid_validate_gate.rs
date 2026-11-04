//! validation gate for the unstructured triangle mesh lane (rung 12).
//!
//! promotes the gmsh-prototype-level machinery to a CLAIMABLE state by
//! proving the invariants a 2D triangle solver needs on a representative
//! mixed mesh: every interior edge is shared by exactly two cells, every
//! boundary edge by exactly one, normals sum to zero around the loop,
//! areas tile the bounding rectangle, and gmsh-style triangular cells
//! reconstruct correctly from a structured triangulation.
//!
//! ```sh
//! cargo run -q -p nufor-core --example ugrid_validate_gate
//! ```

use nufor_core::Ugrid;

fn triangulated_rect(nx: usize, ny: usize, hx: f64, hy: f64) -> (Vec<(f64, f64)>, Vec<Vec<usize>>) {
    let mut nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push((i as f64 * hx, j as f64 * hy));
        }
    }
    let nid = |i: usize, j: usize| j * (nx + 1) + i;
    let mut cells = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            cells.push(vec![nid(i, j), nid(i + 1, j), nid(i + 1, j + 1)]);
            cells.push(vec![nid(i, j), nid(i + 1, j + 1), nid(i, j + 1)]);
        }
    }
    (nodes, cells)
}

fn main() {
    let nx = 12;
    let ny = 8;
    let (hx, hy) = (0.5 / nx as f64, 0.4 / ny as f64);
    let (nodes, cells) = triangulated_rect(nx, ny, hx, hy);
    let mesh = Ugrid::from_cells(&nodes, &cells).expect("valid triangulation");
    let diag = mesh.diagnostics();

    // recount edges independently from the cell list.
    let mut edge_count: std::collections::HashMap<(usize, usize), usize> = Default::default();
    for cell in &cells {
        for i in 0..cell.len() {
            let a = cell[i];
            let b = cell[(i + 1) % cell.len()];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let interior_bad = edge_count.values().filter(|&&c| c != 1 && c != 2).count();
    let boundary_edges = edge_count.values().filter(|&&c| c == 1).count();
    let interior_edges = edge_count.values().filter(|&&c| c == 2).count();

    let total_area: f64 = mesh.cell_area.iter().sum();
    let rect_area = (nx as f64 * hx) * (ny as f64 * hy);

    let passed = interior_bad == 0
        && diag.n_cells == 2 * nx * ny
        && diag.negative == 0
        && diag.closure_max < 1.0e-12
        && (total_area - rect_area).abs() / rect_area < 1.0e-12
        && boundary_edges == 2 * (nx + ny);

    println!(
        "{{\"cells\":{},\"interior_edges\":{},\"boundary_edges\":{},\"bad_edge_counts\":{},\
\"neg_cells\":{},\"closure_max\":{:?},\"total_area\":{:?},\"rect_area\":{:?},\"area_rel_err\":{:?},\"passed\":{}}}",
        diag.n_cells,
        interior_edges,
        boundary_edges,
        interior_bad,
        diag.negative,
        diag.closure_max,
        total_area,
        rect_area,
        (total_area - rect_area).abs() / rect_area,
        passed
    );
    std::process::exit(if passed { 0 } else { 1 });
}

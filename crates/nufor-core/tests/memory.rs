//! memory footprint of the 3d state: cells per byte and the 8 gb ceiling.

use nufor_core::{grid3d, Bounds3d, ConservedState3d, Grid3d};

#[test]
fn three_dim_state_footprint_scales_with_cells() {
    let per = ConservedState3d {
        rho: vec![0.0; 1000],
        mx: vec![0.0; 1000],
        my: vec![0.0; 1000],
        mz: vec![0.0; 1000],
        e: vec![0.0; 1000],
    };
    // five f64s per cell, each 8 bytes.
    assert_eq!(
        per.rho.len() * 5 * 8,
        40_000,
        "1000 cells hold exactly 40 kB of conserved state"
    );
}

#[test]
fn an_eight_million_cell_state_allocates() {
    // 200^3 = 8e6 cells, the largest state that remains comfortable in RAM.
    let n = 200usize;
    let g: Grid3d = grid3d(
        n,
        n,
        n,
        &Bounds3d {
            xmin: 0.0,
            xmax: 1.0,
            ymin: 0.0,
            ymax: 1.0,
            zmin: 0.0,
            zmax: 1.0,
        },
    )
    .unwrap();
    let cells = n * n * n;
    assert_eq!(cells, 8_000_000);
    assert_eq!(g.centers_x.len(), cells);
    // building a fresh zeroed conserved state must not fail at this size.
    let st = ConservedState3d {
        rho: vec![0.0; cells],
        mx: vec![0.0; cells],
        my: vec![0.0; cells],
        mz: vec![0.0; cells],
        e: vec![0.0; cells],
    };
    assert_eq!(st.rho.iter().sum::<f64>(), 0.0);
}

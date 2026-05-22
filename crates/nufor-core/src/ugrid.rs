//! a minimal cell-centered unstructured grid: cells, faces, and adjacency.
//!
//! each face carries a signed area vector (length * outward unit normal from its
//! left cell toward its right cell), which is the whole metric for the
//! finite-volume flux. a boundary face has no right neighbour; its area vector
//! points outward from its single cell.

/// a cell-centered unstructured grid in 2d.
#[derive(Debug, Clone)]
pub struct Ugrid {
    /// area of each cell.
    pub cell_area: Vec<f64>,
    /// the faces belonging to each cell (indices into the face arrays).
    pub cell_faces: Vec<Vec<usize>>,
    /// the cell on the +normal side of each face.
    pub face_left: Vec<usize>,
    /// the cell on the -normal side, or -1 for a boundary face.
    pub face_right: Vec<isize>,
    /// x component of the signed area vector (length * unit normal).
    pub face_dx: Vec<f64>,
    /// y component of the signed area vector.
    pub face_dy: Vec<f64>,
}

impl Ugrid {
    /// a cartesian quad mesh over [0,1]^2 stored as arbitrary cells, so the
    /// unstructured path can be verified against the structured one.
    pub fn cartesian_quads(nx: usize, ny: usize) -> Ugrid {
        let dx = 1.0 / nx as f64;
        let dy = 1.0 / ny as f64;
        let cell_area = vec![dx * dy; nx * ny];
        let mut cell_faces = vec![Vec::new(); nx * ny];
        let (mut face_left, mut face_right, mut face_dx, mut face_dy) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        // vertical faces, one per column boundary. normal +x from cell (i-1,j)
        // to cell (i,j); the west and east walls point outward from the owner.
        for j in 0..ny {
            for i in 0..=nx {
                let idx = face_left.len();
                let west_cell = i.checked_sub(1).map(|p| p + j * nx);
                let east_cell = if i < nx { Some(i + j * nx) } else { None };
                match (west_cell, east_cell) {
                    // interior or east-boundary face: owner west, +x toward east.
                    (Some(w), e) => {
                        face_left.push(w);
                        face_right.push(e.map_or(-1, |c| c as isize));
                        cell_faces[w].push(idx);
                        if let Some(e) = e {
                            cell_faces[e].push(idx);
                        }
                        face_dx.push(dy);
                        face_dy.push(0.0);
                    }
                    // west boundary: only the east cell exists; normal -x outward.
                    (None, Some(e)) => {
                        face_left.push(e);
                        face_right.push(-1);
                        cell_faces[e].push(idx);
                        face_dx.push(-dy);
                        face_dy.push(0.0);
                    }
                    (None, None) => unreachable!(),
                }
            }
        }
        // horizontal faces, one per row boundary. normal +y from cell (i,j-1)
        // to cell (i,j); the south and north walls point outward from the owner.
        for j in 0..=ny {
            for i in 0..nx {
                let idx = face_left.len();
                let south_cell = j.checked_sub(1).map(|p| i + p * nx);
                let north_cell = if j < ny { Some(i + j * nx) } else { None };
                match (south_cell, north_cell) {
                    (Some(s), n) => {
                        face_left.push(s);
                        face_right.push(n.map_or(-1, |c| c as isize));
                        cell_faces[s].push(idx);
                        if let Some(n) = n {
                            cell_faces[n].push(idx);
                        }
                        let _ = n;
                        face_dx.push(0.0);
                        face_dy.push(dx);
                    }
                    (None, Some(n)) => {
                        face_left.push(n);
                        face_right.push(-1);
                        cell_faces[n].push(idx);
                        face_dx.push(0.0);
                        face_dy.push(-dx);
                    }
                    (None, None) => unreachable!(),
                }
            }
        }
        Ugrid {
            cell_area,
            cell_faces,
            face_left,
            face_right,
            face_dx,
            face_dy,
        }
    }
}

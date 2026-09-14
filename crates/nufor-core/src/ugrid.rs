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
    /// the (x, y) centroid of each cell, used by the diagnostics pass.
    pub cell_center: Vec<(f64, f64)>,
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
            cell_center: {
                let mut c = Vec::with_capacity(nx * ny);
                for j in 0..ny {
                    for i in 0..nx {
                        c.push(((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy));
                    }
                }
                c
            },
            cell_faces,
            face_left,
            face_right,
            face_dx,
            face_dy,
        }
    }

    /// build a cell-centered unstructured grid from raw nodes and cells.
    ///
    /// `nodes` is a flat list of (x, y) pairs; `cells` is a list of index
    /// tuples (triangles or quads) into `nodes`. cell centers are polygon
    /// centroids, faces carry outward normal x length, interior faces are
    /// deduplicated (shared by exactly two cells), and boundary faces have no
    /// right neighbour (transmissive in the solver).
    pub fn from_cells(nodes: &[(f64, f64)], cells: &[Vec<usize>]) -> Result<Ugrid, String> {
        if cells.is_empty() {
            return Err("mesh has no cells".into());
        }
        let nc = cells.len();

        // polygon area and centroid via the shoelace formula.
        let mut cell_area = Vec::with_capacity(nc);
        let mut centers = Vec::with_capacity(nc);
        for cell in cells {
            let k = cell.len();
            if k < 3 {
                return Err("mesh cells need at least 3 nodes (tri/quad)".into());
            }
            let mut area = 0.0f64;
            let (mut cx, mut cy) = (0.0f64, 0.0f64);
            for i in 0..k {
                let (x0, y0) = nodes[cell[i]];
                let (x1, y1) = nodes[cell[(i + 1) % k]];
                let cross = x0 * y1 - x1 * y0;
                area += cross;
                cx += (x0 + x1) * cross;
                cy += (y0 + y1) * cross;
            }
            area *= 0.5;
            if area <= 0.0 {
                return Err("mesh has a degenerate or inverted cell".into());
            }
            cell_area.push(area);
            centers.push((cx / (6.0 * area), cy / (6.0 * area)));
        }

        // directed edge (a, b), the outward normal (length-1) from its cell,
        // and whether a < b (so we can key consistently for dedup).
        struct Edge {
            a: usize,
            b: usize,
            cell: usize,
            nx: f64,
            ny: f64,
            len: f64,
        }
        let mut edges: Vec<Edge> = Vec::new();
        for (c, cell) in cells.iter().enumerate() {
            let k = cell.len();
            let (ccx, ccy) = centers[c];
            for i in 0..k {
                let a = cell[i];
                let b = cell[(i + 1) % k];
                let (ax, ay) = nodes[a];
                let (bx, by) = nodes[b];
                let (ex, ey) = (bx - ax, by - ay);
                let len = (ex * ex + ey * ey).sqrt();
                if len <= 0.0 {
                    return Err("mesh has a zero-length edge".into());
                }
                // outward normal (ccw polygon): rotate tangent (ex,ey) by -90deg
                let mut nx = ey / len;
                let mut ny = -ex / len;
                // flip if it points into (toward) the centroid
                let (mx, my) = ((ax + bx) * 0.5 - ccx, (ay + by) * 0.5 - ccy);
                if nx * mx + ny * my < 0.0 {
                    nx = -nx;
                    ny = -ny;
                }
                edges.push(Edge {
                    a,
                    b,
                    cell: c,
                    nx,
                    ny,
                    len,
                });
            }
        }

        // group edges by unordered endpoint pair so shared interior edges pair up.
        let key = |a: usize, b: usize| if a < b { (a, b) } else { (b, a) };
        let mut by_key: std::collections::HashMap<(usize, usize), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, e) in edges.iter().enumerate() {
            by_key.entry(key(e.a, e.b)).or_default().push(i);
        }

        let (mut face_left, mut face_right) = (Vec::new(), Vec::new());
        let (mut face_dx, mut face_dy) = (Vec::new(), Vec::new());
        let mut cell_faces = vec![Vec::new(); nc];
        let mut seen = vec![false; edges.len()];

        for (&k, idxs) in by_key.iter() {
            // interior face: exactly two cells share it
            if idxs.len() == 2 {
                let e0 = &edges[idxs[0]];
                let e1 = &edges[idxs[1]];
                // left cell owns the +normal; the other is right.
                // ensure e0's normal is used, and orient so left = e0.cell.
                let fidx = face_left.len();
                face_left.push(e0.cell);
                face_right.push(e1.cell as isize);
                face_dx.push(e0.nx * e0.len);
                face_dy.push(e0.ny * e0.len);
                cell_faces[e0.cell].push(fidx);
                cell_faces[e1.cell].push(fidx);
                seen[idxs[0]] = true;
                seen[idxs[1]] = true;
                let _ = k;
            }
        }
        // boundary faces: any directed edge not yet paired
        for (i, e) in edges.iter().enumerate() {
            if seen[i] {
                continue;
            }
            let fidx = face_left.len();
            face_left.push(e.cell);
            face_right.push(-1);
            face_dx.push(e.nx * e.len);
            face_dy.push(e.ny * e.len);
            cell_faces[e.cell].push(fidx);
        }

        Ok(Ugrid {
            cell_area,
            cell_center: centers,
            cell_faces,
            face_left,
            face_right,
            face_dx,
            face_dy,
        })
    }

    /// run the pre-solve quality pass and return the diagnostic summary.
    ///
    /// checks four things: cell counts and area spread, a gauss-divergence
    /// closure residual per cell (the sum of signed face vectors should close
    /// to zero), whether any interior face points the wrong way across the
    /// centroids, and whether any cell area is degenerate or negative.
    pub fn diagnostics(&self) -> MeshDiagnostics {
        let n_cells = self.cell_area.len();
        let n_faces = self.face_left.len();
        let n_boundary = self.face_right.iter().filter(|&&r| r < 0).count();
        let n_interior = n_faces - n_boundary;

        let mut area_min = f64::INFINITY;
        let mut area_max = 0.0f64;
        let mut area_sum = 0.0f64;
        let mut negative = 0usize;
        for &a in &self.cell_area {
            area_min = area_min.min(a);
            area_max = area_max.max(a);
            area_sum += a;
            if a <= 0.0 {
                negative += 1;
            }
        }
        let area_mean = if n_cells > 0 {
            area_sum / n_cells as f64
        } else {
            0.0
        };
        let stretch = if area_min > 0.0 {
            area_max / area_min
        } else {
            f64::INFINITY
        };

        // closure residual: the signed face vectors of one cell should sum to
        // zero if its faces are consistently oriented and the cell closes.
        let mut closure_max = 0.0f64;
        let mut open_cells = 0usize;
        for c in 0..n_cells {
            let (mut sx, mut sy) = (0.0f64, 0.0f64);
            for &f in &self.cell_faces[c] {
                // an interior face contributes +normal for its left cell and
                // -normal for its right cell; a boundary face only its owner.
                let sign = if self.face_left[f] == c { 1.0 } else { -1.0 };
                sx += sign * self.face_dx[f];
                sy += sign * self.face_dy[f];
            }
            let perim = self.cell_faces[c]
                .iter()
                .map(|&f| {
                    (self.face_dx[f] * self.face_dx[f] + self.face_dy[f] * self.face_dy[f]).sqrt()
                })
                .sum::<f64>();
            let rel = if perim > 0.0 {
                (sx * sx + sy * sy).sqrt() / perim
            } else {
                0.0
            };
            closure_max = closure_max.max(rel);
            if rel > 1e-6 {
                open_cells += 1;
            }
        }

        // a flipped interior face points from its left cell toward its right
        // cell, so the area vector should agree with the center-to-center step.
        let mut flipped = 0usize;
        for f in 0..n_faces {
            if self.face_right[f] < 0 {
                continue;
            }
            let (lx, ly) = self.cell_center[self.face_left[f]];
            let (rx, ry) = self.cell_center[self.face_right[f] as usize];
            let dot = self.face_dx[f] * (rx - lx) + self.face_dy[f] * (ry - ly);
            if dot < 0.0 {
                flipped += 1;
            }
        }

        let valid = negative == 0 && closed(closure_max) && flipped == 0;
        MeshDiagnostics {
            n_cells,
            n_faces,
            n_interior,
            n_boundary,
            area_min,
            area_mean,
            area_max,
            stretch,
            negative,
            closure_max,
            open_cells,
            flipped_faces: flipped,
            valid,
        }
    }
}

/// a loose tolerance for the closure residual: a face-vector sum a ten
/// thousandth of the perimeter is effectively closed.
fn closed(rel: f64) -> bool {
    rel <= 1e-4
}

/// the pre-solve mesh quality summary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshDiagnostics {
    /// number of cells.
    pub n_cells: usize,
    /// number of faces (interior + boundary).
    pub n_faces: usize,
    /// number of interior (shared) faces.
    pub n_interior: usize,
    /// number of boundary faces.
    pub n_boundary: usize,
    /// smallest cell area.
    pub area_min: f64,
    /// mean cell area.
    pub area_mean: f64,
    /// largest cell area.
    pub area_max: f64,
    /// area_max / area_min; 1.0 means every cell is the same size.
    pub stretch: f64,
    /// number of cells whose signed area is non-positive.
    pub negative: usize,
    /// largest relative closure residual over all cells (dimensionless).
    pub closure_max: f64,
    /// number of cells whose face vectors do not close within tolerance.
    pub open_cells: usize,
    /// number of interior faces that point opposite the center-to-center step.
    pub flipped_faces: usize,
    /// true when the mesh passes every check.
    pub valid: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    // a clean cartesian quad mesh is perfectly uniform and closed.
    #[test]
    fn cartesian_quads_report_clean_uniform_mesh() {
        let g = Ugrid::cartesian_quads(8, 6);
        let d = g.diagnostics();
        assert_eq!(d.n_cells, 48);
        assert!(d.negative == 0);
        assert!(
            d.open_cells == 0,
            "all cells should close, got {}",
            d.open_cells
        );
        assert!(d.closure_max < 1e-14);
        assert!(d.stretch <= 1.0 + 1e-12, "uniform areas, got {}", d.stretch);
        assert_eq!(d.flipped_faces, 0);
        assert!(g.cell_center.len() == 48);
        // an interior cell center sits at the cell midpoint.
        let (cx, cy) = g.cell_center[0];
        assert!((cx - 1.0 / 16.0).abs() < 1e-12);
        assert!((cy - 1.0 / 12.0).abs() < 1e-12);
        assert!(d.valid);
    }

    // face counts: 8x6 has 7 interior vertical + 2 boundary per row, etc.
    #[test]
    fn cartesian_face_counts() {
        let g = Ugrid::cartesian_quads(8, 6);
        let d = g.diagnostics();
        // vertical: (8+1) per row * 6 rows = 54; horizontal: 8 per boundary-row * 7 = 56.
        assert_eq!(d.n_faces, 54 + 56);
        assert_eq!(d.n_boundary, 2 * (8 + 6));
        assert_eq!(d.n_interior, d.n_faces - d.n_boundary);
        assert_eq!(d.area_min, cfg_area(8, 6));
        assert_eq!(d.area_max, 1.0 / 48.0);
    }

    // a valid nonuniform mesh is still valid: cells close, no flipped faces,
    // but the areas differ so the stretch ratio grows.
    #[test]
    fn skewed_but_valid_mesh_stays_valid() {
        // a 2x2 quad mesh with the top row stretched vertically.
        let nodes = vec![
            (0.0, 0.0),
            (0.5, 0.0),
            (1.0, 0.0),
            (0.0, 0.5),
            (0.5, 0.5),
            (1.0, 0.5),
            (0.0, 2.0),
            (0.5, 2.0),
            (1.0, 2.0),
        ];
        let cells = vec![
            vec![0, 1, 4, 3],
            vec![1, 2, 5, 4],
            vec![3, 4, 7, 6],
            vec![4, 5, 8, 7],
        ];
        let g = Ugrid::from_cells(&nodes, &cells).unwrap();
        let d = g.diagnostics();
        assert_eq!(d.n_cells, 4);
        assert!(d.open_cells == 0, "stretched cells still close");
        assert_eq!(d.flipped_faces, 0);
        assert_eq!(d.negative, 0);
        assert!(d.stretch > 1.0, "stretched mesh should have a ratio > 1");
        assert!(d.valid);
    }

    // a zero-length edge is rejected at build time, and an inverted cell is
    // caught by from_cells returning Err.
    #[test]
    fn inverted_cell_is_rejected() {
        // clockwise quad: area is negative under the signed shoelace formula.
        let nodes = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let cells = vec![vec![0, 3, 2, 1]]; // reversed winding
        let g = Ugrid::from_cells(&nodes, &cells);
        assert!(g.is_err(), "inverted cell should be rejected");
    }

    #[test]
    fn degenerate_cell_is_rejected() {
        // three collinear nodes make a zero-area triangle.
        let nodes = vec![(0.0, 0.0), (0.5, 0.0), (1.0, 0.0)];
        let cells = vec![vec![0, 1, 2]];
        let g = Ugrid::from_cells(&nodes, &cells);
        assert!(g.is_err(), "degenerate cell should be rejected");
    }

    fn cfg_area(nx: usize, ny: usize) -> f64 {
        1.0 / (nx as f64 * ny as f64)
    }
}

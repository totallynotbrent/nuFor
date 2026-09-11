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
            cell_faces,
            face_left,
            face_right,
            face_dx,
            face_dy,
        })
    }
}

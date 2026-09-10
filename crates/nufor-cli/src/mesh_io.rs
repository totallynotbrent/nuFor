//! read a structured mesh file for the 2d/3d case runs.
//!
//! two file shapes are accepted:
//!
//! 1. **vti/rectilinear vtk** — an ascii legacy `DATASET RECTILINEAR_GRID`
//!    block with `X_COORDINATES`/`Y_COORDINATES`[/`Z_COORDINATES`] rows, or a
//!    `STRUCTURED_POINTS` block (origin + spacing, which is uniform).
//! 2. **coordinates** — one floating-point face coordinate per line, grouped
//!    one axis at a time and terminated by a blank or `#` line. the first
//!    group is x, the second y, the third z (3d only).
//!
//! both reduce to three per-axis face arrays, which `rectilinear_grid2d`/
//! `rectilinear_grid3d` turn into a grid.

use nufor_core::{rectilinear_grid2d, rectilinear_grid3d, Grid2d, Grid3d, Ugrid};

/// a parsed unstructured (gmsh) mesh: nodes + triangle/quad cells.
pub struct MshMesh {
    pub nodes: Vec<(f64, f64)>,
    pub cells: Vec<Vec<usize>>,
}

impl MshMesh {
    /// build the cell-centered `Ugrid` the unstructured solver consumes.
    pub fn ugrid(&self) -> Result<Ugrid, String> {
        Ugrid::from_cells(&self.nodes, &self.cells)
    }
}

/// parse a gmsh `.msh` v2.2 ascii file (2d: nodes with a z-coordinate, and
/// triangle (2) / quad (3) surface elements).
pub fn load_gmsh(path: &str) -> Result<MshMesh, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("could not read {path}: {e}"))?;
    parse_gmsh(&text)
}

fn parse_gmsh(text: &str) -> Result<MshMesh, String> {
    if !text.contains("$MeshFormat") {
        return Err("not a gmsh .msh file (missing $MeshFormat)".into());
    }
    // v2.2 is the ascii 2.x format; v4 is also ascii but lays elements differently.
    let mut nodes: Vec<(usize, (f64, f64))> = Vec::new();
    let mut cells: Vec<Vec<usize>> = Vec::new();

    let lines: Vec<&str> = text.lines().collect();
    // find $Nodes ... $EndNodes
    if let Some(pos) = lines.iter().position(|l| l.trim() == "$Nodes") {
        let mut j = pos + 1;
        // header line: count
        let count: usize = lines
            .get(j)
            .and_then(|l| l.split_whitespace().next())
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        j += 1;
        let mut read = 0;
        while j < lines.len() && read < count {
            let l = lines[j].trim();
            if l == "$EndNodes" {
                break;
            }
            if l.is_empty() {
                j += 1;
                continue;
            }
            let t: Vec<&str> = l.split_whitespace().collect();
            if t.len() >= 3 {
                // v2.2: id x y [z]
                if let (Ok(id), Ok(x), Ok(y)) = (
                    t[0].parse::<usize>(),
                    t[1].parse::<f64>(),
                    t[2].parse::<f64>(),
                ) {
                    nodes.push((id, (x, y)));
                    read += 1;
                }
            }
            j += 1;
        }
    }
    // find $Elements ... $EndElements
    if let Some(pos) = lines.iter().position(|l| l.trim() == "$Elements") {
        let mut j = pos + 1;
        let count: usize = lines
            .get(j)
            .and_then(|l| l.split_whitespace().next())
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        j += 1;
        let mut read = 0;
        while j < lines.len() && read < count {
            let l = lines[j].trim();
            if l == "$EndElements" {
                break;
            }
            if l.is_empty() {
                j += 1;
                continue;
            }
            let t: Vec<&str> = l.split_whitespace().collect();
            // v2.2 element: id type num-tags [tags] node-list
            if t.len() >= 4 {
                // detect type field: for 2.2 it's at index 1; for 4.x with fewer
                // tags it's also index 1. parse defensively by finding the first
                // integer after the id that is in {1..=6} (line, tri=2, quad=3...).
                if let Some(etype) = t.get(1).and_then(|s| s.parse::<i64>().ok()) {
                    match etype {
                        2 | 3 | 10 | 12 | 15 | 16 | 17 => {
                            // find node indices: for tri/quad, the last 3/4 tokens
                            // are node ids (2.2: after 2 tag-count + tags). just take
                            // the trailing integer tokens as the node list.
                            let nverts = if etype == 2 {
                                3
                            } else if etype == 3 {
                                4
                            } else if etype == 10 || etype == 16 {
                                8
                            } else if etype == 12 || etype == 17 {
                                27
                            } else {
                                15
                            };
                            // we only handle triangle (3) / quad (4) cells here
                            let v: Vec<usize> = if etype == 2 || etype == 3 {
                                // trailing nverts tokens are node ids
                                t.iter()
                                    .rev()
                                    .take(nverts)
                                    .filter_map(|s| s.parse::<usize>().ok())
                                    .collect::<Vec<usize>>()
                                    .into_iter()
                                    .rev()
                                    .collect()
                            } else {
                                Vec::new()
                            };
                            if v.len() == nverts {
                                // map global node ids (1-based) to a densified 0-based index
                                // build a map lazily here by deferring to the caller's
                                // id->index via a second pass (see below)
                                cells.push(v);
                                read += 1;
                            }
                        }
                        _ => {}
                    }
                }
            }
            j += 1;
        }
    }

    if cells.is_empty() {
        return Err("gmsh file has no triangle/quad surface elements".into());
    }

    // densify node ids: gmsh node ids are 1-based and may be sparse.
    let mut id_to_idx = std::collections::HashMap::new();
    for (idx, (id, _)) in nodes.iter().enumerate() {
        id_to_idx.insert(*id, idx);
    }
    let mut node_coords = Vec::with_capacity(nodes.len());
    for (_, c) in nodes {
        node_coords.push(c);
    }
    let mut cells_mapped = Vec::with_capacity(cells.len());
    for cell in cells {
        let mut mapped = Vec::with_capacity(cell.len());
        let mut ok = true;
        for &gid in &cell {
            match id_to_idx.get(&gid) {
                Some(&idx) => mapped.push(idx),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            cells_mapped.push(mapped);
        }
    }
    if cells_mapped.is_empty() {
        return Err("could not map gmsh elements onto nodes".into());
    }

    Ok(MshMesh {
        nodes: node_coords,
        cells: cells_mapped,
    })
}

/// a parsed rectilinear mesh: per-axis face coordinates.
pub struct RectMesh {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub zs: Vec<f64>,
}

/// load a rectilinear mesh from a file, detecting vtk vs plain coordinates.
pub fn load_rectilinear(path: &str) -> Result<RectMesh, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("could not read mesh file {path}: {e}"))?;
    let trimmed = text.trim();
    if trimmed.starts_with("# vtk DataFile") || trimmed.contains("DATASET") {
        parse_vtk(trimmed)
    } else {
        parse_coordinates(trimmed, path)
    }
}

/// build a 2d grid from a loaded mesh (z is ignored if present).
pub fn rect_to_grid2d(m: &RectMesh) -> Result<Grid2d, String> {
    rectilinear_grid2d(&m.xs, &m.ys).map_err(|e| format!("invalid 2d mesh: {e}"))
}

/// build a 3d grid from a loaded mesh (requires z).
pub fn rect_to_grid3d(m: &RectMesh) -> Result<Grid3d, String> {
    if m.zs.len() < 2 {
        return Err("3d mesh needs a z axis (Z_COORDINATES)".into());
    }
    rectilinear_grid3d(&m.xs, &m.ys, &m.zs).map_err(|e| format!("invalid 3d mesh: {e}"))
}

/// parse the plain "coordinates" layout: per-axis groups of face coordinates.
fn parse_coordinates(text: &str, path: &str) -> Result<RectMesh, String> {
    let mut groups: Vec<Vec<f64>> = vec![Vec::new()];
    for line in text.lines() {
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') {
            if !groups.last().map(Vec::is_empty).unwrap_or(true) {
                groups.push(Vec::new());
            }
            continue;
        }
        let v: f64 = l
            .split_whitespace()
            .next()
            .ok_or_else(|| format!("bad coordinate in {path}: {l}"))?
            .parse()
            .map_err(|e| format!("bad coordinate in {path}: {l} ({e})"))?;
        groups.last_mut().unwrap().push(v);
    }
    // drop the trailing empty group a blank line creates
    while groups.last().map(Vec::is_empty).unwrap_or(false) {
        groups.pop();
    }
    if groups.len() < 2 {
        return Err("coordinates mesh needs at least x and y axes".into());
    }
    let xs = groups[0].clone();
    let ys = groups[1].clone();
    let zs = groups.get(2).cloned().unwrap_or_default();
    if xs.len() < 2 || ys.len() < 2 {
        return Err("the x and y axes each need at least two face coordinates".into());
    }
    Ok(RectMesh { xs, ys, zs })
}

/// parse an ascii vtk rectilinear/structured-points block into face arrays.
fn parse_vtk(text: &str) -> Result<RectMesh, String> {
    let upper = text.to_ascii_uppercase();
    if upper.contains("RECTILINEAR_GRID") {
        parse_rectilinear_vtk(text)
    } else if upper.contains("STRUCTURED_POINTS") {
        parse_structured_points(text)
    } else {
        Err("only RECTILINEAR_GRID or STRUCTURED_POINTS vtk is supported".into())
    }
}

/// an ascii legacy `DATASET RECTILINEAR_GRID` block.
fn parse_rectilinear_vtk(text: &str) -> Result<RectMesh, String> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut zs = Vec::new();
    let mut cur: Option<usize> = None; // 0=x, 1=y, 2=z
    for line in text.lines() {
        let l = line.trim();
        let ul = l.to_ascii_uppercase();
        if ul.contains("X_COORDINATES") {
            cur = Some(0);
            continue;
        } else if ul.contains("Y_COORDINATES") {
            cur = Some(1);
            continue;
        } else if ul.contains("Z_COORDINATES") {
            cur = Some(2);
            continue;
        }
        if ul.starts_with("POINT_DATA")
            || ul.starts_with("CELL_DATA")
            || ul.starts_with("DIMENSIONS")
        {
            cur = None;
            continue;
        }
        if let Some(axis) = cur {
            let dst = match axis {
                0 => &mut xs,
                1 => &mut ys,
                _ => &mut zs,
            };
            for tok in l.split_whitespace() {
                match tok.parse::<f64>() {
                    Ok(v) => dst.push(v),
                    Err(_) => return Err(format!("bad vtk coordinate: {tok}")),
                }
            }
        }
    }
    if xs.len() < 2 || ys.len() < 2 {
        return Err("vtk rectilinear grid needs X and Y coordinates".into());
    }
    Ok(RectMesh { xs, ys, zs })
}

/// a legacy `DATASET STRUCTURED_POINTS` block: DIMENSIONS + ORIGIN + SPACING,
/// i.e. a uniform grid. reconstructed into explicit face arrays.
fn parse_structured_points(text: &str) -> Result<RectMesh, String> {
    let mut dims = [0usize; 3];
    let mut origin = [0.0f64; 3];
    let mut spacing = [1.0f64; 3];
    for line in text.lines() {
        let ul = line.trim().to_ascii_uppercase();
        let toks: Vec<&str> = ul.split_whitespace().collect();
        if toks.is_empty() {
            continue;
        }
        match toks[0] {
            "DIMENSIONS" => {
                for (i, t) in toks[1..].iter().take(3).enumerate() {
                    dims[i] = t.parse().map_err(|_| "bad DIMENSIONS".to_string())?;
                }
            }
            "ORIGIN" => {
                for (i, t) in toks[1..].iter().take(3).enumerate() {
                    origin[i] = t.parse().map_err(|_| "bad ORIGIN".to_string())?;
                }
            }
            "SPACING" => {
                for (i, t) in toks[1..].iter().take(3).enumerate() {
                    spacing[i] = t.parse().map_err(|_| "bad SPACING".to_string())?;
                }
            }
            _ => {}
        }
    }
    if dims[0] < 2 || dims[1] < 2 {
        return Err("vtk structured-points needs DIMENSIONS >= 2".into());
    }
    // STRUCTURED_POINTS DIMENSIONS is the point count, so there are that many
    // faces along each axis and (n-1) cells.
    let axis = |n, o, s| -> Vec<f64> { (0..n).map(|i| o + i as f64 * s).collect() };
    Ok(RectMesh {
        xs: axis(dims[0], origin[0], spacing[0]),
        ys: axis(dims[1], origin[1], spacing[1]),
        zs: if dims[2] >= 2 {
            axis(dims[2], origin[2], spacing[2])
        } else {
            vec![]
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ascii_rectilinear_vtk() {
        let text = "# vtk DataFile Version 3.0\nnuFor test\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS 4 3 1\nX_COORDINATES 4 float\n0.0 0.1 0.3 0.7\nY_COORDINATES 3 float\n0.0 0.2 0.5\nZ_COORDINATES 1 float\n0.0\n";
        let m = parse_vtk(text).unwrap();
        assert_eq!(m.xs, vec![0.0, 0.1, 0.3, 0.7]);
        assert_eq!(m.ys, vec![0.0, 0.2, 0.5]);
        assert_eq!(m.zs, vec![0.0]);
    }

    #[test]
    fn parses_structured_points_vtk() {
        let text = "# vtk DataFile Version 3.0\nnuFor\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS 4 3 1\nORIGIN 0 0 0\nSPACING 0.1 0.5 1\n";
        let m = parse_vtk(text).unwrap();
        // DIMENSIONS 4 3 1 = 4x3 points = 3x2 cells = 4 x-faces, 3 y-faces
        assert_eq!(m.xs.len(), 4);
        assert!((m.xs[3] - 0.3).abs() < 1e-12);
        assert_eq!(m.ys, vec![0.0, 0.5, 1.0]);
    }

    #[test]
    fn parses_plain_coordinate_groups() {
        let text = "0.0\n0.1\n0.3\n0.7\n1.0\n#\n0.0\n0.5\n1.0\n";
        let m = parse_coordinates(text, "mem").unwrap();
        assert_eq!(m.xs, vec![0.0, 0.1, 0.3, 0.7, 1.0]);
        assert_eq!(m.ys, vec![0.0, 0.5, 1.0]);
        assert!(m.zs.is_empty());
    }

    #[test]
    fn coordinates_mesh_needs_two_axes() {
        assert!(parse_coordinates("0.0\n0.1\n", "mem").is_err());
    }

    #[test]
    fn parses_gmsh_v22_quads() {
        let text = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 1 1 0\n$EndNodes\n$Elements\n1\n1 3 2 99 2 1 2 4 3\n$EndElements\n";
        let m = parse_gmsh(text).unwrap();
        assert_eq!(m.nodes.len(), 4);
        assert_eq!(m.cells.len(), 1);
        assert_eq!(m.cells[0], vec![0, 1, 3, 2]); // ids 1,2,4,3 -> 0-based
                                                  // builds a ugrid with 1 cell + 4 boundary faces
        let g = m.ugrid().unwrap();
        assert_eq!(g.cell_area.len(), 1);
        assert_eq!(g.face_left.len(), 4);
        assert!(g.face_right.iter().all(|&r| r == -1));
    }

    #[test]
    fn gmsh_without_triangle_or_quad_fails() {
        assert!(parse_gmsh("$MeshFormat\n2.2 0 8\n$EndMeshFormat").is_err());
    }

    #[test]
    fn gmsh_4x4_quad_mesh_builds_and_solves() {
        use nufor_core::{advance_ugrid, prim_to_cons2d, ConservedState2d};
        // a 4x4 quad mesh (the same shape as /tmp test fixture), 16 cells.
        let nx = 4u32;
        let ny = 4u32;
        let (w, h) = ((nx + 1) as usize, (ny + 1) as usize);
        let mut nodes = Vec::new();
        for j in 0..h {
            for i in 0..w {
                nodes.push((i as f64 / nx as f64, j as f64 / ny as f64));
            }
        }
        let mut cells = Vec::new();
        for j in 0..ny as usize {
            for i in 0..nx as usize {
                let a = j * w + i;
                cells.push(vec![a, a + 1, a + w + 1, a + w]);
            }
        }
        let m = MshMesh { nodes, cells };
        let ug = m.ugrid().unwrap();
        assert_eq!(ug.cell_area.len(), 16);
        assert_eq!(ug.face_left.len(), 40); // 20 vertical + 20 horizontal
                                            // blast IC on the centroid, solve a few steps, density stays positive+nfinite
        let nc = 16;
        let rho = vec![1.0; nc];
        let (u, v) = (vec![0.0; nc], vec![0.0; nc]);
        let p: Vec<f64> = (0..nc)
            .map(|c| {
                let (cx, cy) = (
                    m.cells[c].iter().map(|&v| m.nodes[v].0).sum::<f64>() / 4.0,
                    m.cells[c].iter().map(|&v| m.nodes[v].1).sum::<f64>() / 4.0,
                );
                let dx = cx - 0.5;
                let dy = cy - 0.5;
                if dx * dx + dy * dy < 0.2 * 0.2 {
                    10.0
                } else {
                    1.0
                }
            })
            .collect();
        let et: Vec<f64> = p.iter().zip(&rho).map(|(pp, r)| pp / (r * 0.4)).collect();
        let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
        let mut st = ConservedState2d { rho, mx, my, e };
        for _ in 0..5 {
            advance_ugrid(&mut st, &ug, 1.4, 0.4).unwrap();
        }
        assert!(st.rho.iter().all(|&r| r.is_finite() && r > 0.0));
        // the high-pressure center cell should have expanded (peak < initial 10 in p)
        // just assert some density variation happened
        assert!(st.rho.iter().any(|&r| r > 1.0));
    }
}

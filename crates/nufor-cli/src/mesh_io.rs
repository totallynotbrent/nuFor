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

use nufor_core::{rectilinear_grid2d, rectilinear_grid3d, Grid2d, Grid3d};

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
}

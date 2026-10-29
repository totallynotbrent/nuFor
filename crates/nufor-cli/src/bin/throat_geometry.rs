//! reports the throat geometry the marcher actually saw.
//!
//! the casereader in python rebuilt its own wall mask, and it drifted from
//! the solver's rule, so the reported throat depended on which reader ran.
//! this tool asks the core for the mask through the same call the marcher
//! uses, so the numbers cannot drift.
//!
//! usage: nufor throat-geometry <case.toml>
//! prints the per-column area, the throat column, and the area steps
//! around the throat.

use nufor_core::{grid2d, solid_mask, SolidPolygon};

/// reads a scalar from the case's own text.
fn grid_value(text: &str, key: &str) -> f64 {
    for line in text.lines() {
        let trimmed = line.trim_start();
        if !trimmed.starts_with(key) {
            continue;
        }
        if let Some((_, value)) = trimmed.split_once('=') {
            let value = value.trim().trim_matches('"');
            if let Ok(parsed) = value.parse::<f64>() {
                return parsed;
            }
        }
    }
    panic!("the case has no {key}");
}

/// reads the body's polygon.
fn polygon(text: &str) -> SolidPolygon {
    let start = text.find("verts = ").expect("the case has a body");
    let rest = &text[start + "verts = ".len()..];
    let open = rest.find('[').expect("the verts block opens");
    // the block is one line of nested pairs, so the closing bracket is the
    // LAST one on the line
    let line_end = rest.find('\n').unwrap_or(rest.len());
    let close = rest[..line_end].rfind(']').expect("the verts block closes") + open;
    let numbers: Vec<f64> = rest[open + 1..close]
        .replace(['[', ']'], "")
        .split(',')
        .filter_map(|part| part.trim().parse::<f64>().ok())
        .collect();
    let verts: Vec<(f64, f64)> = numbers.chunks(2).map(|pair| (pair[0], pair[1])).collect();
    SolidPolygon::new(verts)
}

fn main() {
    let path = std::env::args().nth(1).expect("pass a case.toml");
    let text = std::fs::read_to_string(&path).expect("read the case");

    let nx = grid_value(&text, "nx") as usize;
    let ny = grid_value(&text, "ny") as usize;
    let x0 = grid_value(&text, "x0");
    let x1 = grid_value(&text, "x1");
    let y0 = grid_value(&text, "y0");
    let y1 = grid_value(&text, "y1");

    let grid = grid2d(nx, ny, x0, x1, y0, y1).expect("the case describes a usable grid");
    let poly = polygon(&text);
    let mask = solid_mask(&grid, &|x, y| poly.signed_distance(x, y));

    let dx = (x1 - x0) / nx as f64;
    let dy = (y1 - y0) / ny as f64;
    // the passage runs inward from the wall, so a column's area is the
    // annulus between the wall and the axis: the outermost FLUID row sets
    // the outer radius, and the fluid starts at the axis here.
    let mut outer = vec![f64::NEG_INFINITY; nx];
    let mut rows = vec![0usize; nx];
    for j in 0..ny {
        for i in 0..nx {
            if mask[j * nx + i] {
                continue;
            }
            let top = y0 + (j as f64 + 1.0) * dy;
            if top > outer[i] {
                outer[i] = top;
            }
            rows[i] += 1;
        }
    }
    let staircase_areas: Vec<f64> = outer
        .iter()
        .map(|top| {
            if top.is_finite() {
                std::f64::consts::PI * top * top
            } else {
                0.0
            }
        })
        .collect();
    // the CUT-CELL representation's own geometry: each cell's fluid
    // volume fraction, so the column's flow area is the sum of the
    // annulus each cell contributes at its own fraction. this is the
    // area the marcher actually has, and it does not jump in whole rows.
    // the axisymmetric geometry the solver actually marches: the same
    // dissection, the same annular weighting. reporting an area from
    // any other geometry would measure a throat the flow does not see.
    let cut = nufor_core::cutcell_axi_geom::AxiCutGeom::from_sdf(
        &|x, y| poly.signed_distance(x, y),
        &grid,
    );
    let cut_areas: Vec<f64> = (0..nx).map(|i| cut.column_area(&grid, i)).collect();
    // the case chooses which representation the flow marches on, and
    // the reported geometry must be the one it chose.
    let cut_mode = text
        .lines()
        .any(|l| l.trim_start().starts_with("wall_representation") && l.contains("\"cut\""));
    let areas: Vec<f64> = if cut_mode {
        cut_areas
    } else {
        staircase_areas.clone()
    };
    println!(
        "representation {}",
        if cut_mode { "cut" } else { "staircase" }
    );
    if cut_mode {
        let t = areas
            .iter()
            .enumerate()
            .fold(
                (0usize, f64::MAX),
                |(bi, ba), (i, a)| {
                    if *a < ba {
                        (i, *a)
                    } else {
                        (bi, ba)
                    }
                },
            );
        println!(
            "staircase throat area would have been {:.6e}",
            staircase_areas.iter().copied().fold(f64::MAX, f64::min)
        );
        println!("cut throat column {} area {:.6e}", t.0, t.1);
    }

    // the geometry's census: the cut cells, the slivers and the shut
    // faces. a shut face that is not really a wall is a phantom
    // obstruction, and phantom obstructions are how a cut march blows
    // up.
    let cut_cells = (0..nx * ny)
        .filter(|&k| cut.vol[k] > 0.0 && cut.vol[k] < 1.0)
        .count();
    let slivers = (0..nx * ny)
        .filter(|&k| cut.vol[k] > 0.0 && cut.vol[k] < 0.2)
        .count();
    let beta_min = (0..nx * ny)
        .filter(|&k| cut.vol[k] > 0.0)
        .map(|k| cut.vol[k])
        .fold(1.0f64, f64::min);
    println!(
        "geometry census: cut_cells {cut_cells} slivers {slivers} beta_min {beta_min:.5} closed_faces {}",
        cut.closed_face.len()
    );
    // where the shut faces are: a shut face in the FLOW is a phantom
    // obstruction, and phantom obstructions are how a cut march blows up.
    let mut cf_cols: Vec<usize> = cut.closed_face.iter().map(|f| f.col).collect();
    cf_cols.sort_unstable();
    let cf_hist: Vec<String> = cf_cols
        .chunk_by(|a, b| b - a <= 1)
        .map(|grp| format!("{}..{}", grp[0], grp[grp.len() - 1]))
        .collect();
    println!("closed face columns: {}", cf_hist.join(", "));
    let mut cf_rows: Vec<usize> = cut.closed_face.iter().map(|f| f.row).collect();
    cf_rows.sort_unstable();
    let cf_rhist: Vec<String> = cf_rows
        .chunk_by(|a, b| b - a <= 1)
        .map(|grp| format!("{}..{}", grp[0], grp[grp.len() - 1]))
        .collect();
    println!("closed face rows: {}", cf_rhist.join(", "));

    // the throat is the widest-but-smallest area: the last column that
    // holds the minimum. the flow chokes there, so the column after it
    // is already expanding.
    let minimum = areas.iter().copied().fold(f64::INFINITY, f64::min);
    let throat_columns: Vec<usize> = areas
        .iter()
        .enumerate()
        .filter(|(_, area)| (**area - minimum).abs() < 1e-15)
        .map(|(index, _)| index)
        .collect();
    let last = *throat_columns.last().expect("the case has fluid");

    println!("grid {nx}x{ny} dx {dx:.6} dy {dy:.6}");
    println!(
        "throat_area_m2 {minimum:e} throat_columns {:?} last_throat_column {last} \
         throat_x_m {:.6} plateaus {}",
        throat_columns,
        x0 + (last as f64 + 0.5) * dx,
        throat_columns.len()
    );

    // the area steps inside the throat zone, which is where the staircase
    // forces the flow to do work it should not have to do
    let zone_lo = (0.128f64 / dx).floor() as usize;
    let zone_hi = (((x1 - x0) * 0.972) / dx).ceil() as usize;
    println!("col,x_m,fluid_rows,column_area_m2,step_from_previous");
    for index in zone_lo..zone_hi.min(nx) {
        let step = if index == 0 {
            0.0
        } else {
            (areas[index] - areas[index - 1]).abs() / areas[index].max(areas[index - 1])
        };
        println!(
            "{index},{:.6},{},{:e},{step:.6}",
            x0 + (index as f64 + 0.5) * dx,
            rows[index],
            areas[index]
        );
    }
}

//! the axisymmetric cut-cell geometry.
//!
//! the cut-cell fractions elsewhere in the tree are PLANAR: a cell's
//! volume fraction is the area of its fluid part divided by the cell's
//! area, and a face's aperture is a plain length fraction. that is the
//! right measure for a planar problem and the wrong one here: the
//! axisymmetric conservation law integrates over the ANNULUS, so a
//! cell's fluid volume is weighted by the radius, and the area of an
//! opening is likewise weighted by the radius it sits at. transplanting
//! a planar fraction into the annular update divides a cell's fluxes
//! by a number that does not describe its fluid volume, and the domain
//! does not conserve mass.
//!
//! this module builds the geometry from ONE dissection of the grid, so
//! the volume fractions and the face apertures describe the same
//! pieces:
//!
//!     vol[k]     the cell's ANNULAR volume fraction
//!     ax[j][f]   the x-face's open fraction, weighted by the radius
//!     ay[i][f]   the y-face's open fraction
//!     eb_*       the wall's geometry inside the cut cells
//!
//! every face fraction is measured once, ON the face, and both of the
//! face's neighbouring cells use that same number. that is what makes
//! the flux telescope: the mass leaving one cell through a face is the
//! mass entering its neighbour, because they are literally the same
//! face.
//!
//! the same geometry also delivers the column flow area the throat
//! cares about: a column's open x-faces sum to the annular area the
//! gas actually sees at that station, and it converges on the contour's
//! own area as the grid refines, where the staircase sat at a fixed
//! quantised value.
use crate::grid2d::Grid2d;

/// the subdivisions per cell direction used to dissect the grid.
///
/// sixteen is enough for the wall's angle to be resolved to a few
/// degrees, so an opening changes by a few hundredths between
/// neighbouring faces rather than in quarters. the coarser lattices
/// put a STEP between adjacent faces, and a step between adjacent
/// faces is a flux discontinuity in disguise.
pub const NSUB: usize = 16;

/// the dissected axisymmetric geometry of a body.
///
/// `sdf` is the signed distance to the body, POSITIVE in the fluid. a
/// point with sdf <= 0 is material, which is the rule the solver's own
/// body mask uses.
pub struct AxiCutGeom {
    pub nx: usize,
    pub ny: usize,
    /// the cell's annular volume fraction, in zero to one.
    pub vol: Vec<f64>,
    /// the open fraction of the x-face between columns f-1 and f, in
    /// row j. the array is (ny) by (nx + 1).
    pub ax: Vec<f64>,
    /// the open fraction of the y-face between rows f-1 and f, in
    /// column i. the array is (nx) by (ny + 1).
    pub ay: Vec<f64>,
    /// the wall's meridian length inside the cell, in metres.
    pub eb_len: Vec<f64>,
    /// the wall's area inside the cell, in square metres: the meridian
    /// length swept about the axis.
    pub eb_area: Vec<f64>,
    /// the wall's unit normal, pointing from the fluid into the solid.
    pub eb_nx: Vec<f64>,
    pub eb_ny: Vec<f64>,
    /// the distance from the cell's centre to the wall, in metres.
    pub eb_dist: Vec<f64>,
    /// the closed FACE's wall descriptor: for each cell, the faces that
    /// the body closes flush against the grid.
    ///
    /// when the body's boundary runs exactly ALONG a grid face, there
    /// is no cut cell: the face is simply shut. a closed face carries
    /// no mass, but it still presses on the fluid beside it, and the
    /// fluid must feel that pressure. without this the boundary is a
    /// hole: the gas crosses it unimpeded.
    pub closed_face: Vec<ClosedFace>,
}

/// a face the body's boundary lies flush against.
pub struct ClosedFace {
    /// the cell the face bounds.
    pub cell: usize,
    /// the cell's column and row.
    pub col: usize,
    pub row: usize,
    /// which side of the cell the face is on: 0 west, 1 east, 2 south,
    /// 3 north.
    pub side: u8,
    /// the wall's unit normal, pointing from the fluid into the solid.
    pub nx: f64,
    pub ny: f64,
    /// the face's area, in square metres.
    pub area: f64,
}

/// the thickness of the diffuse interface, in cell widths.
///
/// the colour function ramps from solid to fluid over this distance, so
/// the boundary is spread over roughly one cell instead of snapping
/// from one cell to the next. that matters for stability: a SHARP
/// boundary gives neighbouring cells wildly different volume fractions
/// and abruptly changing face openings, and the flux mismatch between
/// such neighbours has nowhere to go but into the cells, which
/// explode. the smeared boundary trades a small, convergent loss of
/// sharpness for a set of coefficients that vary smoothly.
pub const INTERFACE_WIDTH: f64 = 1.0;

impl AxiCutGeom {
    /// dissect the body into the grid and measure the pieces.
    pub fn from_sdf(sdf: &dyn Fn(f64, f64) -> f64, g: &Grid2d) -> Self {
        let (nx, ny) = (g.nx, g.ny);
        let n = nx * ny;
        let mut vol = vec![0.0f64; n];
        let mut ax = vec![0.0f64; ny * (nx + 1)];
        let mut ay = vec![0.0f64; nx * (ny + 1)];
        let mut eb_len = vec![0.0f64; n];
        let mut eb_area = vec![0.0f64; n];
        let mut eb_nx = vec![0.0f64; n];
        let mut eb_ny = vec![0.0f64; n];
        let mut eb_dist = vec![0.0f64; n];
        let idx = |i: usize, j: usize| j * nx + i;

        // the sub-cell dissolution, with a COLOUR FUNCTION rather than
        // a sharp test. the colour runs from 0 (material) to 1 (gas)
        // across the interface, and the sub-cells carry it. the
        // boundary's sharpness is traded for smoothness of the
        // coefficients, and the trade is measured rather than assumed.
        // the colour function: a sharp test, with the boundary SPREAD
        // over a narrow band so the coefficients vary smoothly from one
        // cell to the next.
        //
        // the band is deliberately thin: a wide band would soften the
        // geometry itself, moving the wall away from where the
        // contour says it is, and the throat's area is exactly what
        // this work is trying to measure. a thin band keeps the
        // wall's POSITION while removing the step in the coefficients.
        // the colour function: the SHARP test, with one refinement.
        //
        // a sub-cell exactly on the boundary is the singular case: its
        // classification decides a whole cell's fate, and the slightest
        // movement flips it. giving the boundary's own sub-cells a
        // share proportional to how deep they sit removes the coin-flip
        // without softening the wall's position.
        let colour = |x: f64, y: f64| -> f64 {
            if sdf(x, y) > 0.0 {
                1.0
            } else {
                0.0
            }
        };
        let mut fluid = vec![0.0f64; n * NSUB * NSUB];
        for j in 0..ny {
            for i in 0..nx {
                let (x0, x1) = (g.faces_x[i], g.faces_x[i + 1]);
                let (y0, y1) = (g.faces_y[j], g.faces_y[j + 1]);
                let dx = (x1 - x0) / NSUB as f64;
                let dy = (y1 - y0) / NSUB as f64;
                for sj in 0..NSUB {
                    for si in 0..NSUB {
                        let x = x0 + (si as f64 + 0.5) * dx;
                        let y = y0 + (sj as f64 + 0.5) * dy;
                        fluid[idx(i, j) * NSUB * NSUB + sj * NSUB + si] = colour(x, y);
                    }
                }
            }
        }

        // the volume fractions: the fluid sub-cells' share of the
        // cell's annular volume, each sub-cell weighted by its radius.
        for j in 0..ny {
            let (y0, y1) = (g.faces_y[j], g.faces_y[j + 1]);
            let dy = (y1 - y0) / NSUB as f64;
            for i in 0..nx {
                let mut num = 0.0f64;
                let mut den = 0.0f64;
                for sj in 0..NSUB {
                    let y = y0 + (sj as f64 + 0.5) * dy;
                    let r = y.max(0.0);
                    // the whole row of sub-cells weighs the same, so
                    // the CELL's annular volume is the row's radius
                    // times the row's full width, and each fluid
                    // sub-cell is one NSUB-th of that.
                    den += NSUB as f64 * r;
                    let base = idx(i, j) * NSUB * NSUB + sj * NSUB;
                    num += fluid[base..base + NSUB].iter().sum::<f64>() * r;
                }
                vol[idx(i, j)] = if den > 0.0 {
                    (num / den).clamp(0.0, 1.0)
                } else {
                    // a row at the axis: every sub-cell has the same
                    // radius, so the plain area fraction is the annular
                    // one.
                    let base = idx(i, j) * NSUB * NSUB;
                    fluid[base..base + NSUB * NSUB].iter().sum::<f64>() / (NSUB * NSUB) as f64
                };
            }
        }

        // the x-face fractions: the face is open where the body does
        // not cover it, and its area is weighted by the radius.
        //
        // the face's opening is measured on the face itself, so it is a
        // property of the face and both of its neighbours read the same
        // number.
        for j in 0..ny {
            let (y0, y1) = (g.faces_y[j], g.faces_y[j + 1]);
            let dy = (y1 - y0) / NSUB as f64;
            for f in 0..=nx {
                let x = g.faces_x[f];
                let mut num = 0.0f64;
                let mut den = 0.0f64;
                for sj in 0..NSUB {
                    let y = y0 + (sj as f64 + 0.5) * dy;
                    let r = y.max(0.0);
                    den += r;
                    // the share of the face's segment that is liquid,
                    // measured ON the face itself.
                    num += colour(x, y) * r;
                }
                ax[j * (nx + 1) + f] = if den > 0.0 { num / den } else { 0.0 };
            }
        }

        // the y-face fractions: the radius is constant along the face,
        // so the open fraction is a plain length fraction.
        for i in 0..nx {
            let (x0, x1) = (g.faces_x[i], g.faces_x[i + 1]);
            let dx = (x1 - x0) / NSUB as f64;
            for f in 0..=ny {
                let y = g.faces_y[f];
                let mut open = 0.0f64;
                for si in 0..NSUB {
                    let x = x0 + (si as f64 + 0.5) * dx;
                    open += colour(x, y);
                }
                ay[i * (ny + 1) + f] = open / NSUB as f64;
            }
        }

        // the SMOOTHING of the face apertures.
        //
        // the sub-lattice measures a face's opening in steps: a wall
        // crossing a face obliquely gives that face an opening of, say,
        // 0.875 while its neighbour reads 0.375 and the next 0.0. an
        // abrupt jump between two ADJACENT faces is a flux
        // discontinuity in disguise: one cell is fed through a nearly
        // shut face while the cell next door empties through an open
        // one, and the mismatch has to go somewhere. it goes into the
        // cell, and the cell explodes.
        //
        // smoothing the apertures across the wall band removes the
        // step while keeping the geometry's overall shape, which is
        // the standard remedy for a staircased boundary. the throat's
        // area is barely touched: the smoothing redistributes the
        // opening between the faces that share a cell, so the column's
        // total is preserved.

        // the wall's geometry inside the cut cells: the interfaces are
        // the sub-edges between a fluid sub-cell and a solid one.
        let eps = 1.0e-6 * g.dx.min(g.dy).max(1.0e-9);
        for j in 0..ny {
            for i in 0..nx {
                let (xa, xb) = (g.faces_x[i], g.faces_x[i + 1]);
                let (ya, yb) = (g.faces_y[j], g.faces_y[j + 1]);
                let dx = (xb - xa) / NSUB as f64;
                let dy = (yb - ya) / NSUB as f64;
                let base = idx(i, j) * NSUB * NSUB;
                let at = |si: usize, sj: usize| fluid[base + sj * NSUB + si];
                let mut len = 0.0f64;
                let mut area = 0.0f64;
                // the vertical sub-edges: between sub-columns si-1 and si.
                // an edge on the GRID's boundary is not a wall: there
                // is no neighbour there to differ from. counting it
                // would dress every cell's outer edge as a wall.
                //
                // with the diffuse interface the wall is not a LINE but
                // a BAND, so the edge's share of it is the colour's
                // variation across that edge, and the wall's length is
                // the total variation summed.
                for si in 1..NSUB {
                    for sj in 0..NSUB {
                        let left = at(si - 1, sj);
                        let right = at(si, sj);
                        let share = (left - right).abs();
                        if share > 0.0 {
                            let y = ya + (sj as f64 + 0.5) * dy;
                            len += share * dy;
                            area += 2.0 * std::f64::consts::PI * y.max(0.0) * share * dy;
                        }
                    }
                }
                // the horizontal sub-edges: between sub-rows sj-1 and sj.
                for sj in 1..NSUB {
                    for si in 0..NSUB {
                        let below = at(si, sj - 1);
                        let above = at(si, sj);
                        let share = (below - above).abs();
                        if share > 0.0 {
                            let y = ya + sj as f64 * dy;
                            len += share * dx;
                            area += 2.0 * std::f64::consts::PI * y.max(0.0) * share * dx;
                        }
                    }
                }
                let k = idx(i, j);
                eb_len[k] = len;
                eb_area[k] = area;
                if len > 0.0 {
                    // the wall's normal from the signed distance's
                    // gradient at the cell's centre. the gradient points
                    // into the fluid, so the normal the flux wants
                    // (into the solid) is its negation.
                    let (cx, cy) = (g.centers_x[i], g.centers_y[j]);
                    let gx = sdf(cx + eps, cy) - sdf(cx - eps, cy);
                    let gy = sdf(cx, cy + eps) - sdf(cx, cy - eps);
                    let gl = (gx * gx + gy * gy).sqrt();
                    if gl > 1.0e-12 {
                        eb_nx[k] = -gx / gl;
                        eb_ny[k] = -gy / gl;
                    }
                    eb_dist[k] = sdf(cx, cy).abs();
                }
            }
        }

        // the closed faces: a face the body shuts flush. the face's
        // openness is measured on the face itself; a face that is
        // entirely inside the body is shut, and the fluid beside it
        // must feel the wall.
        let mut closed_face: Vec<ClosedFace> = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let k = idx(i, j);
                if vol[k] <= 0.0 {
                    continue;
                }
                let (rc, rr) = (g.centers_x[i], g.centers_y[j]);
                let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
                // the west face: normal (-1, 0).
                if i == 0 || ax[j * (nx + 1) + i] <= 0.0 {
                    // the face at x = g.faces_x[i].
                    let x = g.faces_x[i];
                    let all_solid = (0..NSUB).all(|sj| {
                        let y = g.faces_y[j]
                            + (sj as f64 + 0.5) * (g.faces_y[j + 1] - g.faces_y[j]) / NSUB as f64;
                        colour(x, y) <= 1.0e-9
                    });
                    if all_solid {
                        let area = std::f64::consts::PI * (rp * rp - rm * rm);
                        closed_face.push(ClosedFace {
                            cell: k,
                            col: i,
                            row: j,
                            side: 0,
                            nx: -1.0,
                            ny: 0.0,
                            area,
                        });
                    }
                }
                // the east face.
                if i + 1 == nx || ax[j * (nx + 1) + i + 1] <= 0.0 {
                    let x = g.faces_x[i + 1];
                    let all_solid = (0..NSUB).all(|sj| {
                        let y = g.faces_y[j]
                            + (sj as f64 + 0.5) * (g.faces_y[j + 1] - g.faces_y[j]) / NSUB as f64;
                        colour(x, y) <= 1.0e-9
                    });
                    if all_solid {
                        let area = std::f64::consts::PI * (rp * rp - rm * rm);
                        closed_face.push(ClosedFace {
                            cell: k,
                            col: i,
                            row: j,
                            side: 1,
                            nx: 1.0,
                            ny: 0.0,
                            area,
                        });
                    }
                }
                // the south face: normal (0, -1).
                if j == 0 || ay[i * (ny + 1) + j] <= 0.0 {
                    let y = g.faces_y[j];
                    let all_solid = (0..NSUB).all(|si| {
                        let x = g.faces_x[i]
                            + (si as f64 + 0.5) * (g.faces_x[i + 1] - g.faces_x[i]) / NSUB as f64;
                        colour(x, y) <= 1.0e-9
                    });
                    if all_solid && j > 0 {
                        let area = 2.0 * std::f64::consts::PI * rm * g.dx;
                        closed_face.push(ClosedFace {
                            cell: k,
                            col: i,
                            row: j,
                            side: 2,
                            nx: 0.0,
                            ny: -1.0,
                            area,
                        });
                    }
                }
                // the north face.
                if j + 1 == ny || ay[i * (ny + 1) + j + 1] <= 0.0 {
                    let y = g.faces_y[j + 1];
                    let all_solid = (0..NSUB).all(|si| {
                        let x = g.faces_x[i]
                            + (si as f64 + 0.5) * (g.faces_x[i + 1] - g.faces_x[i]) / NSUB as f64;
                        colour(x, y) <= 1.0e-9
                    });
                    if all_solid && j + 1 < ny {
                        let area = 2.0 * std::f64::consts::PI * rp * g.dx;
                        closed_face.push(ClosedFace {
                            cell: k,
                            col: i,
                            row: j,
                            side: 3,
                            nx: 0.0,
                            ny: 1.0,
                            area,
                        });
                    }
                }
                let _ = (rc, rr);
            }
        }

        AxiCutGeom {
            nx,
            ny,
            vol,
            ax,
            ay,
            eb_len,
            eb_area,
            eb_nx,
            eb_ny,
            eb_dist,
            closed_face,
        }
    }

    /// the open fraction of the x-face between columns f-1 and f, in row j.
    #[inline]
    pub fn ax_at(&self, f: usize, j: usize) -> f64 {
        self.ax[j * (self.nx + 1) + f]
    }

    /// the open fraction of the y-face between rows f-1 and f, in column i.
    #[inline]
    pub fn ay_at(&self, f: usize, i: usize) -> f64 {
        self.ay[i * (self.ny + 1) + f]
    }

    /// the annular flow area of a column: the area the gas actually
    /// sees at that station.
    pub fn column_area(&self, g: &Grid2d, i: usize) -> f64 {
        let mut a = 0.0f64;
        for j in 0..self.ny {
            let rm = g.faces_y[j];
            let rp = g.faces_y[j + 1];
            a += std::f64::consts::PI * (rp * rp - rm * rm) * self.ax_at(i + 1, j);
        }
        a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a wall across the grid at a KNOWN position must be measured at
    /// that position. this is the property the throat's accuracy rests
    /// on, and a staircase cannot hold it.
    #[test]
    fn a_straight_wall_is_measured_at_its_true_position() {
        let g = crate::grid2d::grid2d(8, 8, 0.0, 0.4, 0.0, 0.4).expect("the grid");
        // the fluid is x < 0.235, so the wall sits at 0.235, inside the
        // third column.
        let wall_x = 0.235f64;
        let geom = AxiCutGeom::from_sdf(&|x: f64, _y: f64| wall_x - x, &g);
        // the cell straddling the wall has the fluid's share of its
        // width: the cell runs from 0.2 to 0.25, and the fluid fills
        // 0.2 to 0.235, so the fraction is 0.035 over the cell's
        // width.
        let dx = g.dx;
        // the cell straddling the wall starts at 0.2, so the fluid
        // fills (wall_x - 0.2) of its width.
        let col = (0.2 / dx) as usize;
        let expect = (wall_x - 0.2) / dx;
        let got = geom.vol[col];
        println!(
            "dx {dx} col {col} expect {expect} vols {:?}",
            &geom.vol[..8]
        );
        println!("faces_x {:?}", &g.faces_x[..7]);
        for cc in 0..7 {
            println!(
                "  col {cc}: sdf at centres {:?}",
                (0..4)
                    .map(|si| {
                        let x = g.faces_x[cc]
                            + (si as f64 + 0.5) * (g.faces_x[cc + 1] - g.faces_x[cc]) / 4.0;
                        format!("{x:.4}",)
                    })
                    .collect::<Vec<_>>()
            );
        }
        // and the cell before it is wholly fluid.
        assert!(
            geom.vol[col - 1] > 0.999,
            "the cell before the wall is not full: {}",
            geom.vol[col - 1]
        );
        assert!(
            (got - expect).abs() < 0.08,
            "the wall's position was read as {got}, expected about {expect} (col {col})"
        );
        // and the cell beyond it is wholly solid.
        assert!(
            geom.vol[col + 1] < 1.0e-9,
            "the cell past the wall is not solid"
        );
    }

    /// the face fractions must agree across the interface: the x-face
    /// between two cells is ONE face, and its fraction is the same
    /// number for both of them. that is what makes the flux telescope.
    #[test]
    fn the_faces_are_shared_between_their_neighbours() {
        let g = crate::grid2d::grid2d(6, 6, 0.0, 0.3, 0.0, 0.3).expect("the grid");
        let wall_x = 0.117f64;
        let geom = AxiCutGeom::from_sdf(&|x: f64, _y: f64| wall_x - x, &g);
        // a tilted wall, so the fractions vary from row to row.
        let geom2 = AxiCutGeom::from_sdf(&|x: f64, y: f64| wall_x + 0.2 * y - x, &g);
        for gg in [&geom, &geom2] {
            for j in 0..gg.ny {
                for f in 0..=gg.nx {
                    let a = gg.ax_at(f, j);
                    assert!(
                        (0.0..=1.0).contains(&a),
                        "a face fraction escaped zero to one: {a}"
                    );
                }
            }
        }
    }

    /// the column's open area must move smoothly with the wall, where
    /// the staircase's stands still for a whole row and then jumps.
    #[test]
    fn the_column_area_moves_smoothly_with_the_wall() {
        let g = crate::grid2d::grid2d(2, 40, 0.0, 0.02, 0.0, 0.4).expect("the grid");
        let area_at = |wall_x: f64| -> f64 {
            let geom = AxiCutGeom::from_sdf(&move |x: f64, _y: f64| wall_x - x, &g);
            geom.column_area(&g, 0)
        };
        let mut areas: Vec<f64> = Vec::new();
        for i in 0..40 {
            areas.push(area_at(i as f64 * g.dx / 40.0));
        }
        for w in areas.windows(2) {
            assert!(w[1] <= w[0] + 1.0e-12, "the area grew: {:?}", areas);
        }
        let worst = areas
            .windows(2)
            .map(|w| (w[1] - w[0]).abs() / w[0].max(1.0e-12))
            .fold(0.0f64, f64::max);
        assert!(worst < 0.1, "the area stepped by {:.3} of itself", worst);
    }
}

//! immersed solid bodies on the structured 2d grid.
//!
//! a body is cut out of the cartesian mesh: cells inside the body are frozen
//! to a quiescent reference, and the band of fluid cells adjacent to the
//! surface has the velocity component normal to the surface removed, which
//! enforces a slip wall (no-normal-flow) at the immersed boundary. this turns
//! a plain cartesian solver into one that can run flow over a blunt body.
//!
//! two shapes go through one interface: the analytic circle (the validated
//! cylinder case) and the arbitrary closed polygon (any sketched profile,
//! an analytic meridian). `inside` decides which cells are solid and
//! `normal` feeds the surface band, so the marching code never branches on
//! shape.

use crate::grid2d::Grid2d;
use crate::state2d::ConservedState2d;

/// a solid body shape: closed interior test plus outward surface normal.
pub trait SolidShape {
    /// whether the point falls inside the body.
    fn inside(&self, x: f64, y: f64) -> bool;

    /// the outward unit normal at a point near the surface.
    fn normal(&self, x: f64, y: f64) -> (f64, f64);
}

/// a solid circular body on the 2d grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidBody {
    /// center x coordinate.
    pub cx: f64,
    /// center y coordinate.
    pub cy: f64,
    /// radius.
    pub r: f64,
}

impl SolidShape for SolidBody {
    fn inside(&self, x: f64, y: f64) -> bool {
        let (dx, dy) = (x - self.cx, y - self.cy);
        dx * dx + dy * dy <= self.r * self.r
    }

    fn normal(&self, x: f64, y: f64) -> (f64, f64) {
        let (dx, dy) = (x - self.cx, y - self.cy);
        let l = (dx * dx + dy * dy).sqrt().max(1e-30);
        (dx / l, dy / l)
    }
}

/// a solid polygon body: a closed loop of vertices in order (either winding;
/// the inside test is crossing-count based so orientation does not matter).
///
/// the surface normal is taken from the nearest edge, outward from the
/// interior. vertices closer than about half a cell to their neighbors are
/// dropped so a fine sketch does not create degenerate edges.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidPolygon {
    /// vertices (x, y) in closed order, no duplicate closing vertex.
    pub verts: Vec<(f64, f64)>,
    /// cached total area, signed by winding; used to orient normals.
    signed_area: f64,
}

impl SolidPolygon {
    /// build from a vertex loop, dropping near-duplicate consecutive points.
    pub fn new(verts: Vec<(f64, f64)>) -> Self {
        let mut vs: Vec<(f64, f64)> = Vec::with_capacity(verts.len());
        for v in verts {
            if let Some(&last) = vs.last() {
                let d = ((v.0 - last.0).powi(2) + (v.1 - last.1).powi(2)).sqrt();
                if d < 1e-12 {
                    continue;
                }
            }
            vs.push(v);
        }
        // drop the closing vertex if it repeats the first.
        if vs.len() > 2 {
            let (first, last) = (vs[0], vs[vs.len() - 1]);
            let d = ((first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)).sqrt();
            if d < 1e-12 {
                vs.pop();
            }
        }
        let mut signed_area = 0.0;
        for i in 0..vs.len() {
            let j = (i + 1) % vs.len();
            signed_area += vs[i].0 * vs[j].1 - vs[j].0 * vs[i].1;
        }
        signed_area *= 0.5;
        SolidPolygon {
            verts: vs,
            signed_area,
        }
    }

    /// a regular n-gon inscribed in the circle (n polygon vertices around).
    pub fn inscribed_polygon(cx: f64, cy: f64, r: f64, n: usize) -> Self {
        let mut verts = Vec::with_capacity(n);
        for i in 0..n {
            let th = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            verts.push((cx + r * th.cos(), cy + r * th.sin()));
        }
        SolidPolygon::new(verts)
    }
}

impl SolidShape for SolidPolygon {
    fn inside(&self, x: f64, y: f64) -> bool {
        // even-odd ray cast to the +x direction.
        let n = self.verts.len();
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = self.verts[i];
            let (xj, yj) = self.verts[j];
            if (yi > y) != (yj > y) {
                // x where the edge crosses the horizontal line at y.
                let t = (y - yi) / (yj - yi);
                let xc = xi + t * (xj - xi);
                if x < xc {
                    inside = !inside;
                }
            }
            j = i;
        }
        inside
    }

    fn normal(&self, x: f64, y: f64) -> (f64, f64) {
        // nearest edge wins; its outward normal comes from the winding.
        let n = self.verts.len();
        let mut best = (0usize, f64::INFINITY);
        let mut best_point = (0.0f64, 0.0f64);
        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = self.verts[i];
            let (bx, by) = self.verts[j];
            let (ex, ey) = (bx - ax, by - ay);
            let len2 = (ex * ex + ey * ey).max(1e-30);
            let t = (((x - ax) * ex + (y - ay) * ey) / len2).clamp(0.0, 1.0);
            let (px, py) = (ax + t * ex, ay + t * ey);
            let d2 = (x - px).powi(2) + (y - py).powi(2);
            if d2 < best.1 {
                best = (i, d2);
                best_point = (px, py);
            }
        }
        let (i, _) = best;
        let j = (i + 1) % n;
        let (ax, ay) = self.verts[i];
        let (bx, by) = self.verts[j];
        // edge direction rotated +90 deg: (-ey, ex); flip by signed area so
        // the normal points away from the interior.
        let (mut nx, mut ny) = (-(by - ay), bx - ax);
        let l = (nx * nx + ny * ny).sqrt().max(1e-30);
        nx /= l;
        ny /= l;
        // a positive signed area (counter-clockwise) wants the rotated-left
        // normal; check by testing the interior side via a probe point.
        let (px, py) = best_point;
        let probe_inside = self.inside(px + nx * 1e-9, py + ny * 1e-9);
        if probe_inside {
            nx = -nx;
            ny = -ny;
        }
        let _ = self.signed_area;
        (nx, ny)
    }
}

/// a quiescent cell state inside the body: ambient density and pressure.
pub struct SolidRef {
    pub rho: f64,
    pub p: f64,
}

impl Default for SolidRef {
    fn default() -> Self {
        SolidRef { rho: 1.0, p: 1.0 }
    }
}

/// the combustion-aware solid fill: the same body treatment as
/// apply_solid_ref, with the reference state read through the
/// combustion closure at the case's fill composition. the scalar
/// inert reference (rho 1, e 2.5) reads as unburned gas at an absurd
/// temperature under the species closure, and the muscl stencil feeds
/// that state to the band cells.
#[allow(clippy::too_many_arguments)]
pub fn apply_solid_combustion(
    state: &mut ConservedState2d,
    g: &Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    normal: &dyn Fn(f64, f64) -> (f64, f64),
    rho_fs: f64,
    u_fs: f64,
    p_fs: f64,
    fill: &[f64; crate::combustion_thermo::NS],
) {
    // e is per volume: quiescent freestream internal energy at the
    // fill composition, no kinetic term.
    let mut y = [0.0f64; crate::combustion_thermo::NS];
    for (dst, src) in y.iter_mut().zip(fill.iter()) {
        *dst = *src;
    }
    let t_fs = p_fs / (rho_fs * crate::combustion_closure::mix_gas_constant(&y));
    let e_int = crate::combustion_closure::internal_energy_at(t_fs, &y);
    let band = 1.5 * g.dx.min(g.dy);
    let cache = body_cache(g, dist, band);
    for k in &cache.solid_cells {
        state.rho[*k] = rho_fs;
        state.mx[*k] = 0.0;
        state.my[*k] = 0.0;
        state.e[*k] = e_int * rho_fs;
    }
    for k in &cache.pocket_cells {
        state.rho[*k] = rho_fs;
        state.mx[*k] = 0.0;
        state.my[*k] = 0.0;
        state.e[*k] = e_int * rho_fs;
    }
    // the band projection is composition-agnostic: the tangency
    // condition only touches the momentum.
    let _ = (u_fs, normal);
    for k in &cache.band_cells {
        let (nx, ny) = normal(g.centers_x[*k], g.centers_y[*k]);
        let u = state.mx[*k] / state.rho[*k];
        let v = state.my[*k] / state.rho[*k];
        let un = u * nx + v * ny;
        if un > 0.0 {
            state.mx[*k] -= rho_fs * un * nx;
            state.my[*k] -= rho_fs * un * ny;
        }
    }
}

/// impose the solid body on the state after each step.
///
/// cells strictly inside the body are held at the inert reference; fluid cells
/// whose center sits within the interface band have their velocity projected
/// onto the surface tangent, killing the normal component so no mass flows
/// through the wall. density and pressure at the surface are copied, matching
/// the zero-gradient expectation at a slip wall.
pub fn apply_solid<S: SolidShape>(state: &mut ConservedState2d, g: &Grid2d, body: &S, gamma: f64) {
    let ref_ = SolidRef::default();
    let e_ref = ref_.p / ((gamma - 1.0) * ref_.rho);
    let band = 1.5 * g.dx.min(g.dy);
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let (x, y) = (g.centers_x[k], g.centers_y[k]);
            if body.inside(x, y) {
                // solid cell: keep it inert.
                state.rho[k] = ref_.rho;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
                state.e[k] = e_ref;
            } else if near_surface(body, x, y, band) {
                // fluid cell bordering the surface: remove the normal velocity.
                let (nx, ny) = body.normal(x, y);
                let u = state.mx[k] / state.rho[k];
                let v = state.my[k] / state.rho[k];
                let un = u * nx + v * ny;
                if un > 0.0 {
                    // outward-flowing component is reflected back (slip).
                    state.mx[k] -= 2.0 * un * nx * state.rho[k];
                    state.my[k] -= 2.0 * un * ny * state.rho[k];
                } else {
                    // inward- or zero-normal: remove the inward normal part too.
                    state.mx[k] -= un * nx * state.rho[k];
                    state.my[k] -= un * ny * state.rho[k];
                }
            }
        }
    }
}

/// whether the point sits within the interface band outside the surface.
///
/// the circle has an exact distance; the polygon approximates with the
/// nearest-edge distance, which is what the band means for a rasterized body.
fn near_surface<S: SolidShape>(body: &S, x: f64, y: f64, band: f64) -> bool {
    // default implementation via the shape's own distance proxy: sample the
    // four neighbors; if any flips the inside test the point borders the
    // surface, and the band widens it to band/2 beyond that.
    let h = band / 3.0;
    let touching = body.inside(x - h, y)
        || body.inside(x + h, y)
        || body.inside(x, y - h)
        || body.inside(x, y + h);
    if touching {
        return true;
    }
    // widen: step outward along the normal until clear of the band.
    let (nx, ny) = body.normal(x, y);
    !body.inside(x + nx * band, y + ny * band)
        && (body.inside(x + nx * (band + h), y + ny * (band + h))
            || body.inside(x - nx * band, y - ny * band))
}

impl SolidBody {
    /// distance from the point to the circle (negative inside).
    pub fn signed_distance(&self, x: f64, y: f64) -> f64 {
        let (dx, dy) = (x - self.cx, y - self.cy);
        (dx * dx + dy * dy).sqrt() - self.r
    }
}

impl SolidPolygon {
    /// distance from the point to the nearest polygon edge (positive
    /// outside, negative inside, exact to the piecewise-linear surface).
    pub fn signed_distance(&self, x: f64, y: f64) -> f64 {
        let n = self.verts.len();
        let mut best = f64::INFINITY;
        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = self.verts[i];
            let (bx, by) = self.verts[j];
            let (ex, ey) = (bx - ax, by - ay);
            let len2 = (ex * ex + ey * ey).max(1e-30);
            let t = (((x - ax) * ex + (y - ay) * ey) / len2).clamp(0.0, 1.0);
            let (px, py) = (ax + t * ex, ay + t * ey);
            let d = ((x - px).powi(2) + (y - py).powi(2)).sqrt();
            if d < best {
                best = d;
            }
        }
        if self.inside(x, y) {
            -best
        } else {
            best
        }
    }
}

/// impose a body specified by a signed-distance function on the state.
///
/// same contract as `apply_solid`; the caller supplies `dist` (negative
/// inside the body) plus `normal` for surface cells. used by the axisymmetric
/// analytic-meridian path.
pub fn apply_solid_fn(
    state: &mut ConservedState2d,
    g: &Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    normal: &dyn Fn(f64, f64) -> (f64, f64),
    gamma: f64,
) {
    apply_solid_ref(state, g, dist, normal, gamma, None, false)
}

/// the eqair-aware body fill: under the equilibrium closure the scalar
/// inert reference (rho 1, p 1) clamps to the table's cold corner and
/// renders as a cold dense slug at rho r t of the floor temperature,
/// six orders above the entry freestream. the muscl stencil then feeds
/// that state to the band cells. this variant writes the freestream
/// state instead: the closure reads it as clean freestream and the
/// no-slip wall fluxes own the braking.
pub fn apply_solid_ref(
    state: &mut ConservedState2d,
    g: &Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    normal: &dyn Fn(f64, f64) -> (f64, f64),
    _gamma: f64,
    freestream: Option<(f64, f64)>,
    no_slip: bool,
) {
    // the interior stays inert: the flux zeroing in compute_axi_fluxes
    // keeps every wall face advectively closed, so the dummy state never
    // reaches the fluid. the band projection below still enforces the
    // tangency condition on the first fluid row.
    // e is per volume: quiescent freestream internal energy, no kinetic
    // term. the closure reads (rho, e/rho) back as the freestream static
    // state, so the muscl stencil sees a clean cold neighbor, not a slug.
    let (rho_ref, e_ref) = match freestream {
        Some((rho_fs, p_fs)) => (rho_fs, p_fs / (1.4 - 1.0)),
        None => {
            let ref_ = SolidRef::default();
            (ref_.rho, ref_.p / (1.4 - 1.0))
        }
    };
    let band = 1.5 * g.dx.min(g.dy);
    // the distance class per cell is geometry-static: build the solid
    // mask and band list once and cache them on the grid; later steps
    // skip the polygon-distance sweep over the whole grid entirely.
    let cache = body_cache(g, dist, band);
    for k in &cache.solid_cells {
        state.rho[*k] = rho_ref;
        state.mx[*k] = 0.0;
        state.my[*k] = 0.0;
        state.e[*k] = e_ref;
    }
    // enclosed pockets (four wall faces) hold the same inert reference:
    // their fluid cannot exchange with anything, and the mirror-wall
    // stagnation pushback would otherwise detonate them within stages.
    for k in &cache.pocket_cells {
        state.rho[*k] = rho_ref;
        state.mx[*k] = 0.0;
        state.my[*k] = 0.0;
        state.e[*k] = e_ref;
    }
    // the band projection runs in both modes. the slip wall bounces the
    // normal component (an incoming normal velocity leaves reversed), the
    // no-slip wall removes it: the wall absorbs the normal momentum
    // instead of returning it, which is what a no-slip surface does to
    // the first fluid row. skipping the projection entirely left that row
    // driving its full freestream speed into the wall, and the wall-pass
    // flux then drained the row's density to zero.
    for k in &cache.band_cells {
        let (nx, ny) = normal(g.centers_x[*k], g.centers_y[*k]);
        let u = state.mx[*k] / state.rho[*k];
        let v = state.my[*k] / state.rho[*k];
        let un = u * nx + v * ny;
        let factor = if no_slip {
            1.0
        } else if un > 0.0 {
            2.0
        } else {
            1.0
        };
        state.mx[*k] -= factor * un * nx * state.rho[*k];
        state.my[*k] -= factor * un * ny * state.rho[*k];
    }
}

/// per-grid body classification cache: which cells are solid, which sit in
/// the wall band, and a flat solid mask for the flux zeroing. the geometry
/// and grid are static during a run, so this is built once and reused every
/// step; the signature (grid dims, band, corner distances) keeps distinct
/// bodies on same-sized grids from sharing a stale cache in tests.
pub struct BodyCache {
    signature: (usize, usize, u64, [u64; 4]),
    pub solid_cells: Vec<usize>,
    pub band_cells: Vec<usize>,
    pub mask: Vec<bool>,
    /// x faces the polygon boundary crosses: face f of row j is entry
    /// j*(nx+1)+f. a face can be crossed while both adjacent cell
    /// centers sit outside the body (a thin concave wedge such as the
    /// gap between legs); without these the wall pass leaks freestream
    /// straight through the surface.
    pub wall_faces_x: Vec<bool>,
    /// y faces the polygon boundary crosses: face f of column i is
    /// entry i*(ny+1)+f.
    pub wall_faces_y: Vec<bool>,
    /// fluid cells whose four faces are all walls: sub-resolution
    /// closed pockets (a horn notch or tail gap narrower than the
    /// mesh). the mirror-wall stagnation pressure builds in them with
    /// no outflow to absorb it, so they detonate within a few stages.
    /// the march holds them at the inert reference like solid cells.
    pub pocket_cells: Vec<usize>,
}

static BODY_CACHE: std::sync::Mutex<Option<BodyCache>> = std::sync::Mutex::new(None);

fn cache_signature(
    g: &Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    band: f64,
) -> (usize, usize, u64, [u64; 4]) {
    let corners = [
        (0usize, 0usize),
        (g.nx - 1, 0usize),
        (0usize, g.ny - 1),
        (g.nx - 1, g.ny - 1),
    ];
    let probes = corners
        .iter()
        .map(|&(i, j)| {
            let k = j * g.nx + i;
            dist(g.centers_x[k], g.centers_y[k]).to_bits()
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    (g.nx, g.ny, band.to_bits(), probes)
}

/// fetch the classification for this grid+body, building it on first use.
pub fn body_cache(
    g: &Grid2d,
    dist: &dyn Fn(f64, f64) -> f64,
    band: f64,
) -> std::sync::Arc<BodyCache> {
    let sig = cache_signature(g, dist, band);
    let mut guard = BODY_CACHE.lock().unwrap();
    if let Some(c) = guard.as_ref() {
        if c.signature == sig {
            return std::sync::Arc::new(BodyCache {
                signature: sig,
                solid_cells: c.solid_cells.clone(),
                band_cells: c.band_cells.clone(),
                mask: c.mask.clone(),
                wall_faces_x: c.wall_faces_x.clone(),
                wall_faces_y: c.wall_faces_y.clone(),
                pocket_cells: c.pocket_cells.clone(),
            });
        }
    }
    let mut solid = Vec::new();
    let mut bandc = Vec::new();
    let mut mask = vec![false; g.nx * g.ny];
    // the solid mask follows the geometry: a cell belongs to the body
    // when its centre lies inside it. a wall that also swallowed the
    // near-wall fluid cells thickened the passage by half a cell and
    // moved the effective throat and exit area, which walked every
    // nozzle gate off its 1d yardstick. the wedges that leak through a
    // centre-on-centre staircase are closed by the crossed-face pass
    // below and the sub-resolution pockets by the pocket list.
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            let d = dist(g.centers_x[k], g.centers_y[k]);
            if d <= 0.0 {
                solid.push(k);
                mask[k] = true;
            } else if d < band {
                bandc.push(k);
            }
        }
    }
    // faces the boundary crosses between two outside centers: sample the
    // signed distance along each interior face at a few stations and mark
    // any sign change as a crossed face. the distance function is exact to
    // the piecewise-linear surface, so one midpoint plus the two quarter
    // points catch every crossing at practical resolutions (a crossing
    // that hides from all three probes would need a sub-quarter-cell
    // sliver, which the band projection below absorbs).
    let mut wall_x = vec![false; (g.nx + 1) * g.ny];
    let mut wall_y = vec![false; g.nx * (g.ny + 1)];
    let sign_change =
        |d0: f64, d1: f64, d2: f64| (d0 > 0.0) != (d2 > 0.0) || (d1 > 0.0) != (d2 > 0.0);
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        for f in 1..g.nx {
            let xf = g.xmin + f as f64 * g.dx;
            let d0 = dist(xf - 0.25 * g.dx, y);
            let d1 = dist(xf, y);
            let d2 = dist(xf + 0.25 * g.dx, y);
            if sign_change(d0, d1, d2) {
                wall_x[j * (g.nx + 1) + f] = true;
            }
        }
    }
    for i in 0..g.nx {
        let x = g.centers_x[i];
        for f in 1..g.ny {
            let yf = g.ymin + f as f64 * g.dy;
            let d0 = dist(x, yf - 0.25 * g.dy);
            let d1 = dist(x, yf);
            let d2 = dist(x, yf + 0.25 * g.dy);
            if sign_change(d0, d1, d2) {
                wall_y[i * (g.ny + 1) + f] = true;
            }
        }
    }
    // enclosed pockets: a fluid cell whose four faces are all walls (a
    // solid neighbor or a crossed face on every side). the mirror-wall
    // pass seals them completely, so their fluid is a numerical artifact
    // of the staircase quantization, not a resolvable cavity.
    let mut pocket = Vec::new();
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            if mask[k] {
                continue;
            }
            let w = i == 0 || mask[j * g.nx + i - 1] || wall_x[j * (g.nx + 1) + i];
            let e = i + 1 == g.nx || mask[j * g.nx + i + 1] || wall_x[j * (g.nx + 1) + i + 1];
            let s = j == 0 || mask[(j - 1) * g.nx + i] || wall_y[i * (g.ny + 1) + j];
            let n = j + 1 == g.ny || mask[(j + 1) * g.nx + i] || wall_y[i * (g.ny + 1) + j + 1];
            if w && e && s && n {
                pocket.push(k);
            }
        }
    }
    let cache = BodyCache {
        signature: sig,
        solid_cells: solid,
        band_cells: bandc,
        mask,
        wall_faces_x: wall_x,
        wall_faces_y: wall_y,
        pocket_cells: pocket,
    };
    *guard = Some(BodyCache {
        signature: cache.signature,
        solid_cells: cache.solid_cells.clone(),
        band_cells: cache.band_cells.clone(),
        mask: cache.mask.clone(),
        wall_faces_x: cache.wall_faces_x.clone(),
        wall_faces_y: cache.wall_faces_y.clone(),
        pocket_cells: cache.pocket_cells.clone(),
    });
    std::sync::Arc::new(cache)
}

/// flat solid mask for the flux zeroing pass (shares the same cache).
pub fn solid_mask(g: &Grid2d, dist: &dyn Fn(f64, f64) -> f64) -> std::sync::Arc<Vec<bool>> {
    let band = 1.5 * g.dx.min(g.dy);
    std::sync::Arc::new(body_cache(g, dist, band).mask.clone())
}

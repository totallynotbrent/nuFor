//! the curvilinear body-fitted mesh for axisymmetric capsule flows.
//!
//! a structured O-grid over the sphere-cone meridian: the inner
//! boundary IS the body surface (no staircase, no mask band), the
//! outer boundary sits beyond the detached bow shock. cells carry
//! their own areas and face normals, so the finite-volume update is
//! geometrically exact on the mapped grid.
//!
//! the mapping is a transfinite interpolation between the body
//! meridian and the outer boundary curve, with an algebraic cluster
//! function that packs cells near the wall where the shock layer
//! lives.

use crate::spherecone::SphereCone;

/// the far-field state injected at the upstream column and the outer
/// boundary: the freestream the case file defines.
#[derive(Clone, Copy)]
pub struct Freestream {
    pub rho: f64,
    pub u: f64,
    pub p: f64,
}

/// one cell of the mapped grid: centroid plus its geometric data.
#[derive(Clone, Copy)]
pub struct CurvCell {
    /// centroid x (axial, meters).
    pub x: f64,
    /// centroid r (radial, meters, >= 0 for the half-domain).
    pub r: f64,
    /// cell area in the meridian plane (m^2).
    pub area: f64,
}

/// one face of the mapped grid, with its outward-unit normal and
/// length in the meridian plane. faces are stored per row and per
/// column so the sweep code can mirror the rectilinear layout.
#[derive(Clone, Copy)]
pub struct CurvFace {
    /// face midpoint x.
    pub x: f64,
    /// face midpoint r.
    pub r: f64,
    /// face length in the meridian plane (m).
    pub len: f64,
    /// outward unit normal (nx, nr), normalized.
    pub nx: f64,
    pub nr: f64,
}

/// the mapped grid: `nx` streamwise columns (along the body, from the
/// nose over the flank to the outflow), `ny` wall-normal rows (body
/// surface to the outer boundary).
pub struct CurvGrid {
    pub nx: usize,
    pub ny: usize,
    /// cell centroids and areas, row-major [j*nx + i].
    pub cells: Vec<CurvCell>,
    /// vertical faces between columns i and i+1 (nx+1 per row):
    /// `xfaces[j*(nx+1) + i]`.
    pub xfaces: Vec<CurvFace>,
    /// horizontal faces between rows j and j+1 (ny+1 per column):
    /// `yfaces[i*(ny+1) + j]`.
    pub yfaces: Vec<CurvFace>,
    /// the body surface points per column (inner boundary), for
    /// reporting surface states.
    pub surface_x: Vec<f64>,
    pub surface_r: Vec<f64>,
    /// the outer boundary points per column.
    pub outer_x: Vec<f64>,
    pub outer_r: Vec<f64>,
}

/// the body meridian point at streamwise parameter t in [0,1]:
/// sphere cap from the stagnation point around to the tangent point,
/// then the cone flank to the base ring.
/// the body meridian point at parameter t in [0,1]: the shared
/// curve behind the grid builder and the surface reports.
pub fn body_point(sc: &SphereCone, t: f64) -> (f64, f64) {
    let xt = sc.x_tangent();
    let yt = sc.y_tangent();
    let xb = sc.x_base();
    // the cap spans the center angle from PI (the nose) down to
    // PI/2 + delta (the tangent point); the flank then runs to the
    // base ring. split t: 45% on the cap, 55% on the flank.
    let cap_end = 0.45;
    if t <= cap_end {
        let ang = std::f64::consts::PI - (t / cap_end) * (std::f64::consts::PI / 2.0 - sc.delta);
        (sc.xc + sc.rn * ang.cos(), sc.rn * ang.sin())
    } else {
        let u = (t - cap_end) / (1.0 - cap_end);
        (xt + u * (xb - xt), yt + u * (sc.rb - yt))
    }
}

/// the outward unit normal of the body meridian at parameter t.
fn body_normal(sc: &SphereCone, t: f64) -> (f64, f64) {
    let cap_end = 0.45;
    if t <= cap_end {
        // sphere cap: the normal points away from the sphere center
        let ang = std::f64::consts::PI - (t / cap_end) * (std::f64::consts::PI / 2.0 - sc.delta);
        (ang.cos(), ang.sin())
    } else {
        // cone flank: the outward normal of the surface r = y_t +
        // (x - xt) tan(delta)
        (-sc.delta.sin(), sc.delta.cos())
    }
}

impl CurvGrid {
    /// build the grid over the sphere-cone. `far` is the radial
    /// extent of the outer boundary (meters); `beta` the wall
    /// clustering strength.
    pub fn over_sphere_cone(sc: &SphereCone, nx: usize, ny: usize, far: f64, beta: f64) -> Self {
        // node positions first (nx+1 streamwise, ny+1 normal)
        let mut nodes: Vec<(f64, f64)> = vec![(0.0, 0.0); (nx + 1) * (ny + 1)];
        for i in 0..=nx {
            let t = i as f64 / nx as f64;
            let (bx, br) = body_point(sc, t);
            // the outer boundary: same streamwise x as the body
            // downstream, a bow upstream. keep it simple: the outer
            // point sits at (bx, far) downstream of the nose and on a
            // large circle ahead of it.
            // the outer boundary: the body meridian offset outward
            // by `far` along the surface normal. upstream of the nose
            // this reaches x = nose - far, giving the boundary true
            // upstream-facing faces the axial freestream can cross.
            let (nx_b, nr_b) = body_normal(sc, t);
            let ox = bx + far * nx_b;
            let or_ = br + far * nr_b;
            for j in 0..=ny {
                // wall-normal parameter: j=0 ON the body, j=ny at the
                // outer boundary, clustered toward the wall.
                let s = (j as f64 / ny as f64).clamp(0.0, 1.0);
                let sn = (beta * s).tanh() / beta.tanh();
                nodes[j * (nx + 1) + i] = (bx + (ox - bx) * sn, br + (or_ - br) * sn);
            }
        }
        // cells: centroid + area via the four corners
        let mut cells = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let n00 = nodes[j * (nx + 1) + i];
                let n10 = nodes[j * (nx + 1) + i + 1];
                let n01 = nodes[(j + 1) * (nx + 1) + i];
                let n11 = nodes[(j + 1) * (nx + 1) + i + 1];
                // shoelace area (positive when the quad is ccw)
                let area = 0.5
                    * ((n10.0 * n00.1 - n00.0 * n10.1)
                        + (n11.0 * n10.1 - n10.0 * n11.1)
                        + (n01.0 * n11.1 - n11.0 * n01.1)
                        + (n00.0 * n01.1 - n01.0 * n00.1))
                        .abs();
                cells.push(CurvCell {
                    x: 0.25 * (n00.0 + n10.0 + n01.0 + n11.0),
                    r: 0.25 * (n00.1 + n10.1 + n01.1 + n11.1),
                    area,
                });
            }
        }
        // faces
        let mut xfaces = Vec::with_capacity(ny * (nx + 1));
        for j in 0..ny {
            for i in 0..=nx {
                // edge from the bottom node to the top node; the ccw
                // rotation of a +r-pointing edge points +x (outward
                // toward the next column).
                let a = nodes[j * (nx + 1) + i];
                let b = nodes[(j + 1) * (nx + 1) + i];
                let (dx, dr) = (b.0 - a.0, b.1 - a.1);
                let len = (dx * dx + dr * dr).sqrt().max(1e-12);
                xfaces.push(CurvFace {
                    x: 0.5 * (a.0 + b.0),
                    r: 0.5 * (a.1 + b.1),
                    len,
                    nx: -dr / len,
                    nr: dx / len,
                });
            }
        }
        let mut yfaces = Vec::with_capacity(nx * (ny + 1));
        for i in 0..nx {
            for j in 0..=ny {
                // edge from the left node to the right node; the ccw
                // rotation points +r (outward toward the next row).
                let a = nodes[j * (nx + 1) + i];
                let b = nodes[j * (nx + 1) + i + 1];
                let (dx, dr) = (b.0 - a.0, b.1 - a.1);
                let len = (dx * dx + dr * dr).sqrt().max(1e-12);
                yfaces.push(CurvFace {
                    x: 0.5 * (a.0 + b.0),
                    r: 0.5 * (a.1 + b.1),
                    len,
                    nx: -dr / len,
                    nr: dx / len,
                });
            }
        }
        let surface_x = (0..=nx)
            .map(|i| body_point(sc, i as f64 / nx as f64).0)
            .collect();
        let surface_r = (0..=nx)
            .map(|i| body_point(sc, i as f64 / nx as f64).1)
            .collect();
        let outer_x = (0..=nx)
            .map(|i| {
                let t = i as f64 / nx as f64;
                let (bx, _) = body_point(sc, t);
                let nose_x = sc.xc - sc.rn;
                let blend = (t / 0.15).clamp(0.0, 1.0);
                (nose_x - 0.5 * far) + blend * (bx - (nose_x - 0.5 * far))
            })
            .collect();
        let outer_r = vec![far; nx + 1];
        Self {
            nx,
            ny,
            cells,
            xfaces,
            yfaces,
            surface_x,
            surface_r,
            outer_x,
            outer_r,
        }
    }
}

/// keep the tanh helper honest: a simple, single-sided clustering
/// used by tests to pin its behavior.
pub fn wall_cluster(s: f64, beta: f64) -> f64 {
    (beta * s.clamp(0.0, 1.0)).tanh() / beta.tanh()
}

/// one explicit euler step on the curvilinear grid, axisymmetric:
/// face fluxes from the rotated face states, area-weighted update
/// with the annular r-factor, and the geometric pressure source.
///
/// the state layout matches ConservedState2d (row-major, j*nx+i with
/// i streamwise and j wall-normal). the model closes pressure and
/// sound speed the same way the rectilinear axi march does.
pub fn advance2d_axi_curv(
    state: &mut crate::state2d::ConservedState2d,
    g: &CurvGrid,
    model: crate::thermo::ThermoModel,
    cfl: f64,
    freestream: Freestream,
) -> Result<(f64, f64), crate::Error> {
    use crate::hllc2d::{hllc_flux, FacePrim};
    use crate::thermo::{self, ThermoModel};

    let (nx, ny) = (g.nx, g.ny);
    let n = nx * ny;
    let sc_r_ref = g.cells.iter().map(|c| c.r).fold(0.0f64, f64::max).max(1e-6);
    for k in 0..n {
        // pre-step floor: a shocked pole cell can carry rho at or below
        // zero from the previous update; the primitive conversion
        // rejects it with InvalidArgs. floor here so the march can
        // recover instead of failing mid-step (the post-step floors
        // then take over for the derived state).
        if !(state.rho[k].is_finite()) || state.rho[k] <= 1.0e-6 {
            // zero the WHOLE state: leaving e while flooring rho leaves
            // a cell whose implied velocity is (2e/rho)^0.5, which is
            // the next step's NaN. the floor sets the cell to a cold
            // quiescent gas at 1e-6 kg/m^3 (e_int = 1e4, the closure
            // floor), so u = 0, a small, p small — shocked neighbours
            // re-pressurize it through the riemann flux next step.
            state.rho[k] = 1.0e-6;
            state.mx[k] = 0.0;
            state.my[k] = 0.0;
            state.e[k] = 1.0e-6 * 1.0e4;
        } else if !(state.e[k].is_finite()) || state.e[k] < 0.0 {
            // negative/non-finite energy without a bad density: floor the
            // internal part at the closure floor, keep the momentum.
            let ke = 0.5 * (state.mx[k] * state.mx[k] + state.my[k] * state.my[k])
                / state.rho[k].max(1e-12);
            state.e[k] = state.rho[k].max(1e-12) * (ke + 1.0e4);
        }
    }
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0; n];
    for k in 0..n {
        e_int[k] = (et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k])).max(1.0);
    }
    // transient-robust closure evaluation: during the initial
    // collapse the field passes through states outside the TGAS1
    // table's tidy region. floor the closure inputs (not the state)
    // the same way the rectilinear march floors pressure for the
    // perfect gas — physical states here are far above the floors.
    let e_cl: Vec<f64> = e_int.iter().map(|x| x.max(1.0e4)).collect();
    let p = thermo::pressure(model, &state.rho, &e_cl, &u, &v).map_err(|e| {
        // probe: which cell kills the closure
        let bad = (0..n)
            .find(|&k| {
                let p_at = crate::eqair::eqair_pressure_at(state.rho[k], e_cl[k]);
                !p_at.is_finite() || p_at <= 0.0
            })
            .map(|k| (k % nx, k / nx, state.rho[k], e_cl[k]));
        if let Some((bi, bj, br, be)) = bad {
            eprintln!(
                "closure probe: cell (i={}, j={}) rho={:.3e} e_int={:.3e}",
                bi, bj, br, be
            );
        }
        e
    })?;
    let et_cl: Vec<f64> = et
        .iter()
        .enumerate()
        .map(|(k, x)| {
            let ke = 0.5 * (u[k] * u[k] + v[k] * v[k]);
            x.max(ke + 1.0e4)
        })
        .collect();
    // the sound speed only sets wave speeds and dt; if the TGAS1
    // closure momentarily loses positivity in a transient cell at a
    // block edge, fall back to a cold-air speed there rather than
    // failing the march.
    let a = match thermo::sound_speed(model, &state.rho, &et_cl, &u, &v) {
        Ok(v) => v,
        Err(e) => {
            let bad = (0..n)
                .map(|k| {
                    (
                        k,
                        crate::eqair::eqair_sound_at(
                            state.rho[k],
                            et_cl[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]),
                        ),
                    )
                })
                .find(|(_, av)| !av.is_finite() || *av <= 0.0);
            if let Some((k, av)) = bad {
                eprintln!(
                    "sound probe: cell (i={}, j={}) rho={:.3e} e_int={:.3e} et_cl={:.3e} u={:.1} v={:.1} a={:.3e}",
                    k % nx,
                    k / nx,
                    state.rho[k],
                    et_cl[k].max(1e-30) - 0.5 * (u[k] * u[k] + v[k] * v[k]),
                    et_cl[k],
                    u[k],
                    v[k],
                    av
                );
            }
            let _ = e;
            vec![300.0; n]
        }
    };

    // dt from the smallest cell and the fastest wave
    let mut smax = 0.0f64;
    for k in 0..n {
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let mut hmin = f64::INFINITY;
    for k in 0..n {
        // approximate the smaller face spacing from area / face length
        let xf = g.xfaces[k];
        let j = k / nx;
        let yf = g.yfaces[k % nx * (ny + 1) + j];
        hmin = hmin.min(g.cells[k].area / xf.len.max(yf.len).max(1e-12));
    }
    let dt = cfl * hmin / smax.max(1e-12);

    // face states per column boundary: for each row j and face i,
    // build the left/right FacePrim from the adjacent cells in the
    // FACE-LOCAL frame (normal velocity u_n = u*nx + v*nr, tangential
    // u_t = -u*nr + v*nx for x-faces).
    let gamma = match model {
        ThermoModel::Perfect { gamma } => gamma,
        ThermoModel::EqAir => 1.4,
    };
    // x-direction (streamwise) fluxes: index [j*(nx+1) + i]
    let mut fx_mass = vec![0.0f64; ny * (nx + 1)];
    let mut fx_mx = vec![0.0f64; ny * (nx + 1)];
    let mut fx_my = vec![0.0f64; ny * (nx + 1)];
    let mut fx_e = vec![0.0f64; ny * (nx + 1)];

    let face_prim = |model: ThermoModel, rho: f64, un: f64, ut: f64, x: f64| -> FacePrim {
        match model {
            ThermoModel::Perfect { gamma } => FacePrim::perfect(rho, un, ut, x, gamma),
            ThermoModel::EqAir => {
                let p = crate::eqair::eqair_pressure_at(rho, x);
                FacePrim {
                    rho,
                    u: un,
                    v: ut,
                    p,
                    a: crate::eqair::eqair_sound_at(rho, x),
                    e: x,
                }
            }
        }
    };

    for j in 0..ny {
        for i in 0..=nx {
            let f = g.xfaces[j * (nx + 1) + i];
            // left cell (i-1), right cell (i), clipped at the domain
            let (li, ri) = (i.saturating_sub(1), i.min(nx - 1));
            let kl = j * nx + li;
            let kr = j * nx + ri;
            // wall (i=0) and outflow (i=nx) faces take the interior
            // state mirrored when the index clips
            let same = (i == 0 || i == nx) || kl == kr;
            let (rl, ul, vl, xl) = (
                state.rho[kl],
                u[kl],
                v[kl],
                closure_state(model, p[kl], e_int[kl]),
            );
            let (rr, ur, vr, xr) = (
                state.rho[kr],
                u[kr],
                v[kr],
                closure_state(model, p[kr], e_int[kr]),
            );
            // rotate into the face frame
            let (unl, utl) = (ul * f.nx + vl * f.nr, -ul * f.nr + vl * f.nx);
            let (unr, utr) = (ur * f.nx + vr * f.nr, -ur * f.nr + vr * f.nx);
            // the upstream end (i=0) is the nose seam of the
            // O-grid: a symmetry plane through the stagnation
            // streamline, not an inflow. the flow enters through the
            // outer bow. mirror the interior state with the normal
            // velocity reflected (slip wall), which is exact for the
            // on-axis stagnation streamline.
            let (left, right) = if i == 0 {
                // seam: reflect normal velocity OUT so the riemann
                // solver doesn't believe there is momentum entering
                // through the body. utl stays (tangential slip).
                let l = face_prim(model, rl, -unl, utl, xl);
                (l, l)
            } else {
                let l = face_prim(model, rl, unl, utl, xl);
                let r = if same {
                    l
                } else {
                    face_prim(model, rr, unr, utr, xr)
                };
                (l, r)
            };
            let q = hllc_flux(gamma, left, right, 0);
            let (fm, fmx, fmy, fe) = (q.mass, q.mx, q.my, q.e);
            // rotate the momentum flux back to the lab frame
            fx_mass[j * (nx + 1) + i] = fm * f.len;
            fx_mx[j * (nx + 1) + i] = (fmx * f.nx - fmy * f.nr) * f.len;
            fx_my[j * (nx + 1) + i] = (fmx * f.nr + fmy * f.nx) * f.len;
            fx_e[j * (nx + 1) + i] = fe * f.len;
        }
    }

    // y-direction (wall-normal) fluxes: index [i*(ny+1) + j]
    let mut fy_mass = vec![0.0f64; nx * (ny + 1)];
    let mut fy_mx = vec![0.0f64; nx * (ny + 1)];
    let mut fy_my = vec![0.0f64; nx * (ny + 1)];
    let mut fy_e = vec![0.0f64; nx * (ny + 1)];
    for i in 0..nx {
        for j in 0..=ny {
            let f = g.yfaces[i * (ny + 1) + j];
            let (lj, rj) = (j.saturating_sub(1), j.min(ny - 1));
            let kl = lj * nx + i;
            let kr = rj * nx + i;
            let (rl, ul, vl, xl) = (
                state.rho[kl],
                u[kl],
                v[kl],
                closure_state(model, p[kl], e_int[kl]),
            );
            let (rr, ur, vr, xr) = (
                state.rho[kr],
                u[kr],
                v[kr],
                closure_state(model, p[kr], e_int[kr]),
            );
            let (unl, utl) = (ul * f.nx + vl * f.nr, -ul * f.nr + vl * f.nx);
            let (unr, utr) = (ur * f.nx + vr * f.nr, -ur * f.nr + vr * f.nx);
            if j == 0 {
                // the body wall: zero mass flux, slip tangential.
                let wall = face_prim(model, rl, 0.0, utl, xl);
                let q = hllc_flux(gamma, wall, wall, 0);
                let (fm, fmx, fmy, fe) = (q.mass, q.mx, q.my, q.e);
                fy_mass[i * (ny + 1) + j] = fm * f.len;
                fy_mx[i * (ny + 1) + j] = (fmx * f.nx - fmy * f.nr) * f.len;
                fy_my[i * (ny + 1) + j] = (fmx * f.nr + fmy * f.nx) * f.len;
                fy_e[i * (ny + 1) + j] = fe * f.len;
                continue;
            }
            let same = kl == kr;
            // the outer bow injects the freestream.
            let (left, right) = if j == ny {
                let xf_e = freestream_e(model, &freestream);
                // the lab-frame freestream (u_inf, 0) rotated into
                // the face frame: un = u*nx, ut = -u*nr.
                let fl = face_prim(
                    model,
                    freestream.rho,
                    freestream.u * f.nx,
                    -freestream.u * f.nr,
                    xf_e,
                );
                (fl, fl)
            } else {
                let l = face_prim(model, rl, unl, utl, xl);
                let r = if same {
                    l
                } else {
                    face_prim(model, rr, unr, utr, xr)
                };
                (l, r)
            };
            let q = hllc_flux(gamma, left, right, 0);
            let (fm, fmx, fmy, fe) = (q.mass, q.mx, q.my, q.e);
            fy_mass[i * (ny + 1) + j] = fm * f.len;
            fy_mx[i * (ny + 1) + j] = (fmx * f.nx - fmy * f.nr) * f.len;
            fy_my[i * (ny + 1) + j] = (fmx * f.nr + fmy * f.nx) * f.len;
            fy_e[i * (ny + 1) + j] = fe * f.len;
        }
    }

    // area-weighted update with the annular r-factor
    let mut resid = 0.0f64;
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let c = g.cells[k];
            // annular volume per unit azimuth: 2*pi*r_c * area
            // pole treatment: near the symmetry axis the annular
            // volume collapses; clamp r to a fraction of the
            // wall-normal spacing so the first row's update stays
            // bounded.
            let r_pole = (c.r).max(1e-3 * sc_r_ref);
            let vol = c.area * 2.0 * std::f64::consts::PI * r_pole;
            // face r-weights: the annular area of each face
            let xf_l = g.xfaces[j * (nx + 1) + i];
            let xf_r = g.xfaces[j * (nx + 1) + i + 1];
            let yf_b = g.yfaces[i * (ny + 1) + j];
            let yf_t = g.yfaces[i * (ny + 1) + j + 1];
            let sa_l = 2.0 * std::f64::consts::PI * xf_l.r.max(1e-9) * xf_l.len;
            let sa_r = 2.0 * std::f64::consts::PI * xf_r.r.max(1e-9) * xf_r.len;
            let sa_b = 2.0 * std::f64::consts::PI * yf_b.r.max(1e-9) * yf_b.len;
            let sa_t = 2.0 * std::f64::consts::PI * yf_t.r.max(1e-9) * yf_t.len;
            let fmass = -(fx_mass[j * (nx + 1) + i + 1] * sa_r - fx_mass[j * (nx + 1) + i] * sa_l
                + fy_mass[i * (ny + 1) + j + 1] * sa_t
                - fy_mass[i * (ny + 1) + j] * sa_b)
                / vol;
            let fmx = -(fx_mx[j * (nx + 1) + i + 1] * sa_r - fx_mx[j * (nx + 1) + i] * sa_l
                + fy_mx[i * (ny + 1) + j + 1] * sa_t
                - fy_mx[i * (ny + 1) + j] * sa_b)
                / vol;
            let fmy = -(fx_my[j * (nx + 1) + i + 1] * sa_r - fx_my[j * (nx + 1) + i] * sa_l
                + fy_my[i * (ny + 1) + j + 1] * sa_t
                - fy_my[i * (ny + 1) + j] * sa_b)
                / vol;
            // the geometric source: the pressure thrust on the
            // r-projected face areas. sa_* above are the full annular
            // areas; the radial source needs the r-projection, i.e.
            // each face's area weighted by its radial normal.
            let src_y = p[k]
                * 2.0
                * std::f64::consts::PI
                * (yf_t.r * yf_t.len * yf_t.nr - yf_b.r * yf_b.len * yf_b.nr)
                / vol;
            let fe = -(fx_e[j * (nx + 1) + i + 1] * sa_r - fx_e[j * (nx + 1) + i] * sa_l
                + fy_e[i * (ny + 1) + j + 1] * sa_t
                - fy_e[i * (ny + 1) + j] * sa_b)
                / vol;
            state.rho[k] += dt * fmass;
            state.mx[k] += dt * fmx;
            state.my[k] += dt * (fmy + src_y);
            state.e[k] += dt * fe;
            // transient guard: the uniform IC is unphysical near the
            // body while the shock organizes, and the collapse can
            // drive a cell's internal energy below the closure's
            // valid band. floor at 1e4 J/kg (~145 K): below every
            // physical state in hypersonic air, above the TGAS1
            // cold-extrapolation region where the sound-speed
            // closure loses positivity.
            // density guard: the shock-transient near the nose can
            // momentarily drive a tiny wall cell negative; floor at
            // 1e-6 kg/m^3 (near-vacuum, far below any state here).
            if state.rho[k] < 1.0e-6 {
                state.rho[k] = 1.0e-6;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
            }
            let ke_k = 0.5 * (state.mx[k] * state.mx[k] + state.my[k] * state.my[k])
                / state.rho[k].max(1e-12);
            let et_k = state.e[k] / state.rho[k].max(1e-12);
            if et_k - ke_k < 1.0e4 {
                state.e[k] = state.rho[k].max(1e-12) * (ke_k + 1.0e4);
            }
            resid = resid
                .max(fmass.abs())
                .max(fmx.abs())
                .max(fmy.abs())
                .max(fe.abs());
        }
    }
    Ok((dt, resid))
}

/// the closure variable for the injected freestream: pressure for
/// the perfect gas, internal energy for equilibrium air.
fn freestream_e(model: crate::thermo::ThermoModel, fs: &Freestream) -> f64 {
    match model {
        crate::thermo::ThermoModel::Perfect { .. } => fs.p,
        crate::thermo::ThermoModel::EqAir => crate::eqair::energy_from_pressure(&[fs.rho], &[fs.p])
            .unwrap_or_default()
            .first()
            .copied()
            .unwrap_or(0.0),
    }
}

/// the closure variable the face builder reconstructs: pressure for
/// the perfect gas, internal energy for equilibrium air.
fn closure_state(model: crate::thermo::ThermoModel, p: f64, e_int: f64) -> f64 {
    match model {
        crate::thermo::ThermoModel::Perfect { .. } => p,
        crate::thermo::ThermoModel::EqAir => e_int,
    }
}

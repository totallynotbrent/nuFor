//! the axisymmetric viscous combustion wall on a CUT-CELL geometry.
//!
//! WHY THIS EXISTS
//!
//! the engine wall bracket needs two things at once: a wall that takes
//! measurable heat out of the gas, and a throat whose area is not
//! quantised by the grid. the staircase representation cannot give both
//! at once: its throat area snaps to whole radial rows, and its wall
//! stencil is a three-point extrapolation across a jagged boundary.
//!
//! in a cut cell the wall's position is KNOWN: the cell is cut against
//! the body's exact signed distance function, so the distance from the
//! cell's centre to the wall is the signed distance itself. the
//! one-sided derivative to the wall value is then a linear ramp over
//! that exact distance, not an extrapolation across a staircase.
//!
//! the wall here is:
//
//!   - isothermal, with the wall's temperature given
//!   - no slip: the tangential velocity vanishes at the wall
//!   - impermeable: no mass crosses it
//!   - non-reactive: the wall neither creates nor absorbs species, so
//!     the wall's composition gradient carries no flux
//!
//! the diffusion is molecular: Wilke's mixture rule for the viscosity
//! and the conductivity, and the effective diffusivities for the
//! species, the same coefficients the planar viscous lane and the
//! bridge lane use.

use crate::combustion_closure as closure9;
use crate::combustion_march::CombustionState;
use crate::combustion_thermo as ct;
use crate::cutcell_axi_geom::AxiCutGeom;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::state2d::{cons_to_prim2d, ConservedState2d};

/// the ratio of the case's absolute viscosity to the gas's MOLECULAR
/// viscosity at its hottest cell.
///
/// the case names an absolute viscosity, while the closure gives the
/// molecule's own. the scale bridges the two, so the cut-cell lane
/// transports the gas with the viscosity the case asked for, keeping
/// the molecular diffusivities in proportion. a case that asks for the
/// molecular viscosity gets a scale of one and the molecular transport
/// stands as it is.
pub fn viscosity_scale(state: &ConservedState2d, species: &CombustionState, mu_case: f64) -> f64 {
    if mu_case <= 0.0 {
        return 1.0;
    }
    let (u, v, et) = match cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e) {
        Ok(v) => v,
        Err(_) => return 1.0,
    };
    let mut mol_mu = 0.0f64;
    let mut t_hot = f64::MIN;
    for k in 0..state.rho.len() {
        let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] = species.y[s][k];
        }
        let t = closure9::temperature_at(e_int, &y);
        if t > t_hot {
            t_hot = t;
            mol_mu = crate::combustion_transport::mixture_viscosity(&y, t);
        }
    }
    if mol_mu > 1.0e-12 {
        mu_case / mol_mu
    } else {
        1.0
    }
}

/// the transport properties of the combustion set, per cell.
pub(crate) struct TransportProps {
    pub t: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub d_s: Vec<[f64; ct::NS]>,
}

/// the per-cell transport properties.
pub(crate) fn transport_props(
    state: &ConservedState2d,
    species: &CombustionState,
    g: &Grid2d,
    mu_scale: f64,
) -> Result<TransportProps, Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut tr = TransportProps {
        t: vec![0.0f64; n],
        mu: vec![0.0f64; n],
        k: vec![0.0f64; n],
        d_s: vec![[0.0f64; ct::NS]; n],
    };
    for k in 0..n {
        let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] = species.y[s][k];
        }
        tr.t[k] = closure9::temperature_at(e_int, &y);
        let p = closure9::pressure_at(state.rho[k], e_int, &y);
        tr.mu[k] = crate::combustion_transport::mixture_viscosity(&y, tr.t[k]) * mu_scale;
        let mut cps = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            cps[s] = ct::cp_mass(s, tr.t[k]);
        }
        tr.k[k] = crate::combustion_transport::mixture_conductivity(&y, tr.t[k], &cps);
        tr.d_s[k] = crate::combustion_transport::effective_diffusivity(&y, tr.t[k], p.max(1.0));
    }
    Ok(tr)
}

/// one viscous application on a cut-cell grid, with the wall on the
/// immersed boundary.
///
/// the diffusion runs on the OPEN faces between cells, weighted by
/// those faces' apertures, so the composition and the momentum obey the
/// same wall. the wall's own flux is applied at the immersed boundary,
/// where the wall's exact distance from the cell centre is known.
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_cut_combustion_viscous(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cut: &AxiCutGeom,
    tr: &TransportProps,
    dt: f64,
    wall_temperature: f64,
    sdf: &dyn Fn(f64, f64) -> f64,
) -> Result<(), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let n = nx * ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;

    // the diffused fluxes on the open faces.
    let mut tx_x = vec![0.0f64; (nx + 1) * ny];
    let mut ty_x = vec![0.0f64; (nx + 1) * ny];
    let mut te_x = vec![0.0f64; (nx + 1) * ny];
    let mut tx_y = vec![0.0f64; nx * (ny + 1)];
    let mut ty_y = vec![0.0f64; nx * (ny + 1)];
    let mut te_y = vec![0.0f64; nx * (ny + 1)];
    let mut ts_x: Vec<Vec<f64>> = vec![vec![0.0f64; (nx + 1) * ny]; ct::NS];
    let mut ts_y: Vec<Vec<f64>> = vec![vec![0.0f64; nx * (ny + 1)]; ct::NS];

    // the axial faces
    for j in 0..ny {
        for f in 1..nx {
            let a = idx(f - 1, j);
            let b = idx(f, j);
            // the shared face's opening: ONE number for the face,
            // measured on the face itself, so the two cells cannot
            // disagree about how much of it is open.
            let ap = cut.ax_at(f, j);
            if ap <= 0.0 {
                continue;
            }
            let k = f * ny + j;
            let mu_f = 0.5 * (tr.mu[a] + tr.mu[b]);
            let k_f = 0.5 * (tr.k[a] + tr.k[b]);
            let dx = g.centers_x[b] - g.centers_x[a];
            let dudx = (u[b] - u[a]) / dx;
            let dvdx = (v[b] - v[a]) / dx;
            let dtdx = (tr.t[b] - tr.t[a]) / dx;
            let rc = g.centers_y[a];
            let div = dudx + v[a] / rc.max(1.0e-12) + (v[b] - v[a]) / g.dys[j];
            let txx = mu_f * (2.0 * dudx - 2.0 / 3.0 * div);
            let txy = mu_f * dvdx;
            let qx = k_f * dtdx;
            let rho_f = 0.5 * (state.rho[a] + state.rho[b]);
            let mut q_enth = 0.0f64;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let dydx = (species.y[s][b] - species.y[s][a]) / dx;
                let d_f = 0.5 * (tr.d_s[a][s] + tr.d_s[b][s]);
                let js = -rho_f * d_f * dydx;
                ts_x[s][k] = js;
                q_enth += ct::h_mass(s, 0.5 * (tr.t[a] + tr.t[b])) * js;
            }
            te_x[k] = 0.5 * (u[a] + u[b]) * txx + 0.5 * (v[a] + v[b]) * txy + qx - q_enth;
            tx_x[k] = txx;
            ty_x[k] = txy;
        }
    }

    // the radial faces
    for i in 0..nx {
        for f in 1..ny {
            let a = idx(i, f - 1);
            let b = idx(i, f);
            let ap = cut.ay_at(f, i);
            if ap <= 0.0 {
                continue;
            }
            let k = i * (ny + 1) + f;
            let mu_f = 0.5 * (tr.mu[a] + tr.mu[b]);
            let k_f = 0.5 * (tr.k[a] + tr.k[b]);
            let dy = g.centers_y[b] - g.centers_y[a];
            let dudy = (u[b] - u[a]) / dy;
            let dvdy = (v[b] - v[a]) / dy;
            let dtdy = (tr.t[b] - tr.t[a]) / dy;
            let rc = 0.5 * (g.centers_y[a] + g.centers_y[b]);
            let div = (u[b] - u[a]) / g.dxs[i] + v[a] / rc.max(1.0e-12) + dvdy;
            let txy = mu_f * dudy;
            let tyy = mu_f * (2.0 * dvdy - 2.0 / 3.0 * div);
            let qy = k_f * dtdy;
            let rho_f = 0.5 * (state.rho[a] + state.rho[b]);
            let mut q_enth = 0.0f64;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let dydy = (species.y[s][b] - species.y[s][a]) / dy;
                let d_f = 0.5 * (tr.d_s[a][s] + tr.d_s[b][s]);
                let js = -rho_f * d_f * dydy;
                ts_y[s][k] = js;
                q_enth += ct::h_mass(s, 0.5 * (tr.t[a] + tr.t[b])) * js;
            }
            te_y[k] = 0.5 * (u[a] + u[b]) * txy + 0.5 * (v[a] + v[b]) * tyy + qy - q_enth;
            tx_y[k] = txy;
            ty_y[k] = tyy;
        }
    }

    // the IMMERSED BOUNDARY: the wall's flux in each cut cell.
    //
    // the wall's distance from the cell's centre is the signed distance
    // function's value there (positive inside the fluid), so the
    // one-sided derivative from the cell's centre to the wall is
    // (wall_value - cell_value) / distance. no extrapolation, and no
    // dependence on where the wall happens to fall against the grid.
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            if cut.vol[k] <= 0.0 || cut.vol[k] >= 1.0 || cut.eb_len[k] <= 0.0 {
                continue;
            }
            let (cx, cy) = (g.centers_x[i], g.centers_y[j]);
            let dist = sdf(cx, cy).max(1.0e-9);
            let (enx, eny) = (-cut.eb_nx[k], -cut.eb_ny[k]);
            // the components of the cell's velocity along the wall's
            // tangent and normal
            let _un = u[k] * enx + v[k] * eny;
            let ut = -u[k] * eny + v[k] * enx;
            // the shear: the tangential velocity ramps linearly to zero
            // at the wall over the exact distance.
            let shear = -ut / dist;
            // the conduction: the temperature ramps to the wall's over
            // the same distance.
            //
            // fourier's law: the heat flux along the inward normal is
            // q.n = -k grad(T).n, and the temperature RISES from the
            // wall to the gas when the gas is hotter, so grad(T).n is
            // +(T_cell - T_wall) / dist. the flux INTO the gas is
            // therefore -k (T_cell - T_wall) / dist: NEGATIVE for a
            // cold wall, which is the heat leaving the gas that the
            // bracket measures.
            let flux_into_gas = -tr.k[k] * (tr.t[k] - wall_temperature) / dist;

            // the wall's face: the interface's LENGTH inside this cell,
            // seen obliquely in the axisymmetric plane.
            //
            // the cut cell's interface is the chord where the wall
            // crosses it. for a straight wall crossing a cell at a
            // distance `dist` from the centre, the chord's length is
            // 2 sqrt((half the cell)^2 - dist^2) when the wall is
            // nearly perpendicular to the cell's shorter side. the
            // cell's own dimensions give the transverse extent.
            let (cw, ch) = (g.dxs[i], g.dys[j]);
            // the distance from the cell's centre to its nearer edge
            // along the wall's normal, and the chord's half-length
            let half_normal = if eny.abs() > enx.abs() {
                0.5 * ch
            } else {
                0.5 * cw
            };
            let along = if eny.abs() > enx.abs() { cw } else { ch };
            let half_chord = (half_normal * half_normal - dist * dist).max(0.0).sqrt();
            let iface_len = 2.0 * half_chord;
            // the oblique projection: the wall's face tilts away from
            // the coordinate plane by the angle between its normal and
            // the cell's normal, so the true face is longer than its
            // projected chord.
            let tilt = (along * along + dist * dist).sqrt() / along.max(1.0e-12);
            let r_face = cy.max(1.0e-9);
            // the wall's area: the interface's length times the
            // circumference's sweep. a cylindrical wall's area is
            // 2 pi r times its length along the axis.
            let area = 2.0 * std::f64::consts::PI * r_face * iface_len * tilt;

            // the wall's traction, resolved onto the axes. the shear
            // acts along the tangent; the pressure's normal component
            // is already accounted for by the flow update.
            let tau = tr.mu[k] * shear;
            state.mx[k] -= dt * tau * (-eny) * area / cut.vol[k].max(1.0e-9);
            state.my[k] -= dt * tau * (enx) * area / cut.vol[k].max(1.0e-9);
            // the conductive flux into the gas raises the cell's
            // energy. a cold wall gives a negative flux, so the cell
            // LOSES energy, which is the sign the bracket needs.
            state.e[k] += dt * flux_into_gas * area / cut.vol[k].max(1.0e-9);
        }
    }

    // the divergence update, with the apertures and the cell's own
    // volume, mirroring the cut-cell flow update.
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            if cut.vol[k] <= 0.0 {
                continue;
            }
            let vol = cut.vol[k].max(1.0e-9);
            let dxi = dt / g.dxs[i];
            let kfx_e = (i + 1) * ny + j;
            let kfx_w = i * ny + j;
            let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
            let rc = 0.5 * (rm + rp);
            // the face apertures: the SHARED face's opening, one
            // number for the face, so the two cells cannot disagree.
            let fe = if i + 1 < nx { cut.ax_at(i + 1, j) } else { 0.0 };
            let fw = if i > 0 { cut.ax_at(i, j) } else { 0.0 };
            let fn_ = if j + 1 < ny { cut.ay_at(j + 1, i) } else { 0.0 };
            let fs = if j > 0 { cut.ay_at(j, i) } else { 0.0 };
            // the radial divergence, with the annulus weight AND the
            // radial faces' apertures, exactly as the cut-cell flow
            // update does.
            let radial =
                |arr: &Vec<f64>| (rp * fn_ * arr[j + 1] - rm * fs * arr[j]) / (rc * g.dys[j]);
            let dmx = -(dxi * (fe * tx_x[kfx_e] - fw * tx_x[kfx_w]) + dt * radial(&tx_y)) / vol;
            let dmy = -(dxi * (fe * ty_x[kfx_e] - fw * ty_x[kfx_w]) + dt * radial(&ty_y)) / vol;
            let de = -(dxi * (fe * te_x[kfx_e] - fw * te_x[kfx_w]) + dt * radial(&te_y)) / vol;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let ds = -(dxi * (fe * ts_x[s][kfx_e] - fw * ts_x[s][kfx_w])
                    + dt * radial(&ts_y[s]))
                    / vol;
                species.y[s][k] += ds / state.rho[k].max(1.0e-12);
            }
        }
    }
    for k in 0..n {
        species.renorm(k);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a COLD wall must cool the gas beside it, in a cut cell.
    ///
    /// the cut cell's wall sits at a known distance inside the cell, so
    /// the one-sided derivative is a ramp over that exact distance. the
    /// cell next to the wall must lose energy, and a cell far from the
    /// wall must not gain any.
    #[test]
    fn a_cold_wall_cools_the_cell_next_to_it() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let sdf = |_x: f64, y: f64| y - 0.032;
        let cut = AxiCutGeom::from_sdf(&sdf, &g);
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = transport_props(&st, &sp, &g, 1.0).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        add_cut_combustion_viscous(&mut s1, &mut y1, &g, &cut, &tr, 1.0e-5, 300.0, &sdf).unwrap();
        // the wall at y = 0.032 falls inside row 3 (0.030..0.040), so
        // row 3 is the cut cell with the wall in it.
        let near = 3 * nx;
        let far = (ny - 1) * nx;
        assert!(
            s1.e[near] < st.e[near],
            "the cell next to a cold wall gained energy: {} against {}",
            s1.e[near],
            st.e[near]
        );
        let far_delta = s1.e[far] - st.e[far];
        assert!(
            far_delta.abs() < 1.0e-9,
            "the cell far from the wall changed: {far_delta}"
        );
    }

    /// the wall is IMPERMEABLE: the viscous pass must not move mass.
    /// a wall that leaked mass would change the throat's state, and the
    /// bracket would be measuring the leak rather than the cooling.
    #[test]
    fn the_viscous_pass_leaves_the_mass_alone() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let sdf = |_x: f64, y: f64| y - 0.032;
        let cut = AxiCutGeom::from_sdf(&sdf, &g);
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = transport_props(&st, &sp, &g, 1.0).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        add_cut_combustion_viscous(&mut s1, &mut y1, &g, &cut, &tr, 1.0e-5, 300.0, &sdf).unwrap();
        let before: f64 = st.rho.iter().sum();
        let after: f64 = s1.rho.iter().sum();
        assert!(
            (after - before).abs() < 1.0e-12 * before.abs().max(1.0),
            "the viscous pass changed the mass: {before} against {after}"
        );
    }

    /// the case's ABSOLUTE viscosity must map onto the molecular one by
    /// their ratio. asking for the molecular viscosity must give a scale
    /// of one, so the molecular transport stands unchanged.
    #[test]
    fn the_viscosity_scale_maps_absolute_onto_molecular() {
        let (nx, ny) = (4, 4);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        let tr = transport_props(&st, &sp, &g, 1.0).unwrap();
        // the hottest cell's molecular viscosity
        // the gas is uniform in this fixture, so one cell speaks for
        // all of them.
        let mut y = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] = sp.y[s][0];
        }
        let mol_mu = crate::combustion_transport::mixture_viscosity(&y, tr.t[0]);
        let scale = viscosity_scale(&st, &sp, mol_mu);
        assert!(
            (scale - 1.0).abs() < 1.0e-9,
            "asking for the molecular viscosity gave a scale of {scale}"
        );
        // and twice the molecular viscosity must give exactly two
        let scale2 = viscosity_scale(&st, &sp, 2.0 * mol_mu);
        assert!(
            (scale2 - 2.0).abs() < 1.0e-9,
            "doubling the viscosity gave a scale of {scale2}"
        );
    }
}

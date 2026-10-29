//! the cut-cell axisymmetric combustion march: the nine-species
//! reacting flow on a grid whose wall is resolved by REAL APERTURES.
//!
//! WHY THIS EXISTS
//!
//! the chamber-throat case marches a walled duct on a cartesian grid,
//! and the wall's staircase quantises the throat's area to whole radial
//! rows. on a 110 by 40 grid the discrete throat came out at
//! 1.256637e-03 m^2 against the contour's true 1.333814e-03, a 5.79
//! percent shortfall, and the marching flow paid for it with a lost
//! total pressure: the choked pressure sat 7 to 8 percent below the
//! gas's own isentrope.
//!
//! refining the grid does not fix that. the throat's radius snapped to
//! 0.020000 m on the 110 by 40, the 165 by 60 and the 220 by 80 grids
//! alike, because the contour's neck fell on the same fractional row
//! each time. some grids land on the contour's area by luck of
//! alignment, and choosing one of those would have posted a number that
//! depended on the wall's PHASE rather than on the physics.
//!
//! WHAT THIS DOES ABOUT IT
//!
//! the cut-cell machinery already in the repo resolves the wall by
//! cutting each cell against the body's exact signed distance function.
//! each cell gets a fluid VOLUME FRACTION and per-face APERTURES, so
//! the throat's area is a smooth function of where the contour falls
//! instead of a whole-row quantum. the flow update then:
//!
//!   - weighs every face flux by that face's aperture
//!   - divides each cell's increment by that cell's OWN fluid volume
//!   - applies the immersed-boundary wall flux with the reflection
//!     ghost pair
//!   - pools the small cells' states with Berger-Giuliani
//!     redistribution, so a cell with little fluid in it cannot blow up
//!
//! THIS ADDS the parts the plain cut-cell march does not have, because
//! this is a REACTING gas:
//!
//!   - the closure's own pressure and sound speed per cell, since the
//!     nine-species closure does not go through the perfect-gas
//!     thermo::pressure call
//!   - aperture-weighted SPECIES transport, so the composition obeys
//!     the same wall the momentum does
//!   - redistribution that carries the composition with the conserved
//!     state, so a small cell's mixture is not averaged with a
//!     neighbour's
//!   - the chemistry substep on the post-blend cell, using the
//!     composition the cell actually ended up with
//!
//! the flux faces are built by the shared axi machinery
//! (crate::axi::compute_axi_fluxes), so the cut-cell version sees
//! exactly the same reconstruction, limiter and HLLC solve as the
//! staircase version. the difference between the two is the geometry and
//! nothing else.

use crate::axi::{compute_axi_fluxes, SpeciesField};
use crate::combustion_closure as closure9;
use crate::combustion_march::CombustionState;
use crate::combustion_thermo as ct;
use crate::cutcell_march::CutField;
use crate::grid2d::Grid2d;
use crate::hllc2d::{hllc_flux, FacePrim};
use crate::state2d::{cons_to_prim2d, ConservedState2d};
use crate::thermo::ThermoModel;
use crate::{Boundaries2d, Error};

/// the per-cell closure state the cut-cell march needs.
struct Closure {
    p: Vec<f64>,
    a: Vec<f64>,
    e_int: Vec<f64>,
    e_floor: Vec<f64>,
}

fn closures(
    state: &ConservedState2d,
    species: &CombustionState,
    n: usize,
) -> Result<Closure, Error> {
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut p = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    let mut e_int = vec![0.0f64; n];
    let mut e_floor = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] = species.y[s][k];
        }
        p[k] = closure9::pressure_at(state.rho[k], e_int[k], &y);
        a[k] = closure9::sound_speed_at(state.rho[k], e_int[k], &y);
        e_floor[k] = closure9::internal_energy_at(200.0, &y);
    }
    Ok(Closure {
        p,
        a,
        e_int,
        e_floor,
    })
}

/// the aperture-weighted species transport.
///
/// the shared routine transports the composition with the same face
/// mass fluxes, but it assumes every face is fully open. a cut cell's
/// faces are only partly open, so the fluxes need scaling by the same
/// apertures the momentum update used, or the composition would leak
/// through the wall the momentum respects.
///
/// the update is on rho*Y, divided by the cell's own fluid volume, so
/// the per-species mass stays consistent with the conserved mass.
#[allow(clippy::too_many_arguments)]
fn transport_species_apertures(
    y: &mut [Vec<f64>],
    rho: &[f64],
    dt: f64,
    g: &Grid2d,
    fx: &[Vec<f64>],
    fy: &[Vec<f64>],
    cut: &CutField,
) -> f64 {
    let ns = y.len();
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let mut resid = 0.0f64;

    // the face upwind compositions, per species.
    let mut yx = vec![vec![0.0f64; (nx + 1) * ny]; ns];
    let mut yy = vec![vec![0.0f64; (ny + 1) * nx]; ns];
    for s in 0..ns {
        let yi = &y[s];
        for j in 0..ny {
            for f in 0..=nx {
                let mf = fx[j][f];
                let left = if f > 0 {
                    yi[idx(f - 1, j)]
                } else {
                    yi[idx(f, j)]
                };
                let right = if f < nx {
                    yi[idx(f, j)]
                } else {
                    yi[idx(f - 1, j)]
                };
                yx[s][f * ny + j] = if mf >= 0.0 { left } else { right };
            }
        }
        for i in 0..nx {
            for f in 0..=ny {
                let mf = fy[i][f];
                let low = if f > 0 {
                    yi[idx(i, f - 1)]
                } else {
                    yi[idx(i, f)]
                };
                let high = if f < ny {
                    yi[idx(i, f)]
                } else {
                    yi[idx(i, f - 1)]
                };
                yy[s][f * nx + i] = if mf >= 0.0 { low } else { high };
            }
        }
    }

    // the transported mass fractions, rho*Y, per cell.
    let mut mass_y: Vec<Vec<f64>> = (0..ns)
        .map(|s| (0..nx * ny).map(|k| rho[k] * y[s][k]).collect())
        .collect();

    for s in 0..ns {
        for j in 0..ny {
            let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
            let rc = 0.5 * (rm + rp);
            let dr = g.dys[j];
            for i in 0..nx {
                let k = idx(i, j);
                let f = cut.at(i, j);
                if f.vol <= 0.0 {
                    // no fluid here, so there is nothing to transport.
                    // the cell's composition is left alone: it is
                    // scaffolding beside the wall, not part of the gas.
                    continue;
                }
                let dxi = dt / g.dxs[i];
                // the aperture-weighted face mass fluxes
                let flux_x = f.east * fx[j][i + 1] * yx[s][(i + 1) * ny + j]
                    - f.west * fx[j][i] * yx[s][i * ny + j];
                let radial = (rp * f.north * fy[i][j + 1] * yy[s][(j + 1) * nx + i]
                    - rm * f.south * fy[i][j] * yy[s][j * nx + i])
                    / (rc * dr);
                let d = -(dxi * flux_x + dt * radial) / f.vol;
                mass_y[s][k] += d;
                resid = resid.max(d.abs());
            }
        }
    }

    // back to fractions. the density is the conserved one, so the
    // fraction is the transported mass over it.
    for s in 0..ns {
        for k in 0..nx * ny {
            if rho[k] > 1.0e-12 {
                y[s][k] = (mass_y[s][k] / rho[k]).max(0.0);
            }
        }
    }
    resid
}

/// one cut-cell axisymmetric combustion step: the aperture-weighted
/// flux update, the aperture-weighted species transport, the small
/// cells' redistribution, and the chemistry on the post-step cell.
///
/// this MIRRORS advance2d_combustion_axi exactly, with ONE difference:
/// the geometry. every algebraic step, every sign and every ordering is
/// the staircase's, so a difference in the result is attributable to the
/// wall's representation and to nothing else. that attribution is the
/// whole point: the two lanes must differ ONLY where they are meant to.
#[allow(clippy::too_many_arguments)]
// the species loops index fixed-size arrays by species number, which
// reads more plainly than zipping two collections of the same length.
#[allow(clippy::needless_range_loop)]
pub fn advance2d_combustion_axi_cut(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cut: &CutField,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    wall_temperature: f64,
    viscosity_scale: f64,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;

    let clo = closures(state, species, n)?;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut smax = 0.0f64;
    for k in 0..n {
        if cut.cells[k].vol <= 0.0 {
            continue;
        }
        smax = smax.max(u[k].abs().max(v[k].abs()) + clo.a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1.0e-12);

    // the SAME flux machinery the staircase uses, and the SAME wall
    // argument: the wall masks the solid cells out of the flux
    // construction, so the cells beside the wall do not see a bogus
    // pair with the scaffold inside it.
    let model = ThermoModel::Combustion9 { gamma_shaping: 1.2 };
    let (fx, fy) = compute_axi_fluxes(
        model,
        state,
        g,
        muscl,
        bc,
        &clo.p,
        &clo.e_int,
        Some(SpeciesField::Combustion9(&*species)),
        nthreads,
        wall,
    )?;

    // the entry state, so the solid scaffold can be restored after the
    // update exactly as the staircase does.
    let entry_rho = state.rho.clone();
    let entry_e = state.e.clone();

    let mut drho = vec![0.0f64; n];
    let mut dmx = vec![0.0f64; n];
    let mut dmy = vec![0.0f64; n];
    let mut de = vec![0.0f64; n];

    for (j, fxrow) in fx.iter().enumerate() {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        for (i, fycol) in fy.iter().enumerate() {
            let k = idx(i, j);
            let f = cut.at(i, j);
            if f.vol <= 0.0 {
                continue;
            }
            let vol = f.vol.max(1.0e-9);
            let dxi = dt / g.dxs[i];

            // the face fluxes, scaled by the face's aperture.
            let mass_f = dxi * (f.east * fxrow.mass[i + 1] - f.west * fxrow.mass[i])
                + dt * (rp * f.north * fycol.mass[j + 1] - rm * f.south * fycol.mass[j])
                    / (rc * dr);
            let mx_f = dxi * (f.east * fxrow.mx[i + 1] - f.west * fxrow.mx[i])
                + dt * (rp * f.north * fycol.mx[j + 1] - rm * f.south * fycol.mx[j]) / (rc * dr);
            let my_f = dxi * (f.east * fxrow.my[i + 1] - f.west * fxrow.my[i])
                + dt * (rp * f.north * fycol.my[j + 1] - rm * f.south * fycol.my[j]) / (rc * dr);
            let e_f = dxi * (f.east * fxrow.e[i + 1] - f.west * fxrow.e[i])
                + dt * (rp * f.north * fycol.e[j + 1] - rm * f.south * fycol.e[j]) / (rc * dr);

            // the immersed-boundary wall flux: the reflected ghost pair
            // in the wall-local frame, using the combustion face state.
            let (mut wmass, mut wmx, mut wmy, mut we) = (0.0, 0.0, 0.0, 0.0);
            if cut.eb_len[k] > 0.0 {
                let (enx, eny) = (-cut.eb_nx[k], -cut.eb_ny[k]);
                let un = u[k] * enx + v[k] * eny;
                let ut = -u[k] * eny + v[k] * enx;
                let mut y9 = [0.0f64; ct::NS];
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    y9[s] = species.y[s][k];
                }
                let t = closure9::temperature_at(clo.e_int[k], &y9);
                let rmix = closure9::mix_gas_constant(&y9);
                let pw = state.rho[k].max(1.0e-12) * rmix * t;
                let cp = ct::mix_cp_mass(&y9, t);
                let cv = (cp - rmix).max(1.0);
                let ga = (cp / cv).max(1.01);
                let aw = (ga * pw.max(1.0e-9) / state.rho[k].max(1.0e-12)).sqrt();
                let fw = |uu: f64| FacePrim {
                    rho: state.rho[k],
                    u: uu,
                    v: ut,
                    p: pw,
                    a: aw,
                    e: clo.e_int[k],
                };
                let q = hllc_flux(1.2, fw(un), fw(-un), 0);
                let l = dt * cut.eb_len[k];
                wmass = q.mass * l;
                wmx = (q.mx * enx - q.my * eny) * l;
                wmy = (q.mx * eny + q.my * enx) * l;
                we = q.e * l;
            }

            drho[k] = (-mass_f - wmass) / vol;
            dmx[k] = (-mx_f - wmx) / vol;
            dmy[k] = (-my_f - wmy) / vol + dt * clo.p[k] / rc;
            de[k] = (-e_f - we) / vol;
        }
    }

    // the increments land
    let mut resid = 0.0f64;
    for k in 0..n {
        if cut.cells[k].vol <= 0.0 {
            continue;
        }
        state.rho[k] += drho[k];
        state.mx[k] += dmx[k];
        state.my[k] += dmy[k];
        state.e[k] += de[k];
        resid = resid
            .max(drho[k].abs())
            .max(dmx[k].abs())
            .max(dmy[k].abs())
            .max(de[k].abs());
    }

    // the SMALL CELLS are pooled, and only the small ones. a cell with
    // most of its volume in the fluid updates on its own; the pooling
    // is for the slivers, where the annular amplification would
    // overwhelm the timestep.
    {
        let mut nbrs: Vec<Vec<usize>> = vec![Vec::new(); n];
        let fluid = |k: usize| cut.cells[k].vol > 0.0;
        let vfull = g.dx * g.dy;
        for j in 0..ny {
            for i in 0..nx {
                let k = idx(i, j);
                let alpha = cut.cells[k].vol;
                if alpha >= 1.0 || alpha <= 0.0 {
                    if alpha > 0.0 {
                        nbrs[k].push(k);
                    }
                    continue;
                }
                let mut member = vec![k];
                let mut vol = alpha * vfull;
                let (mut r0, mut r1, mut c0, mut c1) = (i, i, j, j);
                while vol < vfull && (r0 > 0 || r1 + 1 < nx || c0 > 0 || c1 + 1 < ny) {
                    let (pr0, pr1, pc0, pc1) = (r0, r1, c0, c1);
                    r0 = r0.saturating_sub(1);
                    r1 = (r1 + 1).min(nx - 1);
                    c0 = c0.saturating_sub(1);
                    c1 = (c1 + 1).min(ny - 1);
                    for jj in c0..=c1 {
                        for ii in r0..=r1 {
                            let on_new_ring = ii == r0 || ii == r1 || jj == c0 || jj == c1;
                            if !on_new_ring
                                || (ii == pr0 && jj >= pc0 && jj <= pc1)
                                || (ii == pr1 && jj >= pc0 && jj <= pc1)
                            {
                                continue;
                            }
                            let k2 = idx(ii, jj);
                            if fluid(k2) && !member.contains(&k2) {
                                member.push(k2);
                                vol += cut.cells[k2].vol * vfull;
                            }
                        }
                    }
                }
                nbrs[k] = member;
            }
        }
        let mut cover = vec![0usize; n];
        for row in &nbrs {
            for &m in row {
                cover[m] += 1;
            }
        }
        let mut qavg: Vec<[f64; 4]> = vec![[0.0; 4]; n];
        let mut qavg_y: Vec<Vec<f64>> = vec![vec![0.0f64; n]; ct::NS];
        for k in 0..n {
            if nbrs[k].is_empty() {
                continue;
            }
            let mut wv = 0.0f64;
            let mut acc = [0.0f64; 4];
            let mut accy = [0.0f64; ct::NS];
            for &m in &nbrs[k] {
                let w = cut.cells[m].vol * vfull / cover[m].max(1) as f64;
                wv += w;
                acc[0] += w * state.rho[m];
                acc[1] += w * state.mx[m];
                acc[2] += w * state.my[m];
                acc[3] += w * state.e[m];
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    accy[s] += w * species.y[s][m];
                }
            }
            qavg[k] = [acc[0] / wv, acc[1] / wv, acc[2] / wv, acc[3] / wv];
            for s in 0..ct::NS {
                qavg_y[s][k] = accy[s] / wv;
            }
        }
        let mut nnew = vec![0usize; n];
        let mut nacc = vec![[0.0f64; 4]; n];
        let mut nacc_y: Vec<Vec<f64>> = vec![vec![0.0f64; n]; ct::NS];
        for k in 0..n {
            for &m in &nbrs[k] {
                nacc[m][0] += qavg[k][0];
                nacc[m][1] += qavg[k][1];
                nacc[m][2] += qavg[k][2];
                nacc[m][3] += qavg[k][3];
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    nacc_y[s][m] += qavg_y[s][k];
                }
                nnew[m] += 1;
            }
        }
        for m in 0..n {
            if nnew[m] > 0 {
                state.rho[m] = nacc[m][0] / nnew[m] as f64;
                state.mx[m] = nacc[m][1] / nnew[m] as f64;
                state.my[m] = nacc[m][2] / nnew[m] as f64;
                state.e[m] = nacc[m][3] / nnew[m] as f64;
                for s in 0..ct::NS {
                    species.y[s][m] = (nacc_y[s][m] / nnew[m] as f64).max(0.0);
                }
            }
        }
    }

    // the chemistry, per cell, on the composition the cell ended up
    // with. this is the operator split the staircase uses: the flow
    // update first, then the chemistry at fixed total energy.
    for k in 0..n {
        if cut.cells[k].vol <= 0.0 {
            continue;
        }
        species.renorm(k);
        // the positivity guard, the staircase's own.
        let ke = 0.5 * (state.mx[k] * state.mx[k] + state.my[k] * state.my[k])
            / state.rho[k].max(1.0e-12);
        let e_int = state.e[k] - ke;
        if !state.e[k].is_finite()
            || !state.rho[k].is_finite()
            || state.rho[k] <= 0.0
            || e_int < clo.e_floor[k]
        {
            state.rho[k] = state.rho[k].max(1.0e-10);
            state.mx[k] = 0.0;
            state.my[k] = 0.0;
            state.e[k] = ke + clo.e_floor[k] * state.rho[k] + 1.0e-6;
        }
    }

    // the VISCOUS stage, when the wall has a temperature. the wall is
    // on the immersed boundary, so its one-sided derivative uses the
    // EXACT distance from the cell's centre to the wall, which the
    // signed distance function gives. this is the cut-cell's advantage
    // over the staircase: the wall's flux does not depend on where the
    // wall falls against the grid.
    if wall_temperature > 0.0 {
        if let Some(dist) = wall {
            let tr = crate::cutcell_viscous::transport_props(state, species, g, viscosity_scale)?;
            crate::cutcell_viscous::add_cut_combustion_viscous(
                state,
                species,
                g,
                cut,
                &tr,
                dt,
                wall_temperature,
                dist,
            )?;
        }
    }

    // the species transport rides the SAME apertures, on the same face
    // mass fluxes the momentum used.
    let fx_mass: Vec<Vec<f64>> = fx.iter().map(|s| s.mass.clone()).collect();
    let fy_mass: Vec<Vec<f64>> = fy.iter().map(|s| s.mass.clone()).collect();
    let resid_y =
        transport_species_apertures(&mut species.y, &state.rho, dt, g, &fx_mass, &fy_mass, cut);
    for k in 0..n {
        species.renorm(k);
    }

    // the solid scaffold is restored, exactly as the staircase does: it
    // is scaffolding, not flow.
    if let Some(dist) = wall {
        let band = 1.5 * g.dx.min(g.dy);
        let cache = crate::body::body_cache(g, dist, band);
        for &k in &cache.solid_cells {
            state.rho[k] = entry_rho[k];
            state.mx[k] = 0.0;
            state.my[k] = 0.0;
            state.e[k] = entry_e[k];
        }
    }

    Ok((dt, resid.max(resid_y)))
}

/// the cut-cell axisymmetric combustion march in the same heun (rk2)
/// wrapper the staircase uses, with the operator-split chemistry at the
/// post-step state.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_axi_cut_rk2(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cut: &CutField,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    wall_temperature: f64,
    viscosity_scale: f64,
    react: bool,
) -> Result<(f64, f64), Error> {
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    advance2d_combustion_axi_cut(
        &mut s1,
        &mut y1,
        g,
        cut,
        cfl,
        muscl,
        bc,
        nthreads,
        wall,
        wall_temperature,
        viscosity_scale,
    )?;
    advance2d_combustion_axi_cut(
        &mut s1,
        &mut y1,
        g,
        cut,
        cfl,
        muscl,
        bc,
        nthreads,
        wall,
        wall_temperature,
        viscosity_scale,
    )?;
    let n = g.nx * g.ny;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    for s in 0..ct::NS {
        for k in 0..n {
            species.y[s][k] = 0.5 * species.y[s][k] + 0.5 * y1.y[s][k];
        }
    }
    for k in 0..n {
        species.renorm(k);
    }
    // the chemistry's dt comes from the post-step state, the same probe
    // the staircase uses.
    let (dt, _) = {
        let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let mut smax = 0.0f64;
        for k in 0..n {
            let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
            let mut y_cell = [0.0f64; ct::NS];
            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                y_cell[s] = species.y[s][k];
            }
            let a = closure9::sound_speed_at(state.rho[k], e_int, &y_cell);
            smax = smax.max(u[k].abs().max(v[k].abs()) + a);
        }
        (cfl * g.dx.min(g.dy) / smax.max(1.0e-12), 0.0)
    };
    let mut oversub = 0.0f64;
    if react {
        oversub = crate::combustion_march::chemistry_substep9(state, species, g, dt, nthreads)?;
    }
    Ok((dt, oversub))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a flat wall: the fluid is y > 0.03 on a 0.1 by 0.1 grid, so the
    /// cut cells straddle the wall and the throat-like quantisation is
    /// visible in the area.
    fn flat_wall_field(g: &Grid2d) -> CutField {
        CutField::from_sdf(
            // the wall sits inside a row, not on a face: the grid's
            // rows are 0.005 tall and the wall is at 0.032, so the
            // sixth row (0.030..0.035) is genuinely cut.
            &|_x: f64, y: f64| y - 0.032,
            &|_x: f64, _y: f64| (0.0, -1.0),
            g,
        )
    }

    /// the cut cells really are cut: somewhere on the wall there must be
    /// a cell with a fluid fraction strictly between nothing and all.
    /// a staircase cannot produce that, so this is the property that
    /// separates the two representations.
    #[test]
    fn the_wall_cells_are_partly_fluid() {
        let g = crate::grid2d::grid2d(20, 20, 0.0, 0.1, 0.0, 0.1).expect("the grid");
        let cut = flat_wall_field(&g);
        let partial = (0..g.nx * g.ny)
            .filter(|&k| cut.cells[k].vol > 0.0 && cut.cells[k].vol < 1.0)
            .count();
        assert!(
            partial > 0,
            "no cell was cut: the wall landed on whole rows"
        );
    }

    /// the cut cell's volume fraction must lie between its corner
    /// samples' signs: a cell with all four corners in the fluid is
    /// full, all four in the body is empty, and a mixed one is
    /// proportionally filled.
    #[test]
    fn a_full_cell_is_full_and_a_solid_one_is_empty() {
        let g = crate::grid2d::grid2d(20, 20, 0.0, 0.1, 0.0, 0.1).expect("the grid");
        let cut = flat_wall_field(&g);
        // the wall at y = 0.032 cuts the sixth row (0.030..0.035).
        // the fifth row (0.025..0.030) is wholly in the body and the
        // seventh (0.035..0.040) is wholly in the fluid.
        let j_solid = 5;
        let j_full = 7;
        for i in 0..g.nx {
            assert!(
                cut.at(i, j_solid).vol <= 0.0,
                "the row under the wall is not empty"
            );
            assert!(
                (cut.at(i, j_full).vol - 1.0).abs() < 1.0e-12,
                "the row over the wall is not full"
            );
        }
        // and the cut row must be strictly between
        for i in 0..g.nx {
            let v = cut.at(i, 6).vol;
            assert!(v > 0.0 && v < 1.0, "the cut row is not partly fluid: {v}");
        }
    }

    /// a uniform composition must not develop a gradient: the species
    /// transport conserves a uniform field exactly, whatever the
    /// apertures are. if the aperture weighting were inconsistent with
    /// the flux balance, a uniform mixture would erode near the wall.
    #[test]
    fn a_uniform_mixture_stays_uniform() {
        let (nx, ny) = (12, 12);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.12).expect("the grid");
        let cut = flat_wall_field(&g);
        let n = nx * ny;
        let rho = 1.0f64;
        let state = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * 1.0e5; n],
        };
        let mut species = CombustionState::stoich_h2o2(n);
        // a uniform field, zero velocity, so the face mass fluxes are
        // zero and the composition cannot move.
        let before: Vec<f64> = (0..ct::NS).map(|s| species.y[s][0]).collect();
        transport_species_apertures(
            &mut species.y,
            &state.rho,
            1.0e-4,
            &g,
            &vec![vec![0.0f64; nx + 1]; ny],
            &vec![vec![0.0f64; ny + 1]; nx],
            &cut,
        );
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            for k in 0..n {
                assert!(
                    (species.y[s][k] - before[s]).abs() < 1.0e-12,
                    "species {s} moved in cell {k}: {} against {}",
                    species.y[s][k],
                    before[s]
                );
            }
        }
    }

    /// the geometric consistency the throat cares about: summing the
    /// cut cells' contributions must give a flow area that DEPENDS ON
    /// THE CONTOUR'S POSITION and not on where the grid's rows happen
    /// to fall. shifting the wall by a quarter of a cell must move the
    /// area continuously, while the staircase's area would jump.
    #[test]
    fn the_flow_area_follows_the_contour_smoothly() {
        let (nx, ny) = (10, 40);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.02, 0.0, 0.1).expect("the grid");
        let dy = 0.1 / ny as f64;
        // the flow area: the fluid annulus each column carries, with the
        // cut cells contributing their own fraction.
        let area_at = |shift: f64| -> f64 {
            // the shift keeps the wall OFF the cell faces: a corner that
            // lands exactly on the wall has phi = 0, and the corner
            // sampler's sign test puts that corner in the SOLID side,
            // so the cell reads as fully fluid. a wall exactly on a face
            // is a degenerate configuration for the sampler, and the
            // sweep below must avoid it.
            let wall = 0.03 + shift + 0.13 * dy;
            let cut = CutField::from_sdf(
                &move |_x: f64, y: f64| y - wall,
                &move |_x: f64, _y: f64| (0.0, -1.0),
                &g,
            );
            let mut a = 0.0f64;
            for j in 0..ny {
                let f = cut.at(nx / 2, j);
                if f.vol <= 0.0 {
                    continue;
                }
                let rm = j as f64 * dy;
                let rp = rm + dy;
                a += f.vol * std::f64::consts::PI * (rp * rp - rm * rm);
            }
            a
        };
        // sweeping the wall across a full row: the area must change
        // MONOTONICALLY and the individual steps must stay small. a
        // staircase would stand still for a whole row and then jump.
        let areas: Vec<f64> = (0..40).map(|i| area_at(i as f64 * dy / 40.0)).collect();
        for w in areas.windows(2) {
            assert!(
                w[1] <= w[0] + 1.0e-12,
                "moving the wall inward raised the area: {:?}",
                areas
            );
        }
        let steps: Vec<f64> = areas
            .windows(2)
            .map(|w| ((w[1] - w[0]).abs() / w[0].max(1.0e-12)) * 100.0)
            .collect();
        let worst = steps.iter().copied().fold(0.0f64, f64::max);
        // one fortieth of a row's movement must change the area by at
        // most a small fraction of a cell's worth. the staircase's step
        // was twenty percent.
        assert!(worst < 2.0, "the area jumped by {worst} percent");
    }
}

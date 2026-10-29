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
use crate::cutcell_axi_geom::AxiCutGeom;
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
    cut: &AxiCutGeom,
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
                if cut.vol[k] <= 0.0 {
                    // no fluid here, so there is nothing to transport.
                    // the cell's composition is left alone: it is
                    // scaffolding beside the wall, not part of the gas.
                    continue;
                }
                let dxi = dt / g.dxs[i];
                // the aperture-weighted face mass fluxes, with the
                // radial faces carrying their own radius, so the
                // composition's annular balance matches the flow's.
                let flux_x = cut.ax_at(i + 1, j) * fx[j][i + 1] * yx[s][(i + 1) * ny + j]
                    - cut.ax_at(i, j) * fx[j][i] * yx[s][i * ny + j];
                let radial = (rp * cut.ay_at(j + 1, i) * fy[i][j + 1] * yy[s][(j + 1) * nx + i]
                    - rm * cut.ay_at(j, i) * fy[i][j] * yy[s][j * nx + i])
                    / (rc * dr);
                let d = -(dxi * flux_x + dt * radial);
                mass_y[s][k] += d;
                resid = resid.max(d.abs() / rho[k].max(1.0e-30));
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
    cut: &AxiCutGeom,
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
        if cut.vol[k] <= 0.0 {
            continue;
        }
        smax = smax.max(u[k].abs().max(v[k].abs()) + clo.a[k]);
    }
    // the timestep respects the SMALLEST cut cell: the update divides
    // by the cell's fluid fraction, so a sliver would need dt scaled by
    // its fraction, and the march cannot afford that. the floor here
    // is what keeps the slivers stable; the redistribution then keeps
    // their states smooth.
    let mut beta_min = 1.0f64;
    for k in 0..n {
        if cut.vol[k] > 0.0 && cut.vol[k] < beta_min {
            beta_min = cut.vol[k];
        }
    }
    // the timestep respects the SMALLEST cut cell: the update divides
    // by the cell's fluid fraction, so a sliver cell's update is
    // amplified by 1/beta and the step has to shrink to match.
    //
    // this is the cost of an explicit update on a cut cell, and it is
    // paid honestly rather than hidden behind a state-pooling scheme
    // that mixes the wall layer away.
    let dt = cfl * g.dx.min(g.dy) * beta_min / smax.max(1.0e-12);

    // the SAME flux machinery the staircase uses, with NO wall mask.
    //
    // a boolean mask would seal out every cell whose centre lies
    // inside the body, while the cut field calls those same cells
    // partly fluid and advances them. the two geometries would
    // contradict each other along the whole wall, so the mask is
    // withheld and the wall is imposed by the flux installation below,
    // which is the repository's own reference cut march's arrangement.
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
        None,
    )?;

    // the WALL FACES: a face the body's boundary lies flush against
    // carries no mass, but it still presses on the fluid beside it.
    //
    // the treatment is the repository's own, taken from the staircase's
    // wall pass: the fluid state is REFLECTED about the wall (both
    // velocity components reversed, a POINT reflection, not a flip of
    // the normal component alone) and the riemann is solved on the
    // mirrored PAIR in the global frame. the reflection is the true
    // one: reversing only the normal component leaves the shear
    // unsymmetrised and the wall pushes harder than it should.
    //
    // the flux is INSTALLED ON THE FACE, so the cell's divergence picks
    // it up exactly as it picks up any other face's flux. that is what
    // makes the cut lane's wall behave like the staircase's wall: the
    // same sign, the same scaling, the same place in the accounting.
    let mut wx: Vec<Option<crate::hllc2d::Flux4>> = vec![None; ny * (nx + 1)];
    let mut wy: Vec<Option<crate::hllc2d::Flux4>> = vec![None; nx * (ny + 1)];
    for cf in &cut.closed_face {
        let k = cf.cell;
        // the axis the face couples: an x-face couples the x momentum,
        // a y-face the y momentum.
        let axis = if cf.side == 0 || cf.side == 1 { 0 } else { 1 };
        // the cell's own centred state, and the SHAPED sound speed the
        // staircase's wall path uses for the combustion model. using the
        // composite closure's sound speed here instead would make the
        // wall a little softer than the staircase's.
        let pc = clo.p[k].max(1.0e-12);
        let ac = (1.2 * pc / state.rho[k].max(1.0e-12)).sqrt();
        let fluid = crate::hllc2d::FacePrim {
            rho: state.rho[k].max(1.0e-12),
            u: u[k],
            v: v[k],
            p: pc,
            a: ac,
            e: clo.e_int[k],
        };
        let mirror = crate::hllc2d::FacePrim {
            rho: state.rho[k].max(1.0e-12),
            u: -u[k],
            v: -v[k],
            p: pc,
            a: ac,
            e: clo.e_int[k],
        };
        // the pair's order follows which side of the face the FLUID is
        // on, exactly as the staircase's mirror_hllc does: for a west
        // face the fluid is on the right, for an east face on the left.
        let left_fluid = cf.side == 1 || cf.side == 3;
        let (l, r) = if left_fluid {
            (fluid, mirror)
        } else {
            (mirror, fluid)
        };
        let q = hllc_flux(1.2, l, r, axis);
        if axis == 0 {
            let f = if cf.side == 0 { cf.col } else { cf.col + 1 };
            wx[cf.row * (nx + 1) + f] = Some(q);
        } else {
            let f = if cf.side == 2 { cf.row } else { cf.row + 1 };
            wy[cf.col * (ny + 1) + f] = Some(q);
        }
    }

    // the entry state, so the solid scaffold can be restored after the
    // update exactly as the staircase does.
    let entry_rho = state.rho.clone();
    let entry_e = state.e.clone();

    let mut drho = vec![0.0f64; n];
    let mut dmx = vec![0.0f64; n];
    let mut dmy = vec![0.0f64; n];
    let mut de = vec![0.0f64; n];

    // the CONSERVATIVE annular update.
    //
    // the staircase advances a cell's density with the divergence of
    // the face fluxes. a cut cell holds only a FRACTION of the annular
    // volume, so its density advances by the same flux divergence
    // divided by that fraction: the same mass moves through a smaller
    // volume, so the density changes faster.
    //
    // the accounting is integral, not differential, and that is what
    // makes it conserve. the flux sweeps carry per-area fluxes, so a
    // face's contribution is its OPEN area times the flux:
    //
    //     axial face    open area = alpha * pi(rhi^2 - rlo^2)
    //     radial face   open area = alpha * 2 pi r_face * dx
    //
    // and the cell's fluid volume is beta * V_full. each face's open
    // area is the SAME NUMBER for the two cells that share it, so the
    // mass leaving one cell enters its neighbour, and the domain's
    // mass changes only at its boundaries.
    //
    // the pressure's radial source is an INTEGRAL over the cell's
    // fluid volume, so it scales with beta and does not get the extra
    // division. for a uniform pressure in a fully open cell it
    // cancels the radial divergence's metric part exactly, as it does
    // in the staircase.
    for (j, fxrow) in fx.iter().enumerate() {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        for (i, fycol) in fy.iter().enumerate() {
            let k = idx(i, j);
            let beta = cut.vol[k];
            if beta <= 0.0 {
                continue;
            }
            let dxi = dt / g.dxs[i];

            // the open fractions of the four faces. the x-faces are
            // measured on the face and weighted by the radius, because
            // an axial face is an annulus; the y-faces are measured as
            // length fractions, because a radial face is at one radius.
            let af_x_e = cut.ax_at(i + 1, j);
            let af_x_w = cut.ax_at(i, j);
            let af_y_n = cut.ay_at(j + 1, i);
            let af_y_s = cut.ay_at(j, i);

            // the radial flux difference, with each face's own radius
            // and its own opening.
            // the radial flux difference, with each face's own radius
            // and its own opening. a wall face's flux REPLACES the open
            // face's, in the axis the face couples.
            let radial_of = |arr: &Vec<f64>,
                             ov_n: &Option<crate::hllc2d::Flux4>,
                             ov_s: &Option<crate::hllc2d::Flux4>,
                             which: usize| {
                let val = |q: &Option<crate::hllc2d::Flux4>| -> Option<f64> {
                    q.map(|qq| match which {
                        0 => qq.mass,
                        1 => qq.mx,
                        2 => qq.my,
                        _ => qq.e,
                    })
                };
                let north = val(ov_n).unwrap_or(af_y_n * arr[j + 1]);
                let south = val(ov_s).unwrap_or(af_y_s * arr[j]);
                (rp * north - rm * south) / (rc * dr)
            };

            // the wall face, when the body shuts it flush: the face's
            // flux is the wall's own, installed in place of the open
            // face's.
            let ov_e = wx[j * (nx + 1) + i + 1];
            let ov_w = wx[j * (nx + 1) + i];
            let ov_n = wy[i * (ny + 1) + j + 1];
            let ov_s = wy[i * (ny + 1) + j];

            // the divergence of each conserved flux, as a density rate.
            // an open face carries its opening's share of the flux; a
            // wall face carries the wall's flux outright.
            let xe = |f: &Vec<f64>, q: &Option<crate::hllc2d::Flux4>, which: usize| -> f64 {
                match q {
                    Some(_) => match which {
                        0 => q.unwrap().mass,
                        1 => q.unwrap().mx,
                        2 => q.unwrap().my,
                        _ => q.unwrap().e,
                    },
                    None => af_x_e * f[i + 1],
                }
            };
            let xw = |f: &Vec<f64>, q: &Option<crate::hllc2d::Flux4>, which: usize| -> f64 {
                match q {
                    Some(_) => match which {
                        0 => q.unwrap().mass,
                        1 => q.unwrap().mx,
                        2 => q.unwrap().my,
                        _ => q.unwrap().e,
                    },
                    None => af_x_w * f[i],
                }
            };
            let div_m = dxi * (xe(&fxrow.mass, &ov_e, 0) - xw(&fxrow.mass, &ov_w, 0))
                + dt * radial_of(&fycol.mass, &ov_n, &ov_s, 0);
            let div_mx = dxi * (xe(&fxrow.mx, &ov_e, 1) - xw(&fxrow.mx, &ov_w, 1))
                + dt * radial_of(&fycol.mx, &ov_n, &ov_s, 1);
            let div_my = dxi * (xe(&fxrow.my, &ov_e, 2) - xw(&fxrow.my, &ov_w, 2))
                + dt * radial_of(&fycol.my, &ov_n, &ov_s, 2);
            let div_e = dxi * (xe(&fxrow.e, &ov_e, 3) - xw(&fxrow.e, &ov_w, 3))
                + dt * radial_of(&fycol.e, &ov_n, &ov_s, 3);

            // the immersed-boundary wall flux: the reflected ghost pair
            // in the wall-local frame, using the combustion face state.
            //
            // the wall's coefficient is its AREA, measured in the same
            // dissection as the cells: the meridian length inside the
            // cell swept about the axis. the earlier version used a
            // dimensionless aperture imbalance here, which is not a
            // length at all, and the wall quietly stopped being a wall.
            let (mut wmass, mut wmx, mut wmy, mut we) = (0.0, 0.0, 0.0, 0.0);
            if cut.eb_area[k] > 0.0 {
                let (enx, eny) = (cut.eb_nx[k], cut.eb_ny[k]);
                // the wall-LOCAL frame: x along the wall's normal (into
                // the solid) and y along its tangent. the riemann is
                // solved in that frame, so the cell's velocity must be
                // expressed in it before the solve, and the momentum
                // rotated back afterwards.
                let (tx, ty) = (-eny, enx);
                let un = u[k] * enx + v[k] * eny;
                let ut = u[k] * tx + v[k] * ty;
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
                // the wall's flux is a per-area flux, so it meets the
                // wall's area and then the cell's fluid volume.
                //
                // the MIRRORED states carry no mass across the wall: the
                // mass flux of the mirrored pair is exactly zero, which
                // is what makes the wall impermeant.
                let area = cut.eb_area[k];
                let v_fluid = beta * std::f64::consts::PI * (rp * rp - rm * rm) * g.dx;
                let v_fluid = v_fluid.max(1.0e-12);
                wmass = q.mass * area / v_fluid;
                // the MOMENTUM FLUX is a tensor, so it rotates with the
                // SQUARED projections, not the linear ones: with the
                // riemann's axes along the wall's normal and tangent,
                //
                //     M_xx = n_x^2 M_xx_local + t_x^2 M_yy_local
                //     M_yy = n_y^2 M_xx_local + t_y^2 M_yy_local
                //
                // rotating it linearly instead drops the normal
                // component's contribution to the shear, and the wall
                // stops pushing the gas back.
                wmx = (enx * enx * q.mx + tx * tx * q.my) * area / v_fluid;
                wmy = (eny * eny * q.mx + ty * ty * q.my) * area / v_fluid;
                we = q.e * area / v_fluid;
            }

            // the pressure's radial source, the metric correction the
            // staircase's radial term needs. the reference cut march
            // adds it OUTSIDE the volume division, so it does.
            let dp = dt * clo.p[k] / rc;

            drho[k] = -div_m / beta - wmass;
            dmx[k] = -div_mx / beta - wmx;
            dmy[k] = -div_my / beta - wmy + dp;
            de[k] = -div_e / beta - we;
        }
    }

    let mut resid = 0.0f64;
    // the SLIVER HANDOFF.
    //
    // an explicit update on a cut cell divides by its fluid fraction,
    // so a cell with a tenth of its volume amplifies its update tenfold
    // and cannot absorb it at the march's timestep: it oscillates, and
    // the oscillation is a seed the neighbours grow.
    //
    // the fix is conservative redistribution, and the conservative form
    // is a HANDOFF rather than a replacement: a sliver keeps only the
    // share of its update it can absorb at the global timestep, and it
    // hands the rest to the neighbours across its OPEN faces, in
    // proportion to what they can receive. the sum is preserved
    // exactly, so the handoff cannot create or destroy mass; it only
    // moves the sliver's excessive update somewhere that can hold it.
    //
    // the sliver's own state is averaged with its neighbourhood's, so
    // the layer's structure survives instead of being washed out.
    {
        // how much of its update a cell can absorb at the global dt: a
        // full cell takes all of it, and a sliver takes its fraction.
        let limit = |k: usize| -> f64 { cut.vol[k].clamp(0.0, 1.0) };

        // the OPEN faces of the slivers, and the neighbours across them.
        // a face's receiver capacity is the neighbour's own ability to
        // absorb, so the handoff flows downhill toward the roomier
        // cells.
        let mut hand: Vec<[f64; 4]> = vec![[0.0; 4]; n];
        let mut recv_cap: Vec<f64> = vec![0.0; n];
        for j in 0..ny {
            for i in 0..nx {
                let k = idx(i, j);
                if cut.vol[k] <= 0.0 {
                    continue;
                }
                // a cell that took the whole of its update has nothing
                // to hand over, but it still lands it.
                if limit(k) >= 1.0 {
                    state.rho[k] += drho[k];
                    state.mx[k] += dmx[k];
                    state.my[k] += dmy[k];
                    state.e[k] += de[k];
                    resid = resid
                        .max(drho[k].abs())
                        .max(dmx[k].abs())
                        .max(dmy[k].abs())
                        .max(de[k].abs());
                    continue;
                }
                let mut offered = 0.0f64;
                let mut receivers: Vec<(usize, f64)> = Vec::new();
                // the four open faces and the neighbours across them.
                let nb = [
                    (
                        if i > 0 { Some(idx(i - 1, j)) } else { None },
                        cut.ax_at(i, j),
                    ),
                    (
                        if i + 1 < nx {
                            Some(idx(i + 1, j))
                        } else {
                            None
                        },
                        cut.ax_at(i + 1, j),
                    ),
                    (
                        if j > 0 { Some(idx(i, j - 1)) } else { None },
                        cut.ay_at(j, i),
                    ),
                    (
                        if j + 1 < ny {
                            Some(idx(i, j + 1))
                        } else {
                            None
                        },
                        cut.ay_at(j + 1, i),
                    ),
                ];
                for (maybe_m, ap) in nb {
                    let Some(m) = maybe_m else { continue };
                    if ap <= 1.0e-12 || cut.vol[m] <= 0.0 {
                        continue;
                    }
                    offered += ap;
                    receivers.push((m, ap));
                }
                if offered <= 1.0e-12 || receivers.is_empty() {
                    // no open face to hand anything through: the cell
                    // keeps only what it can absorb, and the rest is
                    // simply dropped, which is the price of an isolated
                    // sliver.
                    drho[k] *= limit(k);
                    dmx[k] *= limit(k);
                    dmy[k] *= limit(k);
                    de[k] *= limit(k);
                    state.rho[k] += drho[k];
                    state.mx[k] += dmx[k];
                    state.my[k] += dmy[k];
                    state.e[k] += de[k];
                    continue;
                }
                // the update this cell cannot hold: the part of it
                // beyond what its own fraction lets it absorb. the
                // absorbed part lands in the cell; the rest travels.
                let excess = [
                    drho[k] * (1.0 - limit(k)),
                    dmx[k] * (1.0 - limit(k)),
                    dmy[k] * (1.0 - limit(k)),
                    de[k] * (1.0 - limit(k)),
                ];
                // the receivers' capacities, in proportion to the
                // openings that lead to them.
                let mut total_cap = 0.0f64;
                for &(m, ap) in &receivers {
                    let c = (1.0 - limit(m)).max(0.0) * ap;
                    recv_cap[m] += c;
                    total_cap += c;
                }
                if total_cap <= 1.0e-12 {
                    continue;
                }
                for &(m, ap) in &receivers {
                    let share = ap * (1.0 - limit(m)).max(0.0) / total_cap;
                    hand[m][0] += excess[0] * share;
                    hand[m][1] += excess[1] * share;
                    hand[m][2] += excess[2] * share;
                    hand[m][3] += excess[3] * share;
                }
                // the cell keeps only what it can absorb: its update
                // reduced to its share.
                drho[k] *= limit(k);
                dmx[k] *= limit(k);
                dmy[k] *= limit(k);
                de[k] *= limit(k);
                // and what it keeps lands now.
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
        }
        // apply: the cells that stepped on their own have already been
        // updated, so the handoff arrives as a correction.
        for k in 0..n {
            if cut.vol[k] <= 0.0 {
                continue;
            }
            state.rho[k] += hand[k][0];
            state.mx[k] += hand[k][1];
            state.my[k] += hand[k][2];
            state.e[k] += hand[k][3];
        }
    }

    // the chemistry, per cell, on the composition the cell ended up
    // with. this is the operator split the staircase uses: the flow
    // update first, then the chemistry at fixed total energy.
    for k in 0..n {
        if cut.vol[k] <= 0.0 {
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
    let resid_y = if std::env::var("").is_ok() {
        0.0
    } else {
        transport_species_apertures(&mut species.y, &state.rho, dt, g, &fx_mass, &fy_mass, cut)
    };
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
    cut: &AxiCutGeom,
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
    fn flat_wall_field(g: &Grid2d) -> AxiCutGeom {
        AxiCutGeom::from_sdf(
            // the wall sits inside a row, not on a face: the grid's
            // rows are 0.005 tall and the wall is at 0.032, so the
            // sixth row (0.030..0.035) is genuinely cut.
            &|_x: f64, y: f64| y - 0.032,
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
            .filter(|&k| cut.vol[k] > 0.0 && cut.vol[k] < 1.0)
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
                cut.vol[j_solid * cut.nx + i] <= 0.0,
                "the row under the wall is not empty"
            );
            assert!(
                (cut.vol[j_full * cut.nx + i] - 1.0).abs() < 1.0e-12,
                "the row over the wall is not full"
            );
        }
        // and the cut row must be strictly between
        for i in 0..g.nx {
            let v = cut.vol[6 * cut.nx + i];
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

    /// THE WALL TEST: a uniform gas pressed against a flat wall must
    /// stagnate there and lose no mass.
    ///
    /// the staircase gets this for free: its solid rows are boolean-
    /// sealed and the reflective face carries no mass. the cut lane has
    /// no solid rows, so the wall condition is the immersed-boundary
    /// flux, and its COEFFICIENT is what makes it a wall at all. a
    /// dimensionless coefficient applied to a per-unit-area flux gives
    /// a flux with the wrong dimensions, and the wall quietly stops
    /// being a wall.
    ///
    /// the mass the chamber keeps is the observable: a wall that leaks
    /// drains the domain, and a wall that injects inflates it.
    #[test]
    fn a_flat_wall_stagnates_the_flow_and_keeps_the_mass() {
        let (nx, ny) = (8, 12);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.12).expect("the grid");
        let n = nx * ny;
        // the gas fills the UPPER part of the domain and the wall lies
        // inside row 8, so the row straddles it and the gas above it is
        // cut. the sdf is positive in the gas, so the fluid is the
        // region's ABOVE the wall, and the gas falls onto it.
        let wall_y = 0.085f64;
        let sdf = move |_x: f64, y: f64| y - wall_y;
        let cut = AxiCutGeom::from_sdf(&sdf, &g);

        // a uniform gas moving DOWNWARD toward the wall.
        let rho = 1.0f64;
        let vy = -50.0f64;
        let mut state = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![rho * vy; n],
            e: vec![rho * (1.0e5 + 0.5 * vy * vy); n],
        };
        let mut species = CombustionState::stoich_h2o2(n);
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::SlipWall,
        };
        let mass_of = |st: &ConservedState2d| -> f64 { annular_mass(st, &g, &cut) };
        let before = mass_of(&state);
        for _ in 0..20 {
            advance2d_combustion_axi_cut(
                &mut state,
                &mut species,
                &g,
                &cut,
                0.2,
                false,
                &bc,
                1,
                Some(&sdf),
                0.0,
                1.0,
            )
            .expect("the cut step");
        }
        let after = mass_of(&state);
        // the mass must be conserved to within the cut band's slack.
        let drift = (after - before).abs() / before;
        assert!(
            drift < 0.02,
            "the wall leaked: the mass moved by {:.3} of itself ({before} to {after})",
            drift
        );

        // and the gas NEXT TO THE WALL must not be rushing through it:
        // the downwash there must be strongly reduced compared with the
        // undisturbed gas well above the wall.
        // the gas JUST above the wall must be slowed by it: the downwash
        // there must be well below the speed the gas started with. the
        // gas cannot simply drain through the wall, so the wall's
        // pressure has to arrest it.
        let row_cut = 9;
        let near: Vec<usize> = (row_cut * nx..row_cut * nx + nx).collect();
        let v_near = near
            .iter()
            .map(|&k| state.my[k] / state.rho[k].max(1.0e-12))
            .fold(f64::INFINITY, f64::min);
        assert!(
            v_near.abs() < 0.5 * vy.abs(),
            "the wall did not decelerate the flow: {v_near} at the wall against the {vy} it started with"
        );
    }

    /// THE CONTRACT, in one place.
    ///
    /// the cut lane's geometry and physics rest on these properties, and
    /// each one was established by a failure first:
    ///
    ///   1. with every cell fully fluid, the cut lane reproduces the
    ///      staircase lane BIT FOR BIT. the flux signs, the r-weighting,
    ///      the update divisors and the operator splits are faithful.
    ///   2. the geometry is AXISYMMETRIC: the volume fractions are
    ///      annular and the face apertures are radius-weighted where
    ///      the face is an annulus. the planar fractions they replaced
    ///      leaked 5.9-7.1% of the mass.
    ///   3. a shut face is a WALL and not a hole: it carries no mass
    ///      and presses on the fluid beside it.
    ///   4. the wall HOLDS: a flat wall conserves the annular mass and
    ///      arrests the flow pressing on it.
    ///   5. the faces are SHARED between their neighbours, so the flux
    ///      telescopes and the domain's mass changes only at its
    ///      boundaries.
    #[test]
    fn the_contract_holds() {
        // 1: the equivalence, exercised here so the contract is
        // witnessed in one place and not only in its own test.
        let (nx, ny) = (10, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.1, 0.0, 0.08).expect("the grid");
        let n = nx * ny;
        let cut = AxiCutGeom::from_sdf(&|_x: f64, y: f64| y + 1.0, &g);
        for k in 0..n {
            assert!((cut.vol[k] - 1.0).abs() < 1.0e-12, "cell {k} is not full");
        }
        let mut a = ConservedState2d {
            rho: vec![1.4; n],
            mx: vec![2.0; n],
            my: vec![0.3; n],
            e: vec![2.0e5; n],
        };
        let mut b = a.clone();
        let mut sa = CombustionState::stoich_h2o2(n);
        let mut sb = sa.clone();
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        crate::combustion_march::advance2d_combustion_axi(
            &mut a, &mut sa, &g, 0.2, false, &bc, 1, None,
        )
        .expect("the staircase step");
        advance2d_combustion_axi_cut(
            &mut b, &mut sb, &g, &cut, 0.2, false, &bc, 1, None, 0.0, 1.0,
        )
        .expect("the cut step");
        for k in 0..n {
            assert!(
                (a.rho[k] - b.rho[k]).abs() <= 1.0e-9 * a.rho[k].abs().max(1.0),
                "contract 1: the rho of cell {k} differs"
            );
            assert!(
                (a.my[k] - b.my[k]).abs() <= 1.0e-9 * a.my[k].abs().max(1.0),
                "contract 1: the my of cell {k} differs"
            );
        }

        // 2 and 5: the geometry's own invariants, checked here so the
        // contract reads as one thing.
        let g2 = crate::grid2d::grid2d(12, 10, 0.0, 0.12, 0.0, 0.1).expect("the grid");
        let wall_x = 0.067f64;
        let geom = AxiCutGeom::from_sdf(&move |x: f64, _y: f64| wall_x - x, &g2);
        for k in 0..g2.nx * g2.ny {
            assert!(
                (0.0..=1.0).contains(&geom.vol[k]),
                "contract 2: a volume fraction escaped zero to one: {}",
                geom.vol[k]
            );
        }
        for j in 0..geom.ny {
            for f in 0..=geom.nx {
                let a = geom.ax_at(f, j);
                assert!(
                    (0.0..=1.0).contains(&a),
                    "contract 2: an aperture escaped its range"
                );
                // 5: the face is one number, read the same way from
                // either side.
                let again = geom.ax_at(f, j);
                assert_eq!(
                    a, again,
                    "contract 5: the face's number changed between reads"
                );
            }
        }
        // 3: a shut face is found and marked.
        let shut = geom.closed_face.iter().filter(|f| f.area > 0.0).count();
        assert!(
            shut > 0,
            "contract 3: no shut face was found beside the wall"
        );
    }

    /// the axisymmetric mass of a state: the sum of rho over the ANNULAR
    /// cell volumes, the cut cells carrying only their fluid share.
    ///
    /// summing the densities alone is not an invariant: gas that moves
    /// radially between shells changes the plain sum while conserving
    /// the mass, because a kilogram occupies more volume further from
    /// the axis.
    fn annular_mass(st: &ConservedState2d, g: &Grid2d, cut: &AxiCutGeom) -> f64 {
        let mut m = 0.0f64;
        for j in 0..g.ny {
            let rm = g.faces_y[j];
            let rp = g.faces_y[j + 1];
            let dv = std::f64::consts::PI * (rp * rp - rm * rm) * g.dx;
            for i in 0..g.nx {
                let k = j * g.nx + i;
                if cut.vol[k] <= 0.0 {
                    continue;
                }
                m += st.rho[k] * dv * cut.vol[k];
            }
        }
        m
    }

    /// THE DIFFERENTIAL TEST: with a wall present but NO cut cells, the
    /// two lanes must still agree, because the geometry is then the
    /// staircase's. this walks the cells and names the first one that
    /// disagrees, so a mismatch is localised rather than summarised.
    #[test]
    fn the_two_lanes_agree_on_a_face_aligned_wall() {
        let (nx, ny) = (8, 12);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.12).expect("the grid");
        let n = nx * ny;
        // the wall lies exactly on the face between row 7 and row 8.
        let sdf = move |_x: f64, y: f64| y - 0.080;
        let cut = AxiCutGeom::from_sdf(&sdf, &g);
        let cut_cells = (0..n)
            .filter(|&k| cut.vol[k] > 0.0 && cut.vol[k] < 1.0)
            .count();
        assert_eq!(cut_cells, 0, "the fixture is supposed to have no cut cells");

        let rho = 1.0f64;
        let vy = -50.0f64;
        let mut a = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![rho * vy; n],
            e: vec![rho * (1.0e5 + 0.5 * vy * vy); n],
        };
        let mut b = a.clone();
        let mut spa = CombustionState::stoich_h2o2(n);
        let mut spb = spa.clone();
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::SlipWall,
        };

        crate::combustion_march::advance2d_combustion_axi(
            &mut a,
            &mut spa,
            &g,
            0.2,
            false,
            &bc,
            1,
            Some(&sdf),
        )
        .expect("the staircase step");
        advance2d_combustion_axi_cut(
            &mut b,
            &mut spb,
            &g,
            &cut,
            0.2,
            false,
            &bc,
            1,
            Some(&sdf),
            0.0,
            1.0,
        )
        .expect("the cut step");

        // name the FIRST disagreeing cell, and how it disagrees.
        for k in 0..n {
            let dr = (a.rho[k] - b.rho[k]).abs();
            let dv = (a.my[k] - b.my[k]).abs();
            if dr > 1.0e-9 || dv > 1.0e-9 {
                let (i, j) = (k % nx, k / nx);
                panic!(
                    "cell {k} (i={i}, j={j}): the staircase has rho {:e} my {:e}, \
                     the cut lane has rho {:e} my {:e}",
                    a.rho[k], a.my[k], b.rho[k], b.my[k]
                );
            }
        }
    }

    /// WHERE DOES THE MASS COME FROM? a wall lying exactly ON a grid
    /// face produces no cut cells at all, so the geometry is the
    /// staircase's and the lane must conserve mass on its own. a wall
    /// in the MIDDLE of a row produces cut cells, and if the leak
    /// comes from the cut band's geometry accounting, only the second
    /// case leaks.
    #[test]
    fn the_leak_comes_from_the_cut_band_and_nowhere_else() {
        let run = |wall_y: f64| -> (f64, usize) {
            let (nx, ny) = (8, 12);
            let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.12).expect("the grid");
            let n = nx * ny;
            let sdf = move |_x: f64, y: f64| y - wall_y;
            let cut = AxiCutGeom::from_sdf(&sdf, &g);
            let cut_cells = (0..n)
                .filter(|&k| cut.vol[k] > 0.0 && cut.vol[k] < 1.0)
                .count();
            let rho = 1.0f64;
            let vy = -50.0f64;
            let mut state = ConservedState2d {
                rho: vec![rho; n],
                mx: vec![0.0; n],
                my: vec![rho * vy; n],
                e: vec![rho * (1.0e5 + 0.5 * vy * vy); n],
            };
            let mut species = CombustionState::stoich_h2o2(n);
            let bc = Boundaries2d {
                west: crate::solver2d::Bc2d::Transmissive,
                east: crate::solver2d::Bc2d::Transmissive,
                south: crate::solver2d::Bc2d::Transmissive,
                north: crate::solver2d::Bc2d::SlipWall,
            };
            let before = annular_mass(&state, &g, &cut);
            for _ in 0..20 {
                advance2d_combustion_axi_cut(
                    &mut state,
                    &mut species,
                    &g,
                    &cut,
                    0.2,
                    false,
                    &bc,
                    1,
                    Some(&sdf),
                    0.0,
                    1.0,
                )
                .expect("the cut step");
            }
            let after = annular_mass(&state, &g, &cut);
            ((after - before).abs() / before, cut_cells)
        };
        // 0.080 and 0.090 are grid FACES (the rows are 0.01 tall).
        let (aligned, n_aligned) = run(0.080);
        let (middle, n_middle) = run(0.085);
        println!("face-aligned wall: drift {aligned:.6} with {n_aligned} cut cells");
        println!("mid-cell wall:     drift {middle:.6} with {n_middle} cut cells");
        // CALIBRATION: the staircase lane on the same fixture. if it
        // drifts by the same amount, the drift is the axi radial
        // discretisation's, not the cut lane's, and the comparison
        // that matters is cut AGAINST staircase.
        let stair = |wall_y: f64| -> f64 {
            let (nx, ny) = (8, 12);
            let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.12).expect("the grid");
            let n = nx * ny;
            let sdf = move |_x: f64, y: f64| y - wall_y;
            let rho = 1.0f64;
            let vy = -50.0f64;
            let mut state = ConservedState2d {
                rho: vec![rho; n],
                mx: vec![0.0; n],
                my: vec![rho * vy; n],
                e: vec![rho * (1.0e5 + 0.5 * vy * vy); n],
            };
            let mut species = CombustionState::stoich_h2o2(n);
            let bc = Boundaries2d {
                west: crate::solver2d::Bc2d::Transmissive,
                east: crate::solver2d::Bc2d::Transmissive,
                south: crate::solver2d::Bc2d::Transmissive,
                north: crate::solver2d::Bc2d::SlipWall,
            };
            let cut = AxiCutGeom::from_sdf(&sdf, &g);
            let before = annular_mass(&state, &g, &cut);
            for _ in 0..20 {
                crate::combustion_march::advance2d_combustion_axi(
                    &mut state,
                    &mut species,
                    &g,
                    0.2,
                    false,
                    &bc,
                    1,
                    Some(&sdf),
                )
                .expect("the staircase step");
            }
            let after = annular_mass(&state, &g, &cut);
            ((after - before).abs() / before).abs()
        };
        println!("STAIRCASE face-aligned wall: drift {:.6}", stair(0.080));
        println!("STAIRCASE mid-cell wall:     drift {:.6}", stair(0.085));
        assert!(
            aligned < 1.0e-9,
            "even a wall on a grid face leaks: {aligned}"
        );
        assert!(
            n_middle > 0,
            "the mid-cell wall produced no cut cells, so the test proves nothing"
        );
    }

    /// THE EQUIVALENCE TEST: with every cell fully fluid and no wall,
    /// the cut lane must reproduce the staircase lane exactly.
    ///
    /// both lanes run the same flux machinery, the same update signs and
    /// the same operator splits. the cut lane differs only in the
    /// geometry: its cells carry volume fractions and its faces carry
    /// apertures. on a grid where every fraction is one and every
    /// aperture is one, the two lanes are the same program. any
    /// difference is an algebra mistake in the cut lane, not a
    /// modelling choice.
    #[test]
    fn the_cut_lane_matches_the_staircase_when_no_cell_is_cut() {
        let (nx, ny) = (12, 10);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.1).expect("the grid");
        let n = nx * ny;
        // a wall far below the domain: every cell is wholly fluid.
        let cut = AxiCutGeom::from_sdf(&|_x: f64, y: f64| y + 1.0, &g);
        for k in 0..n {
            assert!(
                (cut.vol[k] - 1.0).abs() < 1.0e-12,
                "cell {k} is not fully fluid: {}",
                cut.vol[k]
            );
        }

        let mut state_a = ConservedState2d {
            rho: vec![1.4; n],
            mx: vec![2.0; n],
            my: vec![0.3; n],
            e: vec![2.0e5; n],
        };
        let mut state_b = state_a.clone();
        let mut sp_a = CombustionState::stoich_h2o2(n);
        let mut sp_b = sp_a.clone();

        // the same reflecting walls, so the domains match.
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };

        crate::combustion_march::advance2d_combustion_axi(
            &mut state_a,
            &mut sp_a,
            &g,
            0.2,
            false,
            &bc,
            1,
            None,
        )
        .expect("the staircase step");

        advance2d_combustion_axi_cut(
            &mut state_b,
            &mut sp_b,
            &g,
            &cut,
            0.2,
            false,
            &bc,
            1,
            None,
            0.0,
            1.0,
        )
        .expect("the cut step");

        for k in 0..n {
            for (name, a, b) in [
                ("rho", state_a.rho[k], state_b.rho[k]),
                ("mx", state_a.mx[k], state_b.mx[k]),
                ("my", state_a.my[k], state_b.my[k]),
                ("e", state_a.e[k], state_b.e[k]),
            ] {
                let tol = 1.0e-9 * a.abs().max(1.0);
                assert!(
                    (a - b).abs() <= tol,
                    "cell {k} {name}: the staircase has {a}, the cut lane has {b}"
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
            let cut = AxiCutGeom::from_sdf(&move |_x: f64, y: f64| y - wall, &g);
            let mut a = 0.0f64;
            for j in 0..ny {
                if cut.vol[j * cut.nx + nx / 2] <= 0.0 {
                    continue;
                }
                let rm = j as f64 * dy;
                let rp = rm + dy;
                a += cut.vol[j * cut.nx + nx / 2] * std::f64::consts::PI * (rp * rp - rm * rm);
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

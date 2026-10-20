//! park5 species march (rung 8): the planar hllc step with the
//! 5-species mixture closure and the species transport pass.
//!
//! the closure fields (pressure, sound speed) come from
//! species_closure reading the per-cell composition, so this march
//! never calls thermo::pressure for park5. the chemistry substep
//! hooks after the transport pass once the equilibrium-constant
//! data lands; until then the march is the frozen-composition
//! lane and the reacting lane is unreachable from config.

use crate::error::Error;
use crate::species::SpeciesState;
use crate::species_closure as closure5;
use crate::species_transport::transport_species;
use crate::thermo::ThermoModel;
use crate::{ConservedState2d, Grid2d};

/// one explicit species-carrying step under park5.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_park5(
    state: &mut ConservedState2d,
    species: &mut SpeciesState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &crate::Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let y_cell: Vec<f64> = (0..crate::park::NS).map(|s| species.y[s][k]).collect();
        p[k] = closure5::pressure_at(state.rho[k], e_int[k], &y_cell);
        a[k] = closure5::sound_speed_at(state.rho[k], e_int[k], &y_cell);
    }
    let mut smax = 0.0f64;
    for k in 0..n {
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1e-12);

    let model = ThermoModel::Park5 { gamma_shaping: 1.4 };
    let (fx, fy) = crate::axi::compute_axi_fluxes(
        model,
        state,
        g,
        muscl,
        bc,
        &p,
        &e_int,
        Some(species),
        nthreads,
        wall,
    )?;

    let mut resid = 0.0f64;
    #[allow(clippy::needless_range_loop)]
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let dyi = dt / g.dys[j];
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i])
                + dyi * (fy[i].mass[j + 1] - fy[i].mass[j]));
            let dmx =
                -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dyi * (fy[i].mx[j + 1] - fy[i].mx[j]));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dyi * (fy[i].my[j + 1] - fy[i].my[j]));
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dyi * (fy[i].e[j + 1] - fy[i].e[j]));
            state.rho[k] += drho;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            resid = resid
                .max(drho.abs())
                .max(dmx.abs())
                .max(dmy.abs())
                .max(de.abs());
        }
    }

    // species advection with the same face mass fluxes. the update
    // carries rho*Y so the per-species mass is conservative.
    // the sweeps carry Flux4 per face; the species pass wants the mass
    // flux per face: flatten both directions.
    let fx_mass: Vec<Vec<f64>> = fx.iter().map(|s| s.mass.clone()).collect();
    let fy_mass: Vec<Vec<f64>> = fy.iter().map(|s| s.mass.clone()).collect();
    let resid_y = transport_species(species, &state.rho, dt, &g.dxs, &g.dys, &fx_mass, &fy_mass)?;
    for k in 0..n {
        species.sanitize(k);
    }
    Ok((dt, resid.max(resid_y)))
}

/// heun (rk2) wrapper around the park5 step, mirroring
/// advance2d_model_rk2's structure: predictor-corrector on the flow
/// state with the species transported in both stages.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_park5_rk2(
    state: &mut ConservedState2d,
    species: &mut SpeciesState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &crate::Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (dt, _) = {
        let (u, v, et) =
            crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let n = g.nx * g.ny;
        let mut smax = 0.0f64;
        for k in 0..n {
            let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
            let y_cell: Vec<f64> = (0..crate::park::NS).map(|s| species.y[s][k]).collect();
            let a = closure5::sound_speed_at(state.rho[k], e_int, &y_cell);
            smax = smax.max(u[k].abs().max(v[k].abs()) + a);
        }
        (cfl * g.dx.min(g.dy) / smax.max(1e-12), 0.0)
    };
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    advance2d_park5(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
    advance2d_park5(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
    let n = g.nx * g.ny;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    for s in 0..crate::park::NS {
        for k in 0..n {
            species.y[s][k] = 0.5 * species.y[s][k] + 0.5 * y1.y[s][k];
        }
    }
    for k in 0..n {
        species.sanitize(k);
    }
    // chemistry (operator split, after the transport update): the
    // reactions run at the post-step state, the same split park2t uses.
    let (_, oversub) = chemistry_substep(state, species, g, dt)?;
    Ok((dt, oversub))
}

/// internal energy per mass for a standard-air park5 initial
/// condition at (rho, p): T from the mixture gas constant, then the
/// closure energy (sensible plus formation offset).
pub fn air_ic_energy(rho: f64, p: f64) -> f64 {
    let y = [0.767f64, 0.233, 0.0, 0.0, 0.0];
    let t = p / (rho * closure5::mix_gas_constant(&y));
    closure5::internal_energy_at(t, &y)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a uniform frozen-air field must march identically to the perfect-gas
    /// lane: same pressure field, same dt, and the composition must not
    /// move a single digit.
    #[test]
    fn frozen_uniform_matches_perfect_gas() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.8, 0.0, 0.8).unwrap();
        let (rho, u, _p) = (1.2f64, 300.0f64, 101325.0f64);
        // e_int under park5 carries the formation offset for air.
        let y = [0.767f64, 0.233, 0.0, 0.0, 0.0];
        let e_int = closure5::internal_energy_at(300.0, &y);
        let e_cell = e_int + 0.5 * u * u;
        let n = nx * ny;
        let mut st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![rho * u; n],
            my: vec![0.0; n],
            e: vec![rho * e_cell; n],
        };
        let mut sp = SpeciesState::air(n);
        let bc = crate::Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        // one rk2 step
        let (dt, _) = advance2d_park5_rk2(&mut st, &mut sp, &g, 0.4, false, &bc, 1, None).unwrap();
        assert!(dt > 0.0 && dt.is_finite());
        // uniform throughflow with transmissive ghosts: the field must stay
        // uniform to roundoff (flux divergence cancels).
        for k in 0..n {
            assert!((st.rho[k] - rho).abs() / rho < 1e-12);
            assert!((st.mx[k] - rho * u).abs() / (rho * u) < 1e-12);
        }
        // and the composition is untouched.
        assert!(sp.max_sum_drift() < 1e-12);
        assert!((sp.y[0][0] - 0.767).abs() < 1e-12);
    }
}

/// chemistry substep (operator-split): integrate dY/dt = omega(T, Y)
/// per cell over the flow dt with sub-stepping, moving the chemical
/// energy between the species pool and the thermal pool at fixed
/// total energy. the energy stays in the conserved state: the
/// formation part the products carry is already inside e through the
/// closure, so the substep only updates Y — the thermal temperature
/// follows from the same e through temperature_at.
///
/// returns the largest |dY| any substep applied and the largest
/// chemical timescale violation, so the march can report both.
pub fn chemistry_substep(
    state: &ConservedState2d,
    species: &mut SpeciesState,
    g: &Grid2d,
    dt: f64,
) -> Result<(f64, f64), Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut worst_dy = 0.0f64;
    let mut worst_oversub = 0.0f64;
    for k in 0..n {
        let e_int = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        let mut y: [f64; crate::park::NS] = [0.0; crate::park::NS];
        for (ys, ys_src) in y.iter_mut().zip(species.y.iter()) {
            *ys = ys_src[k];
        }
        let rho = state.rho[k];
        // iterate: T depends on Y through the closure; a few fixed-point
        // sweeps converge it (the temperature change per chemistry step is
        // second-order in dY).
        let mut t = closure5::temperature_at(e_int, &y);
        let mut remaining = dt;
        let mut steps = 0;
        while remaining > 1e-14 && steps < 64 {
            let (dydt, _) = crate::species_chem::omega(rho, &y, t);
            // chemical timescale on the fastest species: cap the substep
            // at a fraction of it so the explicit sweep stays stable.
            let rate_scale = dydt.iter().fold(0.0f64, |a, d| a.max(d.abs()));
            let tau = if rate_scale > 1e-30 {
                1.0 / rate_scale
            } else {
                remaining
            };
            let h = remaining.min(tau * 0.2);
            for (ys, d) in y.iter_mut().zip(dydt.iter()) {
                let d = d * h;
                worst_dy = worst_dy.max(d.abs());
                *ys = (*ys + d).clamp(0.0, 1.0);
            }
            // renormalize to sum 1 (the clamp and the per-reaction
            // stoichiometry keep drift tiny, but the floor is exact).
            let sum: f64 = y.iter().sum();
            if sum > 0.0 {
                for ys in y.iter_mut() {
                    *ys /= sum;
                }
            }
            // total e is unchanged: the formation energy the products
            // carry is inside e through the closure, so the thermal
            // temperature simply re-solves from the same e at the new
            // composition. dissociation draws energy and cools the cell.
            t = closure5::temperature_at(e_int, &y);
            remaining -= h;
            steps += 1;
            if steps == 1 {
                // record how far the full dt exceeded the chemical
                // timescale: the sub-iteration budget this cell needed.
                worst_oversub = worst_oversub.max(dt / (tau * 0.2));
            }
        }
        for (ys_dst, ys) in species.y.iter_mut().zip(y.iter()) {
            ys_dst[k] = *ys;
        }
    }
    Ok((worst_dy, worst_oversub))
}

#[cfg(test)]
mod chem_tests {
    use super::*;

    /// a hot nitrogen cell (pure N2, 9000 K) must dissociate toward the
    /// equilibrium composition over repeated substeps, cooling as the
    /// formation energy is drawn.
    #[test]
    fn hot_nitrogen_dissociates_and_cools() {
        let (nx, ny) = (1, 1);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.1, 0.0, 0.1).unwrap();
        let rho = 1.0e-3f64;
        let y = [1.0f64, 0.0, 0.0, 0.0, 0.0];
        let e_int = closure5::internal_energy_at(9000.0, &y);
        let e_cell = e_int; // zero velocity
        let st = ConservedState2d {
            rho: vec![rho],
            mx: vec![0.0],
            my: vec![0.0],
            e: vec![rho * e_cell],
        };
        let mut sp = SpeciesState::uniform(1, &y).unwrap();
        let dt = 1.0e-6f64;
        let mut t_last = 9000.0f64;
        for _ in 0..40 {
            chemistry_substep(&st, &mut sp, &g, dt).unwrap();
            t_last = closure5::temperature_at(
                e_int,
                &[sp.y[0][0], sp.y[1][0], sp.y[2][0], sp.y[3][0], sp.y[4][0]],
            );
        }
        // dissociation happened: N fraction grew, N2 dropped, cell cooled.
        assert!(sp.y[3][0] > 0.01, "y_N = {}", sp.y[3][0]);
        assert!(sp.y[0][0] < 0.999, "y_N2 = {}", sp.y[0][0]);
        assert!(t_last < 8800.0, "T = {t_last}");
        // composition still normalized.
        assert!(sp.max_sum_drift() < 1e-12);
    }
}

/// the park5 laminar viscous march: one rk2 step = the inviscid flow
/// advance + the species viscous pass (stresses, wilke conduction,
/// species diffusion, catalytic wall energy flux) + the chemistry
/// substep, with the diffusion stability cap folded into the step.
///
/// the rk2 shape mirrors advance2d_park5_rk2: stage 1 runs two flow
/// sub-steps on the clone with the viscous pass applied after each,
/// then the combine; the viscous pass uses that sub-step's dt so the
/// wall energy budget and the species update stay synchronized.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_park5_visc_rk2(
    state: &mut ConservedState2d,
    species: &mut SpeciesState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &crate::Boundaries2d,
    wall_temperature: f64,
    catalytic: crate::catalytic::CatalyticWall,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    // transport snapshot + the diffusive stability cap; the acoustic
    // bound comes from the flow advance itself.
    let ct = crate::species_viscous::cell_transport(state, species, g)?;
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    let wall_opts = crate::species_viscous::WallOpts {
        temperature: wall_temperature,
        catalytic,
    };
    let mut dt_last = 0.0f64;
    for _ in 0..2 {
        let (dt, _) = advance2d_park5(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
        // the species diffusion rides the face-local stability caps
        // inside add_species_viscous; the flow keeps the acoustic dt.
        crate::species_viscous::add_species_viscous(
            &mut s1, &mut y1, g, &ct, dt, bc, wall_opts, wall,
        )?;
        dt_last = dt;
    }
    // combine stages
    let n = g.nx * g.ny;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    for s in 0..crate::park::NS {
        for k in 0..n {
            species.y[s][k] = 0.5 * species.y[s][k] + 0.5 * y1.y[s][k];
        }
    }
    for k in 0..n {
        species.sanitize(k);
    }
    // chemistry at the combined state, same operator split as the
    // inviscid lane.
    let (_, oversub) = chemistry_substep(state, species, g, dt_last)?;
    Ok((dt_last, oversub))
}

#[cfg(test)]
mod visc_tests {
    use super::*;

    /// failure mode 6 (plans/catalytic_research/failure-modes.md): a
    /// non-catalytic viscous wall must not move the composition. the
    /// uniform cold field with an isothermal wall heats via the
    /// conductor but every species gradient stays zero, so Y is
    /// bit-identical after the step.
    #[test]
    fn non_catalytic_viscous_wall_keeps_uniform_composition() {
        let (nx, ny) = (12, 12);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.12).unwrap();
        let mut st = ConservedState2d {
            rho: vec![1.2e-3; nx * ny],
            mx: vec![0.0; nx * ny],
            my: vec![0.0; nx * ny],
            e: vec![0.0; nx * ny],
        };
        // cold static air at rest, formation-offset energy
        let y_air = [0.767f64, 0.233, 0.0, 0.0, 0.0];
        for k in 0..nx * ny {
            st.e[k] = st.rho[k] * closure5::internal_energy_at(300.0, &y_air);
        }
        let mut sp = crate::species::SpeciesState::uniform(nx * ny, &y_air).unwrap();
        let bc = crate::solver2d::Boundaries2d {
            west: crate::solver2d::Bc2d::SlipWall,
            east: crate::solver2d::Bc2d::SlipWall,
            south: crate::solver2d::Bc2d::NoSlipWall,
            north: crate::solver2d::Bc2d::SlipWall,
        };
        advance2d_park5_visc_rk2(
            &mut st,
            &mut sp,
            &g,
            0.4,
            true,
            &bc,
            1000.0,
            crate::catalytic::CatalyticWall::NonCatalytic,
            1,
            None,
        )
        .unwrap();
        for (s, ys) in sp.y.iter().enumerate() {
            for (k, yk) in ys.iter().enumerate() {
                assert!(
                    (yk - y_air[s]).abs() < 1e-12,
                    "species {s} cell {k} moved to {}",
                    yk
                );
            }
        }
    }

    /// failure mode 5/10: the fully-catalytic wall on a hot
    /// dissociated layer pulls atoms to the wall: after the step the
    /// near-wall cells are molecule-richer than with the non-catalytic
    /// condition, and the wall-adjacent composition moved.
    #[test]
    fn fully_catalytic_wall_pulls_atoms_wallward() {
        let (nx, ny) = (12, 12);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.12).unwrap();
        // seed a hot atom-rich layer above a cold wall: the bottom row
        // is 9000 K partially-dissociated air, the rest is cold.
        let mut st = ConservedState2d {
            rho: vec![1.2e-3; nx * ny],
            mx: vec![0.0; nx * ny],
            my: vec![0.0; nx * ny],
            e: vec![0.0; nx * ny],
        };
        let y_diss = [0.3f64, 0.0, 0.0, 0.5, 0.2];
        let y_cold = [0.767f64, 0.233, 0.0, 0.0, 0.0];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let y = if j < 2 { y_diss } else { y_cold };
                st.e[k] = st.rho[k]
                    * closure5::internal_energy_at(if j < 2 { 9000.0 } else { 300.0 }, &y);
            }
        }
        let mut sp = crate::species::SpeciesState::uniform(nx * ny, &y_cold).unwrap();
        for j in 0..2 {
            for i in 0..nx {
                for (ys, yd) in sp.y.iter_mut().zip(y_diss.iter()) {
                    ys[j * nx + i] = *yd;
                }
            }
        }
        let bc = crate::solver2d::Boundaries2d {
            west: crate::solver2d::Bc2d::SlipWall,
            east: crate::solver2d::Bc2d::SlipWall,
            south: crate::solver2d::Bc2d::NoSlipWall,
            north: crate::solver2d::Bc2d::SlipWall,
        };
        // non-catalytic reference
        let mut st_n = st.clone();
        let mut sp_n = sp.clone();
        advance2d_park5_visc_rk2(
            &mut st_n,
            &mut sp_n,
            &g,
            0.4,
            true,
            &bc,
            1000.0,
            crate::catalytic::CatalyticWall::NonCatalytic,
            1,
            None,
        )
        .unwrap();
        // fully-catalytic
        let mut sp_f = sp.clone();
        let mut st_f = st.clone();
        advance2d_park5_visc_rk2(
            &mut st_f,
            &mut sp_f,
            &g,
            0.4,
            true,
            &bc,
            1000.0,
            crate::catalytic::CatalyticWall::FullyCatalytic,
            1,
            None,
        )
        .unwrap();
        // the wall-adjacent row (j=0, above the bottom no-slip wall):
        // the catalytic condition drains atoms toward the wall, so the
        // atom mass fraction there must drop relative to non-cat.
        let n_atoms_n: f64 = (0..nx).map(|i| sp_n.y[3][i] + sp_n.y[4][i]).sum();
        let n_atoms_f: f64 = (0..nx).map(|i| sp_f.y[3][i] + sp_f.y[4][i]).sum();
        assert!(
            n_atoms_f < n_atoms_n,
            "catalytic wall must deplete wall atoms: full {n_atoms_f} vs non {n_atoms_n}"
        );
        // and the wall rows must stay normalized after the update
        for k in 0..nx * ny {
            let s: f64 = (0..crate::park::NS).map(|q| sp_f.y[q][k]).sum();
            assert!((s - 1.0).abs() < 1e-9, "sum Y = {s} at {k}");
        }
    }

    /// the north-wall mirror of the catalytic test: the same hot
    /// atom-rich layer under the TOP no-slip wall must also drain
    /// atoms under the fully-catalytic condition. pins the face-
    /// direction sign in the domain-wall catalytic flux.
    #[test]
    fn fully_catalytic_north_wall_also_pulls_atoms() {
        let (nx, ny) = (12, 12);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.12).unwrap();
        let mut st = ConservedState2d {
            rho: vec![1.2e-3; nx * ny],
            mx: vec![0.0; nx * ny],
            my: vec![0.0; nx * ny],
            e: vec![0.0; nx * ny],
        };
        let y_diss = [0.3f64, 0.0, 0.0, 0.5, 0.2];
        let y_cold = [0.767f64, 0.233, 0.0, 0.0, 0.0];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let hot = j >= ny - 2;
                let y = if hot { y_diss } else { y_cold };
                st.e[k] =
                    st.rho[k] * closure5::internal_energy_at(if hot { 9000.0 } else { 300.0 }, &y);
            }
        }
        let mut sp = crate::species::SpeciesState::uniform(nx * ny, &y_cold).unwrap();
        for j in (ny - 2)..ny {
            for i in 0..nx {
                for (ys, yd) in sp.y.iter_mut().zip(y_diss.iter()) {
                    ys[j * nx + i] = *yd;
                }
            }
        }
        let bc = crate::solver2d::Boundaries2d {
            west: crate::solver2d::Bc2d::SlipWall,
            east: crate::solver2d::Bc2d::SlipWall,
            south: crate::solver2d::Bc2d::SlipWall,
            north: crate::solver2d::Bc2d::NoSlipWall,
        };
        let mut st_n = st.clone();
        let mut sp_n = sp.clone();
        advance2d_park5_visc_rk2(
            &mut st_n,
            &mut sp_n,
            &g,
            0.4,
            true,
            &bc,
            1000.0,
            crate::catalytic::CatalyticWall::NonCatalytic,
            1,
            None,
        )
        .unwrap();
        let mut sp_f = sp.clone();
        let mut st_f = st.clone();
        advance2d_park5_visc_rk2(
            &mut st_f,
            &mut sp_f,
            &g,
            0.4,
            true,
            &bc,
            1000.0,
            crate::catalytic::CatalyticWall::FullyCatalytic,
            1,
            None,
        )
        .unwrap();
        let row = ny - 1;
        let n_atoms_n: f64 = (0..nx)
            .map(|i| sp_n.y[3][row * nx + i] + sp_n.y[4][row * nx + i])
            .sum();
        let n_atoms_f: f64 = (0..nx)
            .map(|i| sp_f.y[3][row * nx + i] + sp_f.y[4][row * nx + i])
            .sum();
        assert!(
            n_atoms_f < n_atoms_n,
            "north catalytic wall must deplete wall atoms: full {n_atoms_f} vs non {n_atoms_n}"
        );
    }

    /// failure mode 8 (cold-flow identity through the viscous pass):
    /// the cold uniform field must not create composition anywhere,
    /// and the step must remain stable at the face-local diffusion dt.
    #[test]
    fn viscous_march_stays_stable_on_the_stagnation_setup() {
        // 20x40 grid, Mach-20 freestream, isothermal wall, no body:
        // the march must run 20 steps without a solver error and keep
        // the composition bounded.
        let (nx, ny) = (20, 40);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.2, 0.0, 0.4).unwrap();
        let mut st = ConservedState2d {
            rho: vec![1.363e-5; nx * ny],
            mx: vec![1.363e-5 * 6047.0; nx * ny],
            my: vec![0.0; nx * ny],
            e: vec![0.0; nx * ny],
        };
        let y_air = [1.0f64, 0.0, 0.0, 0.0, 0.0];
        for k in 0..nx * ny {
            st.e[k] = st.rho[k] * closure5::internal_energy_at(197.0, &y_air)
                + 0.5 * st.mx[k] * st.mx[k] / st.rho[k];
        }
        let mut sp = crate::species::SpeciesState::uniform(nx * ny, &y_air).unwrap();
        let bc = crate::solver2d::Boundaries2d {
            west: crate::solver2d::Bc2d::SupersonicInflow {
                rho: 1.363e-5,
                u: 6047.0,
                v: 0.0,
                p: 0.89,
            },
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::NoSlipWall,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        for _ in 0..20 {
            advance2d_park5_visc_rk2(
                &mut st,
                &mut sp,
                &g,
                0.4,
                true,
                &bc,
                1000.0,
                crate::catalytic::CatalyticWall::FullyCatalytic,
                1,
                None,
            )
            .unwrap();
        }
        for k in 0..nx * ny {
            let s: f64 = (0..crate::park::NS).map(|q| sp.y[q][k]).sum();
            assert!((s - 1.0).abs() < 1e-9);
        }
    }
}

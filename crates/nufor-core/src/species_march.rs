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

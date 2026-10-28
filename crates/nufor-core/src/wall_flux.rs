//! the shared wall heat-flux estimator: one implementation of the wall
//! conduction and catalytic heating the viscous passes apply, so every
//! reported flux comes from the march's own stencil and transport.

use crate::catalytic::CatalyticWall;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::species::SpeciesState;
use crate::state2d::ConservedState2d;
use crate::thermo::ThermoModel;
use crate::viscous2d::ViscCoeffs;

/// one wetted wall face: where it sits and how much heat crosses it.
#[derive(Clone, Copy, Debug)]
pub struct WallFluxFace {
    /// the surface position along x, metres.
    pub x: f64,
    /// the surface position along y, metres.
    pub y: f64,
    /// the magnitude of the wall energy-flux entry the viscous pass applies here, W/m^2.
    pub q_w_m2: f64,
}

/// the wall heat flux over a marched field.
#[derive(Clone, Debug, Default)]
pub struct WallFluxReport {
    /// the wetted faces, ordered by x.
    pub faces: Vec<WallFluxFace>,
    /// the largest face flux, W/m^2.
    pub peak_w_m2: f64,
    /// the peak's x position, metres.
    pub peak_x: f64,
    /// the peak's y position, metres.
    pub peak_y: f64,
}

/// the transport and boundary numbers a wall-flux estimate needs: the
/// same values the viscous pass hands its conductor.
#[derive(Clone, Copy, Debug)]
pub struct WallFluxSetup {
    /// the run's thermodynamic closure.
    pub model: ThermoModel,
    /// the reference viscosity, kg/(m s).
    pub mu0: f64,
    /// sutherland's reference temperature, K.
    pub t0_ref: f64,
    /// sutherland's constant, K.
    pub s_param: f64,
    /// the prandtl number.
    pub pr: f64,
    /// the wall temperature, K (0 leaves the wall adiabatic).
    pub wall_temperature: f64,
    /// the species lane's wall treatment.
    pub catalytic: CatalyticWall,
}

/// the gamma the viscous marches use for a closure.
fn lane_gamma(model: ThermoModel) -> f64 {
    match model {
        ThermoModel::Perfect { gamma } => gamma,
        ThermoModel::EqAir { .. } => 1.4,
        ThermoModel::Park5 { gamma_shaping } | ThermoModel::Combustion9 { gamma_shaping } => {
            gamma_shaping
        }
    }
}

/// the wall heat flux over a marched field, computed with the march's own
/// wall stencil, temperature field and transport coefficients.
///
/// the species lane uses its wilke transport and catalytic wall term; the
/// single-species lanes use the closure temperature and sutherland's law.
/// a face counts as wetted when exactly one of its two cells is solid,
/// the same classification the viscous sweeps use.
pub fn wall_flux(
    state: &ConservedState2d,
    g: &Grid2d,
    setup: &WallFluxSetup,
    species: Option<&SpeciesState>,
    solid: &[bool],
) -> Result<WallFluxReport, Error> {
    if solid.is_empty() {
        return Ok(WallFluxReport::default());
    }
    let n = g.nx * g.ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut faces: Vec<WallFluxFace> = Vec::new();
    match species {
        Some(sp) => {
            let ct = crate::species_viscous::cell_transport(state, sp, g)?;
            let opts = crate::species_viscous::WallOpts {
                temperature: setup.wall_temperature,
                catalytic: setup.catalytic,
            };
            for j in 0..g.ny {
                for f in 1..g.nx {
                    let left = solid[j * g.nx + f - 1];
                    let right = solid[j * g.nx + f];
                    if left == right {
                        continue;
                    }
                    let (_txx, _txy, q_cond, q_cat, _j_wall) = crate::species_viscous::wall_face_x(
                        left, f, j, &u, &v, &ct, sp, &state.rho, g.nx, g, opts,
                    );
                    let k = if left { j * g.nx + f } else { j * g.nx + f - 1 };
                    faces.push(WallFluxFace {
                        x: g.centers_x[k],
                        y: g.centers_y[k],
                        q_w_m2: q_cond.abs() + q_cat.abs(),
                    });
                }
            }
            for i in 0..g.nx {
                for f in 1..g.ny {
                    let below = solid[(f - 1) * g.nx + i];
                    let above = solid[f * g.nx + i];
                    if below == above {
                        continue;
                    }
                    let (_txy, _tyy, q_cond, q_cat, _j_wall) = crate::species_viscous::wall_face_y(
                        below, i, f, &u, &v, &ct, sp, &state.rho, g.nx, g, opts,
                    );
                    let k = if below {
                        f * g.nx + i
                    } else {
                        (f - 1) * g.nx + i
                    };
                    faces.push(WallFluxFace {
                        x: g.centers_x[k],
                        y: g.centers_y[k],
                        q_w_m2: q_cond.abs() + q_cat.abs(),
                    });
                }
            }
        }
        None => {
            // the single-species lanes: the closure temperature the
            // viscous pass conducts with, and sutherland's law on it.
            let mut e_int = vec![0.0; n];
            for k in 0..n {
                e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
            }
            let t = match setup.model {
                ThermoModel::EqAir { .. } => {
                    crate::eqair_cea::fused_state_mu(&state.rho, &e_int, |_| 0.0).1
                }
                ThermoModel::Perfect { .. } => {
                    let p = crate::thermo::pressure(setup.model, &state.rho, &e_int, &u, &v)?;
                    crate::thermo::temperature(setup.model, &state.rho, &p)?
                }
                ThermoModel::Park5 { .. } | ThermoModel::Combustion9 { .. } => {
                    // the species lanes need their composition: without it
                    // the transport is undefined, so refuse rather than
                    // report a number from the wrong closure.
                    return Err(Error::InvalidArgs);
                }
            };
            let mu_cells: Vec<f64> = t
                .iter()
                .map(|&tk| {
                    crate::viscous2d::sutherland_mu(
                        setup.mu0,
                        tk.max(1.0),
                        setup.t0_ref,
                        setup.s_param,
                    )
                })
                .collect();
            let mu_cell = |k: usize| mu_cells[k];
            let gamma = lane_gamma(setup.model);
            let coeffs = ViscCoeffs {
                kappa: crate::viscous2d::conduction_coeff(Some(setup.model), gamma, setup.pr),
                kappa_t: 0.0,
                gamma,
            };
            for j in 0..g.ny {
                for f in 1..g.nx {
                    let left = solid[j * g.nx + f - 1];
                    let right = solid[j * g.nx + f];
                    if left == right {
                        continue;
                    }
                    let (_txx, _txy, q) = crate::viscous2d::body_wall_x(
                        left,
                        f,
                        j,
                        &u,
                        &v,
                        &t,
                        g.nx,
                        g,
                        &mu_cell,
                        &[],
                        &coeffs,
                        setup.wall_temperature,
                    );
                    let k = if left { j * g.nx + f } else { j * g.nx + f - 1 };
                    faces.push(WallFluxFace {
                        x: g.centers_x[k],
                        y: g.centers_y[k],
                        q_w_m2: q.abs(),
                    });
                }
            }
            for i in 0..g.nx {
                for f in 1..g.ny {
                    let below = solid[(f - 1) * g.nx + i];
                    let above = solid[f * g.nx + i];
                    if below == above {
                        continue;
                    }
                    let (_txy, _tyy, q) = crate::viscous2d::body_wall_y(
                        below,
                        i,
                        f,
                        &u,
                        &v,
                        &t,
                        g.nx,
                        g,
                        &mu_cell,
                        &[],
                        &coeffs,
                        setup.wall_temperature,
                    );
                    let k = if below {
                        f * g.nx + i
                    } else {
                        (f - 1) * g.nx + i
                    };
                    faces.push(WallFluxFace {
                        x: g.centers_x[k],
                        y: g.centers_y[k],
                        q_w_m2: q.abs(),
                    });
                }
            }
        }
    }
    faces.sort_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.q_w_m2.total_cmp(&b.q_w_m2))
    });
    let mut report = WallFluxReport {
        faces,
        ..Default::default()
    };
    if let Some(pk) = report
        .faces
        .iter()
        .max_by(|a, b| a.q_w_m2.total_cmp(&b.q_w_m2))
    {
        report.peak_w_m2 = pk.q_w_m2;
        report.peak_x = pk.x;
        report.peak_y = pk.y;
    }
    Ok(report)
}

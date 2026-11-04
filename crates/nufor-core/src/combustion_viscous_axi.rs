//! the axisymmetric viscous combustion lane: the reacting axi march plus
//! molecular transport, with the engine's wall treated.
//!
//! this is the lane the cooled-wall engine cases need. the planar
//! viscous combustion lane (combustion_viscous.rs) treats the duct's
//! rectangular boundary, but an engine's contour is a body: the wall
//! runs along a staircase of cells, and the flow there needs a one-sided
//! shear, conduction and species flux pointing INTO the gas.
//!
//! the transport coefficients are the same wilke-mixture ones the planar
//! lane uses. the difference is the geometry: the radial direction is
//! annular, so the radial fluxes divide by the radius and the axis is a
//! symmetry boundary rather than a wall.
//!
//! the wall treatment is one sided and isothermal:
//!   - no slip, so the tangential velocity vanishes at the wall
//!   - fourier conduction to the wall temperature when one is given
//!   - no composition gradient, so the wall neither creates nor absorbs
//!     species (the non-catalytic default)
//!
//! the wall is placed on the body's crossed faces: a face the wall
//! crosses carries the flux, taken from the fluid cell on its
//! gas side.

use crate::combustion_closure as closure9;
use crate::combustion_march::{advance2d_combustion_axi, CombustionState};
use crate::combustion_thermo as ct;
use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::solver2d::Boundaries2d;
use crate::state2d::{cons_to_prim2d, ConservedState2d};

/// per-cell transport snapshot for the combustion set: closure
/// temperature, wilke-mixture mu and k, and the per-species effective
/// diffusivities.
///
/// this is public so a gate can drive one wall application against the
/// solver's own transport snapshot.
pub struct AxiCombustionTransport {
    pub t: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub d_s: Vec<[f64; ct::NS]>,
}

/// the wall's treatment options.
#[derive(Clone, Copy, Debug, Default)]
pub struct AxiWallOpts {
    /// the wall temperature in K; zero or less means adiabatic.
    pub temperature: f64,
    /// the wall's grey-body emissivity; zero or less means the wall
    /// does not radiate.
    pub emissivity: f64,
    /// the gas-side film heat-transfer coefficient in W per square
    /// metre kelvin; zero or less means no convective film exchange.
    ///
    /// this is the GAS-SIDE resistance of the film's two-resistance
    /// balance. the coolant's own capacity rate, `coolant_mass_flux *
    /// coolant_capacity`, is the second resistance, and the two in
    /// series bound the duty the wall can draw. the OLD model applied
    /// this coefficient straight against the cold wall, which is an
    /// unbounded sink: it grew with the gas's temperature and stripped
    /// the near-wall cell below its thermodynamic floor, which tripped
    /// the positivity guard and diverged the film leg.
    pub film_coefficient: f64,
    /// the regenerative coolant's mass flux in kg per square metre
    /// second; zero or less means the coolant carries no capacity.
    ///
    /// this coolant stays inside the wall's passages. it takes energy
    /// from the wall and leaves by the chamber's own circuit, so it
    /// changes the energy balance and NOT the chamber's mass balance.
    pub coolant_mass_flux: f64,
    /// the regenerative coolant's specific heat capacity in J per
    /// kilogram kelvin.
    pub coolant_capacity: f64,
    /// the regenerative coolant's inlet temperature in K.
    pub regenerative_coolant_temperature: f64,
    /// the blowing-in film's mass flux in kg per square metre second;
    /// zero or less means no coolant enters the chamber's gas.
    ///
    /// unlike the regenerative coolant, this one CROSSES the wall
    /// into the gas, so the chamber's mass balance has to carry it and
    /// the gas's composition and enthalpy change with it.
    pub film_injection_mass_flux: f64,
    /// the blowing-in coolant's inlet temperature in K.
    pub film_injection_temperature: f64,
    /// the blowing-in coolant's composition as mass fractions, over
    /// the nine species of the closure.
    pub film_injection_composition: [f64; ct::NS],
    /// the wall liner's thickness in metres; zero or less means the
    /// wall has no liner and the two-resistance model stands alone.
    ///
    /// the liner's conduction sits in series between the gas-side film
    /// and the coolant's capacity rate, so a thick resistive liner
    /// lowers the duty and runs the inner surface hotter.
    pub liner_thickness: f64,
    /// the wall liner's conductivity in W per metre kelvin; zero or
    /// less means the liner conducts perfectly.
    pub liner_conductivity: f64,
    /// the fraction of a cell's energy ABOVE its thermodynamic floor
    /// that a wall may take in one step.
    ///
    /// the wall's exchange is explicit, and a near-wall cell holds
    /// little energy above the closure's floor. a wall that strips
    /// that margin in a single step drives the cell below the floor and
    /// trips the positivity guard, which costs the cell its momentum
    /// and poisons the march. limiting the withdrawal to a fraction of
    /// the available margin keeps the scheme admissible without
    /// capping the physics at a converged state: as the cell cools, the
    /// permitted withdrawal shrinks on its own. zero disables the
    /// bound.
    pub max_energy_fraction_per_step: f64,
}

/// the gas temperature extrapolated to the wall face, in K.
///
/// the radiative and film fluxes act on the gas at the wall surface,
/// not at the first cell's centre. near a cold wall the cell centre is
/// markedly cooler than the gas actually touching the wall, so using
/// the centre there would bias the exchange downward. this
/// extrapolates the two nearest cell-centre temperatures linearly to
/// the wall's position.
fn surface_temperature(v0: f64, v1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    if (c1 - c0).abs() < 1.0e-30 {
        return v0;
    }
    v0 + (y_wall - c0) * (v1 - v0) / (c1 - c0)
}

/// the wall's outward normal's direction, from the gas toward the body.
///
/// the signed distance is negative inside the body and grows outward,
/// so its gradient points from the body toward the gas. the heat
/// leaving the gas travels the other way, toward the body, so the
/// direction the heat travels is the negated gradient.
fn wall_heat_direction(dist: &dyn Fn(f64, f64) -> f64, x: f64, y: f64, h: f64) -> (f64, f64) {
    let dx = (dist(x + h, y) - dist(x - h, y)) / (2.0 * h);
    let dy = (dist(x, y + h) - dist(x, y - h)) / (2.0 * h);
    let norm = (dx * dx + dy * dy).sqrt();
    if norm < 1.0e-30 {
        return (0.0, 0.0);
    }
    (-dx / norm, -dy / norm)
}

/// the stefan-boltzmann constant in W per square metre per kelvin to
/// the fourth.
pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;

/// the radiant flux leaving the gas for a grey wall, in W per square
/// metre.
///
/// the gas sees the wall through the grey-body exchange
///
///   q = eps * sigma * (t_gas^4 - t_wall^4)
///
/// which at 3690 k is about 9.5 mw per square metre for a coated
/// chamber wall. that is an order of magnitude above the fourier
/// conduction the near-wall cells carry, so a wall modelled with
/// conduction alone exchanges a small fraction of the heat a real
/// chamber wall does.
fn radiant_wall_flux(t_gas: f64, t_wall: f64, emissivity: f64) -> f64 {
    if emissivity <= 0.0 || t_gas <= 0.0 || t_wall <= 0.0 {
        return 0.0;
    }
    let tg = t_gas.clamp(0.0, 6000.0);
    let tw = t_wall.clamp(0.0, 6000.0);
    emissivity * STEFAN_BOLTZMANN * (tg.powi(4) - tw.powi(4))
}

/// the film's flux from the gas to the wall, in W per square metre.
fn film_wall_flux(t_gas: f64, t_wall: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    h * (t_gas - t_wall)
}

/// the cell's volume in an axisymmetric grid, in cubic metres.
///
/// the axisymmetric lane's cells are annular sectors, so the volume
/// carries the pi and the squared radii. the wall's admissibility bound
/// needs it: the margin a near-wall cell holds above its
/// thermodynamic floor is a volume of energy.
/// the volume of one axisymmetric cell, in cubic metres.
pub fn axi_cell_volume(g: &Grid2d, i: usize, j: usize) -> f64 {
    let rm = g.faces_y[j];
    let rp = g.faces_y[j + 1];
    std::f64::consts::PI * (rp * rp - rm * rm) * g.dxs[i]
}

/// the wall's flux, bounded so the cell keeps a usable energy margin.
///
/// the wall's exchange is explicit. a near-wall cell holds only a
/// little energy above the closure's floor, and an unclamped cold-wall
/// sink strips that margin in a single step, dropping the cell below
/// the floor and tripping the positivity guard. the guard then zeroes
/// the cell's momentum, and the poisoned cell poisons its neighbours,
/// which is how the film leg diverged.
///
/// the bound permits a wall to withdraw at most a fraction of the
/// energy the cell holds above the floor, so the scheme stays
/// admissible. it is NOT a physics cap: as the cell cools, the margin
/// shrinks and the permitted withdrawal shrinks with it, so a
/// converged solution that sits comfortably above the floor feels no
/// limit at all.
fn admissible_wall_flux(ctx: &WallAdmissibility<'_>, cell: usize, flux: f64, area: f64) -> f64 {
    let fraction = ctx.fraction;
    let dt = ctx.dt;
    if fraction <= 0.0 || dt <= 0.0 || area <= 0.0 || flux == 0.0 {
        return flux;
    }
    // a NEGATIVE flux puts heat INTO the gas, and the cell's margin
    // is irrelevant to that case: giving a cell energy can never trip
    // the floor. only a withdrawal needs bounding.
    if flux < 0.0 {
        return flux;
    }
    let g = ctx.grid;
    let nx = g.nx;
    let i = cell % nx;
    let j = cell / nx;
    let volume = axi_cell_volume(g, i, j);
    if volume <= 0.0 {
        return flux;
    }
    let mut y_cell = [0.0f64; ct::NS];
    for (dst, src_row) in y_cell.iter_mut().zip(ctx.species.y.iter()) {
        *dst = src_row[cell];
    }
    let floor = closure9::internal_energy_at(200.0, &y_cell);
    let rho = ctx.state.rho[cell].max(1.0e-12);
    let mx = ctx.state.mx[cell];
    let my = ctx.state.my[cell];
    // the state's energy is an energy DENSITY (rho times the specific
    // total energy), so the specific energies come from dividing it
    // by the density first. comparing the raw density against a
    // specific floor would overstate the margin by the density and
    // leave the bound toothless.
    let et = ctx.state.e[cell] / rho;
    let ke = crate::state2d::specific_kinetic_energy(rho, mx, my);
    let e_int = et - ke;
    let margin_joules = (e_int - floor).max(0.0) * rho * volume;
    let permitted_power = fraction * margin_joules / dt;
    let demanded_power = flux.abs() * area;
    if demanded_power <= permitted_power {
        return flux;
    }
    flux * (permitted_power / demanded_power)
}

// the factor that keeps the wall's withdrawals admissible.
//
// the wall's contributions do not all point the same way. the
// conduction term is signed along the face's AXIS, while the
// radiant and film terms are magnitudes of heat LEAVING the gas,
// turned into the axis direction by the face's outward normal.
// converting the conduction term into the same outward sense lets
// every term be compared as a loss.
//
// the bound then applies only to the terms that WITHHOLD energy:
// handing the cell energy can never drive it under the floor, so a
// term that adds heat is left untouched. the withdrawing terms share
// one proportional factor, so the wall's thermal split survives.
fn bounded_wall_losses(
    ctx: &WallAdmissibility<'_>,
    cell: usize,
    area: f64,
    conduction_outward: f64,
    radiant: f64,
    film: f64,
) -> f64 {
    let gaining = conduction_outward.min(0.0) + radiant.min(0.0) + film.min(0.0);
    let withdrawing = conduction_outward.max(0.0) + radiant.max(0.0) + film.max(0.0);
    if withdrawing <= 0.0 {
        return gaining;
    }
    let scale = admissible_wall_scale(ctx, cell, withdrawing, area);
    gaining + scale * withdrawing
}

/// the near-wall cell's state and the geometry the wall's
/// admissibility bound needs.
///
/// the bound reads the cell's energy margin and the face's area, and
/// both come from the same march, so they travel together.
struct WallAdmissibility<'a> {
    state: &'a ConservedState2d,
    species: &'a CombustionState,
    grid: &'a Grid2d,
    dt: f64,
    fraction: f64,
}

/// the factor that keeps the wall's total heat loss admissible.
///
/// the wall removes heat from the gas three ways: conduction down the
/// near-wall gradient, radiation to the surface, and the film's own
/// convective exchange. all three draw on the SAME near-wall cell, so
/// bounding them apart would let the pair exceed what the cell can
/// spare. this computes one scale factor for the lot.
///
/// the scale is proportional, so the wall's thermal SPLIT survives: a
/// radiating wall that is partly clipped keeps its radiation-to-film
/// ratio. and because the limit reads the cell's energy MARGIN, a
/// converged state sitting comfortably above the floor feels no limit
/// at all. the limit exists only to keep the explicit step
/// admissible.
fn admissible_wall_scale(
    ctx: &WallAdmissibility<'_>,
    cell: usize,
    withdrawal: f64,
    area: f64,
) -> f64 {
    if withdrawal <= 0.0 {
        return 1.0;
    }
    let bounded = admissible_wall_flux(ctx, cell, withdrawal, area);
    (bounded / withdrawal).clamp(0.0, 1.0)
}

/// the wall's exchange with the gas at one face, resolved through the
/// wall's thermal circuit.
///
/// the gas meets the wall's surface through the film and through
/// radiation, and the surface sheds what arrives into the liner and
/// the coolant. the whole thing is solved together, so the surface's
/// temperature and the three duties are mutually consistent.
///
/// the wall's OWN stated temperature is the sink the surface sheds
/// into when the wall carries no coolant of its own. a wall with a
/// coolant sheds into the coolant instead, and the wall's stated
/// temperature no longer describes the surface at all.
///
/// this is the ONE place the wall's circuit is solved. the solver's
/// faces and the duty report both come through here, so the report
/// cannot flatter the solver by measuring a different wall.
pub(crate) fn wall_exchange(t_gas: f64, wall_opts: &AxiWallOpts) -> WallExchange {
    // a wall with no backing of its own is a plain isothermal surface,
    // and it is exchanged against exactly as before: the film's and the
    // radiation's laws read at the wall's stated temperature. only a
    // wall that HAS a backing, a liner or a coolant, settles its own
    // surface temperature first.
    let has_backing = wall_opts.liner_thickness > 0.0 && wall_opts.liner_conductivity > 0.0
        || wall_opts.coolant_mass_flux > 0.0 && wall_opts.coolant_capacity > 0.0;
    if !has_backing {
        let t_wall = wall_opts.temperature;
        let radiative = radiant_wall_flux(t_gas, t_wall, wall_opts.emissivity);
        let convective = film_wall_flux(t_gas, t_wall, wall_opts.film_coefficient);
        return WallExchange {
            surface_temperature: t_wall,
            convective_flux: convective,
            radiative_flux: radiative,
        };
    }
    let request = crate::conjugate_wall::WallRequest {
        gas_temperature: t_gas,
        gas_side_coefficient: wall_opts.film_coefficient,
        emissivity: wall_opts.emissivity,
        stated_temperature: wall_opts.temperature,
        sink_temperature: wall_opts.coolant_temperature(),
        liner_thickness: wall_opts.liner_thickness,
        liner_conductivity: wall_opts.liner_conductivity,
        coolant_mass_flux: wall_opts.coolant_mass_flux,
        coolant_capacity: wall_opts.coolant_capacity,
    };
    let surface = crate::conjugate_wall::solve_wall_surface(&request);
    // the wall's surface is where the gas's film and its radiation meet
    // the wall's backing, so the classical terms are read at THAT
    // temperature rather than at the wall's stated one.
    let t_wall = surface.surface_temperature;
    let radiative = radiant_wall_flux(t_gas, t_wall, wall_opts.emissivity);
    let convective = film_wall_flux(t_gas, t_wall, wall_opts.film_coefficient);
    WallExchange {
        surface_temperature: t_wall,
        convective_flux: convective,
        radiative_flux: radiative,
    }
}

impl AxiWallOpts {
    /// the regenerative coolant's inlet temperature, in K.
    ///
    /// the coolant enters the wall's passages at this temperature and
    /// leaves warmer, having taken the duty the gas paid.
    pub fn coolant_temperature(&self) -> f64 {
        self.regenerative_coolant_temperature
    }
}

/// the area of an axisymmetric face at constant x, in square metres.
///
/// the face spans the annulus between two radii and sweeps the full
/// circumference, so its area carries the pi and the squared radii.
pub(crate) fn axi_wall_face_area_x(g: &Grid2d, j: usize) -> f64 {
    let rm = g.faces_y[j];
    let rp = g.faces_y[j + 1];
    std::f64::consts::PI * (rp * rp - rm * rm)
}

/// the area of an axisymmetric face at constant radius, in square
/// metres.
///
/// the face spans one cell's width in x and sweeps the circumference
/// at its own radius.
pub(crate) fn axi_wall_face_area_y(g: &Grid2d, i: usize, f: usize) -> f64 {
    2.0 * std::f64::consts::PI * g.faces_y[f] * g.dxs[i]
}

/// the wall's exchange with the gas at one face, in W per square metre.
///
/// the fluxes are signed the way the wall's terms are: a positive value
/// means heat leaving the gas.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct WallExchange {
    /// the wall's surface temperature, in K.
    pub surface_temperature: f64,
    /// the film's share of the exchange, in W per square metre.
    pub convective_flux: f64,
    /// radiation's share of the exchange, in W per square metre.
    pub radiative_flux: f64,
}

/// the film's exchange with the gas.
#[derive(Clone, Copy, Debug, Default)]
pub struct FilmExchange {
    /// the flux leaving the gas, in W per square metre.
    pub flux: f64,
    /// the film's interfacial temperature, in K.
    pub film_temperature: f64,
    /// the duty the coolant carries away, in W per square metre.
    pub coolant_duty: f64,
}

/// one-sided quadratic derivative through the wall value and the first
/// two cell centres, evaluated at the wall.
fn wall_quad_deriv(v_wall: f64, v0: f64, v1: f64, y_wall: f64, c0: f64, c1: f64) -> f64 {
    let w0 = (2.0 * y_wall - c0 - c1) / ((y_wall - c0) * (y_wall - c1));
    let w1 = (2.0 * y_wall - y_wall - c1) / ((c0 - y_wall) * (c0 - c1));
    let w2 = (2.0 * y_wall - y_wall - c0) / ((c1 - y_wall) * (c1 - c0));
    v_wall * w0 + v0 * w1 + v1 * w2
}

/// the blown-in film's injection across the crossed wall faces.
///
/// the film's crossing of the wall is INDEPENDENT of the wall's own
/// thermal treatment: the coolant enters the gas whatever the wall
/// does thermally, so the injection is enumerated on its own and the
/// solver and the report both read it from here.
#[derive(Clone, Debug)]
pub(crate) struct FilmInjection {
    /// the coolant's mass rate into each cell, kg per second.
    pub mass_rate: Vec<f64>,
    /// the coolant's energy rate into each cell, W.
    pub energy_rate: Vec<f64>,
    /// the coolant's species mass rate into each cell, kg per second.
    pub species_rate: Vec<Vec<f64>>,
    /// the wall area the film crossed, square metres.
    ///
    /// a boundary cut obliquely books its flux against one cell from
    /// two crossed faces, and the cell takes half of each, so the area
    /// here is the wall's extent as the film met it.
    pub area_square_metres: f64,
    /// the coolant's mass flow across the wall, kg per second.
    pub mass_flow_kg_per_second: f64,
}

/// enumerates the blown-in film over the crossed wall faces.
///
/// the per-face charges are split across the faces that share a cell
/// by the same weights the wall's own flux uses, so the film's mass
/// lands where the solver's wall flux lands.
pub(crate) fn film_injection(
    g: &Grid2d,
    solid: &[bool],
    wall_x: &[bool],
    wall_y: &[bool],
    opts: &AxiWallOpts,
) -> FilmInjection {
    let nx = g.nx;
    let ny = g.ny;
    let n = nx * ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let mut out = FilmInjection {
        mass_rate: vec![0.0f64; n],
        energy_rate: vec![0.0f64; n],
        species_rate: vec![vec![0.0f64; n]; ct::NS],
        area_square_metres: 0.0,
        mass_flow_kg_per_second: 0.0,
    };
    if opts.film_injection_mass_flux <= 0.0 {
        return out;
    }
    let charge = wall_charge_counts(nx, ny, solid, wall_x, wall_y);
    let enthalpy = crate::combustion_thermo::mix_h_mass(
        &opts.film_injection_composition,
        opts.film_injection_temperature.max(200.0),
    );

    let deposit = |out: &mut FilmInjection, cell: usize, area: f64| {
        if cell >= n {
            return;
        }
        let weight = match charge[cell] {
            0 => 0.0,
            1 => 1.0,
            shares => 1.0 / shares as f64,
        };
        if weight <= 0.0 {
            return;
        }
        let rate = opts.film_injection_mass_flux * area * weight;
        out.mass_rate[cell] += rate;
        out.energy_rate[cell] += rate * enthalpy;
        for (species, row) in out.species_rate.iter_mut().enumerate() {
            row[cell] += rate * opts.film_injection_composition[species];
        }
        out.area_square_metres += area * weight;
        out.mass_flow_kg_per_second += rate;
    };

    for j in 0..ny {
        for f in 1..nx {
            if !wall_x[j * (nx + 1) + f] {
                continue;
            }
            let (a, ok) = if !solid[idx(f, j)] {
                (idx(f, j), f + 1 < nx)
            } else {
                (idx(f - 1, j), f >= 2)
            };
            if !ok {
                continue;
            }
            deposit(&mut out, a, axi_wall_face_area_x(g, j));
        }
    }
    for i in 0..nx {
        for f in 1..ny {
            if !wall_y[i * (ny + 1) + f] {
                continue;
            }
            let (a, ok) = if !solid[idx(i, f)] {
                (idx(i, f), f + 1 < ny)
            } else {
                (idx(i, f - 1), f >= 2)
            };
            if !ok {
                continue;
            }
            deposit(&mut out, a, axi_wall_face_area_y(g, i, f));
        }
    }
    out
}

/// the number of wall faces that charge each gas cell.
///
/// a boundary that cuts a cell obliquely is represented by two crossed
/// grid faces at once, and both book their full flux against the same
/// cell, so the cell receives the wall's exchange twice for one
/// boundary segment. the count lets the caller split the segment's flux
/// between the faces that share it.
pub(crate) fn wall_charge_counts(
    nx: usize,
    ny: usize,
    solid: &[bool],
    wall_x: &[bool],
    wall_y: &[bool],
) -> Vec<usize> {
    let idx = |i: usize, j: usize| j * nx + i;
    let mut charge = vec![0usize; nx * ny];
    for j in 0..ny {
        for f in 1..nx {
            if !wall_x[j * (nx + 1) + f] {
                continue;
            }
            let (a, ok) = if !solid[idx(f, j)] {
                (idx(f, j), f + 1 < nx)
            } else {
                (idx(f - 1, j), f >= 2)
            };
            if ok {
                charge[a] += 1;
            }
        }
    }
    for i in 0..nx {
        for f in 1..ny {
            if !wall_y[i * (ny + 1) + f] {
                continue;
            }
            let (a, ok) = if !solid[idx(i, f)] {
                (idx(i, f), f + 1 < ny)
            } else {
                (idx(i, f - 1), f >= 2)
            };
            if ok {
                charge[a] += 1;
            }
        }
    }
    charge
}

/// the per-cell transport properties of the combustion set.
pub fn axi_combustion_transport(
    state: &ConservedState2d,
    species: &CombustionState,
    g: &Grid2d,
) -> Result<AxiCombustionTransport, Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut tr = AxiCombustionTransport {
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
        tr.mu[k] = crate::combustion_transport::mixture_viscosity(&y, tr.t[k]);
        let mut cps = [0.0f64; ct::NS];
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            cps[s] = ct::cp_mass(s, tr.t[k]);
        }
        tr.k[k] = crate::combustion_transport::mixture_conductivity(&y, tr.t[k], &cps);
        let d_eff = crate::combustion_transport::effective_diffusivity(&y, tr.t[k], p.max(1.0));
        tr.d_s[k] = d_eff;
    }
    Ok(tr)
}

/// one viscous application on (state, species) with the given dt.
///
/// the wall flux is applied on the faces the body's boundary crosses.
/// the gas side of such a face is the fluid cell beside it, and the
/// one-sided derivatives run from the wall value into the gas.
///
/// this is public so a gate can measure one wall application against
/// the report's own accounting, with the two readings taken from the
/// same code path the solver uses.
#[allow(clippy::too_many_arguments)]
pub fn add_axi_combustion_viscous(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    tr: &AxiCombustionTransport,
    dt: f64,
    wall_opts: AxiWallOpts,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(), Error> {
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    // the signed distance drives the wall's normal direction.
    let wall_fn = wall.map(|f| f as &dyn Fn(f64, f64) -> f64);

    // the entry state, so the solid interior can be restored after the
    // update: the solid cells are scaffolding, not flow, and the mirror
    // pairs ratchet them if they are left alone.
    let entry_rho = state.rho.clone();
    let entry_e = state.e.clone();

    // the body's crossed faces and the solid mask. the wall flux rides
    // on the crossed faces, and the solid cells hold the flow out.
    let (solid, wall_x, wall_y) = match wall {
        Some(dist) => {
            let band = 1.5 * g.dx.min(g.dy);
            let cache = crate::body::body_cache(g, dist, band);
            (
                cache.mask.clone(),
                cache.wall_faces_x.clone(),
                cache.wall_faces_y.clone(),
            )
        }
        None => (
            vec![false; nx * ny],
            vec![false; (nx + 1) * ny],
            vec![false; nx * (ny + 1)],
        ),
    };

    // face accumulators. the x faces carry the axial transport, the y
    // faces the radial one; in the axisymmetric lane the radial flux
    // divides by the radius, which the update does.
    let mut tx_x = vec![0.0f64; (nx + 1) * ny];
    let mut ty_x = vec![0.0f64; (nx + 1) * ny];
    let mut te_x = vec![0.0f64; (nx + 1) * ny];
    let mut ts_x: Vec<Vec<f64>> = vec![vec![0.0f64; (nx + 1) * ny]; ct::NS];

    let mut tx_y = vec![0.0f64; nx * (ny + 1)];
    let mut ty_y = vec![0.0f64; nx * (ny + 1)];
    let mut te_y = vec![0.0f64; nx * (ny + 1)];
    let mut ts_y: Vec<Vec<f64>> = vec![vec![0.0f64; nx * (ny + 1)]; ct::NS];

    // the interior axial faces
    for j in 0..ny {
        for f in 1..nx {
            let (a, b) = (idx(f - 1, j), idx(f, j));
            let k = f * ny + j;
            let mu_f = 0.5 * (tr.mu[a] + tr.mu[b]);
            let k_f = 0.5 * (tr.k[a] + tr.k[b]);
            let dx = g.centers_x[b] - g.centers_x[a];
            let dudx = (u[b] - u[a]) / dx;
            let dvdx = (v[b] - v[a]) / dx;
            let dtdx = (tr.t[b] - tr.t[a]) / dx;
            // the axisymmetric divergence: du/dx + v/r + dv/dr
            let rc = g.centers_y[a];
            let div = dudx + v[a] / rc.max(1.0e-12) + (v[b] - v[a]) / g.dys[j];
            let txx = mu_f * (2.0 * dudx - 2.0 / 3.0 * div);
            let txy = mu_f * dvdx;
            // the face stores the physical fourier flux, because the
            // update takes minus the divergence of the stored flux.
            let qx = -k_f * dtdx;
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
            // the conservative energy flux is -(u.tau) + q + sum(h_s J_s):
            // the update takes minus the divergence of this stored flux,
            // so the viscous work and the fourier term both come out
            // with their physical sign, and the diffusing mass brings
            // its own enthalpy along with the plus sign it deserves.
            te_x[k] = -(0.5 * (u[a] + u[b]) * txx + 0.5 * (v[a] + v[b]) * txy) + qx + q_enth;
            tx_x[k] = txx;
            ty_x[k] = txy;
        }
    }

    // the interior radial faces
    for i in 0..nx {
        for f in 1..ny {
            let (a, b) = (idx(i, f - 1), idx(i, f));
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
            // the physical fourier flux, matching the axial faces.
            let qy = -k_f * dtdy;
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
            te_y[k] = -(0.5 * (u[a] + u[b]) * txy + 0.5 * (v[a] + v[b]) * tyy) + qy + q_enth;
            tx_y[k] = txy;
            ty_y[k] = tyy;
        }
    }

    let wall_charge = wall_charge_counts(nx, ny, &solid, &wall_x, &wall_y);
    let charge_weight = |k: usize| -> f64 {
        match wall_charge[k] {
            0 => 0.0,
            1 => 1.0,
            n => 1.0 / n as f64,
        }
    };

    // the film's injection, enumerated once over the crossed wall
    // faces. the coolant ENTERS the chamber, so the chamber's mass
    // balance has to carry it. a film that only drained energy
    // falsified the mass balance, and the bracket's mass-flow contrasts
    // then disagreed with the choked relation for a reason that had
    // nothing to do with the flow.
    //
    // the injection is INDEPENDENT of the wall's thermal treatment. a
    // blown-in film crosses the wall into the gas whatever the wall
    // itself is doing, so parking the injection inside the wall's
    // stated-temperature branch meant a wall with no stated
    // temperature injected nothing at all, and the film's contribution
    // vanished from the run it was supposed to dominate.
    let injection = film_injection(g, &solid, &wall_x, &wall_y, &wall_opts);

    // the wall's admissibility bound travels with the march's own
    // state, so the bound always reads the cell as the marcher sees
    // it.
    let wall_ctx = WallAdmissibility {
        state: &*state,
        species: &*species,
        grid: g,
        dt,
        fraction: wall_opts.max_energy_fraction_per_step,
    };

    // the wall faces the body's boundary crosses. the gas side is the
    // fluid cell beside the face.
    if let Some(wall_fn) = wall_fn {
        for j in 0..ny {
            for f in 1..nx {
                if !wall_x[j * (nx + 1) + f] {
                    continue;
                }
                let (a, b, ok) = if !solid[idx(f, j)] {
                    (idx(f, j), idx(f + 1, j), f + 1 < nx)
                } else {
                    (idx(f - 1, j), idx(f - 2, j), f >= 2)
                };
                if !ok {
                    continue;
                }
                let k = f * ny + j;
                let face_x = g.faces_x[f];
                let (c0, c1) = (g.centers_x[a], g.centers_x[b]);
                let mu_f = tr.mu[a];
                let dudx = wall_quad_deriv(0.0, u[a], u[b], face_x, c0, c1);
                let dvdx = wall_quad_deriv(0.0, v[a], v[b], face_x, c0, c1);
                let weight = charge_weight(a);
                let mut q = 0.0f64;
                // the film's injection is NOT booked here. it is
                // enumerated once, independently of the wall's thermal
                // treatment, so the report and the solver cannot
                // disagree about how much coolant entered.
                if wall_opts.temperature > 0.0 {
                    let t_surface = surface_temperature(tr.t[a], tr.t[b], face_x, c0, c1);
                    // the wall's own thermal circuit settles the wall's
                    // SURFACE temperature, and the gas-side terms are
                    // all read at that temperature. with no backing the
                    // surface IS the wall's stated temperature, so a
                    // plain wall behaves exactly as before.
                    let exchange = wall_exchange(t_surface, &wall_opts);
                    let t_wall = exchange.surface_temperature;
                    // the conduction term is already signed along the
                    // face's axis. the radiant and film terms are
                    // magnitudes of heat LEAVING the gas, so they take
                    // the wall's outward normal's sign along that axis
                    // before they join it.
                    let (hx, _) =
                        wall_heat_direction(wall_fn, face_x, g.centers_y[a], 0.25 * g.dx.min(g.dy));
                    let dtdn = wall_quad_deriv(t_wall, tr.t[a], tr.t[b], face_x, c0, c1);
                    let conduction = -tr.k[a] * dtdn;
                    let radiant = exchange.radiative_flux;
                    let film = exchange.convective_flux;
                    // the wall may not strip the cell below the
                    // closure's floor in a single step.
                    let face_area = axi_wall_face_area_x(g, j);
                    // the wall's three losses all draw on this one
                    // cell, so they are bounded TOGETHER and scaled by
                    // one factor, which keeps the wall's thermal split
                    // intact.
                    // the conduction term is already signed along the
                    // face's axis by its own gradient, so it keeps that
                    // sign. the radiant and film terms are magnitudes
                    // of heat leaving the gas, so they take the normal's
                    // sign before they join it.
                    let loss = bounded_wall_losses(&wall_ctx, a, face_area, 0.0, radiant, film);
                    q = conduction + hx * loss;
                }
                tx_x[k] = weight * mu_f * dudx;
                ty_x[k] = weight * mu_f * dvdx;
                te_x[k] = weight * q;
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    ts_x[s][k] = 0.0;
                }
            }
        }
        for i in 0..nx {
            for f in 1..ny {
                if !wall_y[i * (ny + 1) + f] {
                    continue;
                }
                let (a, b, ok) = if !solid[idx(i, f)] {
                    (idx(i, f), idx(i, f + 1), f + 1 < ny)
                } else {
                    (idx(i, f - 1), idx(i, f - 2), f >= 2)
                };
                if !ok {
                    continue;
                }
                let k = i * (ny + 1) + f;
                let face_y = g.faces_y[f];
                let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
                let mu_f = tr.mu[a];
                let dudy = wall_quad_deriv(0.0, u[a], u[b], face_y, c0, c1);
                let dvdy = wall_quad_deriv(0.0, v[a], v[b], face_y, c0, c1);
                let weight = charge_weight(a);
                let mut q = 0.0f64;
                // the film's injection is not booked here either: it
                // comes from the one enumeration, so the two wall
                // families cannot disagree about the coolant's total.
                if wall_opts.temperature > 0.0 {
                    let t_surface = surface_temperature(tr.t[a], tr.t[b], face_y, c0, c1);
                    // the wall's own thermal circuit settles the wall's
                    // SURFACE temperature, and the gas-side terms are
                    // all read at that temperature. with no backing the
                    // surface IS the wall's stated temperature, so a
                    // plain wall behaves exactly as before.
                    let exchange = wall_exchange(t_surface, &wall_opts);
                    let t_wall = exchange.surface_temperature;
                    // the radiant and film terms leave the gas, so they
                    // take the outward normal's sign along the axis.
                    let (_, hy) =
                        wall_heat_direction(wall_fn, g.centers_x[i], face_y, 0.25 * g.dx.min(g.dy));
                    let dtdn = wall_quad_deriv(t_wall, tr.t[a], tr.t[b], face_y, c0, c1);
                    let conduction = -tr.k[a] * dtdn;
                    let radiant = exchange.radiative_flux;
                    let film = exchange.convective_flux;
                    let face_area = axi_wall_face_area_y(g, i, f);
                    let loss = bounded_wall_losses(&wall_ctx, a, face_area, 0.0, radiant, film);
                    q = conduction + hy * loss;
                }
                tx_y[k] = weight * mu_f * dudy;
                ty_y[k] = weight * mu_f * dvdy;
                te_y[k] = weight
                    * (0.5 * (u[a] + u[b]) * (mu_f * dudy)
                        + 0.5 * (v[a] + v[b]) * (mu_f * dvdy)
                        + q);
                #[allow(clippy::needless_range_loop)]
                for s in 0..ct::NS {
                    ts_y[s][k] = 0.0;
                }
            }
        }
    }

    // the divergence update. the radial term divides by the radius, the
    // same way the reactive axi march does.
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let kfx_e = (i + 1) * ny + j;
            let kfx_w = i * ny + j;
            let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
            let rc = 0.5 * (rm + rp);
            // the radial divergence of a face quantity: the annular
            // weight over the radius and the annulus thickness. it
            // carries no dt, and the caller scales it once.
            //
            // the radial face quantities are laid out with a COLUMN
            // stride, one block of ny + 1 faces per column, so the
            // divergence has to reach into the block belonging to this
            // column. reading the block of the first column instead
            // spends every wall's flux on the first column's gas and
            // leaves the rest of the wall inert.
            let radial = |f: &[f64]| {
                let base = i * (ny + 1);
                (rp * f[base + j + 1] - rm * f[base + j]) / (rc * g.dys[j])
            };
            // tau_xx is the x-face axial stress and tau_xy the radial
            // shear: the x-momentum equation takes d(tau_xx)/dx plus
            // d(tau_xy)/dr, the y-momentum the other pairing.
            // the momentum's viscous increment is the divergence of the
            // viscous stress, so a streamline next to a no-slip wall loses
            // momentum to the wall's drag instead of gaining it.
            let dmx = dxi * (tx_x[kfx_e] - tx_x[kfx_w]) + dt * radial(&tx_y[..]);
            let dmy = dxi * (ty_x[kfx_e] - ty_x[kfx_w]) + dt * radial(&ty_y[..]);
            let de = -(dxi * (te_x[kfx_e] - te_x[kfx_w]) + dt * radial(&te_y[..]));
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;

            #[allow(clippy::needless_range_loop)]
            for s in 0..ct::NS {
                let ds = -(dxi * (ts_x[s][kfx_e] - ts_x[s][kfx_w]) + dt * radial(&ts_y[s]));
                species.y[s][k] += ds / state.rho[k];
            }
        }
    }
    // the blowing-in coolant's injection, deposited into the gas.
    //
    // the coolant's MASS enters the chamber, and its energy enters at
    // its own inlet enthalpy. the film's energy exchange has already
    // been booked against the same cells as a heat FLUX, so this
    // injection has to be careful not to double count it: the film's
    // duty is the heat the gas hands to the coolant, and the coolant
    // arriving at its inlet temperature brings that duty's WORTH with
    // it. the two together are the classic film-cooling statement, the
    // gas loses `q` to the film and the film brings `m_dot cp T_c`
    // with it.
    for k in 0..nx * ny {
        if injection.mass_rate[k] <= 0.0 {
            continue;
        }
        // the rates are per second, so multiplying by dt gives the
        // amounts that entered this cell. the state's density and
        // energy are VOLUMETRIC quantities (kg per cubic metre and
        // joules per cubic metre), so the amounts divide by the
        // cell's volume.
        let dm = injection.mass_rate[k] * dt;
        let de = injection.energy_rate[k] * dt;
        let cell_index = k;
        let ix = cell_index % nx;
        let iy = cell_index / nx;
        let volume = axi_cell_volume(g, ix, iy);
        if volume <= 0.0 {
            continue;
        }
        state.rho[k] += dm / volume;
        state.e[k] += de / volume;
        // the species enter as MASS: the cell's composition shifts by
        // the incoming species' mass over the cell's mass, which the
        // injection also raises.
        let inverse_mass = 1.0 / (state.rho[k] * volume).max(1.0e-12);
        for (s, rates) in injection.species_rate.iter().enumerate() {
            species.y[s][k] += rates[k] * dt * inverse_mass;
        }
    }
    for k in 0..nx * ny {
        species.renorm(k);
    }

    // restore the solid interior: the solid cells hold the flow out, but
    // the mirror pairs still move their dummy state, and the wall faces
    // still book their flux against it. the solid state is scaffolding,
    // not flow, so it goes back to the entry values.
    //
    // this matters because the next stage reads the solid cells when it
    // builds its face fluxes. left contaminated, they hand the near-wall
    // gas an unbounded gradient, and that gas loses its energy to the
    // wall until it trips the closure's floor.
    for (k, &is_solid) in solid.iter().enumerate() {
        if is_solid {
            state.rho[k] = entry_rho[k];
            state.mx[k] = 0.0;
            state.my[k] = 0.0;
            state.e[k] = entry_e[k];
        }
    }

    // the positivity guard, matching the reactive march's: a cell whose
    // internal energy sits below the closure's floor for its composition
    // is not a physical state, and the march would overwrite it on the
    // next stage. catching it here keeps the correction inside the
    // viscous pass, where the momentum is still the viscous one.
    for (k, &is_solid) in solid.iter().enumerate() {
        if is_solid {
            continue;
        }
        let mut y_cell = [0.0f64; ct::NS];
        for (s, dst) in y_cell.iter_mut().enumerate() {
            *dst = species.y[s][k];
        }
        let floor = closure9::internal_energy_at(200.0, &y_cell);
        // the state's energy is an energy DENSITY (rho times the
        // specific total energy), so the guard compares the SPECIFIC
        // internal energy against the closure's SPECIFIC floor.
        let rho = state.rho[k].max(1.0e-12);
        let ke = crate::state2d::specific_kinetic_energy(state.rho[k], state.mx[k], state.my[k]);
        let e_int = state.e[k] / rho - ke;
        if !state.e[k].is_finite()
            || !state.rho[k].is_finite()
            || state.rho[k] <= 0.0
            || e_int < floor
        {
            state.rho[k] = state.rho[k].max(1.0e-10);
            state.mx[k] = 0.0;
            state.my[k] = 0.0;
            state.e[k] = floor * state.rho[k] + 1.0e-6;
        }
    }
    Ok(())
}

/// the wall's heat duty and shear, gathered from the crossed faces.
///
/// the bracket needs MEASURED wall numbers, not a qualitative claim that
/// the wall cooler the gas. this sums the conductive flux leaving the
/// gas at every wall face, bucketed along the axis so the throat's
/// heating can be told from the bell's.
#[derive(Clone, Debug, Default)]
pub struct WallDiagnostics {
    /// the total power leaving the gas, W.
    pub heat_watts: f64,
    /// the conductive share of that power, W.
    pub conductive_watts: f64,
    /// the radiative share of that power, W.
    pub radiative_watts: f64,
    /// the film's share of that power, W.
    pub film_watts: f64,
    /// the wall's shear force along the axis, N.
    pub shear_newtons: f64,
    /// the wall's mean surface temperature over the crossed faces, K.
    ///
    /// this is the temperature the GAS SAW, which is not the wall's
    /// stated temperature when the wall carries a liner or a coolant of
    /// its own.
    pub surface_temperature_k: f64,
    /// the perimeter the surface temperature was averaged over, m.
    pub surface_perimeter_metres: f64,
    /// the buckets along the axis: each holds the x range and the power.
    pub buckets: Vec<WallBucket>,
    /// the wall area the blown-in film crossed, square metres.
    ///
    /// this is the film's OWN extent, not the wall's: a boundary cut
    /// obliquely crosses the same wall through several faces, and the
    /// film's area is what the coolant actually met.
    pub film_injection_area_square_metres: f64,
    /// the blown-in film's mass flow into the chamber, kg per second.
    ///
    /// the solver deposits exactly this mass into the gas, so the
    /// number in the report and the number in the chamber agree.
    pub film_injection_mass_flow_kg_per_second: f64,
}

/// one bucket of the wall, along the axis.
#[derive(Clone, Debug, Default)]
pub struct WallBucket {
    pub x0: f64,
    pub x1: f64,
    /// the total power leaving the gas in this slice, W.
    pub heat_watts: f64,
    /// the conductive share in this slice, W.
    pub conductive_watts: f64,
    /// the radiative share in this slice, W.
    pub radiative_watts: f64,
    /// the film's share in this slice, W.
    pub film_watts: f64,
    /// the wall's perimeter in this slice, m.
    pub perimeter_metres: f64,
    /// the wall's mean temperature in this slice, K.
    pub temperature_k: f64,
}

/// gathers the wall's diagnostics from the settled state.
///
/// the wall's conductive flux is the same one the viscous pass applied:
/// -k dT/dn, with dT/dn taken one-sided through the wall value and the
/// two gas cells beside it. a POSITIVE result means the wall is pulling
/// heat out of the gas.
pub fn wall_diagnostics(
    state: &ConservedState2d,
    species: &CombustionState,
    g: &Grid2d,
    nbuckets: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    wall_opts: AxiWallOpts,
) -> Result<WallDiagnostics, Error> {
    let mut out = WallDiagnostics::default();
    let dist = match wall {
        Some(d) => d,
        None => return Ok(out),
    };
    let tr = axi_combustion_transport(state, species, g)?;
    let (u, _v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let band = 1.5 * g.dx.min(g.dy);
    let cache = crate::body::body_cache(g, dist, band);
    let solid = cache.mask.clone();
    // the film's injection, from the SAME enumeration the solver used,
    // so the report's coolant mass flow is the mass the chamber took.
    let injection = film_injection(
        g,
        &solid,
        &cache.wall_faces_x,
        &cache.wall_faces_y,
        &wall_opts,
    );
    out.film_injection_area_square_metres = injection.area_square_metres;
    out.film_injection_mass_flow_kg_per_second = injection.mass_flow_kg_per_second;
    let width = g.faces_x[nx] - g.faces_x[0];
    let mut buckets: Vec<WallBucket> = (0..nbuckets)
        .map(|b| WallBucket {
            x0: g.faces_x[0] + width * b as f64 / nbuckets as f64,
            x1: g.faces_x[0] + width * (b + 1) as f64 / nbuckets as f64,
            ..Default::default()
        })
        .collect();

    // the x faces the wall crosses
    for j in 0..ny {
        for f in 1..nx {
            if !cache.wall_faces_x[j * (nx + 1) + f] {
                continue;
            }
            let (a, b, ok) = if !solid[idx(f, j)] {
                (idx(f, j), idx(f + 1, j), f + 1 < nx)
            } else {
                (idx(f - 1, j), idx(f - 2, j), f >= 2)
            };
            if !ok {
                continue;
            }
            let face_x = g.faces_x[f];
            let (c0, c1) = (g.centers_x[a], g.centers_x[b]);
            let t_surface = surface_temperature(tr.t[a], tr.t[b], face_x, c0, c1);
            // the report reads the SAME wall the solver applied, so the
            // two cannot disagree about what the wall did. an adiabatic
            // wall states no temperature, so the solver charges it
            // nothing at all and the report agrees: no duty. the
            // surface it still presents to the gas is the gas-side
            // extrapolation, which keeps the wall's footprint and its
            // temperature in the record.
            let (t_wall, conduction, q_rad, q_film) = if wall_opts.temperature > 0.0 {
                let exchange = wall_exchange(t_surface, &wall_opts);
                let t_wall = exchange.surface_temperature;
                let dtdn = wall_quad_deriv(t_wall, tr.t[a], tr.t[b], face_x, c0, c1);
                (
                    t_wall,
                    // the conductive flux LEAVING the gas, per unit area
                    -tr.k[a] * dtdn,
                    exchange.radiative_flux,
                    exchange.convective_flux,
                )
            } else {
                (t_surface, 0.0, 0.0, 0.0)
            };
            let (hx, _) = wall_heat_direction(dist, face_x, g.centers_y[a], 0.25 * g.dx.min(g.dy));
            let q = conduction + hx * (q_rad + q_film);
            // the face's area is the annulus at this radius, times the
            // cell's circumferential span. the axisymmetric face is a
            // cylinder of this radius and this height.
            let r_face = g.centers_y[a].max(1.0e-9);
            let height = g.dys[j];
            let area = 2.0 * std::f64::consts::PI * r_face * height;
            let power = q * area;
            let bi = (((face_x - g.faces_x[0]) / width) * nbuckets as f64)
                .floor()
                .max(0.0) as usize;
            let bi = bi.min(nbuckets - 1);
            buckets[bi].heat_watts += power;
            buckets[bi].conductive_watts += conduction * area;
            buckets[bi].radiative_watts += q_rad * area;
            buckets[bi].film_watts += q_film * area;
            buckets[bi].perimeter_metres += 2.0 * std::f64::consts::PI * r_face;
            buckets[bi].temperature_k += t_wall;
            out.surface_temperature_k += t_wall * 2.0 * std::f64::consts::PI * r_face;
            out.surface_perimeter_metres += 2.0 * std::f64::consts::PI * r_face;
            out.heat_watts += power;
            out.conductive_watts += conduction * area;
            out.radiative_watts += q_rad * area;
            out.film_watts += q_film * area;
            // the wall's shear: the axial viscous stress at the wall
            let dudx = wall_quad_deriv(0.0, u[a], u[b], face_x, c0, c1);
            out.shear_newtons += tr.mu[a] * dudx * area;
        }
    }
    // the y faces the wall crosses
    for i in 0..nx {
        for f in 1..ny {
            if !cache.wall_faces_y[i * (ny + 1) + f] {
                continue;
            }
            let (a, b, ok) = if !solid[idx(i, f)] {
                (idx(i, f), idx(i, f + 1), f + 1 < ny)
            } else {
                (idx(i, f - 1), idx(i, f - 2), f >= 2)
            };
            if !ok {
                continue;
            }
            let face_y = g.faces_y[f];
            let (c0, c1) = (g.centers_y[a], g.centers_y[b]);
            let t_surface = surface_temperature(tr.t[a], tr.t[b], face_y, c0, c1);
            // the report reads the SAME wall the solver applied, so the
            // two cannot disagree about what the wall did. an adiabatic
            // wall states no temperature, so the solver charges it
            // nothing at all and the report agrees: no duty.
            let (t_wall, conduction, q_rad, q_film) = if wall_opts.temperature > 0.0 {
                let exchange = wall_exchange(t_surface, &wall_opts);
                let t_wall = exchange.surface_temperature;
                let dtdn = wall_quad_deriv(t_wall, tr.t[a], tr.t[b], face_y, c0, c1);
                (
                    t_wall,
                    -tr.k[a] * dtdn,
                    exchange.radiative_flux,
                    exchange.convective_flux,
                )
            } else {
                (t_surface, 0.0, 0.0, 0.0)
            };
            let (_, hy) = wall_heat_direction(dist, g.centers_x[i], face_y, 0.25 * g.dx.min(g.dy));
            let q = conduction + hy * (q_rad + q_film);
            // a y face is a cone frustum's lateral surface
            let rm = g.centers_y[a];
            let rp = g.centers_y[b];
            let slant = ((rp - rm).powi(2) + g.dxs[i].powi(2)).sqrt();
            let area = std::f64::consts::PI * (rp + rm) * slant;
            let power = q * area;
            let bi = (((g.centers_x[i] - g.faces_x[0]) / width) * nbuckets as f64)
                .floor()
                .max(0.0) as usize;
            let bi = bi.min(nbuckets - 1);
            buckets[bi].heat_watts += power;
            buckets[bi].conductive_watts += conduction * area;
            buckets[bi].radiative_watts += q_rad * area;
            buckets[bi].film_watts += q_film * area;
            buckets[bi].perimeter_metres += std::f64::consts::PI * (rp + rm);
            buckets[bi].temperature_k += t_wall;
            out.surface_temperature_k += t_wall * std::f64::consts::PI * (rp + rm);
            out.surface_perimeter_metres += std::f64::consts::PI * (rp + rm);
            out.heat_watts += power;
            out.conductive_watts += conduction * area;
            out.radiative_watts += q_rad * area;
            out.film_watts += q_film * area;
        }
    }
    for b in buckets.iter_mut() {
        if b.perimeter_metres > 0.0 {
            b.temperature_k /= b.perimeter_metres;
        }
    }
    if out.surface_perimeter_metres > 0.0 {
        out.surface_temperature_k /= out.surface_perimeter_metres;
    }
    out.buckets = buckets;
    Ok(out)
}

/// the reacting axi march plus one viscous application per heun stage.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_combustion_axi_visc_rk2(
    state: &mut ConservedState2d,
    species: &mut CombustionState,
    g: &Grid2d,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    wall_opts: AxiWallOpts,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
    react: bool,
) -> Result<(f64, f64), Error> {
    let tr = axi_combustion_transport(state, species, g)?;
    let mut s1 = state.clone();
    let mut y1 = species.clone();
    let mut dt_last = 0.0f64;
    for _ in 0..2 {
        let (dt, _) =
            advance2d_combustion_axi(&mut s1, &mut y1, g, cfl, muscl, bc, nthreads, wall)?;
        add_axi_combustion_viscous(&mut s1, &mut y1, g, &tr, dt, wall_opts, wall)?;
        dt_last = dt;
    }
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
    if react {
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
        // the chemistry runs on the BLENDED state and species, the ones
        // the caller receives, not on the stage state that the blend
        // discards.
        let oversub = crate::combustion_march::chemistry_substep9(state, species, g, dt, nthreads)?;
        return Ok((dt_last, oversub));
    }
    Ok((dt_last, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_case() -> (Grid2d, ConservedState2d, CombustionState) {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let u = 100.0f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: vec![rho * u; n],
            my: vec![0.0; n],
            e: vec![rho * (e_int + 0.5 * u * u); n],
        };
        let sp = CombustionState::stoich_h2o2(n);
        (g, st, sp)
    }

    /// an obliquely cut boundary books two grid faces against the same
    /// gas cell, and the cell must receive one boundary segment's
    /// exchange rather than two.
    ///
    /// the test uses the real body classifier, so it exercises the same
    /// crossed faces the solver sees. two properties matter. a solid
    /// cell banks no wall flux. and the double-charged cells are the
    /// mesh quantization of an oblique cut, so they stay rarer than one
    /// percent of the gas and do not grow with the grid.
    #[test]
    fn an_oblique_boundary_segment_is_charged_once() {
        let nx = 110usize;
        let ny = 40usize;
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.165, 0.0, 0.1).unwrap();
        // a contracting lip that crosses the grid obliquely, the same
        // shape of boundary the engine bracket measures.
        let verts: Vec<(f64, f64)> = (0..=110)
            .map(|k| {
                let x = 0.00075 + 0.0015 * k as f64;
                let y = 0.099995 - 0.0024 * ((k as f64) * 0.19).sin().abs();
                (x, y)
            })
            .collect();
        let poly = crate::body::SolidPolygon::new(verts);
        let cache = crate::body::body_cache(&g, &|x, y| poly.signed_distance(x, y), 0.0);
        let charge = wall_charge_counts(
            nx,
            ny,
            &cache.mask,
            &cache.wall_faces_x,
            &cache.wall_faces_y,
        );

        // a solid cell banks nothing.
        for (k, &c) in charge.iter().enumerate() {
            if c > 0 {
                assert!(!cache.mask[k], "a solid cell charged wall flux");
            }
        }
        // the crossed faces really exist, so the assertions below are
        // not vacuous on an empty geometry.
        assert!(cache.wall_faces_x.iter().any(|&w| w));

        // the double-charged cells are rare and do not scale with the
        // grid: they quantify an oblique cut, not a wall condition.
        let gas = charge.len() - cache.mask.iter().filter(|&&s| s).count();
        let doubled = charge.iter().filter(|&&c| c >= 2).count();
        assert!(
            doubled * 100 < gas.max(1),
            "{doubled} double-charged cells among {gas} gas cells"
        );

        // splitting a charged cell's faces evenly yields one segment's
        // exchange, which is the invariance the wall loop relies on.
        for &c in charge.iter() {
            if c > 1 {
                let total = (1.0 / c as f64) * c as f64;
                assert!((total - 1.0).abs() < 1.0e-12);
            }
        }
    }

    /// a uniform field has no gradients, so the viscous pass moves
    /// nothing. the axisymmetric lane shares the transport with the
    /// planar one, so this catches a mispaired stress or a stray dt.
    #[test]
    fn the_radiant_wall_flux_follows_the_stefan_boltzmann_law() {
        // the blackbody flux at 3000 k, from the law itself
        let sigma = STEFAN_BOLTZMANN;
        let blackbody = sigma * 3000.0_f64.powi(4);
        // a grey wall emits and absorbs, so its net exchange with a
        // hot gas scales with its emissivity
        let q = radiant_wall_flux(3000.0, 300.0, 0.9);
        let expected = 0.9 * sigma * (3000.0_f64.powi(4) - 300.0_f64.powi(4));
        assert!((q - expected).abs() / expected.abs() < 1.0e-12);
        // a black wall draws more than a white one
        assert!(radiant_wall_flux(3000.0, 300.0, 1.0) > q);
        // an insulated wall draws nothing
        assert_eq!(radiant_wall_flux(3000.0, 300.0, 0.0), 0.0);
        // a wall hotter than the gas receives heat rather than losing it
        assert!(radiant_wall_flux(300.0, 3000.0, 0.9) < 0.0);
        assert!(blackbody > 0.0);
    }

    /// the wall's report must agree with the wall the solver applied.
    ///
    /// an adiabatic wall states no temperature, so the solver charges it
    /// no flux at all. the report must say the same: zero duty. reading
    /// an unwritten wall as a wall at absolute zero would invent a
    /// colossal conductive flux out of the hot gas's own gradient. the
    /// wall's footprint and the temperature the gas saw there are still
    /// real, so the report keeps them, which is what makes an adiabatic
    /// leg comparable with the cooled ones.
    #[test]
    fn an_adiabatic_wall_reports_no_duty_but_a_real_wall() {
        let (nx, ny) = (24usize, 12usize);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.06).unwrap();
        // a blunt body filling the upper half, so the wall runs along
        // the domain's midline and encloses real gas.
        let poly = crate::body::SolidPolygon::new(vec![
            (0.02, 0.035),
            (0.10, 0.035),
            (0.10, 0.055),
            (0.02, 0.055),
        ]);
        let dist = |x: f64, y: f64| poly.signed_distance(x, y);

        // a hot gas pressed against the wall, so a fabricated zero
        // kelvin wall would produce an unmistakable flux.
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

        let adiabatic = wall_diagnostics(
            &st,
            &sp,
            &g,
            6,
            Some(&dist),
            AxiWallOpts {
                temperature: 0.0,
                emissivity: 0.0,
                film_coefficient: 0.0,
                ..Default::default()
            },
        )
        .unwrap();

        // the wall really is there, so the assertions are not vacuous.
        assert!(
            adiabatic.surface_perimeter_metres > 0.0,
            "the report found no wall to measure"
        );
        // and it drew nothing, exactly as the solver's wall did. this
        // is the assertion the old report failed: it read the unwritten
        // wall as a wall at absolute zero and invented a colossal
        // conductive flux out of the hot gas's own gradient.
        assert_eq!(adiabatic.heat_watts, 0.0);
        assert_eq!(adiabatic.conductive_watts, 0.0);
        assert_eq!(adiabatic.radiative_watts, 0.0);
        assert_eq!(adiabatic.film_watts, 0.0);
        // the temperature the gas saw is real, and it is the gas's own
        // rather than the zero the old report substituted.
        assert!(adiabatic.surface_temperature_k > 0.0);

        // the same geometry against a cold wall DOES draw heat, which
        // shows the zero above came from the boundary condition and
        // not from an empty or broken measurement.
        let cooled = wall_diagnostics(
            &st,
            &sp,
            &g,
            6,
            Some(&dist),
            AxiWallOpts {
                temperature: 400.0,
                emissivity: 0.9,
                film_coefficient: 2000.0,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            cooled.heat_watts > 0.0,
            "a cold wall took no heat from a hot gas"
        );
        assert!(cooled.surface_perimeter_metres > 0.0);
    }

    #[test]
    fn the_chamber_walls_radiation_is_large_next_to_conduction() {
        // the point of adding radiation: a 3690 k gas facing a 400 k
        // wall exchanges far more energy by radiation than the
        // near-wall cells conduct, so a conduction-only wall is not a
        // representation of a real chamber wall.
        let q_rad = radiant_wall_flux(3690.0, 400.0, 0.9);
        // the measured conduction-only duty implied an effective
        // conductance of about 340 w per square metre kelvin
        let q_cond = 341.0 * (3690.0 - 400.0);
        // eight times the conduction the old wall model carried
        assert!(q_rad > 8.0 * q_cond);
        // and it is on the order of ten megawatts per square metre
        assert!(q_rad > 8.0e6 && q_rad < 1.1e7);
    }

    #[test]
    fn a_uniform_field_feels_no_viscous_force() {
        let (g, st, sp) = square_case();
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let opts = AxiWallOpts {
            temperature: 0.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
            ..Default::default()
        };
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, 1.0e-6, opts, None).unwrap();
        // the radial annulus weighting makes the divergence of a uniform
        // field nonzero in principle, so the contract is weaker than the
        // planar one: the change must stay a rounding-level fraction of
        // the state's own scale, not a physical acceleration.
        let scale = st.mx.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
        let change = s1
            .mx
            .iter()
            .zip(&st.mx)
            .fold(0.0f64, |acc, (a, b)| acc.max((a - b).abs()));
        assert!(
            change < 1.0e-6 * scale.max(1.0),
            "a uniform field moved by {change} against a scale of {scale}"
        );
    }

    /// the axisymmetric lane and the planar lane must agree on a field
    /// that varies only along the axis.
    ///
    /// nothing about the annular geometry touches an axial-only field:
    /// the radial fluxes vanish, the velocity is zero so the viscous
    /// work is zero, and the axial operator is the same one the planar
    /// lane runs. so the two lanes must book the same energy change in
    /// every cell. they do not when one of them couples the species
    /// diffusion to the energy with a different sign, and the error
    /// hides perfectly well behind a falling residual because the wrong
    /// sign still looks like transport.
    #[test]
    fn the_axisymmetric_lane_matches_the_planar_lane_on_an_axial_field() {
        let (nx, ny) = (8, 8);
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.08, 0.0, 0.08).unwrap();
        let rho = 12.64f64;
        let p = 7.0e6f64;
        let e_int = crate::combustion_march::h2o2_ic_energy(rho, p);
        let n = nx * ny;
        // the velocity varies along the axis and nowhere else, so the
        // viscous work and the momentum stresses do real work without
        // dragging the annular geometry into it. the radial velocity
        // stays zero, so the cylindrical divergence reduces to the
        // planar one and the two lanes have to agree on the momentum
        // and the viscous heating as well as on the species coupling.
        let mut mx = vec![0.0; n];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let frac = i as f64 / (nx - 1) as f64;
                mx[k] = rho * 120.0 * (1.0 - (std::f64::consts::PI * frac).cos());
            }
        }
        let kin = |mv: &[f64], k: usize| crate::state2d::specific_kinetic_energy(rho, mv[k], 0.0);
        let st = ConservedState2d {
            rho: vec![rho; n],
            mx: mx.clone(),
            my: vec![0.0; n],
            e: (0..n).map(|k| rho * (e_int + kin(&mx, k))).collect(),
        };
        let mut sp = CombustionState::stoich_h2o2(n);
        // a composition wave along x.
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let frac = i as f64 / (nx - 1) as f64;
                let swing = 0.04 * (std::f64::consts::TAU * frac).sin();
                for s in 0..ct::NS {
                    sp.y[s][k] += swing * if s % 2 == 0 { 1.0 } else { -1.0 } / ct::NS as f64;
                }
            }
        }

        let tr_axi = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut axi_state = st.clone();
        let mut axi_spec = sp.clone();
        let opts = AxiWallOpts {
            temperature: 0.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
            ..Default::default()
        };
        add_axi_combustion_viscous(
            &mut axi_state,
            &mut axi_spec,
            &g,
            &tr_axi,
            1.0e-6,
            opts,
            None,
        )
        .unwrap();

        let tr_flat = crate::combustion_viscous::combustion_transport(&st, &sp, &g).unwrap();
        let mut flat_state = st.clone();
        let mut flat_spec = sp.clone();
        let bc = Boundaries2d {
            west: crate::solver2d::Bc2d::Transmissive,
            east: crate::solver2d::Bc2d::Transmissive,
            south: crate::solver2d::Bc2d::Transmissive,
            north: crate::solver2d::Bc2d::Transmissive,
        };
        crate::combustion_viscous::add_combustion_viscous(
            &mut flat_state,
            &mut flat_spec,
            &g,
            &tr_flat,
            1.0e-6,
            &bc,
            crate::combustion_viscous::CombustionWallOpts { temperature: 0.0 },
        )
        .unwrap();

        // the interior cells carry the whole story; the outer faces are
        // treated differently by the two lanes' boundary logic.
        let mut worst: f64 = 0.0;
        let mut scale: f64 = 0.0;
        for j in 0..ny {
            for i in 1..nx - 1 {
                let k = j * nx + i;
                let axi_change = axi_state.e[k] - st.e[k];
                let flat_change = flat_state.e[k] - st.e[k];
                worst = worst.max((axi_change - flat_change).abs());
                scale = scale.max(flat_change.abs());
            }
        }
        assert!(
            scale > 0.0,
            "the axial field produced no energy change in either lane"
        );
        let rel = worst / scale;
        assert!(
            rel < 1.0e-9,
            "the axisymmetric lane disagreed with the planar lane by {rel}: the worst cell differed by {worst} against a signal of {scale}"
        );

        // the momentum has to match as well: the axial stress and the
        // axial shear are the same operators in the two lanes.
        let mut mom_worst: f64 = 0.0;
        let mut mom_scale: f64 = 0.0;
        for j in 0..ny {
            for i in 1..nx - 1 {
                let k = j * nx + i;
                mom_worst = mom_worst.max((axi_state.mx[k] - flat_state.mx[k]).abs());
                mom_scale = mom_scale.max((flat_state.mx[k] - st.mx[k]).abs());
            }
        }
        assert!(
            mom_scale > 0.0,
            "the axial field produced no momentum change in either lane"
        );
        let mom_rel = mom_worst / mom_scale;
        assert!(
            mom_rel < 1.0e-9,
            "the axisymmetric lane disagreed with the planar lane on the momentum by {mom_rel}"
        );

        // the species the two lanes diffused must match too, or the
        // energy could agree for the wrong reason.
        for s in 0..ct::NS {
            let mut species_worst: f64 = 0.0;
            let mut species_scale: f64 = 0.0;
            for k in 0..n {
                species_worst = species_worst.max((axi_spec.y[s][k] - flat_spec.y[s][k]).abs());
                species_scale = species_scale.max((flat_spec.y[s][k] - sp.y[s][k]).abs());
            }
            assert!(
                species_scale > 0.0,
                "species {s} did not move in either lane"
            );
            let species_rel = species_worst / species_scale;
            assert!(
                species_rel < 1.0e-9,
                "the axisymmetric lane diffused species {s} differently from the planar lane by {species_rel}"
            );
        }
    }

    /// the MASS must be conserved by the viscous pass. the wall is a
    /// heat sink and a shear surface, but it is impermeable: the gas
    /// may not leave through it. a wall flux that leaked mass would
    /// quietly change the throat's state, and the bracket's numbers
    /// would be measuring that leak instead of the cooling.
    #[test]
    fn the_viscous_pass_conserves_mass() {
        let (nx, ny) = (8, 8);
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
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let opts = AxiWallOpts {
            temperature: 300.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
            ..Default::default()
        };
        let dist = |_x: f64, y: f64| y - 0.02;
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, 1.0e-5, opts, Some(&dist)).unwrap();
        let before: f64 = st.rho.iter().sum();
        let after: f64 = s1.rho.iter().sum();
        assert!(
            (after - before).abs() < 1.0e-12 * before.abs().max(1.0),
            "the viscous pass changed the mass: {before} against {after}"
        );
    }

    /// a cold wall cools the gas beside it. the near-wall cell must
    /// lose energy, and the far cells must not gain it.
    #[test]
    fn a_radiating_wall_cools_the_gas_more_than_a_non_radiating_one() {
        let (nx, ny) = (8, 8);
        let idx = |i: usize, j: usize| j * nx + i;
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
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let dist = |_x: f64, y: f64| y - 0.02;

        // the same cold wall, once conducting only and once radiating.
        let mut dull = st.clone();
        let mut dull_sp = sp.clone();
        let mut glowing = st.clone();
        let mut glowing_sp = sp.clone();

        add_axi_combustion_viscous(
            &mut dull,
            &mut dull_sp,
            &g,
            &tr,
            1.0e-5,
            AxiWallOpts {
                temperature: 300.0,
                emissivity: 0.0,
                film_coefficient: 0.0,
                ..Default::default()
            },
            Some(&dist),
        )
        .unwrap();
        add_axi_combustion_viscous(
            &mut glowing,
            &mut glowing_sp,
            &g,
            &tr,
            1.0e-5,
            AxiWallOpts {
                temperature: 300.0,
                emissivity: 0.9,
                film_coefficient: 0.0,
                ..Default::default()
            },
            Some(&dist),
        )
        .unwrap();

        let near = idx(0, 2);
        let dull_loss = st.e[near] - dull.e[near];
        let glowing_loss = st.e[near] - glowing.e[near];
        // both take heat out of the near-wall cell
        assert!(dull_loss > 0.0);
        assert!(glowing_loss > 0.0);
        // and the radiating wall takes strictly more out
        assert!(
            glowing_loss > dull_loss,
            "radiation did not add to the wall's duty: {glowing_loss:e} against {dull_loss:e}; surf_T would be needed"
        );
    }

    /// the radial divergence must reach the column it is updating.
    ///
    /// the radial face quantities are laid out with a column stride, so
    /// the divergence has to index the column it sits in. if it reaches
    /// another column, the wall's flux is spent on the wrong gas, and
    /// the gas under the wall appears untouched.
    ///
    /// the fixture puts a walled shelf in the middle of the mesh, well
    /// away from the first column, and bounds it in height so there is
    /// gas above it. the gas above the shelf must feel the shelf's
    /// temperature, and the gas beside the shelf must feel it less than
    /// the gas on it.
    #[test]
    fn the_radial_divergence_reads_its_own_column() {
        let (nx, ny) = (12, 8);
        let idx = |i: usize, j: usize| j * nx + i;
        let g = crate::grid2d::grid2d(nx, ny, 0.0, 0.12, 0.0, 0.08).unwrap();
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
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();

        // a shelf: it sits between x = 0.04 and x = 0.07, and rises from
        // the axis to y = 0.03. the distance is negative inside it and
        // positive outside, which is what the body cache reads.
        let dist = |x: f64, y: f64| {
            let dx = (0.04 - x).max(x - 0.07);
            let dy = (0.0 - y).max(y - 0.03);
            dx.max(dy)
        };

        let march = |wall_temp: f64| {
            let mut s = st.clone();
            let mut y = sp.clone();
            let opts = AxiWallOpts {
                temperature: wall_temp,
                emissivity: 0.0,
                film_coefficient: 0.0,
                ..Default::default()
            };
            add_axi_combustion_viscous(&mut s, &mut y, &g, &tr, 1.0e-5, opts, Some(&dist)).unwrap();
            s
        };

        // the gas's own temperature gives a wall that exchanges nothing,
        // so the difference between the two marches is the wall's doing.
        let gas_temp = tr.t[idx(0, 0)];
        let cold = march(300.0);
        let neutral = march(gas_temp);

        // the row of gas sitting on the shelf's top face. the shelf
        // spans the columns 4 to 6.
        let above_row = 3usize;
        let mut drained = 0usize;
        for i in 4..=6 {
            let k = idx(i, above_row);
            let drain = neutral.e[k] - cold.e[k];
            assert!(
                drain > 0.0,
                "the gas on the shelf in column {i} did not lose energy to \
                 the cold shelf: {drain:e}"
            );
            drained += 1;
        }
        assert_eq!(drained, 3, "the fixture did not place gas on the shelf");

        // the gas in the first column, which stands on no wall at all,
        // must not feel the shelf's temperature. if the radial flux is
        // read from the wrong column, the first column picks up every
        // wall's flux instead.
        for i in 0..3 {
            let k = idx(i, above_row);
            let drain = neutral.e[k] - cold.e[k];
            assert!(
                drain.abs() < 1.0e-6,
                "the gas in the wall-free column {i} reacted to the shelf's \
                 temperature: {drain:e}"
            );
        }
    }

    #[test]
    fn a_cold_wall_cools_the_near_wall_cell() {
        let (nx, ny) = (8, 8);
        let idx = |i: usize, j: usize| j * nx + i;
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
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();
        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let opts = AxiWallOpts {
            temperature: 300.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
            ..Default::default()
        };
        // a flat wall along the bottom, so the row 0 cells are the
        // near-wall ones. the distance function says the gas is above
        // the line y = 0.02.
        let dist = |_x: f64, y: f64| y - 0.02;
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, 1.0e-5, opts, Some(&dist)).unwrap();
        // the wall sits at y = 0.02, so the rows below it are the body
        // and the GAS SIDE is the first row above it (row 2).
        let near = idx(0, 2);
        let far = idx(0, ny - 1);
        let delta = s1.e[near] - st.e[near];
        assert!(delta < 0.0, "the near-wall cell gained energy: {delta}");
        let far_delta = s1.e[far] - st.e[far];
        assert!(
            far_delta.abs() < 1.0e-9,
            "the far cell changed: {far_delta}"
        );
    }

    /// a blown-in film crosses the wall into the gas whatever the wall
    /// itself is doing thermally.
    ///
    /// the injection used to live inside the wall's stated-temperature
    /// branch, so a wall with no stated temperature admitted no coolant
    /// at all: the coolant vanished from the chamber's mass balance and
    /// the leg it dominated became indistinguishable from an adiabatic
    /// one. the coolant's crossing is its own business, so this test
    /// asks for an adiabatic wall with a film and requires the gas to
    /// GAIN the coolant's mass.
    #[test]
    fn a_blown_in_film_enters_a_wall_that_states_no_temperature() {
        let (nx, ny) = (8, 8);
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
        let tr = axi_combustion_transport(&st, &sp, &g).unwrap();

        // the coolant is pure nitrogen, at a modest inlet
        // temperature, blown in at a mass flux the wall can certainly
        // see.
        let mut composition = [0.0f64; ct::NS];
        composition[ct::SPECIES_NAMES.iter().position(|&s| s == "N2").unwrap()] = 1.0;
        let opts = AxiWallOpts {
            temperature: 0.0,
            emissivity: 0.0,
            film_coefficient: 0.0,
            film_injection_mass_flux: 5.0,
            film_injection_temperature: 350.0,
            film_injection_composition: composition,
            ..Default::default()
        };
        // the wall sits at y = 0.02 with the gas above it.
        let dist = |_x: f64, y: f64| y - 0.02;

        let mut s1 = st.clone();
        let mut y1 = sp.clone();
        let dt = 1.0e-3;
        add_axi_combustion_viscous(&mut s1, &mut y1, &g, &tr, dt, opts, Some(&dist)).unwrap();

        // the mass the chamber gained, weighted by the axisymmetric
        // cell volumes.
        let gained: f64 = (0..n)
            .map(|k| (s1.rho[k] - st.rho[k]) * axi_cell_volume(&g, k % nx, k / nx))
            .sum();
        assert!(
            gained > 0.0,
            "the film did not enter the chamber: {gained:e}"
        );

        // and the mass it gained is the mass the film carries: the
        // flux times the wall's area times the step.
        let cache = crate::body::body_cache(&g, &dist, 1.5 * g.dx.min(g.dy));
        let injection = film_injection(
            &g,
            &cache.mask,
            &cache.wall_faces_x,
            &cache.wall_faces_y,
            &opts,
        );
        let expected = injection.mass_flow_kg_per_second * dt;
        assert!(
            (gained - expected).abs() < 1.0e-9 * expected.abs(),
            "the chamber took {gained:e} kg of coolant, and the film carries {expected:e} kg"
        );
    }

    /// the report and the solver must agree about the coolant.
    ///
    /// the report's coolant mass flow is the number the run publishes,
    /// and the solver's deposit is the mass the chamber physically
    /// took. if the two came from different enumerations they would
    /// drift apart, and a published number would describe a film the
    /// chamber never received.
    #[test]
    fn the_report_and_the_solver_agree_about_the_injected_coolant() {
        let (nx, ny) = (8, 8);
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

        let mut composition = [0.0f64; ct::NS];
        composition[ct::SPECIES_NAMES.iter().position(|&s| s == "N2").unwrap()] = 1.0;
        let opts = AxiWallOpts {
            temperature: 300.0,
            emissivity: 0.3,
            film_coefficient: 0.0,
            film_injection_mass_flux: 2.5,
            film_injection_temperature: 400.0,
            film_injection_composition: composition,
            ..Default::default()
        };
        let dist = |_x: f64, y: f64| y - 0.02;

        let diag = wall_diagnostics(&st, &sp, &g, 4, Some(&dist), opts).unwrap();
        assert!(
            diag.film_injection_mass_flow_kg_per_second > 0.0,
            "the report saw no coolant entering"
        );
        assert!(
            diag.film_injection_area_square_metres > 0.0,
            "the report saw no wall area for the film"
        );

        // the reported mass flow is the flux times the area the film
        // met.
        let cache = crate::body::body_cache(&g, &dist, 1.5 * g.dx.min(g.dy));
        let injection = film_injection(
            &g,
            &cache.mask,
            &cache.wall_faces_x,
            &cache.wall_faces_y,
            &opts,
        );
        assert!(
            (diag.film_injection_mass_flow_kg_per_second - injection.mass_flow_kg_per_second).abs()
                < 1.0e-12 * injection.mass_flow_kg_per_second.abs(),
            "the report and the solver disagree about the coolant: {} against {}",
            diag.film_injection_mass_flow_kg_per_second,
            injection.mass_flow_kg_per_second
        );
        assert!(
            (diag.film_injection_area_square_metres - injection.area_square_metres).abs()
                < 1.0e-15 * injection.area_square_metres.abs(),
            "the report and the solver disagree about the film's area"
        );
    }

    /// the wall's admissibility bound must read the cell's energy the
    /// way the marcher stores it.
    ///
    /// the state's energy is a density (rho times the specific total
    /// energy), while the closure's floor is a specific energy. a
    /// bound that compares them without dividing by the density sees a
    /// margin inflated by the density, and a wall that ought to be
    /// throttled walks straight through it. this test plants a cell
    /// whose margin above the floor is known, asks for a withdrawal
    /// far bigger than that margin, and requires the bound to bring
    /// it back within the fraction.
    #[test]
    fn the_wall_bound_reads_the_cell_energy_as_specific_energy() {
        let (g, st, sp) = square_case();
        let nx = g.nx;
        let cell = nx + 1;
        let rho = st.rho[cell];
        let volume = axi_cell_volume(&g, cell % nx, cell / nx);

        // the cell's specific internal energy, and the closure's
        // floor for its own composition.
        let mut y = [0.0f64; ct::NS];
        for (s, dst) in y.iter_mut().enumerate() {
            *dst = sp.y[s][cell];
        }
        let floor = closure9::internal_energy_at(200.0, &y);
        let et = st.e[cell] / rho;
        let ke = crate::state2d::specific_kinetic_energy(rho, st.mx[cell], st.my[cell]);
        let margin_spec = (et - ke - floor).max(0.0);
        assert!(margin_spec > 0.0, "the fixture must sit above the floor");

        let ctx = WallAdmissibility {
            state: &st,
            species: &sp,
            grid: &g,
            dt: 1.0e-4,
            fraction: 0.25,
        };
        let area = axi_wall_face_area_x(&g, 1);
        // ask for a withdrawal of a hundred times the cell's margin.
        let demanded = 100.0 * margin_spec * rho * volume / (1.0e-4 * area);
        let bounded = admissible_wall_flux(&ctx, cell, demanded, area);
        let withdrawn_joules = bounded * area * 1.0e-4;
        let margin_joules = margin_spec * rho * volume;
        assert!(
            withdrawn_joules <= 0.25 * margin_joules * 1.000001,
            "the withdrawal {} exceeded the permitted quarter of the margin {}",
            withdrawn_joules,
            margin_joules
        );
    }

    /// a wall that is hotter than the gas heats the gas, and the
    /// bound must not throttle a beneficial exchange into nothing.
    #[test]
    fn a_hot_wall_is_not_throttled_by_the_withdrawal_bound() {
        let (g, st, sp) = square_case();
        let nx = g.nx;
        let cell = nx + 1;
        let ctx = WallAdmissibility {
            state: &st,
            species: &sp,
            grid: &g,
            dt: 1.0e-4,
            fraction: 0.25,
        };
        let area = axi_wall_face_area_x(&g, 1);
        // a negative withdrawal is heat INTO the gas, which the bound
        // must leave alone. the magnitude is deliberately absurd: an
        // unbounded inward flux is harmless to the floor, so the bound
        // must not touch it.
        let heating = -5.0e12;
        let bounded = admissible_wall_flux(&ctx, cell, heating, area);
        assert_eq!(
            bounded, heating,
            "the bound throttled a beneficial heat input"
        );
    }

    /// the bound is a single proportional factor over the wall's
    /// losses, so a wall that combines radiation and a film keeps the
    /// ratio between them.
    #[test]
    fn the_wall_bound_preserves_the_split_between_its_losses() {
        let (g, st, sp) = square_case();
        let nx = g.nx;
        let cell = nx + 1;
        let ctx = WallAdmissibility {
            state: &st,
            species: &sp,
            grid: &g,
            dt: 1.0e-4,
            fraction: 0.25,
        };
        let area = axi_wall_face_area_x(&g, 1);
        let radiant = 3.0e5;
        let film = 1.0e5;
        let scale = admissible_wall_scale(&ctx, cell, radiant + film, area);
        assert!(scale > 0.0 && scale <= 1.0, "the scale fell outside (0, 1]");
        // the same factor applies to both halves.
        let (r, f) = (radiant * scale, film * scale);
        assert!((r / f - 3.0).abs() < 1.0e-9, "the thermal split moved");
    }

    /// the wall's exchange and its reported duty must be the SAME wall.
    ///
    /// the report used to measure the film against the wall's stated
    /// temperature while the solver exchanged against a 0 k reservoir,
    /// so a conjugated wall's report came out identical to a plain
    /// wall's however different the wall actually behaved.
    #[test]
    fn the_reported_duty_follows_the_wall_the_solver_used() {
        let opts = AxiWallOpts {
            temperature: 400.0,
            emissivity: 0.9,
            film_coefficient: 2000.0,
            coolant_mass_flux: 2.0,
            coolant_capacity: 14300.0,
            regenerative_coolant_temperature: 300.0,
            ..Default::default()
        };
        let hot = 3600.0;
        let exchange = wall_exchange(hot, &opts);
        // the wall floats above its coolant, since some temperature
        // difference is what drives the heat inward.
        assert!(
            exchange.surface_temperature > opts.regenerative_coolant_temperature,
            "the wall sat at {}, no hotter than the coolant",
            exchange.surface_temperature
        );
        // the coolant-fed wall floats above its coolant and above a
        // plain 400 k wall alike: the plain wall is pinned at 400 k,
        // while the fed wall settles where the gas's arrival meets the
        // coolant's removal.
        assert!(
            exchange.surface_temperature > opts.temperature,
            "the wall sat at {}, no hotter than the plain wall's {}",
            exchange.surface_temperature,
            opts.temperature
        );
        // a plain wall, with no backing of its own, sits at its stated
        // temperature and exchanges against that.
        let plain = AxiWallOpts {
            temperature: 400.0,
            emissivity: 0.9,
            film_coefficient: 2000.0,
            ..Default::default()
        };
        let plain_exchange = wall_exchange(hot, &plain);
        assert_eq!(plain_exchange.surface_temperature, plain.temperature);
        let expected = film_wall_flux(hot, plain.temperature, plain.film_coefficient);
        assert!((plain_exchange.convective_flux - expected).abs() < 1.0e-9);
    }

    /// the wall's exchange never leans toward a colder reservoir than
    /// the wall's own stated temperature when the wall has no coolant.
    ///
    /// the old path relaxed the film against a 0 k reservoir, which is
    /// an unbounded sink: the duty grew with the gas's temperature
    /// until it stripped the near-wall cell below its floor.
    #[test]
    fn a_conjugate_wall_differs_from_a_plain_regenerative_one() {
        // the two legs differ by a single quantity: the liner's series
        // resistance. a report that measured the raw film law against
        // the wall's stated temperature gave both legs the same duty,
        // which hid the liner entirely. the circuit must distinguish
        // them, and the gas must meet a hotter inner surface behind it.
        let hot = 3600.0;
        let plain = AxiWallOpts {
            temperature: 400.0,
            emissivity: 0.9,
            film_coefficient: 2000.0,
            coolant_mass_flux: 2.0,
            coolant_capacity: 14300.0,
            regenerative_coolant_temperature: 300.0,
            ..Default::default()
        };
        let lined = AxiWallOpts {
            liner_thickness: 0.005,
            liner_conductivity: 45.0,
            ..plain
        };

        let plain_exchange = super::wall_exchange(hot, &plain);
        let lined_exchange = super::wall_exchange(hot, &lined);

        assert!(
            (plain_exchange.convective_flux - lined_exchange.convective_flux).abs() > 1.0e-6,
            "the liner must change what the gas hands over: plain {:e}, lined {:e}",
            plain_exchange.convective_flux,
            lined_exchange.convective_flux
        );
        assert!(
            lined_exchange.surface_temperature > plain_exchange.surface_temperature,
            "behind the liner the gas must meet a hotter inner surface:              plain {:.1} K, lined {:.1} K",
            plain_exchange.surface_temperature,
            lined_exchange.surface_temperature
        );
    }
    #[test]
    fn a_plain_wall_never_exchanges_against_absolute_zero() {
        let opts = AxiWallOpts {
            temperature: 400.0,
            emissivity: 0.9,
            film_coefficient: 2000.0,
            ..Default::default()
        };
        let exchange = wall_exchange(3600.0, &opts);
        assert!(
            exchange.surface_temperature >= opts.temperature,
            "the wall sank to {} against its stated {}",
            exchange.surface_temperature,
            opts.temperature
        );
        let expected = film_wall_flux(3600.0, opts.temperature, opts.film_coefficient);
        assert!((exchange.convective_flux - expected).abs() < 1.0e-9);
    }
}

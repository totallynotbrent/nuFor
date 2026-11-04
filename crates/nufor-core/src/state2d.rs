//! 2d conserved state and its conversions to primitive variables.

use crate::Error;

/// three flat cell arrays returned together from a state conversion.
type Triple = (Vec<f64>, Vec<f64>, Vec<f64>);
/// the four conserved variables over the cells of a 2d grid, flat row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservedState2d {
    /// density rho per cell.
    pub rho: Vec<f64>,
    /// x-momentum rho*u per cell.
    pub mx: Vec<f64>,
    /// y-momentum rho*v per cell.
    pub my: Vec<f64>,
    /// total energy rho*e_t per cell.
    pub e: Vec<f64>,
}

/// the specific kinetic energy of a cell, in joules per kilogram.
///
/// the state's momenta are MOMENTUM DENSITIES, rho times the velocity,
/// so the kinetic energy density is 0.5*(mx^2 + my^2)/rho and dividing
/// that by the density once more gives the specific form. dividing by
/// the density only once leaves the result in newtons, which is not an
/// energy and throws every energy floor that uses it out by a factor
/// of the density.
pub fn specific_kinetic_energy(rho: f64, mx: f64, my: f64) -> f64 {
    let r = rho.max(1.0e-12);
    0.5 * (mx * mx + my * my) / (r * r)
}

/// primitive -> conserved, rejecting non-positive or non-finite density.
pub fn prim_to_cons2d(rho: &[f64], u: &[f64], v: &[f64], et: &[f64]) -> Result<Triple, Error> {
    let n = rho.len();
    if n == 0 || u.len() != n || v.len() != n || et.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut mx = vec![0.0; n];
    let mut my = vec![0.0; n];
    let mut e = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        mx[i] = r * u[i];
        my[i] = r * v[i];
        e[i] = r * et[i];
    }
    Ok((mx, my, e))
}

/// conserved -> primitive (u, v, total specific energy), rejecting non-positive density.
pub fn cons_to_prim2d(rho: &[f64], mx: &[f64], my: &[f64], e: &[f64]) -> Result<Triple, Error> {
    let n = rho.len();
    if n == 0 || mx.len() != n || my.len() != n || e.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut u = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut et = vec![0.0; n];
    for i in 0..n {
        let r = rho[i];
        if !r.is_finite() || r <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        u[i] = mx[i] / r;
        v[i] = my[i] / r;
        et[i] = e[i] / r;
    }
    Ok((u, v, et))
}

/// a physical-validity report over a 2d conserved state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalCheck2d {
    /// true when every cell has finite, positive density and finite conserved values.
    pub ok: bool,
    /// the smallest density seen.
    pub min_rho: f64,
    /// the first bad cell index, if any.
    pub bad_cell: Option<usize>,
}

/// scan a 2d conserved state for non-finite or non-positive density.
///
/// pressure is derived elsewhere (see eos); this catches density and
/// continuity failures, which are the ones that corrupt the momentum/energy split.
pub fn check_physical2d(state: &ConservedState2d) -> PhysicalCheck2d {
    let mut chk = PhysicalCheck2d {
        ok: true,
        min_rho: f64::INFINITY,
        bad_cell: None,
    };
    for i in 0..state.rho.len() {
        let r = state.rho[i];
        chk.min_rho = chk.min_rho.min(r);
        let good = r.is_finite()
            && r > 0.0
            && state.mx[i].is_finite()
            && state.my[i].is_finite()
            && state.e[i].is_finite();
        if !good {
            chk.ok = false;
            chk.bad_cell = chk.bad_cell.or(Some(i));
        }
    }
    chk
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the specific kinetic energy follows from momentum densities.
    ///
    /// with rho = 2 and u = 3, the momentum density is 6 and the
    /// kinetic energy density is 0.5*2*3*3 = 9, so the specific form
    /// is 9/2 = 4.5.
    #[test]
    fn the_specific_kinetic_energy_is_half_velocity_squared() {
        let rho = 2.0;
        let u = 3.0;
        let v = 4.0;
        let ke = specific_kinetic_energy(rho, rho * u, rho * v);
        assert!((ke - 0.5 * (u * u + v * v)).abs() < 1.0e-12);
    }

    /// the kinetic energy a floor subtracts must not carry the density.
    ///
    /// the guards compare a cell's specific internal energy against a
    /// floor in joules per kilogram. the momentum is a momentum
    /// density, so the kinetic energy it implies has to be taken to
    /// the specific form before it is subtracted. subtracting the
    /// kinetic energy density instead leaves a term in newtons, and
    /// the error scales with the density rather than cancelling out.
    #[test]
    fn the_kinetic_energy_the_floor_sees_scales_with_the_density_not_against_it() {
        // the same velocity field at four densities.
        let u = 900.0;
        let v = 0.0;

        for rho in [0.5, 1.0, 4.0, 40.0] {
            let mx = rho * u;
            let my = rho * v;

            // what the faulty form computed: the kinetic energy
            // density, wrongly subtracted from a specific energy.
            let ke_density = 0.5 * (mx * mx + my * my) / rho;
            // what the correct form computes.
            let ke_specific = specific_kinetic_energy(rho, mx, my);

            // the correct form is the kinetic energy of the velocity
            // field, and nothing else.
            assert!(
                (ke_specific - 0.5 * (u * u + v * v)).abs() < 1.0e-9,
                "at rho = {rho} the specific form gave {ke_specific}"
            );

            // the faulty form carries the density with it, so the two
            // differ by exactly the density's contribution.
            assert!(
                (ke_density - ke_specific - 0.5 * (rho - 1.0) * (u * u + v * v)).abs() < 1.0e-6,
                "at rho = {rho} the two forms differed unexpectedly"
            );
            assert!((ke_density - ke_specific).abs() > 1.0e-6 || rho == 1.0);

            // the consequence for the guard: a cell holding an
            // internal energy of 1 MJ/kg must look the same to the
            // floor at every density, because the floor is in joules
            // per kilogram. recovering the internal energy from the
            // conserved energy must therefore reproduce it at every
            // density.
            let internal_specific = 1.0e6;
            let total_density = rho * (internal_specific + ke_specific);
            let recovered = total_density / rho - ke_specific;
            assert!(
                (recovered - internal_specific).abs() < 1.0e-6,
                "the guard sees {recovered} J/kg at rho = {rho}, \
                 not {internal_specific} J/kg"
            );

            // the faulty form misses by exactly the spurious term it
            // subtracted, which at a density of forty is enormous
            // beside the internal energy the guard meant to check.
            let recovered_old = total_density / rho - ke_density;
            let missed = internal_specific - recovered_old;
            assert!((missed - (ke_density - ke_specific)).abs() < 1.0e-6);
        }

        // at the densest case the fault is worth millions of kilojoules
        // per kilogram, which dwarfs the megajoule the guard cares
        // about. this is the magnitude that made the floors fire on
        // perfectly healthy cells.
        let rho = 40.0;
        let ke_specific = specific_kinetic_energy(rho, rho * u, rho * v);
        let ke_density = 0.5 * (rho * u) * (rho * u) / rho;
        assert!(ke_density - ke_specific > 1.0e7);
    }
}

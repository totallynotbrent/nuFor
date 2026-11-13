//! the fortran-owned park 5-species kernels.
//!
//! the fortran core owns the numerical kernels for this lane:
//! `fortran/api/nuforpark5.f90` holds the mixture closure, the park 1990
//! reaction rates, the operator-split chemistry substep and the upwind
//! species transport. this module is the rust side of that ownership:
//! the shipping march reaches the kernels through here, so there is
//! exactly one implementation in the build.
//!
//! the rust copies under `species_closure.rs`, `species_chem.rs` and
//! `species_transport.rs` stay in the tree as the independent reference
//! the parity gate compares against. they are the yardstick, not the
//! shipping code, so the two sides cannot agree trivially.

use crate::error::Error;

use std::os::raw::c_int;

/// the number of species in the park 5-species air set.
pub const NS: usize = 5;

/// the mixture fields for a batch of cells.
pub struct Fields {
    /// the mixture pressure per cell, Pa.
    pub p: Vec<f64>,
    /// the mixture temperature per cell, K.
    pub t: Vec<f64>,
    /// the mixture sound speed per cell, m/s.
    pub a: Vec<f64>,
}

/// the per-cell production rates and the chemical energy deposition.
///
/// the rates are per species, the energy deposition is per mass of mixture.
pub struct Sources {
    /// the per-species production rate, one plane per species.
    pub dydt: Vec<Vec<f64>>,
    /// the chemical energy deposited per mass of mixture, J/kg.
    pub de_chem: Vec<f64>,
}

/// the fortran side reports a malformed call with a nonzero code.
fn check(err: c_int) -> Result<(), Error> {
    if err == 0 {
        Ok(())
    } else {
        Err(Error::InvalidArgs)
    }
}

/// packs the species planes into the flat, species-major layout the
/// fortran side reads: species s occupies `y[(s-1)*n .. s*n]`.
fn flatten(planes: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut flat = vec![0.0f64; NS * n];
    for s in 0..NS {
        flat[s * n..(s + 1) * n].copy_from_slice(&planes[s]);
    }
    flat
}

/// the mixture closure for a batch of cells: pressure, temperature and
/// sound speed from the density, the internal energy per mass and the
/// composition.
pub fn closure(rho: &[f64], e_int: &[f64], planes: &[Vec<f64>]) -> Result<Fields, Error> {
    let n = rho.len();
    if n == 0 || e_int.len() != n || planes.len() < NS {
        return Err(Error::InvalidArgs);
    }
    for plane in planes.iter().take(NS) {
        if plane.len() != n {
            return Err(Error::InvalidArgs);
        }
    }
    let flat = flatten(planes, n);
    let mut p = vec![0.0f64; n];
    let mut t = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    let mut err: c_int = 0;
    unsafe {
        crate::ffi_park5::nfor_park5_closure(
            n as c_int,
            rho.as_ptr(),
            e_int.as_ptr(),
            flat.as_ptr(),
            p.as_mut_ptr(),
            t.as_mut_ptr(),
            a.as_mut_ptr(),
            &mut err,
        );
    }
    check(err)?;
    Ok(Fields { p, t, a })
}

/// the reaction source for a batch of cells at the given temperature.
pub fn omega(rho: &[f64], t: &[f64], planes: &[Vec<f64>]) -> Result<Sources, Error> {
    let n = rho.len();
    if n == 0 || t.len() != n || planes.len() < NS {
        return Err(Error::InvalidArgs);
    }
    for plane in planes.iter().take(NS) {
        if plane.len() != n {
            return Err(Error::InvalidArgs);
        }
    }
    let flat = flatten(planes, n);
    let mut dydt = vec![0.0f64; NS * n];
    let mut de_chem = vec![0.0f64; n];
    let mut err: c_int = 0;
    unsafe {
        crate::ffi_park5::nfor_park5_omega(
            n as c_int,
            rho.as_ptr(),
            t.as_ptr(),
            flat.as_ptr(),
            dydt.as_mut_ptr(),
            de_chem.as_mut_ptr(),
            &mut err,
        );
    }
    check(err)?;
    let mut dydt_planes = vec![vec![0.0f64; n]; NS];
    for s in 0..NS {
        dydt_planes[s].copy_from_slice(&dydt[s * n..(s + 1) * n]);
    }
    Ok(Sources {
        dydt: dydt_planes,
        de_chem,
    })
}

/// the operator-split chemistry substep for a batch of cells.
///
/// the composition advances in place. the two diagnostics are the
/// largest composition excursion and the worst species oversubscription
/// seen across the batch.
pub fn chem_substep(
    rho: &[f64],
    e_int: &[f64],
    dt: f64,
    planes: &mut [Vec<f64>],
) -> Result<(f64, f64), Error> {
    let n = rho.len();
    if n == 0 || e_int.len() != n || planes.len() < NS {
        return Err(Error::InvalidArgs);
    }
    for plane in planes.iter().take(NS) {
        if plane.len() != n {
            return Err(Error::InvalidArgs);
        }
    }
    let mut flat = flatten(planes, n);
    let mut worst_dy: f64 = 0.0;
    let mut worst_oversub: f64 = 0.0;
    let mut err: c_int = 0;
    unsafe {
        crate::ffi_park5::nfor_park5_chem_substep(
            n as c_int,
            rho.as_ptr(),
            e_int.as_ptr(),
            dt,
            flat.as_mut_ptr(),
            &mut worst_dy,
            &mut worst_oversub,
            &mut err,
        );
    }
    check(err)?;
    for s in 0..NS {
        planes[s].copy_from_slice(&flat[s * n..(s + 1) * n]);
    }
    Ok((worst_dy, worst_oversub))
}

/// the face layout the transport kernel reads: the x faces row major
/// (row j, then its nx+1 faces), the y faces column major (column i,
/// then its ny+1 faces).
pub struct Faces<'a> {
    /// the x-direction face mass fluxes.
    pub fx: &'a [Vec<f64>],
    /// the y-direction face mass fluxes.
    pub fy: &'a [Vec<f64>],
}

/// the upwind species transport on the flow's face mass fluxes.
///
/// the update carries rho*Y, so the per-species mass is conservative.
/// the returned residual is the largest fractional change any cell saw.
#[allow(clippy::too_many_arguments)]
pub fn transport(
    nx: usize,
    ny: usize,
    rho: &[f64],
    dt: f64,
    dxs: &[f64],
    dys: &[f64],
    faces: &Faces<'_>,
    planes: &mut [Vec<f64>],
) -> Result<f64, Error> {
    let n = nx.checked_mul(ny).ok_or(Error::InvalidArgs)?;
    if n == 0
        || rho.len() != n
        || dxs.len() != nx
        || dys.len() != ny
        || faces.fx.len() != ny
        || faces.fy.len() != nx
        || planes.len() < NS
    {
        return Err(Error::InvalidArgs);
    }
    for plane in planes.iter().take(NS) {
        if plane.len() != n {
            return Err(Error::InvalidArgs);
        }
    }
    let mut flat = flatten(planes, n);
    let fx_flat: Vec<f64> = faces.fx.iter().flatten().copied().collect();
    let fy_flat: Vec<f64> = faces.fy.iter().flatten().copied().collect();
    let mut resid: f64 = 0.0;
    let mut err: c_int = 0;
    unsafe {
        crate::ffi_park5::nfor_park5_transport(
            nx as c_int,
            ny as c_int,
            rho.as_ptr(),
            dt,
            dxs.as_ptr(),
            dys.as_ptr(),
            fx_flat.as_ptr(),
            fy_flat.as_ptr(),
            flat.as_mut_ptr(),
            &mut resid,
            &mut err,
        );
    }
    check(err)?;
    for s in 0..NS {
        planes[s].copy_from_slice(&flat[s * n..(s + 1) * n]);
    }
    Ok(resid)
}

//! the fortran park 5-species lane, called from rust.
//!
//! this is the safe surface over `fortran/api/nuforpark5.f90`. the
//! fortran module is the authority for the MPI-native phase of the
//! roadmap's rung 14: it carries the mixture closure, the park 1990
//! kinetics, the upwind species transport and the operator-split
//! chemistry substep, and the rust reference lane in
//! `species_closure.rs`, `species_chem.rs`, `species_transport.rs`
//! and `species_march.rs` is what the parity gate compares it
//! against.
//!
//! every entry point speaks the same layout the fortran side declares:
//! the species planes are contiguous, species s occupying
//! `y[(s-1)*n .. s*n]`, and the face mass fluxes arrive row-major for
//! the x faces (row j, then its nx+1 faces) and column-major for the
//! y faces (column i, then its ny+1 faces).

use std::os::raw::c_int;

use crate::error::{codes, Error};
use crate::ffi_park5::{
    nfor_cart_park5_step, nfor_park5_chem_substep, nfor_park5_closure, nfor_park5_omega,
    nfor_park5_transport,
};
use crate::park::NS;

/// flatten the rust species planes into the fortran layout: species s
/// occupies `out[(s-1)*n .. s*n]`.
fn flatten_planes(planes: &[Vec<f64>]) -> Vec<f64> {
    let n = planes.first().map_or(0, Vec::len);
    let mut out = vec![0.0f64; NS * n];
    for (s, plane) in planes.iter().take(NS).enumerate() {
        out[s * n..(s + 1) * n].copy_from_slice(plane);
    }
    out
}

/// scatter the fortran layout back into the rust species planes.
fn scatter_planes(flat: Vec<f64>, n: usize) -> Vec<Vec<f64>> {
    (0..NS).map(|s| flat[s * n..(s + 1) * n].to_vec()).collect()
}

/// translate a fortran error code.
fn check(err: c_int) -> Result<(), Error> {
    if err == codes::OK {
        Ok(())
    } else {
        crate::error::from_code(err)
    }
}

/// the mixture closure from the fortran lane: the pressure, the
/// temperature and the sound speed per cell.
///
/// `e_int` is the internal energy per mass and `planes` the species
/// mass fractions.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureFields {
    /// the mixture pressure per cell.
    pub p: Vec<f64>,
    /// the mixture temperature per cell.
    pub t: Vec<f64>,
    /// the mixture sound speed per cell.
    pub a: Vec<f64>,
}

pub fn closure(rho: &[f64], e_int: &[f64], planes: &[Vec<f64>]) -> Result<MixtureFields, Error> {
    let n = rho.len();
    if n == 0 || e_int.len() != n || planes.len() < NS {
        return Err(Error::InvalidArgs);
    }
    let flat = flatten_planes(planes);
    if flat.len() != NS * n {
        return Err(Error::InvalidArgs);
    }
    let mut p = vec![0.0f64; n];
    let mut t = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    let mut err: c_int = 0;
    unsafe {
        nfor_park5_closure(
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
    Ok(MixtureFields { p, t, a })
}

/// the reaction source from the fortran lane: the per-species
/// production rates and the chemical energy deposition per mass.
pub fn omega(
    rho: &[f64],
    t: &[f64],
    planes: &[Vec<f64>],
) -> Result<(Vec<Vec<f64>>, Vec<f64>), Error> {
    let n = rho.len();
    if n == 0 || t.len() != n || planes.len() < NS {
        return Err(Error::InvalidArgs);
    }
    let flat = flatten_planes(planes);
    if flat.len() != NS * n {
        return Err(Error::InvalidArgs);
    }
    let mut dydt_flat = vec![0.0f64; NS * n];
    let mut de_chem = vec![0.0f64; n];
    let mut err: c_int = 0;
    unsafe {
        nfor_park5_omega(
            n as c_int,
            rho.as_ptr(),
            t.as_ptr(),
            flat.as_ptr(),
            dydt_flat.as_mut_ptr(),
            de_chem.as_mut_ptr(),
            &mut err,
        );
    }
    check(err)?;
    Ok((scatter_planes(dydt_flat, n), de_chem))
}

/// the operator-split chemistry substep from the fortran lane: the
/// composition advances in place, and the two diagnostics come back.
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
    let mut flat = flatten_planes(planes);
    if flat.len() != NS * n {
        return Err(Error::InvalidArgs);
    }
    let mut worst_dy = 0.0f64;
    let mut worst_oversub = 0.0f64;
    let mut err: c_int = 0;
    unsafe {
        nfor_park5_chem_substep(
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
    for (s, plane) in planes.iter_mut().take(NS).enumerate() {
        plane.copy_from_slice(&flat[s * n..(s + 1) * n]);
    }
    Ok((worst_dy, worst_oversub))
}

/// the grid and the face fluxes a species transport runs on.
///
/// the faces follow the rust lane's own convention: `fx` carries the
/// x-face mass fluxes row by row (row j, then its nx+1 faces) and `fy`
/// the y-face fluxes column by column (column i, then its ny+1 faces).
#[derive(Clone, Debug)]
pub struct TransportGrid<'a> {
    /// the cells along x.
    pub nx: usize,
    /// the cells along y.
    pub ny: usize,
    /// the cell widths along x.
    pub dxs: &'a [f64],
    /// the cell heights along y.
    pub dys: &'a [f64],
    /// the x-face mass fluxes, row by row.
    pub fx: &'a [Vec<f64>],
    /// the y-face mass fluxes, column by column.
    pub fy: &'a [Vec<f64>],
}

impl<'a> TransportGrid<'a> {
    /// check the grid shapes against the caller's expectations.
    fn validate(&self) -> Result<usize, Error> {
        let n = self.nx * self.dys.len();
        if self.nx == 0 || self.ny == 0 || self.dxs.len() != self.nx || self.dys.len() != self.ny {
            return Err(Error::InvalidArgs);
        }
        if self.fx.len() != self.ny || self.fy.len() != self.nx {
            return Err(Error::InvalidArgs);
        }
        if self.fx.iter().any(|row| row.len() != self.nx + 1)
            || self.fy.iter().any(|col| col.len() != self.ny + 1)
        {
            return Err(Error::InvalidArgs);
        }
        Ok(n)
    }
}

/// the upwind species transport from the fortran lane, riding the flow's
/// face mass fluxes.
pub fn transport(
    grid: &TransportGrid<'_>,
    rho: &[f64],
    dt: f64,
    planes: &mut [Vec<f64>],
) -> Result<f64, Error> {
    let (nx, ny) = (grid.nx, grid.ny);
    let n = grid.validate()?;
    if rho.len() != n || planes.len() < NS {
        return Err(Error::InvalidArgs);
    }
    let mut flat = flatten_planes(planes);
    if flat.len() != NS * n {
        return Err(Error::InvalidArgs);
    }
    // the flat face arrays follow the fortran declaration: the x faces
    // row-major (row j, then its nx+1 faces), the y faces column-major
    // (column i, then its ny+1 faces).
    let fx_flat: Vec<f64> = grid.fx.iter().flatten().copied().collect();
    let fy_flat: Vec<f64> = grid.fy.iter().flatten().copied().collect();

    let mut resid = 0.0f64;
    let mut err: c_int = 0;
    unsafe {
        nfor_park5_transport(
            nx as c_int,
            ny as c_int,
            rho.as_ptr(),
            dt,
            grid.dxs.as_ptr(),
            grid.dys.as_ptr(),
            fx_flat.as_ptr(),
            fy_flat.as_ptr(),
            flat.as_mut_ptr(),
            &mut resid,
            &mut err,
        );
    }
    check(err)?;
    for (s, plane) in planes.iter_mut().take(NS).enumerate() {
        plane.copy_from_slice(&flat[s * n..(s + 1) * n]);
    }
    Ok(resid)
}

/// the boundary condition codes the cartesian march understands.
///
/// transmissive lets the flow through untouched, inflow imposes the
/// supplied freestream state, outflow copies the edge cell, and the
/// two wall cases mirror the normal velocity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CartBc {
    /// copy the neighbouring cell.
    Transmissive,
    /// impose the supplied state.
    Inflow,
    /// copy the edge cell outward.
    Outflow,
    /// mirror the normal velocity.
    SlipWall,
    /// mirror both velocities.
    NoSlipWall,
}

impl CartBc {
    fn code(self) -> c_int {
        match self {
            Self::Transmissive => 0,
            Self::Inflow => 1,
            Self::Outflow => 2,
            Self::SlipWall => 3,
            Self::NoSlipWall => 4,
        }
    }
}

/// the flow and the composition one cartesian park5 step produces.
#[derive(Clone, Debug, PartialEq)]
pub struct CartStepResult {
    /// the updated density per cell.
    pub rho: Vec<f64>,
    /// the updated x momentum per cell.
    pub mx: Vec<f64>,
    /// the updated y momentum per cell.
    pub my: Vec<f64>,
    /// the updated total energy per cell.
    pub e: Vec<f64>,
    /// the updated species mass fraction planes.
    pub planes: Vec<Vec<f64>>,
    /// the timestep the march chose.
    pub dt: f64,
}

/// one cartesian park5 step from the fortran lane.
///
/// the grid is plain `nx` by `ny` with `k = (j-1)*nx + i`; `bc_vals`
/// carries the imposed state per side as `[density, u, v, pressure,
/// internal energy]` and is ignored for the non-inflow cases.
///
/// the arguments are numerous because the c surface is flat, so the
/// lint for long argument lists is waived here.
#[allow(clippy::too_many_arguments)]
pub fn cart_step(
    nx: usize,
    ny: usize,
    dxs: &[f64],
    dys: &[f64],
    cfl: f64,
    bc_w: CartBc,
    bc_e: CartBc,
    bc_s: CartBc,
    bc_n: CartBc,
    bc_vals: &[[f64; 5]; 4],
    rho: &[f64],
    mx: &[f64],
    my: &[f64],
    e: &[f64],
    planes: &[Vec<f64>],
    y_fill: &[f64],
) -> Result<CartStepResult, Error> {
    let n = nx.checked_mul(ny).ok_or(Error::InvalidArgs)?;
    if n == 0
        || dxs.len() != nx
        || dys.len() != ny
        || rho.len() != n
        || mx.len() != n
        || my.len() != n
        || e.len() != n
        || planes.len() < NS
        || y_fill.len() < NS
    {
        return Err(Error::InvalidArgs);
    }
    let flat = flatten_planes(planes);
    if flat.len() != NS * n {
        return Err(Error::InvalidArgs);
    }
    // the fortran declaration takes the boundary states side-major.
    let mut bc_flat = [0.0f64; 20];
    for (side, vals) in bc_vals.iter().enumerate() {
        bc_flat[side * 5..side * 5 + 5].copy_from_slice(vals);
    }

    let mut rho2 = vec![0.0f64; n];
    let mut mx2 = vec![0.0f64; n];
    let mut my2 = vec![0.0f64; n];
    let mut e2 = vec![0.0f64; n];
    let mut y2 = flat.clone();
    let mut dt = 0.0f64;
    let mut err: c_int = 0;
    unsafe {
        nfor_cart_park5_step(
            nx as c_int,
            ny as c_int,
            rho.as_ptr(),
            mx.as_ptr(),
            my.as_ptr(),
            e.as_ptr(),
            flat.as_ptr(),
            dxs.as_ptr(),
            dys.as_ptr(),
            cfl,
            bc_w.code(),
            bc_e.code(),
            bc_s.code(),
            bc_n.code(),
            bc_flat.as_ptr(),
            y_fill.as_ptr(),
            rho2.as_mut_ptr(),
            mx2.as_mut_ptr(),
            my2.as_mut_ptr(),
            e2.as_mut_ptr(),
            y2.as_mut_ptr(),
            &mut dt,
            &mut err,
        );
    }
    check(err)?;
    Ok(CartStepResult {
        rho: rho2,
        mx: mx2,
        my: my2,
        e: e2,
        planes: scatter_planes(y2, n),
        dt,
    })
}

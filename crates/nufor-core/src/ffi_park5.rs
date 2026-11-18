//! raw C ABI into the statically linked fortran library for the park
//! 5-species lane.
//!
//! these are the bind(c) entry points of `fortran/api/nuforpark5.f90`,
//! which mirrors the rust park5 lane operation for operation:
//! `crates/nufor-core/src/{species_closure,species_chem,
//! species_transport,species_march}.rs`. the fortan side is the
//! authority for the MPI-native phase, the rust side is the reference
//! the parity gate compares against.

use std::os::raw::c_int;

#[link(name = "nuforkernels", kind = "static")]
extern "C" {
    /// the mixture closure: pressure, temperature and sound speed per
    /// cell from the density, the internal energy per mass and the
    /// species planes.
    pub(crate) fn nfor_park5_closure(
        n: c_int,
        rho: *const f64,
        e_int: *const f64,
        y: *const f64,
        p: *mut f64,
        t: *mut f64,
        a: *mut f64,
        err: *mut c_int,
    );
    /// the reaction source: the per-species production rates and the
    /// chemical energy deposition at the given temperature.
    pub(crate) fn nfor_park5_omega(
        n: c_int,
        rho: *const f64,
        t: *const f64,
        y: *const f64,
        dydt: *mut f64,
        de_chem: *mut f64,
        err: *mut c_int,
    );
    /// the operator-split chemistry substep: the composition advances
    /// in place and the two diagnostics come back.
    pub(crate) fn nfor_park5_chem_substep(
        n: c_int,
        rho: *const f64,
        e_int: *const f64,
        dt: f64,
        y: *mut f64,
        worst_dy: *mut f64,
        worst_oversub: *mut f64,
        err: *mut c_int,
    );
    /// the upwind species transport on the flow's face mass fluxes.
    pub(crate) fn nfor_park5_transport(
        nx: c_int,
        ny: c_int,
        rho: *const f64,
        dt: f64,
        dxs: *const f64,
        dys: *const f64,
        fx: *const f64,
        fy: *const f64,
        y: *mut f64,
        resid: *mut f64,
        err: *mut c_int,
    );
    /// the batch transport-property kernel: the wilke mixture viscosity,
    /// the frozen mixture conductivity and the effective binary
    /// diffusivity per cell.
    pub(crate) fn nfor_park5_props(
        n: c_int,
        y: *const f64,
        t: *const f64,
        p: *const f64,
        mu: *mut f64,
        kappa: *mut f64,
        d: *mut f64,
        err: *mut c_int,
    );
    /// the cartesian park5 march: the flow advances under the mixture
    /// closure, the species ride the face mass fluxes and the
    /// operator-split chemistry follows.
    pub(crate) fn nfor_cart_park5_step(
        nx: c_int,
        ny: c_int,
        rho: *const f64,
        mx: *const f64,
        my: *const f64,
        e: *const f64,
        y: *const f64,
        dxs: *const f64,
        dys: *const f64,
        faces_y: *const f64,
        geom: c_int,
        cfl: f64,
        bc_w: c_int,
        bc_e: c_int,
        bc_s: c_int,
        bc_n: c_int,
        bc_vals: *const f64,
        y_fill: *const f64,
        rho2: *mut f64,
        mx2: *mut f64,
        my2: *mut f64,
        e2: *mut f64,
        y2: *mut f64,
        dt_out: *mut f64,
        err: *mut c_int,
    );
    /// the per-species formation enthalpy per mass, J/kg.
    pub(crate) fn nfor_park5_formation_enthalpy(hf: *mut f64, err: *mut c_int);
    /// the internal energy per mass at a stated temperature and
    /// composition.
    /// the equilibrium constants of the five park5 reactions.
    pub(crate) fn nfor_park5_equilibria(n: c_int, t: *const f64, kc: *mut f64, err: *mut c_int);
    pub(crate) fn nfor_park5_internal_energy(
        n: c_int,
        t: *const f64,
        y: *const f64,
        e_out: *mut f64,
        err: *mut c_int,
    );
    /// the seeded temperature and internal energy per mass from a
    /// density, pressure and composition.
    pub(crate) fn nfor_park5_seed_energy(
        n: c_int,
        rho: *const f64,
        p: *const f64,
        y: *const f64,
        t_out: *mut f64,
        e_out: *mut f64,
        err: *mut c_int,
    );
}

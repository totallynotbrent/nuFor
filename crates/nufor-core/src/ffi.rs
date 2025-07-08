//! raw C ABI into the statically linked fortran library (bind(C) names).

use std::ffi::c_char;
use std::os::raw::c_int;

#[link(name = "nuforkernels", kind = "static")]
extern "C" {
    pub(crate) fn nfor_version(ver: *mut c_char, ver_len: c_int, err: *mut c_int);
    pub(crate) fn nfor_saxpy(count: c_int, alpha: f64, x: *const f64, y: *mut f64, err: *mut c_int);
    pub(crate) fn nfor_grid1d_init(
        n: c_int,
        xmin: f64,
        xmax: f64,
        centers: *mut f64,
        faces: *mut f64,
        dx: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_prim_to_cons(
        n: c_int,
        rho: *const f64,
        u: *const f64,
        et: *const f64,
        m: *mut f64,
        e: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_cons_to_prim(
        n: c_int,
        rho: *const f64,
        m: *const f64,
        e: *const f64,
        u: *mut f64,
        et: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_eos_pressure(
        gamma: f64,
        n: c_int,
        rho: *const f64,
        et: *const f64,
        u: *const f64,
        p: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_eos_sound_speed(
        gamma: f64,
        n: c_int,
        rho: *const f64,
        p: *const f64,
        a: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_eos_mach(
        n: c_int,
        u: *const f64,
        a: *const f64,
        mach: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_eos_temperature(
        r: f64,
        n: c_int,
        rho: *const f64,
        p: *const f64,
        t: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_hll_flux(
        gamma: f64,
        n: c_int,
        rho_l: *const f64,
        m_l: *const f64,
        e_l: *const f64,
        rho_r: *const f64,
        m_r: *const f64,
        e_r: *const f64,
        f_rho: *mut f64,
        f_m: *mut f64,
        f_e: *mut f64,
        err: *mut c_int,
    );
    pub(crate) fn nfor_cfl_dt(
        gamma: f64,
        cfl: f64,
        dx: f64,
        n: c_int,
        rho: *const f64,
        m: *const f64,
        e: *const f64,
        s_max: *mut f64,
        dt: *mut f64,
        err: *mut c_int,
    );
}

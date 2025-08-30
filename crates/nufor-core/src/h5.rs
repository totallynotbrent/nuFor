//! hdf5 read/write of the 1d state through direct calls into libhdf5.
//!
//! the system hdf5 is already part of the build toolchain, so the hdf5 crate is
//! avoided; this keeps to the datasets the solver owns and matches the schema
//! the python tooling reads with h5py (centers/rho/m/e plus gamma and time).

use std::ffi::CString;
use std::path::Path;

use crate::output::OutputState;
use crate::Error;

type H5Id = i64;
type H5Herr = i32;
type H5Size = u64;

const H5P_DEFAULT: H5Id = 0;
const H5F_ACC_TRUNC: u32 = 0x0002;
const H5F_ACC_RDONLY: u32 = 0x0000;
const NAMES: [&str; 4] = ["centers", "rho", "m", "e"];

/// the four datasets read back in a snapshot.
type Snapshot = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

#[link(name = "hdf5")]
extern "C" {
    fn H5Fcreate(name: *const i8, flags: u32, fcpl: H5Id, fapl: H5Id) -> H5Id;
    fn H5Fopen(name: *const i8, flags: u32, fapl: H5Id) -> H5Id;
    fn H5Fclose(file: H5Id) -> H5Herr;
    fn H5Screate_simple(rank: i32, dims: *const H5Size, maxdims: *const H5Size) -> H5Id;
    fn H5Sclose(space: H5Id) -> H5Herr;
    fn H5Dcreate2(
        loc: H5Id,
        name: *const i8,
        dtype: H5Id,
        space: H5Id,
        lcpl: H5Id,
        dcpl: H5Id,
        dapl: H5Id,
    ) -> H5Id;
    fn H5Dopen2(loc: H5Id, name: *const i8, dapl: H5Id) -> H5Id;
    fn H5Dwrite(
        dset: H5Id,
        mtype: H5Id,
        mspace: H5Id,
        fspace: H5Id,
        xfer: H5Id,
        buf: *const u8,
    ) -> H5Herr;
    fn H5Dread(
        dset: H5Id,
        mtype: H5Id,
        mspace: H5Id,
        fspace: H5Id,
        xfer: H5Id,
        buf: *mut u8,
    ) -> H5Herr;
    fn H5Dclose(dset: H5Id) -> H5Herr;
    fn H5Dget_space(dset: H5Id) -> H5Id;
    fn H5Sget_simple_extent_npoints(space: H5Id) -> i64;

    /// the predefined double datatype; a library-exported compile-time handle.
    static H5T_NATIVE_DOUBLE_g: H5Id;
}

fn native_double() -> H5Id {
    unsafe { H5T_NATIVE_DOUBLE_g }
}

fn cstr(s: &str) -> Result<CString, Error> {
    CString::new(s).map_err(|_| Error::InvalidArgs)
}

fn cpath(path: &Path) -> Result<CString, Error> {
    cstr(path.to_str().ok_or(Error::InvalidArgs)?)
}

fn checked_len(state: &OutputState) -> Result<usize, Error> {
    let n = state.centers.len();
    if n < 2 || state.rho.len() != n || state.m.len() != n || state.e.len() != n {
        return Err(Error::InvalidArgs);
    }
    Ok(n)
}

/// writes the state to an hdf5 file: datasets centers/rho/m/e.
pub fn write_h5(path: &Path, state: &OutputState, _time: f64) -> Result<(), Error> {
    let n = checked_len(state)?;
    let cpath = cpath(path)?;
    // SAFETY: failed hdf5 calls return a negative handle or herr, checked below.
    let file = unsafe { H5Fcreate(cpath.as_ptr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) };
    if file < 0 {
        return Err(Error::KernelFailure);
    }
    let mut ok = true;
    for (name, data) in [
        ("centers", state.centers),
        ("rho", state.rho),
        ("m", state.m),
        ("e", state.e),
    ] {
        if write_dataset(file, name, data, n).is_err() {
            ok = false;
            break;
        }
    }
    unsafe {
        H5Fclose(file);
    }
    if ok {
        Ok(())
    } else {
        Err(Error::KernelFailure)
    }
}

/// writes a single 1-d f64 dataset.
fn write_dataset(file: H5Id, name: &str, data: &[f64], n: usize) -> Result<(), Error> {
    let cname = cstr(name)?;
    let dim: H5Size = n as H5Size;
    // SAFETY: a single-rank dataspace of n doubles; the pointers are valid for the call.
    let space = unsafe { H5Screate_simple(1, &dim, core::ptr::null()) };
    if space < 0 {
        return Err(Error::KernelFailure);
    }
    let dset = unsafe {
        H5Dcreate2(
            file,
            cname.as_ptr(),
            native_double(),
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT,
        )
    };
    if dset < 0 {
        unsafe {
            H5Sclose(space);
        }
        return Err(Error::KernelFailure);
    }
    let rc = unsafe {
        H5Dwrite(
            dset,
            native_double(),
            space,
            space,
            H5P_DEFAULT,
            data.as_ptr() as *const u8,
        )
    };
    unsafe {
        H5Dclose(dset);
        H5Sclose(space);
    }
    if rc < 0 {
        return Err(Error::KernelFailure);
    }
    Ok(())
}

/// reads back the four datasets written by `write_h5` in the order centers,rho,m,e.
pub fn read_h5(path: &Path) -> Result<Snapshot, Error> {
    let cpath = cpath(path)?;
    // SAFETY: a read-only file handle with every frame validated on use.
    let file = unsafe { H5Fopen(cpath.as_ptr(), H5F_ACC_RDONLY, H5P_DEFAULT) };
    if file < 0 {
        return Err(Error::KernelFailure);
    }
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(NAMES.len());
    let mut ok = true;
    for name in NAMES {
        match read_dataset(file, name) {
            Ok(vec) => out.push(vec),
            Err(_) => {
                ok = false;
                break;
            }
        }
    }
    unsafe {
        H5Fclose(file);
    }
    if !ok || out.len() != 4 {
        return Err(Error::KernelFailure);
    }
    let mut it = out.into_iter();
    Ok((
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
    ))
}

/// reads one dataset into a vector sized by its on-disk extent.
fn read_dataset(file: H5Id, name: &str) -> Result<Vec<f64>, Error> {
    let cname = cstr(name)?;
    let dset = unsafe { H5Dopen2(file, cname.as_ptr(), H5P_DEFAULT) };
    if dset < 0 {
        return Err(Error::KernelFailure);
    }
    let space = unsafe { H5Dget_space(dset) };
    if space < 0 {
        unsafe {
            H5Dclose(dset);
        }
        return Err(Error::KernelFailure);
    }
    let n = unsafe { H5Sget_simple_extent_npoints(space) };
    if n < 0 {
        unsafe {
            H5Sclose(space);
            H5Dclose(dset);
        }
        return Err(Error::KernelFailure);
    }
    let mut buf = vec![0.0f64; n as usize];
    let rc = unsafe {
        H5Dread(
            dset,
            native_double(),
            space,
            space,
            H5P_DEFAULT,
            buf.as_mut_ptr() as *mut u8,
        )
    };
    unsafe {
        H5Sclose(space);
        H5Dclose(dset);
    }
    if rc < 0 {
        return Err(Error::KernelFailure);
    }
    Ok(buf)
}

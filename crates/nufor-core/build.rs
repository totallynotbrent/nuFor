//! build script drives CMake to compile the fortran kernel library (spec 57).

use std::env;
use std::path::PathBuf;
use std::process::Command;

/// emit a rustc link-search for the hdf5 shared library so `-lhdf5` resolves.
///
/// the crate links libhdf5 via #[link(name = "hdf5")]; the library's directory
/// is not always on the linker's default path (ubuntu puts it under
/// <libdir>/hdf5/serial), so find its directory and expose it here.
fn add_hdf5_link_search() {
    let libdir = Command::new("pkg-config")
        .args(["--variable=libdir", "hdf5"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());
    if let Some(dir) = libdir {
        println!("cargo:rustc-link-search=native={dir}");
        return;
    }
    for dir in [
        "/usr/lib/x86_64-linux-gnu/hdf5/serial",
        "/usr/lib64/hdf5/serial",
        "/usr/lib/hdf5/serial",
        "/usr/lib64",
        "/usr/lib",
        "/usr/local/lib",
    ] {
        if PathBuf::from(dir).join("libhdf5.so").exists() {
            println!("cargo:rustc-link-search=native={dir}");
            return;
        }
    }
}

fn main() {
    add_hdf5_link_search();
    add_mpi_link();
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap())
        .join("..")
        .join("..")
        .canonicalize()
        .expect("resolve repository root");
    let build_dir = PathBuf::from(env::var("OUT_DIR").unwrap()).join("nufortran");

    let run = |args: &[&str]| {
        let status = Command::new(args[0])
            .args(&args[1..])
            .current_dir(&root)
            .status()
            .expect("spawn cmake");
        assert!(status.success(), "cmake failed: {args:?}");
    };

    // the fortran build type is overridable for bounds-checked debug runs.
    let fbuild = std::env::var("NF_FORTRAN_BUILD_TYPE").unwrap_or_else(|_| "Release".into());
    run(&[
        "cmake",
        "-S",
        &root.join("fortran").to_string_lossy(),
        "-B",
        &build_dir.to_string_lossy(),
        &format!("-DCMAKE_BUILD_TYPE={fbuild}"),
    ]);
    run(&[
        "cmake",
        "--build",
        &build_dir.to_string_lossy(),
        "--config",
        "Release",
    ]);

    // CMake places libnuforkernels.a at the build-dir root (see fortran/CMakeLists.txt).
    println!(
        "cargo:rustc-link-search=native={}",
        build_dir.to_string_lossy()
    );
    println!("cargo:rustc-link-lib=static=nuforkernels");
    println!(
        "cargo:rerun-if-changed={}",
        root.join("fortran").to_string_lossy()
    );
}

/// emit link-search/link flags for openmpi so the nufor_mpi/bridge modules'
/// mpi_* symbols resolve. pkg-config first (works on both EL and ubuntu),
/// fallback to the known lib dirs.
fn add_mpi_link() {
    let ok = Command::new("pkg-config")
        .args(["--libs", "--cflags", "ompi"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).to_string());
    if let Some(libs) = ok {
        for tok in libs.split_whitespace() {
            if let Some(dir) = tok.strip_prefix("-L") {
                println!("cargo:rustc-link-search=native={dir}");
            } else if let Some(lib) = tok.strip_prefix("-l") {
                println!("cargo:rustc-link-lib=dylib={lib}");
            }
        }
        // the gfortran-compiled `use mpi` module resolves through the
        // mpifh/usempif08 shims on top of libmpi
        println!("cargo:rustc-link-lib=dylib=mpi_mpifh");
        println!("cargo:rustc-link-lib=dylib=mpi_usempif08");
        println!("cargo:rustc-link-lib=dylib=gfortran");
        return;
    }
    for dir in [
        "/usr/lib64/openmpi/lib",
        "/usr/lib/x86_64-linux-gnu/openmpi/lib",
        "/usr/lib/openmpi/lib",
    ] {
        if PathBuf::from(dir).join("libmpi.so").exists() {
            println!("cargo:rustc-link-search=native={dir}");
            println!("cargo:rustc-link-lib=dylib=mpi");
            println!("cargo:rustc-link-lib=dylib=mpi_mpifh");
            println!("cargo:rustc-link-lib=dylib=mpi_usempif08");
            return;
        }
    }
}

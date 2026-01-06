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

    run(&[
        "cmake",
        "-S",
        &root.join("fortran").to_string_lossy(),
        "-B",
        &build_dir.to_string_lossy(),
        "-DCMAKE_BUILD_TYPE=Release",
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

//! Build script: drive CMake to compile the Fortran kernel library (spec 57).
//! Kept on the stdlib rather than the `cmake` crate so the whole native build
//! stays reproducible with no third-party build dependencies.

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
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

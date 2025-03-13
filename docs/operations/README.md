# Operations

Build, install, benchmark, and deployment notes.

## Building the mixed library (PLAN step 2)

`cargo test --workspace` builds and tests the whole stack: `build.rs` drives
CMake to compile the Fortran kernels (`fortran/`), links the static archive,
and runs the FFI integration tests. Prerequisites: gfortran, CMake, cargo.

To build the Fortran kernels standalone (no Rust):

```
cmake -S fortran -B build/fortran -DCMAKE_BUILD_TYPE=Release
cmake --build build/fortran --config Release
```

Flags: Debug is `-O0 -g -Wall -Wextra -fcheck=all -fbacktrace`; Release is
`-O2 -funroll-loops`. Compiler identity is printed by CMake for the benchmark
ledger.
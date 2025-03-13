# CMake coordinates the Fortran numerical layer (spec 57). The Rust layer
# drives this via build.rs; `cmake -S fortran -B <build>` also works standalone.
cmake_minimum_required(VERSION 3.16)
project(nufortran VERSION 0.1.0 LANGUAGES Fortran)

# Debug: strong diagnostics. Release: reproducible optimization (spec 57).
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -Wextra -fcheck=all -fbacktrace")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -funroll-loops")

add_library(nuforkernels STATIC
  api/nuforkernels.f90
)

# Drop the archive at the build-dir root so build.rs has a single known path.
set_target_properties(nuforkernels PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

message(STATUS
  "nuFor Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
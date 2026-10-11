# CMake coordinates the Fortran numerical layer (spec 57). The Rust layer
# drives this via build.rs; `cmake -S fortran -B <build>` also works standalone.
cmake_minimum_required(VERSION 3.16)
project(nufortran VERSION 0.2.0 LANGUAGES Fortran)

# the archive links into rust binaries/tests: position-independent required.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# MPI: find it and compile through the mpif90 wrapper so `use mpi` resolves.
find_package(MPI COMPONENTS Fortran)
if(MPI_Fortran_FOUND)
  message(STATUS "nuFor Fortran MPI: ${MPI_Fortran_LIBRARIES}")
else()
  message(STATUS "nuFor Fortran: MPI not found; nufor_mpi module unavailable")
endif()

# Debug: strong diagnostics. Release: reproducible optimization (spec 57).
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -Wextra -fcheck=all -fbacktrace")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -funroll-loops")

if(MPI_Fortran_FOUND)
  add_library(nuforkernels STATIC
    api/nuforkernels.f90
    api/nufor_mpi.f90
    api/nuforsolver2d.f90
    api/nuforbridge.f90
  )
  target_link_libraries(nuforkernels PUBLIC MPI::MPI_Fortran)
else()
  add_library(nuforkernels STATIC
    api/nuforkernels.f90
    api/nuforsolver2d.f90
    api/nuforbridge.f90
  )
endif()

# Drop the archive at the build-dir root so build.rs has a single known path.
set_target_properties(nuforkernels PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

message(STATUS
  "nuFor Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
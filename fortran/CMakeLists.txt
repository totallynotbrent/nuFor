# CMake coordinates the Fortran numerical layer (spec 57). The Rust layer
# drives this via build.rs; `cmake -S fortran -B <build>` also works standalone.
cmake_minimum_required(VERSION 3.16)
project(nufortran VERSION 0.2.0 LANGUAGES Fortran)

# the archive links into rust binaries/tests: position-independent required.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# MPI: find it and compile through the mpif90 wrapper so `use mpi` resolves.
# the wrapper lives beside the openmpi libs on some layouts and is not on
# PATH, so point cmake at it when a plain find_package misses.
find_package(MPI COMPONENTS Fortran)
if(NOT MPI_Fortran_FOUND)
  find_program(MPIF90_WRAPPER NAMES mpif90 mpifort
    HINTS /usr/lib64/openmpi/bin /usr/lib/x86_64-linux-gnu/openmpi/bin
          /usr/lib/openmpi/bin /usr/local/bin)
  if(MPIF90_WRAPPER)
    message(STATUS "nuFor Fortran: using the mpi wrapper ${MPIF90_WRAPPER}")
    set(CMAKE_Fortran_COMPILER "${MPIF90_WRAPPER}" CACHE FILEPATH "" FORCE)
    find_package(MPI COMPONENTS Fortran)
  endif()
endif()
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
  # the solver always imports the nufor_mpi module, so the file has to
  # be in the build even when cmake cannot find the mpi package: the
  # mpi shim module ships with the openmpi installation and the
  # wrapper-less compiler finds it through the include path below.
  add_library(nuforkernels STATIC
    api/nuforkernels.f90
    api/nufor_mpi.f90
    api/nuforsolver2d.f90
    api/nuforbridge.f90
  )
  find_path(MPI_MOD_INCLUDE_DIR mpi.mod
    HINTS /usr/lib64/openmpi/lib /usr/lib/x86_64-linux-gnu/openmpi/lib
          /usr/lib/openmpi/lib /usr/local/lib
    PATH_SUFFIXES openmpi)
  if(MPI_MOD_INCLUDE_DIR)
    target_include_directories(nuforkernels PUBLIC ${MPI_MOD_INCLUDE_DIR})
  endif()
  target_link_libraries(nuforkernels PUBLIC mpi)
endif()

# Drop the archive at the build-dir root so build.rs has a single known path.
set_target_properties(nuforkernels PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

message(STATUS
  "nuFor Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
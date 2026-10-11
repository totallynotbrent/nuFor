---
title: Parallel
---

How the solver uses the CPU: shared-memory threading, the AVX2 SIMD kernels,
and the MPI-parallel Fortran march.

- [Shared-memory threading](threading.md) - the thread-sharded 2D step and its measured scaling
- [SIMD](simd.md) - the AVX2 kernels
- [MPI](mpi.md) - the multi-rank Fortran march, parity gates, and measured scaling

---
title: SIMD vectorization
---

The inner loop of a finite-volume step is a conservation update: each cell's
conserved variables are changed by the divergence of the face fluxes. That
update is element-wise over a contiguous run of cells, which makes it a clean
target for vector instructions. nuFor uses runtime-dispatched AVX2 kernels
where the host reports them, with a scalar fallback everywhere else.

## The kernel

`vectorize::apply_divergence` computes, for one conserved array, the update

```
a[i] -= dt/dx * (f[i+1] - f[i]) + y[i]
```

over `n` cells, where `f` holds the `n+1` face fluxes and `y` the transverse
contribution. On x86-64 with AVX2 present, four cells are updated per
instruction via 256-bit lanes; otherwise the same loop runs scalar. The two
paths have identical arithmetic inside a lane, so results agree to float
rounding rather than being bit-identical across machines.

The 2D solver reorganizes its update into per-row calls to this kernel: the
x-flux difference is contiguous along a row, and the strided y-flux term is
gathered into a per-row vector first, after which the update is vectorized
along the contiguous direction. The residual is still reduced scalar.

## Policy

- **Capability is reported, not assumed.** `simd_capability` reports the best
  level the host actually advertises (avx2, sse2, or scalar), surfaced in the
  web UI's benchmark panel via `/api/simd`.
- **Numerically equivalent, not identical.** The SIMD policy keeps kernels
  equal within documented floating-point tolerance across architectures; it
  does not demand bitwise identity from different SIMD sets.
- **Fallback is total.** Any host without the vector extension, or any
  non-x86 build, runs the scalar path unchanged.

## Verification

- Unit tests drive the kernel with random data and an odd length (so the
  vector tail is exercised) and assert the AVX2 result matches the scalar
  reference within `1e-12`.
- The full 2D solver suite — oblique shock, shock reflection, channel flow,
  the 2D ring — runs against the vectorized update, so the physics is
  verified on top of the new kernel.

## See also

- [[numerics/threading|Threading]]
- [[numerics/2d-euler|2D Euler solver]]
- [[performance/optimization-ledger|Optimization ledger]]
- [[numerics/memory-layout|Memory layout]]
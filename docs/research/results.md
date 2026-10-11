---
title: Research results
---

The numbers behind the solver features, measured on the shipped cases. Every
figure here regenerates from the case run logs in the repo, so the numbers age
with the code instead of drifting from it.

## Aeroshell station sweep

Four trajectory stations at 22.5 degrees angle of attack, equilibrium air,
viscous with a 300 K wall, 320x176 grid, marched to t = 15 ms. This is the
sweep that anchors the hypersonic 2D path; the 34.6 km station is the peak-q
point.

![March steps and wall time per trajectory station](station-sweep.png)

| station | steps | wall (s) | max density ratio | peak p (Pa) | max T (K) |
|---------|-------|----------|-------------------|-------------|-----------|
| 30.0 km | 1440  | 98.5     | 27.0              | 1.00e6      | 6623      |
| 34.6 km | 1148  | 168.9    | 15.6              | 5.69e5      | 7533      |
| 40.0 km | 2554  | 177.5    | 56.2              | 1.00e6      | 8131      |
| 42.5 km | 2703  | 188.4    | 64.6              | 1.00e6      | 7425      |

The step count grows with altitude because the density ratio across the bow
shock does: thinner shock layers need smaller time steps to stay stable, and
the higher stations more than double the march length. The 30 km station is
the cheapest despite the densest freestream because its shock stands close to
the body with a modest ratio.

Peak pressure saturates at the 1e6 Pa table ceiling on the two highest
stations; the closed-form read of those runs treats peak pressure as
"at or above" the ceiling rather than an exact value.

## The trim anchor

The finest shipped case is the trimmed capsule: 384x192, viscous, a 1600 K
radiation-equilibrium wall, 1750 steps to t = 2 ms in 293.7 s. Peak surface
temperature 27194 K at the nose shoulder, maximum field temperature 28935 K.
This is the case the wake visualizations come from, including the header
render on the repository front page.

## MPI scaling

The multi-rank Fortran march measured on the 34.6 km station case, 320x176,
30 steps, on the 8-core development box at its 35 W package limit:

![MPI strong scaling on the station case](mpi-scaling.png)

| ranks  | wall (s) | speedup |
|--------|----------|---------|
| serial | 1.5      | 1.00    |
| np=1   | 1.5      | 1.00    |
| np=2   | 0.8      | 1.88    |
| np=4   | 0.8      | 1.88    |

np=4 does not beat np=2 on this grid. At four ranks each column slab is 80
cells wide and the halo exchange plus the per-step body gather/broadcast cost
more than the thinner slabs save. The scaling page discusses the mechanism;
the honest summary is that the march scales to two ranks on entry-capsule
grids and needs wider slabs (a bigger grid, or a 3D decomposition) before
more ranks pay.

The parity discipline behind these numbers: every rank-count configuration is
gated against the serial Rust march, exact zero for the body-free cases and
4e-15 worst relative difference for the capsule body case.

## Shared-memory threading (prior result)

The thread-sharded 2D step, measured in the 1.5.0 era on a 384x384 mesh,
40 steps: 1.86x at four physical cores, 2.04x at six threads, backsliding to
1.75x at eight (hyperthread oversubscription). The step is memory-bound, so
the plateau tracks memory bandwidth rather than core count. Full table on the
[threading page](numerics/parallel/threading.md).

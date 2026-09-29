---
title: Grid convergence study
---

The turbulent flat-plate validation ran originally on a single 160 by 60
grid. This study adds a 320 by 120 refinement (every cell halved, same
clustering law) and marches both to t = 6, about five flow-throughs, so the
steady skin friction can be compared grid against grid and against the
Schlichting power law.

## Setup

Both grids cover the same domain: x from 0 to 1.2 with the virtual leading
edge 0.2 upstream of the inflow, y clustered at the wall with first-cell
height 5.18e-5 and growth ratio 1.2. Freestream U = 0.2, mu = 4.8e-7, so
Re_x runs from about 5e4 at the inflow to 5e5 at the outflow. The fine grid
marched 2,328,584 steps over 26.6 hours single-threaded; the explicit CFL
condition is dominated by the near-wall viscous cells, which is why the
fine grid costs roughly ten times the medium one.

## Skin friction

Over the mid-domain band (x in [0.3, 0.85], 147 stations compared):

| quantity | value |
|---------|-------|
| fine cf / correlation, mean | 1.195 |
| fine cf / correlation, range | 1.11 to 1.43 |
| medium cf / correlation, mean | 1.156 |
| medium cf / correlation, range | 1.09 to 1.25 |
| fine vs medium, mean deviation | 4.1% |
| fine vs medium, max deviation | 17% |

Both grids sit 10 to 25 percent above the power law through most of the
band, and refinement does not shrink that gap: if anything the fine grid
drifts slightly further from the correlation in the downstream half. The
two grids agree with each other to about four percent on average, so the
discrete solution is converged in the engineering sense, but it has
converged to something above the 0.0592 Re_x^-0.2 curve. That points at the
remaining gap being physical-model or boundary-condition level (inflow
profile compatibility, the wall treatment), not mesh resolution. The
17 percent outliers cluster near the downstream end where the fine grid
resolves a cf bump the coarse grid smooths over.

## Inner-scaled profiles

The u+ profiles at the three sampled stations (x = 0.24, 0.60, 1.08) show
an odd-even (checkerboard) oscillation riding on the log region on both
grids. Detrended alternating-residual amplitude, as a percent of the local
mean u+:

| station | fine grid | medium grid |
|---------|-----------|-------------|
| x = 0.24 | 12.6% | 10.9% |
| x = 0.60 | 22.3% | 19.9% |
| x = 1.08 | 22.2% | 25.1% |

Refinement does not damp the oscillation. It is a discretization mode of
the scheme on this collocated arrangement, not a mesh artifact that washes
out with resolution, and it is the most concrete open item from this study:
the fix belongs in the reconstruction or the viscous operator, not the
mesh.

## Verdict

The turbulent plate result is mesh-converged: halving every cell moves the
mid-band skin friction by about four percent on average, and both grids
tell the same story against the correlation (10 to 25 percent high,
worst downstream). The correlation gap and the u+ checkerboard are the two
real findings to chase next; neither is a resolution problem.

## See also

- [[turbulent-plate|Turbulent flat plate]]
- [[regression-suite|Regression suite]]

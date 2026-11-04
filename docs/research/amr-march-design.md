# AMR three-level march: design, gates, and the regridding hook

Status: prototype. The structure (`amr2d.rs`: `RefineBlock`, `normalize`,
`cell_gradient_magnitude`, `refine_blocks`) and a compilable proof
(`examples/amr_twolevel_gate.rs`, JSON verdict `"ok":true`) land in this
rung. This doc sketches how `solver2d.rs` turns that structure into a
real two-level time march.

## Where the seam sits today

`solver2d.rs::advance2d_capped_dts(state, grid, gamma, cfl, muscl, bc, tc)`
is one Euler step on a single uniform grid: compute primitives → `dti`
from `cfl * min(dx,dy) / smax` → assemble x/y face fluxes (`Sweep`) →
single divergence update under `TimeControl::Global(dt)` or
`TimeControl::Local(per-cell)`. The viscous wrappers (`viscous2d.rs`)
call this per Heun stage and add the diffusive cap via `max_dt`. That
function is the seam: everything above the flux/update plumbing is reusable
unchanged on either level.

## Proposed consumption pattern (two-level, ratio 2)

1. **Regrid cadence.** Every `regrid_every` steps (and at step 0), rebuild
   blocks: sample conserved `rho` on the coarse grid →
   `cell_gradient_magnitude` → `normalize` → `refine_blocks(threshold)`.
   Hysteresis: keep blocks from the previous regridding whose indicator
   stayed above `threshold_lo`, drop below `threshold_off`, so block
   boundaries don't oscillate frame to frame.

2. **Composite state.** Maintain ONE coarse `ConservedState2d` covering
   the whole domain PLUS one fine `ConservedState2d`-equivalent buffer
   per `RefineBlock` (4 cells each, stored contiguously in a slab keyed
   by `(i, j)`). Coarse cells under a block hold the *restricted*
   (area-averaged) values and are excluded from the coarse update.

3. **Prolongation (fill).** Before each fine step, initialize each block's
   2×2 patch from its parent coarse cell (constant prolongation for the
   prototype; linear limited slopes are the upgrade). Ghost exchange at
   block boundaries uses coarse neighbors where no fine neighbor exists.

4. **Subcycled march.** With coarse `dt` from the coarse CFL estimate:
   - step ALL coarse cells (masking block-covered parents) by `dt`
     through the existing `advance2d_capped_dts` path with a mask flag;
   - step each fine block by `dt/2`, twice, calling the same step body
     on the 2×2 fine grid with `hx = dx/2` (its own CFL recomputed, which
     usually permits slightly larger substeps; clamp to `dt/2` for sync);
   - after the fine substeps, **restrict**: overwrite each parent coarse
     cell with the arithmetic mean of its 4 children, then redo the coarse
     step on the covered cells OR simply adopt the restricted values as
     the synchronized state (conservative by construction either way).

5. **Synchronization order matters.** Standard Berger–Colella ordering:
   coarse step → fine substep ×2 with boundary data interpolated in
   space/time → restrict into covered coarse cells. Flux correction at
   coarse–fine interfaces (replacing the coarse face flux with the sum
   of fine face fluxes) is the conservation requirement for production;
   the prototype gate skips it because the toy operator is diffusive,
   not flux-based.

6. **Where the code goes in solver2d.** Add to `advance2d_capped_dts`'s
   caller (a new `advance2d_amr` wrapper, NOT a change to the existing
   function, keeps the uniform path bit-identical). Build/hold the block
   list, allocate fine slabs, call the step body three times (coarse +
   2 subcycles), restrict. `TimeControl::Local` composes naturally: fine
   cells take half increments inherently via the subcycling. Viscous lane
   unchanged: the wrapper adds its diffusive cap on both levels.

## Verification ladder (what the gate proves, what remains)

- Proven by `amr_twolevel_gate.rs`, indicator, blocks, subcycled step, then
  restrict agrees with the stepped coarse cell (< 1e-2); composite
  prolongation agrees with uniform-fine away from interfaces (< 2e-2).
- Third level landed: `promote_to_level2` nests a 2x2 quarter-spacing patch
  inside promoted fine blocks, subcycled at cfl/4 and restricted back.
  `amr_threelevel_gate.rs` passes: positive on all three levels after 51
  steps, hierarchical mass within 3.7 percent, and measurable level-2
  evolution (spread 1.4e-2 against the fine-parent mean). Prolongation had
  to become bilinear from the coarse neighborhood; constant filling froze
  transmissive fine patches at the initial state.
- Hysteretic regridding landed and is now consumed by the march loop:
  `maybe_regrid` consults `select_blocks_hyst(theta, keep_theta)` every
  `regrid_every` steps, rebuilds the fine level only when the hysteretic
  set actually moves, and re-invalidates so prolongation re-runs.
  `amr_regrid_gate.rs` passes on a translating pulse: 5 block-set
  migrations out of 15 regrid opportunities (steady regions held steady),
  positivity on both levels, hierarchical mass drift 0.0023 over 60
  steps. Standalone frontier churn: 2579 vs plain retriggering's 4482 on a
  sliding noisy plateau (`amr_hyst_gate.rs`).
- Remaining: conservation across coarse-fine faces on the real flux solver
  (block seams are transmissive, hence the 3.7 percent figure); limiter
  behavior at interfaces; load balancing when threading (`shard_slice`
  extends trivially, since the blocks are disjoint); picking
  `regrid_every` per physics case is still empirical.

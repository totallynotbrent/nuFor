---
title: Wall-pump diagnosis, Mach 22 overpressure (since fixed)
---

# Diagnosis: the reflected-ghost wall Riemann pump at Mach 22

Scope: the cut-cell lane's wall Riemann pump overreading near-wall
cells at hypersonic Mach (roadmap rung 13's open defect: stagnation
Cp 3.91–3.99 vs a physical ~1.7–1.9, C_A 1.86–1.91 vs Newtonian 0.603
at 160 cells; field destruction within 700 steps at 640). This note
traces the mechanism in code, ranks candidate fixes, and gives a
reproduction recipe. Resolution: candidate 1 from section (b) shipped as the
supersonic-wall branch in `cutcell_march.rs` (wall-frame normal speed
above the local sound speed switches the EB wall to the exact slip flux
`(0, p n, 0)`), with the h-box donor cap limiting redistribution to the
post-shock layer. Measured after the fix: pump ratio exactly 1.0 at every
Mach from 2 to 22 (`cases/wallpump/gate.json`), and the Mach 22 capsule
peaks at 822 Pa against the Rayleigh-pitot 624 Pa, ratio 1.32
(`cases/capsule-peak-q/gate-mach22-fixed.json`). The reproduction recipe
below remains valid for demonstrating the unfixed branch.

## (a) Mechanism, traced in code

### Which BC constructs the ghost state

The EB wall flux is built entirely in
`crates/nufor-core/src/cutcell_march.rs`, pass 1 of
`advance2d_axi_cut`:

- `cutcell_march.rs:225-271`: for every cut cell with
  `cut.eb_len[k] > 0`, the cell state is rotated into the wall-local
  frame using the EB normal `(enx, eny) = (-cut.eb_nx[k], -cut.eb_ny[k])`
  (line 227): `un = u·en + v·en`, `ut` tangential (lines 228-229).
- The ghost ("right") state mirrors the normal velocity and copies
  everything else (lines 230-257): `wr = FacePrim { rho: rho[k], u: -un,
  v: ut, p: p_wall, a: a_wall, e: e_wall }` against the interior left
  state `wl` with `u: +un`. Same `rho`, same `p`, same `e`; only the
  normal velocity flips sign.
- Line 265 solves this pair with `crate::hllc2d::hllc_flux(gamma, wl,
  wr, 0)`, and lines 266-270 fold the resulting flux into the
  divergence sum multiplied by `dt * cut.eb_len[k]`, with the same
  outward-positive sign as the open faces (line 276: `drho[k] =
  (-mass_f - wmass) / vol`).

So the wall is *not* the exact inviscid slip statement `(0, p·n, 0)`, and
that form was tried and rejected (module docs lines 12-14; comment at
lines 217-223) because it rings and collapses dt. The reflected-ghost
HLLC pair is deliberate: its numerical dissipation is load-bearing.

### What goes wrong at high Mach

In `crates/nufor-core/src/hllc2d.rs:86-88` the PVRS star pressure for
this mirror pair is

```
p_star = 0.5(pl + pr) - 0.5 (un_r - un_l) rho_bar a_bar
       = p - un * rho * a          (since pl = pr = p, un_r = -un_l = -un)
```

With the post-shock near-wall cell at Mach 22 carrying a large
post-normal-shock normal velocity `un ≈ u∞ cosθ` and a strong
compressed state, `rho·a·un` is huge, so `p_star` is pushed far above
`p`. Both `ql`/`qr` (lines 89-98) exceed 1, the estimated wave speeds
`sl, sr` (lines 99-100) widen, but the *contact* `sstar`
(line 102) and the star-region momentum flux route this inflated
`p_star`-scale momentum back into the cell through `wmx/wmy`
(lines 268-269). Because the two states differ only in the sign of
`un`, the Riemann problem is a symmetric collision: physically the
correct answer is `p_star → p_wall` with zero mass/energy flux, but
the HLLC PVRS estimate with a supersonic `|un|/a` yields a star
pressure approaching the *ram* scale `rho u²` rather than the
static pressure. The cell receives `dt · A_eb · (ρ u_n²-scale momentum
flux)` per step, repeated every step at the global dt, this pumps
the near-wall cut cells toward roughly the stagnation/ram pressure,
i.e. about twice the physical stagnation pressure when combined with
the post-shock compression the SRD neighborhood then shares. Hence
the measured Cp ≈ 2× physical.

Contributing amplifier, not the root cause: the divergence divide
`(...)/vol` (line 276-279) multiplies the pumped wall impulse by
`1/α_V` on small cut cells, and pass-2 SRD (lines 303+) then averages
the pumped state across the 9–15-cell neighborhood that straddles the
1–2-cell shock layer, spreading the overpressure along the flank
(>10⁴ negative densities at 640 cells). The h-box donor cap
(`hbox_redist.rs`, max donor distance 2) limits the pool straddle but
does not touch the pump itself.

Viscous lane: `cutcell_viscous.rs:242-295` applies the wall at the IB
via one-sided derivatives over the exact SDF distance (shear
`-ut/dist`, Fourier conduction `-k(T_cell - T_wall)/dist`). These
terms dissipate correctly and are not implicated; the overpressure
originates in the inviscid EB flux of `cutcell_march.rs`.

## (b) Candidate fixes, ranked

1. **Wall-pressure-clamped ghost (recommended).** Keep the
   reflected-ghost pair for its dissipation, but construct the ghost
   with the *wall* pressure obtained from the physical wall condition
   instead of letting HLLC generate a ram-scale star pressure , 
   concretely, after rotating to the wall frame, set the pair's
   normal velocities to a small slip value (|un| ≤ ε·a_wall, e.g.
   clip un to ±0.05·a) before calling `hllc_flux`, or evaluate the
   exact slip flux `(0, p n, 0)` blended with a Rusanov/HLL flux at
   the wall so the star pressure is bounded by `~p`. Bounding `un_r -
   un_l` bounds `p_star` via hllc2d.rs:88 directly at the source of
   the pump. Small change confined to `cutcell_march.rs:230-265`;
   retains stability (still a genuine two-state Riemann solve), kills
   the Mach-dependent overpressure analytically.
2. **Local wall-timestep / impulse limiter on the EB flux.** Cap the
   wall momentum impulse per step at the physical wall pressure
   impulse, `min(|wm|, α_V·Δx·p_max)·sign`, with `p_max` the maximum
   plausible wall pressure (e.g. post-shock stagnation estimate).
   Sturdy fallback if (1) worsens the stagnation-point accuracy, but
   introduces a heuristic clamp that must be tuned per regime.
3. **Shock-layer-aware redistribution (fix the straddle, not the
   pump).** Shrink SRD neighborhoods to stay inside the 1–2 cell
   post-shock layer (h-box with radius 1, or delta-redistribution
   revived for the non-axisymmetric flank, noting the documented
   near-axis vacuum failure at lines 296-302 must be avoided). Does
   not remove the ~2× pump at the wall cell itself; treats the symptom
   downstream. Necessary eventually (rung 13's stated "wall-local
   scheme"), insufficient alone.

## (c) Reproduction recipe

1. Branch off `experimental`; do not modify sources.
2. Case: `cases/capsule-peak-q/case.toml` (Mach 22 perfect gas,
   axisymmetric Euler, sphere-cone, 320x160; copy to scratch and drop
   to 160 cells wide for the quick over-read, keep 640-wide to watch
   the field die within ~700 steps).
3. Run the cut-cell march (`advance2d_axi_cut`) with
   `LD_LIBRARY_PATH=/usr/lib64/openmpi/lib` and
   `timeout 300 cargo run --release -- run <scratch-case>` (foreground
   cargo beyond ~420 s times out; use `timeout 300`).
4. Probe: dump `p` and compute surface Cp
   `(p_wall - p_inf)/(0.5 ρ_inf u_inf²)` at the stagnation cell; the
   defect reads Cp ≈ 3.9–4.0 vs physical ≈ 1.7–1.9, and C_A ≈ 1.86+
   vs Newtonian 0.603. Log `cut.eb_len`, the wall-frame `un`, and
   `p_star` from `hllc_flux` at the stagnation cut cell for one step to
   show `p_star ≫ p_wall` at Mach 22 and `p_star ≈ p_wall` at Mach 2.
5. Optional instrumented micro-repro (scratch only, not in repo): a
   single cut cell with a post-shock state fed through lines 225-271
   at increasing Mach, plotting star pressure vs Mach.

## Outcome

The ranked fix 1 landed in rung 13; see the banner above and the gates
under `cases/wallpump/` and `cases/capsule-peak-q/`. This note stays as
the mechanism record behind that fix.

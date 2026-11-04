# Laminar viscous terms on the unstructured triangle duals

Scope: survey the viscous machinery, name the seam, and lay out the
implementation plan for a laminar viscous flux on `Ugrid`. The plan has
since shipped; see the landed outcome under the gate section below.

## (a) What viscous machinery exists today

### Structured 2D (Cartesian `Grid2d`), the reusable physics
All in `crates/nufor-core/src/viscous2d.rs`:
- `sutherland_mu(mu0, t, t0, s)`, `viscous2d.rs:17`: temperature-dependent laminar law.
- Stress definition `viscous_flux`, `viscous2d.rs:165-177`: Newtonian stress tensor
  `txx = mu*(2 du/dx - 2/3 div)`, `tyy = mu*(2 dv/dy - 2/3 div)`, `txy = mu*(du/dy + dv/dx)`,
  heat flux `q = mu*kappa*grad(T)` with `kappa = conductivity/mu` (`Cp/Pr`), signs chosen
  so a hot peak cools (`viscous2d.rs:170`).
- Eddy/laminar split `face_flux_split`, `viscous2d.rs:196-215`: stresses use
  `mu_lam + mu_t_face`; heat splits `mu_lam*kappa + mu_t_face*kappa_t` (zero `mu_t`
  is the laminar limit).
- Entry points `add_viscous` (`viscous2d.rs:402`) and boundary-aware `add_viscous_bc`
  (`viscous2d.rs:423`): centered cell gradients (`cell_grads`, `:217`), Lagrangian
  3-point derivatives on true positions (`quad_deriv_at`, `:450`).
- Wall treatment: one-sided quadratic wall shear and `q_wall` pinned to the no-slip /
  isothermal wall value (`body_wall_x`-family, `viscous2d.rs:480-588`; reduced
  one-sided form `(9 u0 - u1)/(3 dy)` on uniform grids, `:467`). No-slip ⇒ stress
  work at the wall vanishes (`:533`).

### Cut-cell axisymmetric variant, `cutcell_viscous.rs`
- `viscosity_scale`, `cutcell_viscous.rs:51`: scales a case-named absolute `mu` to
  the molecular Wilke-mixture value.
- `transport_props`, `cutcell_viscous.rs:90-123`: per-cell `T`, `mu` (Wilke
  mixture), `k`, species diffusivities; combustion-set only.
- `add_cut_combustion_viscous`, `cutcell_viscous.rs:133`: aperture-weighted face
  diffusion `mu_f * (phi_b - phi_a)/dx` accumulated per face (`cutcell_viscous.rs:174-180`)
  plus exact-SDF-distance wall fluxes at the immersed boundary (`:125-131`). Pattern to
  copy: per-face scalar flux, applied symmetrically to both cells.

### Turbulence (SA), structured only
- Closures in `sa.rs`: `fv1`, `stilde`, `fw`, `source` (`:34-97`), eddy viscosity
  `mu_t = rho*nu_tilde*fv1(chi)` (`sa.rs:103`).
- March `advance_turb`, `turb2d.rs:162` over `Grid2d` with padded ghosts; wall
  distance from `wall_dist.rs` (structured stencil). Nothing on `Ugrid` yet.

### Unstructured lane status, `solver_u.rs` / `ugrid.rs`
- `advance_ugrid` (`solver_u.rs:78`) and `advance_ugrid_recon` (`solver_u.rs:304`,
  MUSCL/minmod, `:569`): HLLC inviscid flux only, accumulated as `± len * flux` per
  face (`:217-226`, `:528-537`). Slip wall group 1 is the exact `[0, p n, 0]` flux
  (`:148-154`, `:397-401`). Freestream group 2, outflow group 3, transmissive 0.
- Geometry (`ugrid.rs`): faces carry the signed area vector `face_dx/face_dy`
  (= length × unit normal, `ugrid.rs:23-26`), `face_left/face_right` adjacency
  (`:19-22`), `cell_area`/`cell_center` (`:11-14`), `cell_faces` (`:16`).
- Gradient tool already present: `green_gauss_gradient`, `solver_u.rs:275-297`,
  unlimited Green–Gauss `∇φ = (1/A) Σ φ_f S_f`; used today for MUSCL hints.

Net: the stress/energy formulas, Prandtl splitting, Sutherland law, and wall fluxes
exist and are portable; what is missing on `Ugrid` is purely the geometric plumbing
(cell gradients from the dual, face-normal derivatives, wall distance).

## (b) The closest seam

The face loop of `advance_ugrid*` is the seam. Every face already knows its area
vector `S = (face_dx, face_dy)` and its two cells, and the accumulator pattern
`acc[left] -= len*F`, `acc[right] += len*F` (solver_u.rs:217-226) is exactly the
finite-volume form of a viscous flux integral `∮ F_v · n dS`. So viscous terms plug
in as a second per-face flux added to the same accumulators:

```
interior face f:  F_v(f) = [0, tau_f · n̂, (u tau)_f · n̂ + k_f gradT_f · n̂] * len_f
                  acc[L] -= F_v ; acc[R] += F_v          (sign opposite to Euler flux direction convention)
boundary wall (group 1): F_v = [0, tau_wall · n̂, q_wall]  (one-sided, wall-value-pinned)
other boundaries (0/2/3): zero-gradient viscous flux (consistent with current inviscid ghosts)
```

Gradient strategy on the dual (two viable options, see ranking): either reuse
`green_gauss_gradient` per primitive field (already written, solver_u.rs:275) and
average the two cell gradients onto the face, or use the direct two-point normal
derivative `(φ_R − φ_L) / |x_R − x_L|` projected on `n̂`. Both only need
`cell_center`, `face_left/right`, `face_dx/dy`, all present in `Ugrid`.

Wall quantities need one new datum not currently in `Ugrid`: the wall distance of a
boundary cell. Minimal option: perpendicular distance from `cell_center[L]` to the
boundary face (computable from `face_nodes` + node coords, or approximately
`cell_area / (sum of boundary face lengths / 2)`); refine later with a proper
nearest-node wall distance like `wall_dist.rs`.

## (c) Minimal implementation plan (laminar Newtonian on triangles)

New module `crates/nufor-core/src/viscous_u.rs` (mirrors `viscous2d.rs` naming):

1. **`pub fn sutherland_mu` reuse**: call `viscous2d::sutherland_mu` directly, no duplication.
2. **Per-cell scalars**: `T` from the EOS/model (perfect: `T = p/(rho R)`; eqair via
   `eqair_cea`), `mu_c` from Sutherland or constant-case mu, `kappa = Cp/Pr`
   (perfect-gas form documented at `viscous2d.rs:21-24`).
3. **Cell gradients**: one generalized helper calling `green_gauss_gradient`
   (move it from `solver_u.rs:275` into `viscous_u.rs` or make it `pub(crate)`) for
   `u, v, T` with arithmetic face averages, unlimited first pass, limiter optional
   later via the existing `minmod` (`solver_u.rs:569`).
4. **Interior face flux**: `grad_f = 0.5*(grad_L + grad_R)`; decompose into normal/
   tangential parts using `(nx, ny)` from `face_dx/dy` (pattern already in
   `flux_on_normal`, `solver_u.rs:49-71`); assemble the stress exactly as
   `viscous2d.rs:165-177`; contract with `S_f`; add to the same accumulators.
5. **Group-1 wall faces**: no-slip/isothermal one-sided flux, `tau_wall = mu_f *
   (u_L − 0)/d_wall` tangential, `q_wall = mu_f*kappa*dT/d_wall` with `d_wall` the
   face-perpendicular distance (step 0 datum above); zero stress work (velocity is 0
   at the wall, cf. `viscous2d.rs:533`).
6. **Hook**: in `advance_ugrid`/`advance_ugrid_recon`, after the HLLC loop, run the
   viscous face loop with the same `dt` (explicit, dt from the existing `cfl*h/smax`
   computation at `solver_u.rs:130`/`:359`). Guard with a `viscous: bool` +
   `mu, pr` parameter; `advance_ugrid`'s signature grows or a wrapper
   `advance_ugrid_viscous` is added so existing callers/tests are untouched.
7. **Stability**: diffuse explicitly ⇒ additional dt constraint
   `dt ≤ cfl_v * h²/(2ν)` folded into the existing dt calc.

### Ranked approaches
1. **Green–Gauss cell gradients + face-averaged contraction (recommended)**, reuses
   `green_gauss_gradient` verbatim, second-order on smooth meshes, single code path
   for recon and non-recon. Cost: 3 extra gradient passes.
2. **Direct two-point normal derivative `(φ_R−φ_L)·n̂/|Δx|`**, cheaper, exact normal
   derivative on orthogonal meshes, but only first-order accurate on skewed
   triangles; fall back here if (1) shows mesh-dependence errors.
3. **Least-squares gradients from cell neighbors**, most accurate on irregular
   triangles, but requires a neighbor-list (only `cell_faces` exists; derivable) and
   more code; defer until (1)/(2) fail a skewness study.
4. **Reusing `viscous2d::add_viscous` wholesale**, rejected: hard-wired to
   `Grid2d` indexing (`viscous2d.rs:217+`); only the scalar kernels (steps 1–2, 4)
   should be lifted.

Turbulent extension (later rung): `mu_t` from `sa::eddy_viscosity` fed into the same
face flux via the `face_flux_split` pattern (`viscous2d.rs:196`); needs SA march +
wall distance on `Ugrid` first, out of scope for the laminar minimal pass.

### Landed outcome

This plan shipped as `viscous_u` driving Newtonian stress and conduction
from green-gauss gradients on the `Ugrid` faces. The gate results, from
`cases/ugrid-viscous/gate.json`: pairwise-conservative construction
(exactly zero residual on a uniform field), energy sums holding to
1.7e-12 over the decay run, and monotonic diffusion decay of a
sinusoidal energy perturbation. Curved-geometry validation followed in
`cases/ugrid-annulus/gate.json`: a 2304-triangle quarter annulus tiles
its arc-chord area to 1.8e-4 and conserves total energy to 2.8e-15 over
50 viscous steps with positivity. Turbulent extension (SA eddy viscosity
on the duals) remains the noted later step.

### Gate
Verification gate, in order:
1. **Uniform-field preservation**: uniform `(rho,u,v,p,T)` on `Ugrid::cartesian_quads`
   with all-wall or periodic-transmissive boundaries must hold bit-close to machine
   precision (extends the existing `uniform_eqair_field_holds_on_quads`,
   `solver_u.rs:588`).
2. **Parity with the structured lane**: same Cartesian problem run through
   `advance_ugrid_viscous` vs `add_viscous` (`viscous2d.rs:402`) must agree to a
   tight tolerance, the unstructured machinery ships only if it reproduces the
   trusted structured flux on quads.
3. **Analytic laminar case**: plane Couette (linear velocity profile, exact for any
   consistent discretization) or Blasius flat plate vs `blasius.rs` skin-friction,
   mirroring the existing `tests/viscous2d.rs` gates.
Gate to ship: (1) exact, (2) L∞ agreement ≪ 1% on quads, (3) qualitative boundary-layer
profile on a triangle mesh. Plus `cargo clippy -D warnings` with
`LD_LIBRARY_PATH=/usr/lib64/openmpi/lib`, `timeout 300` wrapping.

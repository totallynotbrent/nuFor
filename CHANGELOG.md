# Changelog

## [Unreleased]

### Added
- Live case console for every run: `nufor run` boots a viewer on a
  background thread and prints the link, `nufor watch <case.toml> [port]`
  reopens one later. The console reads the run heartbeat (state, step,
  dt, wall time), renders the finished fields (density, mach, pressure,
  temperature, grayscale schlieren) plus the mesh, and carries a hover
  probe (cell indices, coordinates, local values), force and peak
  blocks, history plots, the log as selectable text, and a case card
  with the freestream and numerics. A standalone HTML report with the
  fields inlined lands beside the outputs when the run finishes. Docs
  in `docs/numerics/tools/web-ui.md`.
- `nufor-core::fay_riddell`: the stagnation-point cold-wall heating
  closed form (perfect-gas Lewis-one, Sutherland transport) with the
  post-shock edge chain, pinned by Rankine-Hugoniot table gates and
  scaling gates. The stagnation-heating anchor; docs in
  `docs/numerics/foundations/fay-riddell.md`.

Two-temperature relaxation: `physics.eos = "park2t"` carries a
live vibrational energy field on the nitrogen Euler march (harmonic
oscillator, Millikan-White plus Park cutoff relaxation, backward-Euler)
and writes the vibrational temperature to `<case>-tv.vtk`. Gates pin the
T_v inversion, the analytic Landau-Teller decay, and coupled-box energy
conservation to the fixed point. The Mach 20 Casseau cylinder anchor
shows the expected shock-layer T/T_v split.

The MPI march gained a laminar viscous lane: equilibrium-air cases with
`physics.mu` set and no turbulence model now run their two-stage Heun
viscous step inside the Fortran core across ranks (Sutherland transport
from the CEA closure temperature), instead of falling back to the serial
Rust march under `mpirun`. Bit-exact against the serial viscous march on
the body-free class at np=1 and np=2, pinned by a new parity gate wired
into CI alongside the existing Euler gates.

The equilibrium-air closure table rebuilt from a wider CEA sweep:
101 pressures from 0.01 Pa to 10.1 MPa (was 100 Pa to 1 MPa) and 30
temperatures from 200 K to 30000 K (was 200 K to 12000 K), with the
energy axis converted to sensible internal energy so it matches the
march's conserved variable. Peak pressures above the old 1 MPa ceiling
no longer clamp, and cold-region temperatures read physically instead
of inflated.

2D restart files: `[time] restart_in` seeds a run from a prior run's
state and `restart_out` writes one at the end, on both the serial and
MPI paths. The NUF2 layout stores rho, mx, my, e with the grid shape,
time, and step count; a restarted run resumes at the file's time and
continues bit-consistently (verified exact against an uninterrupted
run). Also fixes the chunked MPI march's clock, which reset every call
so time-based stop criteria never fired on body cases.

MPI-parallel Fortran march for the equilibrium-air 2D Euler path.
`mpirun -np N nufor run case.toml` splits the grid into column slabs with
two ghost columns per side; the Fortran core owns the decomposition, the
halo exchange, and the RK2 step, with a global CFL reduce for the time
step. Body cases march in single-step chunks: gather on rank 0, solid-cell
projection, broadcast back out. Three parity gates lock the Fortran march
to the serial Rust march: sod 2D and a uniform equilibrium-air field at
exact zero, and a capsule body at 22.5 degrees angle of attack at 4e-15
worst relative difference. Measured on the 34.6 km station case
(320x176, 30 steps): 1.88x at two ranks, with four ranks plateauing at
the same wall time because 80-column slabs are too thin to pay for the
exchange. See [the MPI page](docs/numerics/parallel/mpi.md).

The closure-aware viscous march landed before that: a planar 2D case with
`physics.mu` above zero steps with a viscous march that closes over
whichever EOS the case names, per-cell Sutherland viscosity against the
cell temperature, conductivity following `physics.pr`. The wall thermal
side is the new `physics.wall_temperature` key; zero keeps the adiabatic
wall, any positive value a fixed cold wall. Every case run now writes a
`run-log.txt` summary and appends a `run-log.csv` row next to the case
file: case, grid, step count, wall clock, viscous/eos flags, freestream
state, peak-pressure and peak-temperature extrema with locations, and the
body force coefficients (C_D, C_L, C_m) when the case has a body. The
equilibrium-air closure is a CEA-derived table for cell thermodynamics
plus TGAS1 curve fits at faces, and the four-station aeroshell sweep
(30.0 to 42.5 km, 0 and 22.5 degrees) runs against it with a 1280x640
grid-convergence anchor.

Laminar viscous flow met the thermodynamic closure first: a planar 2D
case with `physics.mu` above zero steps with a viscous march that closes
over whichever closure the case names: per-cell Sutherland viscosity
against the cell temperature, conductivity following `physics.pr`
scaled by the closure's cold-air gas constant. The wall thermal side
is a new `physics.wall_temperature` key (K); zero keeps the adiabatic
wall, any positive value a fixed cold wall.

Axisymmetric Euler solver and masked solid bodies. The 2D planar flux
machinery now drives an annular update for bodies of revolution: radial
faces carry the cell radii and radial momentum gains the geometric p/r
source, with the slip-wall ghost doubling as the symmetry axis
condition. Validated on the Mach-2 sphere bow shock: standoff within
1.2 percent of the Ambrosio-Wortman correlation, post-shock density in
the Rankine-Hugoniot band. The solid-cell mask generalizes to any shape
that answers inside plus normal: the circle keeps its cylinder contract,
a polygon body reproduces the cylinder bow shock, and an analytic
signed-distance sphere-cone reproduces the trade-study dimensions of the
capsule geometry exactly. The sphere-cone runs at the peak-q trajectory
point (Mach 22, perfect gas) form a detached bow shock with
post-shock density near the high-Mach limit, and the surface Cp
extraction feeds an axial-force cross-check.

Navier-Stokes in 3D. The viscous operator works on the uniform structured
3D grid, with no-slip adiabatic walls on any of the six domain faces and
a Heun march capped by the explicit diffusion bound. Two analytic flows
validate it: the Poiseuille parabola held in a 3D slab (worst profile
error about two percent of the peak) and the square-duct series solution
on all four lateral walls (worst error under one percent). A symmetry
contract drains a y-shear exactly like the z-shear, guarding the axis
path the channel case cannot see. The 2D Spalart-Allmaras march also
gained a threaded counterpart, bit-identical to the serial march at
every thread count, worth about 1.4x on fine 384x384 grids.

## [1.5.0] - 2026-09-27

The turbulent flat plate validates. The Spalart-Allmaras model carries a
genuinely turbulent boundary layer at Reynolds numbers where the flat-plate
correlations apply: `nu_tilde` grows to hundreds of molecular viscosity at
mid-chord, and the discrete skin friction sits 12 to 19 percent above the
Schlichting power law over the interior stations, the same sign and scale
as the spread between the NASA CFL3D SA reference and the correlations
themselves. A constant-by-constant audit against the published model fixed
a real defect on the way: the `fw` damping function carried a stray `+1`
in its denominator, which bent the log-layer equilibrium by a quarter
percent everywhere. The Blasius inflow now refuses to anchor at the leading
edge, where the similarity solution is singular, and the shipped plate case
moves its virtual leading edge 0.2 upstream. An experimental
`local_time_stepping` flag advances each cell at its own stability bound;
its fixed-point algebra is verified, but the iterator does not converge on
stiff wall layers, so the time-accurate march remains the validated path.

## [1.4.2] - 2026-09-23

Inflow can carry a boundary layer now. A case side set to `profile_inflow`
prescribes u(y)/v(y) along the inflow plane instead of one uniform state,
from either the built-in Blasius laminar layer (anchored at a virtual leading
edge, the integrated similarity table matching the published values) or a
tabulated `y u v` file. The inviscid ghosts, the reconstruction, and the
viscous ghost gradients all evaluate the same profile at each row's true
height, so the incoming layer is fed consistently from the first cell.

With that, the laminar flat plate validates quantitatively: the analytic
Blasius field laid into a clustered grid and fed by its own profile is held
to about two percent, and the discrete skin friction matches the
0.664/sqrt(Re_x) correlation within a few percent over the interior
stations. The flat-plate case ships with the profile inflow.

Also fixed in passing: the Cf report now differentiates at the true wall
geometry (it was using the uniform-spacing shorthand, wrong on the shipped
clustered mesh), and the SA vorticity uses true neighbor spacing on
stretched grids.

## [1.4.1] - 2026-09-23

Stretched grids solve correctly now. The 2D Euler, viscous, and SA transport
steps are metric-aware on any non-uniform rectilinear mesh: reconstruction
uses the true cell centers and faces, gradients are the exact three-point
Lagrange derivatives (the plain central difference on uniform spacing), the
wall shear uses the true wall geometry, and the divergence divides by each
cell's own width. Uniform grids keep the original fast path bit for bit, so
existing cases are untouched. The stretched Poiseuille hold balances to
machine precision and a contained bump conserves mass to roundoff on a
clustered grid.

Cases can now ask for a clustered mesh directly: `mesh.source = "clustered"`
with `first_cell`, `growth`, and `cluster = "wall" | "channel" | "top"`. The
flat-plate case ships with a wall-clustered grid, the piece the turbulent
wall treatment needs.

Imported non-uniform meshes were silently solved with uniform-spacing
formulas before this; the same fix corrects them.

## [1.4.0] - 2026-09-21

Turbulence: the solver now carries the Spalart-Allmaras model. A case with
`equations = "rans_2d_sa"` solves the 2D Navier-Stokes equations coupled to
the transported SA variable, with the eddy viscosity feeding the momentum
and heat fluxes. The model follows the 1994 baseline (no trip terms) with
the standard guards, closes at `nu_tilde = 0` exactly, and reproduces the
laminar solver bit for bit in that limit.

No-slip walls now shear the flow. The viscous pass carries a
quadratic-consistent one-sided wall flux, so a boundary layer actually
develops against a wall; the Poiseuille hold balances to machine precision
on the discrete stencil. `nufor run` prints skin friction along solid walls
after an SA case, and `cases/flat-plate/case.toml` ships as the runnable
plate example.

Also new: `[boundaries]` gained `top`/`bottom` sides for the 2D solvers, and
`[physics]` gained `mu`, `pr`, and a `[physics.turbulence]` block.

## [1.3.0] - 2026-09-14

Mesh diagnostics: `nufor mesh-check case.msh` validates an imported
unstructured mesh before solving. It reports cell/face counts, area
min/mean/max and the stretch ratio, a per-cell closure residual, flipped
interior faces, and negative or degenerate volumes, then gives a verdict.
The same summary is available at `GET /api/mesh-stats?path=...` for tooling.

`Ugrid` now stores cell centroids, which feeds the orientation and closure
checks.

## [1.2.0] - 2026-09-09

`nufor run <case.toml>` now dispatches on `physics.equations` and solves 1D
(`euler_1d`), 2D (`euler_2d`), or 3D (`euler_3d`) from a single case file.
Results (csv/vtk, plus a PNG for 2D) are written next to the case and a
summary is printed.

- `[mesh]` adds the 2D/3D edges (`ny/y0/y1`, `nz/z0/z1`). `source = "file"`
  plus `path` loads cell-center coordinates from a file. Coordinates must be
  uniformly spaced; non-uniform input is refused with a clear message.
- New `blast` initial condition, an over-pressured fireball in ambient
  surroundings, for 2D/3D. `uniform` and `two_state` stay for 1D. `nufor init`
  writes a runnable template. Validation rejects a non-positive blast radius.
- `write_vtk2d` and `write_vtk3d` emit ParaView-readable `STRUCTURED_POINTS`
  blocks (density/pressure/energy).
- The web Run button drives `/api/run-case` with the configured dimension,
  cell count, final time, and CFL, saves `results/<name>.vtk`, and shows the
  convergence summary (steps, time, wall clock) in the Monitor dock.
- The UI header version now reads from the crate instead of a hardcoded string.

## [1.1.0] - 2026-09-09

### Web UI

- Rebuilt as one workspace like ParaView/Tecplot: a central flow viewport, a
  data/display sidebar on the left, a diagnostics dock below. No per-feature
  tabs.
- Viewport shows the solved field with a toggleable cartesian mesh overlay, a
  pinned colorbar, and play/scrub over the blast time.
- 3D mode: slice planes through a solved 3D sphere-blast, with the axis and
  slice position selectable.
- The sidebar holds simulation controls (case, cells, t, gamma, cfl, run) and
  display controls (dimension, field, layers, mesh density, playback). The
  dock holds the 1D monitor, line probes, the comparison overlay, and run
  history.
- Served on `0.0.0.0` so it is reachable over LAN or Tailscale.

### Running the web UI

```
cargo run --release -- serve 8060
```

then open `http://<host>:8060`. Use the release binary; the flow solves are
orders of magnitude faster than the debug build.

## [1.0.0] - 2026-09-05

First stable public release. The solver covers the 1D, 2D, and 3D Euler
equations and the 2D Navier-Stokes terms in a Rust application shell with a
CPU-first design. It ships with a Quartz documentation site, a web dashboard,
and a verification suite.

### Numerics

- 1D finite-volume Euler solver: HLL and HLLC fluxes, CFL time-step control,
  MUSCL reconstruction with the van Leer limiter, verified against the exact
  Riemann solution (Sod, Lax, and the four-step ladder).
- 2D Euler solver with HLLC flux. Oblique-shock and shock-reflection cases
  verified against shock relations. Positivity-preserving update.
- 2D viscous Navier-Stokes terms, Poiseuille-channel validation.
- 3D structured solver with a 3D HLLC step, verified against the Sod tube, and
  a documented memory budget against the 8 GB target.
- Unstructured cell-centered prototype, verified on a quad mesh. A
  Spalart-Allmaras research note scopes the first turbulence model.
- Supersonic-cylinder immersed body verified against a detached bow shock
  reproducing the normal-shock density ratio.
- Runtime-dispatched AVX2 kernels in the conservation update, with a scalar
  fallback and numerical-equivalence tests.
- Restart format with an explicit schema version, size and physical-sanity
  validation, and legacy reads.

### Interfaces

- Web UI with live 1D plots, 2D field images, line probes, and a comparison
  dashboard.
- CSV, VTK, and HDF5 output plus a bit-exact binary restart.
- CLI: `run`, `inspect`, `export`, `benchmark`, `serve`, `version`.
- Quartz documentation site with a page graph, mermaid, and LaTeX, published
  on GitHub Pages.

### Verification

Cross-platform CI (Linux and Windows) runs the Rust build, the Fortran-linked
FFI integration tests, and a skeleton/LF gate on every change. The full
workspace suite is reproducible from a clean checkout.
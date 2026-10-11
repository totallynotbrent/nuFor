# MPI port plan — audit, decisions, checklist

Decision locked with the user: **Fortran owns the parallel solver core, Rust stays
the front door.** Do NOT base any of this on what exists already (user directive) —
the 1D Fortran path and the Rust threading are precedent for style only. The MPI
march is new code.

## Status (2026-10-11, C1+C2 landed)

G1 parity gates all green: mpi-vs-rust worst rel diff 0.000e0 on sod 2d
AND uniform eqair, at np=1/2/4/8, one-call and chunked (NF_MPI_CHUNKED)
modes; workspace suite + clippy -D warnings clean; CLI E2E serial vs
mpirun np=2/4 vtk outputs bit-identical. Two hard-won fixes of note:
(1) bind(c) scalars MUST carry `value` — nfor_set_tgas1's r0/e0 lacked
it and the fortran dereferenced the f64 bits (SEGV in the CLI, silent
garbage in the tests); (2) the ein y-sweep rework had accidentally
deleted the x-sweep's HLLC pass — freshly allocated fx read as
uninitialized heap, which only surfaced in the chunked CLI loop.
Known box limit: the user systemd slice caps at 11 GiB (user.slice
memory.max) and kincsem sits pinned near it — CLI runs at np>=5 can
get ranks killed by cgroup reclaim during output. Library-level parity
at np=8 is clean; either raise the slice cap or keep CLI runs at
np<=4 on this box.

Two mpi-exit rules baked in after the CI round: every rank must call
MPI_Finalize before process exit or mpirun reports abnormal
termination (ffi_mpi::mpi_shutdown, idempotent, on every CLI exit
path), and the parity tests run their whole mpi session on a spawned
thread — libtest's worker thread + finalize makes the launcher exit 1
even when every rank passes.

## Architecture (the locked shape)

- ONE binary (`nufor`), launched with `mpirun -np N nufor ...`.
- Every rank runs the Rust main. Rust never calls MPI.
- Rank 0: loads case.toml, drives the run, writes VTK/run-log/restart.
- All ranks: call the Fortran core. Fortran (mpif90, `use mpi`) does
  decomposition, halo exchange, and the march.
- Boundary = coarse FFI: state+grid+bc in per phase, snapshots out. No per-step
  Rust<->Fortran chatter.
- n_ranks = 1 runs the SAME MPI code path (Cartesian comm of one rank) and must
  reproduce the serial results BIT-IDENTICALLY. That is the regression anchor
  for every later change.

## Full audit findings (what the port must handle)

- March to parallelize first: `advance2d_model_rk2` (euler_2d, the capsule path),
  called at crates/nufor-cli/src/main.rs:1421/1434. RK2 Heun: full state clone,
  two `advance2d_model` sweeps, global CFL from a max-reduce.
- State layout: ConservedState2d = rho/mx/my/e, flat row-major (j*nx+i) Vecs.
  Grid2d = uniform or per-axis dxs/dys, centers/faces arrays. SAME flat layout
  ships to Fortran; no repacking layer.
- Flux: hllc2d::hllc_flux per face, per-side gamma from (a,p,rho) with floors.
  MUSCL reconstruction on primitive strips, `face_states` builds padded strips.
- BCs (solver2d.rs:146 Bc2d enum): Transmissive, SupersonicInflow,
  SupersonicOutflow, SlipWall, NoSlipWall, ProfileInflow. Capsule cases use
  supersonic_inflow (W+S) / supersonic_outflow (E+N). BCs apply at DOMAIN
  edges only; subdomain edges become halo faces — the Fortran side needs its
  own bc enum mirroring these five.
- Body: SolidPolygon inside/normal via body.rs; solid cells pinned, wall-face
  fluxes zeroed via cached mask (OnceLock cache). In MPI the classification is
  geometry-static: compute ONCE on rank 0 before decomposition, ship the
  per-rank solid/band lists + mask with the decomposition, never recompute
  per step per rank.
- EqAir: CEA table (eqair_cea.rs, OnceLock 40x60) + TGAS1 face fits
  (eqair.rs gamm_and_partials, the load-bearing wake path). Both must exist
  Fortran-side as whole-array kernels; the TGAS1 block dispatch is
  mechanical (12 blocks, grabau sigmoids).
- Run log: per-run csv fields (steps/time/wall/eos/rho_ratio/p_peak/t_max/q_w
  peak/forces C_D C_L C_m) written by the CLI from the final state; force
  integration needs p along the body polygon — rank 0 gathers the wall-face
  band cells after the march.
- Restart: RESTART_VERSION 1, plain binary dump of the 4 arrays. Serial read
  on rank 0, then scatter; gather+write on rank 0. No format change needed.
- Output: VTK snapshots every interval_steps from rank 0 after a gather;
  unchanged file format (gathered full-domain arrays).
- Build: CMake already drives gfortran via build.rs; CI installs gfortran.
  ADD: openmpi-devel + mpif90 in CMake (find_package(MPI COMPONENTS Fortran)),
  openmpi on kincsem (dnf install), apt openmpi in CI. CI runners have
  openmpi available via apt; guard the MPI tests with `mpirun --version`.
- kincsem: 8 cores, no openmpi installed yet. `sudo dnf install openmpi`
  (EL10: module load mpi or /usr/lib64/openmpi/bin on PATH).
- Web UI (:8060 serve path) untouched: it launches the CLI binary; when the
  CLI is mpirun-launched, pass -np from cfg.numerics (new `mpi_ranks` key).

## What moves to Fortran (per the audit)

Phase F1 (parallel euler march): decomposition + halo exchange + one RK2
march for euler_2d, eqair closure whole-array kernels (p_t_at bilinear,
gamm_and_partials TGAS1 for faces), HLLC flux array kernel, MUSCL strip
reconstruction, BC application, body classification shipping, global CFL
(MPI_Allreduce), run-log gather.
Phase F2 (viscous+SA): viscous flux kernels (x+y faces), Sutherland mu,
SA source kernel, the viscous dt cap per cell.
Phase F3 (new physics MPI-first): rung 5 live 2T relaxation kernel,
rung 8 species continuity + VODE-class stiff integrator binding.

## Phase gates (ALL must pass before a phase is "done")

G0 environment: mpirun runs `nufor --version` multi-rank on kincsem + CI.
G1 (F1 correctness): sod 2D + cylinder case, np=1 bit-identical to the
current serial march (all 5 VTK fields, diff == 0.0); then np=2,4,8 with
max rel diff < 1e-12 vs np=1, step counts identical, run-log fields match
to 1e-12.
G2 (capsule anchor): aeroshell-h42.5-a22.5 np=8 runs to t=0.015, L/D and
peak-T match the serial run inside 1e-10; wall-time recorded for the
strong-scaling table.
G3 (strong scaling): np = 1,2,4,8 on the capsule case (320x176), wall times
+ CPU% table in docs/numerics/parallel/scaling.md; honest numbers under the
35 W cap, whatever they are.
G4 (F2): viscous capsule case (aeroshell-capsule-trim) np=8 vs np=1 same gate.
G5 (every later phase): same G1 discipline on the new physics + the phase's
analytic gates (Landau-Teller box for 2T, mass conservation for species).

## File changes (the full checklist)

CHECKLIST (in order; ~2-3 sessions of work on this box):

C1  build surface  [DONE 2026-10-10/11 — see status below]
    [x] fortran/api/nufor_mpi.f90 — new module: decomposition (1D slabs along
        x for v1), halo pack/send/recv (MPI_Isend/Irecv, 4 ghost rows),
        global reductions. mirror bc enum in Fortran.
    [x] fortran/api/nuforsolver2d.f90 — the march: RK2 heun, MUSCL strips,
        HLLC array kernel, eqair whole-array kernels (bilinear table +
        TGAS1 gamm), body mask apply, CFL Allreduce. one `march()` entry
        taking flat arrays + rank info, all state in module variables
        between calls (no SAVE-derived surprises; explicit, reentrant per rank).
    [x] CMakeLists.txt — find_package(MPI COMPONENTS Fortran REQUIRED);
        compile with mpif90 wrappers.
    [x] crates/nufor-core/build.rs — link MPI libs (emit cargo:rustc-link-lib
        for mpi_mpifh or use pkg-config openmpi).
    [x] crates/nufor-core/src/ffi_mpi.rs — Rust side: case-to-flat arrays,
        per-phase calls (init/march-steps/finalize), snapshot gather.
    [x] crates/nufor-config/src/schema.rs — `mpi_ranks: usize` under
        [numerics] (default 1), validation.
    [ ] .github/workflows/ci.yml — apt install openmpi; job step runs the
        np=1 gate + a smoke np=2 sod gate. (STILL OPEN)
    [x] kincsem: sudo dnf install openmpi; PATH fix in the run scripts.

C2  rust-side plumbing  [DONE 2026-10-11 — see status below]
    [x] cli dispatch: when launched under mpirun (OMPI_COMM_WORLD_SIZE > 1)
        OR env NFOR_MPI=1, route eqair inviscid euler_2d marches through
        ffi_mpi (march_euler2d_mpi). np=1 still goes through MPI (the
        bit-identity anchor). all other case types unchanged this phase.
        body cases march in single-step chunks: gather on rank 0, band
        projection (apply_solid_fn), then one packed MPI_Bcast of the 4n
        state buffer so every rank re-carves identical slabs.
    [x] rank-0-only output: the mpi branch returns before output on
        non-root ranks; the serial fallback under mpirun also skips
        non-root ranks (they would clobber the output files).
    [ ] restart read on rank 0 -> scatter; gather -> write. (STILL OPEN)

C3  gates + docs (G0..G3 above)
    [ ] tests/mpi_parity.rs (rust): np=1 sod/cylinder bit-identity vs serial.
    [ ] scratch script: strong-scaling run matrix + table generator.
    [ ] docs/numerics/parallel/mpi-architecture.md — the locked shape above,
        the decision rationale (user directive), halo/ghost discipline.
    [ ] docs/numerics/parallel/scaling.md — measured table.
    [ ] nuFor_README.md sections 71-77 year-4 plan: MPI rung marked in
        progress/done as it lands; roadmap.md gets the MPI rung inserted
        after the threading note (as a parallelism rung, not physics).
    [ ] README.md summary line: "MPI-parallel Fortran core, Rust front end".

C4  F2 viscous+SA port (separate phase, same gates G4)
    [ ] nuforsolver2d.f90 gains the viscous kernels (x+y faces, Sutherland,
        per-cell dt cap, SA source).
    [ ] cli: viscous_2d + sa_2d case types route through ffi_mpi.
    [ ] viscous parity test + trim-case gate.

C5  F3 physics-first (rung 5 and 8 land here, MPI-native)
    [ ] 2T relaxation kernel (Landau-Teller box gate).
    [ ] species continuity (mass conservation gate) + stiff integrator.

## Verification discipline (unchanged doctrine, new thresholds)

- np=1 bit-identity is the anchor; np>1 rel-diff < 1e-12 is the MPI gate.
- analytic gates never skip: sod exact, cylinder standoff, capsule
  Rankine-Hugoniot, modified-Newtonian Cp.
- no cargo test while a long march is live; 420s terminal timeout ->
  background scripts for the scaling matrix; checkpoint the ledger.
- every gate lands as a file on disk (run-log, scaling.md table), not a
  claim in chat.

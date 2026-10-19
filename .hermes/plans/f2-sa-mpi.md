# F2 completion: SA turbulence through the MPI Fortran core

## scope

the rans_2d_sa eqair case class currently runs serial-only (the cli mpi
routing gate excludes turb). the goal: run it through the fortran core
under mpirun with an np=1 bit-parity gate against the serial rust march
`advance2d_sa_model_rk2`, then np=2 drift checks, then the ci smoke.

## design

1. fortran `sa_transport_slab` (nuforsolver2d.f90): advances nu_tilde on
   the slab. exact mirror of turb2d.rs advance_turb on uniform axes:
   - vorticity s = |dv/dx - du/dy|, centered differences, pad2d ghosts
   - first-order upwind advection on face-averaged (u,v,rho), conservative
     in rho*nu_tilde (divide the net flux by rho)
   - central diffusion with (nu_lam + nu_tilde_f)/sigma, plus the
     c_b2/sigma (|grad nu_tilde|^2) cross term on centered grads
   - the 1994 source (production/destruction closures from sa.rs, same
     constants), wall distance d shipped per slab from rust
   - wall ghosts: nu_tilde = 0 at no-slip/slip walls (pad_scalar rule),
     open sides copy the interior
   - positivity floor at the end
   the halo: nu_tilde needs a 2-column halo exchange like the mean flow
   (add halo_exchange_nt or fold into the existing exchange as a 5th array;
   separate routine keeps the laminar path untouched).

2. fortran `viscous_add_slab` gains an optional eddy-viscosity pair
   (mu_t field + kappa_t scalar) mirroring face_flux_split:
   - muf = mu_lam_f + 0.5(mu_t_a + mu_t_b) for the stresses
   - q = (mu_lam_f*kappa + mu_t_f*kappa_t)*grad(T), kappa_t =
     gamma/((gamma-1)*pr_t) (r-free, matching the rust lane exactly,
     warts and all - parity over beauty)
   pass mu_t as an optional array argument; laminar callers pass null.

3. bridge entry `nfor_mpi_step_sa` (nuforbridge.f90): the coupled heun
   exactly as advance2d_sa_model_rk2:
   - stage: euler2d_step_rk2 stage-1 half (the same per-stage call the
     visc path uses), halo, viscous_add_slab with mu_t from the current
     nu_tilde, halo, sa_transport_slab on the post-viscous state
   - stage 2 likewise from the stage-1 state
   - heun average of state and nu_tilde
   - dt cap: dt_sa = 0.25*min(dx,dy)^2 * sigma/(mu/rho_min + nt_max)
     computed on the slab, min-reduced across ranks (allreduce min, same
     pattern as the cfl reduce)
   - returns dt_out, the updated slab state + slab nu_tilde

4. rust side (ffi_mpi.rs): `march_rans2d_mpi` wrapper:
   - ships nu_tilde0 + d per slab (d from wall_distance2d, slab-carved)
   - per call: one step (the cli chunked pattern), returns the gathered
     nu_tilde on rank 0 alongside the state
   - gathers nu_tilde via nfor_mpi_gather pattern (a scalar variant or a
     5th array through the existing gather)

5. cli routing: relax the mpi gate to admit turb cases when
   mpi_launched() && eqair && !visc_body (turb body cases stay serial for
   now - the masked-body re-projection lane is a separate chunked path
   and the sa wall distance for bodies is not wall_distance2d but the
   body distance field; keep the existing serial fallback for those).

6. gate: tests/mpi_parity.rs `mpi_np1_sa_matches_rust`:
   - body-free rans case (no-slip south wall so sa has a real wall):
     uniform eqair freestream, mu=1.716e-5, nu_tilde_inf=3*nu
   - n steps serial rust advance2d_sa_model_rk2 vs march_rans2d_mpi np=1
   - assert rel diff < 1e-6 on rho, mx, my, e, AND nu_tilde
   - ci: add `check 1 mpi_np1_sa_matches_rust` + `check 2 ...` lines

## order of work

fortran sa kernel -> parity smoke locally (np=1) -> bridge entry -> rust
wrapper -> parity test green -> np=2 drift -> cli routing -> ci lines ->
docs (mpi.md F2 section rewrite) -> commit.

## risks

- the rust advance_turb uses padded arrays with true-position lagrange
  grads only on stretched axes; the mpi core is uniform-dx (dx, dy
  scalars), so the uniform branch is the only one to port. verified: all
  mpi-routed grids are uniform (the bridge takes scalar dx, dy).
- bit-parity on the source term needs the exact fp order of sa.rs
  (stilde closure chain, r clip, g, fw powf). f90 and rust both do
  ieee754 doubles; write the f90 in the same operation order.
- the heun dt cap in rust uses nu_tilde_inf as the nt_max floor; mirror
  that (max over field vs the inf value).

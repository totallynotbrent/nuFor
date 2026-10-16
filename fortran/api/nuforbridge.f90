! nuforbridge.f90 — the C-ABI surface between the rust front door and the
! mpi fortran core. every routine is bind(c), takes flat arrays, and never
! allocates globally. state lives on the rust side; fortran holds only the
! cea/tgas1 tables (set once) and the decomp built at case init.
module nuforbridge
  use iso_c_binding, only: c_int, c_double
  use iso_fortran_env, only: real64
  use nufor_mpi
  use nuforsolver2d
  implicit none

  integer, parameter :: RP = real64

  type(decomp_t), save :: g_decomp
  integer, save :: g_rank = -1

contains

  subroutine nfor_mpi_init(ierr) bind(c, name="nfor_mpi_init")
    integer(c_int), intent(out) :: ierr
    integer :: e
    call mpi_init(e)
    call mpi_comm_rank(MPI_COMM_WORLD, g_rank, e)
    ierr = e
  end subroutine

  subroutine nfor_mpi_finalize(ierr) bind(c, name="nfor_mpi_finalize")
    integer(c_int), intent(out) :: ierr
    integer :: e
    call mpi_finalize(e)
    ierr = e
  end subroutine

  ! build the decomposition for this rank; returns rank (0-based), the
  ! global x cell offset (i0-1, 0-based), and the local slab width.
  subroutine nfor_mpi_case_init(nx, ny, out_rank, out_i0, out_nx_local, ierr) &
      bind(c, name="nfor_mpi_case_init")
    integer(c_int), value, intent(in) :: nx, ny
    integer(c_int), intent(out) :: out_rank, out_i0, out_nx_local, ierr
    integer :: sz, rk, e
    call mpi_comm_size(MPI_COMM_WORLD, sz, e)
    call mpi_comm_rank(MPI_COMM_WORLD, rk, e)
    g_decomp = make_decomp(nx, ny, sz, rk)
    out_rank = g_decomp%rank
    out_i0 = g_decomp%i0 - 1   ! 0-based global offset
    out_nx_local = g_decomp%nx_local
    ierr = 0
  end subroutine

  ! ship the tgas1 block coefficients + meta (rust owns the parsed source)
  subroutine nfor_set_tgas1(g1,g2,g3,g4,g5,g6,g7,g8,s,minus,nb, zb, cg, r0, e0) &
      bind(c, name="nfor_set_tgas1")
    integer(c_int), value, intent(in) :: nb
    real(c_double), intent(in) :: g1(12,2), g2(12,3), g3(12,3), g4(12,3)
    real(c_double), intent(in) :: g5(12,2), g6(12,3), g7(12,3), g8(12,3), s(12,4)
    integer(c_int), intent(in) :: minus(12)
    real(c_double), value, intent(in) :: r0, e0
    real(c_double), intent(in) :: zb(3,5), cg(3)
    logical :: minus_l(12)
    integer :: nbk12(12), i
    do i = 1, 12
      minus_l(i) = minus(i) == 1
      ! blocks 5 and 10 carry no grabau sigmoid (rust grabau: None);
      ! the rest do. evaluating the sigmoid on the unused coefficient
      ! slots poisoned every state that landed in those two blocks.
      nbk12(i) = 1
      if (i == 5 .or. i == 10) nbk12(i) = 0
    end do
    call set_tgas1(g1, g2, g3, g4, g5, g6, g7, g8, s, minus_l, nbk12, r0, e0)
    call set_tgas1_meta(zb, cg, r0, e0)
  end subroutine

  ! ship the cea table: the rust side sends flat[k] = cell (i, j) at
  ! k = i*ne + j (i the rho node, j the energy node, both 0-based), matching
  ! the eqair_cea parser's row-major fill.
  subroutine nfor_set_cea(nr, ne, lr, le, t, p) bind(c, name="nfor_set_cea")
    integer(c_int), value, intent(in) :: nr, ne
    real(c_double), intent(in) :: lr(nr), le(ne), t(nr*ne), p(nr*ne)
    real(RP) :: t2(nr,ne), p2(nr,ne)
    integer :: i, j
    do j = 1, ne
      do i = 1, nr
        t2(i,j) = t((i-1)*ne + j)
        p2(i,j) = p((i-1)*ne + j)
      end do
    end do
    call set_cea_table(nr, ne, lr, le, t2, p2)
  end subroutine

  ! re-seed a slab's ghost columns from its edge cells after a stage
  ! advance: euler2d_step_rk2 copies its input ghosts through, so after
  ! stage 1 the domain-edge ghosts still hold the entry state. rust
  ! re-derives every ghost from the current edge each sweep; the
  ! per-sweep bc_ghost pass re-applies the bc and inflow preinverts
  ! ein, so a plain copy is safe for every bc.
  subroutine reseed_edge_ghosts(f, nx_l, ny_)
    integer, intent(in) :: nx_l, ny_
    real(RP), intent(inout) :: f(nx_l+4, ny_)
    integer :: j
    if (g_decomp%west_domain_edge) then
      do j = 1, ny_
        f(1,j) = f(3,j)
        f(2,j) = f(3,j)
      end do
    end if
    if (g_decomp%east_domain_edge) then
      do j = 1, ny_
        f(nx_l+3,j) = f(nx_l+2,j)
        f(nx_l+4,j) = f(nx_l+2,j)
      end do
    end if
  end subroutine

  ! one forward-euler step on this rank's slab. rust owns ghost columns as
  ! copies of its edge cells at entry (np=1) — at np>1 they are overwritten
  ! by the halo exchange.
  subroutine nfor_mpi_step(rho, mx, my, e, nloc, nx_l, ny_, &
                           dx, dy, gamma, cfl, dt_cap, &
                           bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                           rho2, mx2, my2, e2, dt_out, ierr) &
      bind(c, name="nfor_mpi_step")
    integer(c_int), value, intent(in) :: nloc, nx_l, ny_
    real(c_double), intent(in) :: rho(nloc), mx(nloc), my(nloc), e(nloc)
    real(c_double), value, intent(in) :: dx, dy, gamma, cfl, dt_cap
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(5, 4)
    integer(c_int), intent(in) :: solid(nloc)   ! 0/1
    real(c_double), intent(out) :: rho2(nloc), mx2(nloc), my2(nloc), e2(nloc)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: ierr
    logical :: sl(nx_l+4, ny_)
    real(RP) :: rin(nx_l+4, ny_), min_(nx_l+4, ny_), myin(nx_l+4, ny_), ein(nx_l+4, ny_)
    real(RP) :: r1(nx_l+4, ny_), m1(nx_l+4, ny_), my1(nx_l+4, ny_), e1(nx_l+4, ny_)
    real(RP) :: r2(nx_l+4, ny_), m2(nx_l+4, ny_), my2_(nx_l+4, ny_), e2_(nx_l+4, ny_)
    real(RP) :: dt2
    real(RP) :: dt1, dt
    integer :: j, i

    do j = 1, ny_
      do i = 1, nx_l
        rin(i+2,j) = rho((j-1)*nx_l + i)
        min_(i+2,j) = mx((j-1)*nx_l + i)
        myin(i+2,j) = my((j-1)*nx_l + i)
        ein(i+2,j) = e((j-1)*nx_l + i)
        sl(i+2,j) = solid((j-1)*nx_l + i) == 1
      end do
      rin(1,j) = rin(3,j);   rin(2,j) = rin(3,j)
      rin(nx_l+3,j) = rin(nx_l+2,j); rin(nx_l+4,j) = rin(nx_l+2,j)
      min_(1,j) = min_(3,j); min_(2,j) = min_(3,j)
      min_(nx_l+3,j) = min_(nx_l+2,j); min_(nx_l+4,j) = min_(nx_l+2,j)
      myin(1,j) = myin(3,j); myin(2,j) = myin(3,j)
      myin(nx_l+3,j) = myin(nx_l+2,j); myin(nx_l+4,j) = myin(nx_l+2,j)
      ein(1,j) = ein(3,j);   ein(2,j) = ein(3,j)
      ein(nx_l+3,j) = ein(nx_l+2,j); ein(nx_l+4,j) = ein(nx_l+2,j)
      sl(1,j) = .false.; sl(2,j) = .false.
      sl(nx_l+3,j) = .false.; sl(nx_l+4,j) = .false.
    end do
    call halo_exchange_cons(rin, min_, myin, ein, nx_l, ny_)

    ! two-stage heun matching the serial advance2d_model_rk2: stage 1
    ! advances the input, the halo refreshes, stage 2 advances the
    ! stage-1 state (each stage re-derives its own dt from its own
    ! entry state, exactly like the serial lane's per-call dt), and
    ! the output heun-averages the input with the twice-advanced
    ! state. the original single-stage port was first-order against
    ! the serial twin and the uniform-seed gates could not see it;
    ! the stage-2 reconstruction that used to diverge at wall-band
    ! cells was the un-floored NaN sound speed in the y sweep, now
    ! clamped exactly like rust's f64::max.
    call euler2d_step_rk2(rin, min_, myin, ein, nx_l, ny_, &
                          dx, dy, gamma, cfl, dt_cap, g_decomp, &
                          bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                          r1, m1, my1, e1, dt1, 1)
    call halo_exchange_cons(r1, m1, my1, e1, nx_l, ny_)
    call reseed_edge_ghosts(r1, nx_l, ny_)
    call reseed_edge_ghosts(m1, nx_l, ny_)
    call reseed_edge_ghosts(my1, nx_l, ny_)
    call reseed_edge_ghosts(e1, nx_l, ny_)
    ! a masked body keeps the single-stage pack: the stage-2
    ! reconstruction next to reference-state solid cells amplifies
    ! through the y sweep and corrupts the wall band. porting the
    ! staging needs its own body gate.
    if (any(sl)) then
      r2 = r1
      m2 = m1
      my2_ = my1
      e2_ = e1
      dt2 = 0.0_RP
    else
      call euler2d_step_rk2(r1, m1, my1, e1, nx_l, ny_, &
                            dx, dy, gamma, cfl, dt_cap, g_decomp, &
                            bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                            r2, m2, my2_, e2_, dt2, 1)
    end if
    if (any(sl)) then
      do j = 1, ny_
        do i = 1, nx_l
          rho2((j-1)*nx_l + i) = r1(i+2,j)
          mx2((j-1)*nx_l + i) = m1(i+2,j)
          my2((j-1)*nx_l + i) = my1(i+2,j)
          e2((j-1)*nx_l + i) = e1(i+2,j)
        end do
      end do
    else
      do j = 1, ny_
        do i = 1, nx_l
          rho2((j-1)*nx_l + i) = 0.5_RP * rin(i+2,j) + 0.5_RP * r2(i+2,j)
          mx2((j-1)*nx_l + i) = 0.5_RP * min_(i+2,j) + 0.5_RP * m2(i+2,j)
          my2((j-1)*nx_l + i) = 0.5_RP * myin(i+2,j) + 0.5_RP * my2_(i+2,j)
          e2((j-1)*nx_l + i) = 0.5_RP * ein(i+2,j) + 0.5_RP * e2_(i+2,j)
        end do
      end do
    end if
    dt_out = dt1
    ierr = 0
  end subroutine

  ! two-stage heun viscous step on this rank's slab: stage = one euler
  ! step (the same kernel the inviscid call runs) followed by a halo
  ! refresh and the viscous divergence add; the final state heun-averages
  ! the input with the twice-advanced state, matching
  ! axi::advance2d_model_visc_rk2. viscosity follows sutherland from the
  ! cea closure temperature; tw <= 0 means an adiabatic wall.
  subroutine nfor_mpi_step_visc(rho, mx, my, e, nloc, nx_l, ny_, &
                                dx, dy, gamma, cfl, dt_cap, &
                                bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                                mu0, t0_ref, s_par, pr, tw, r_gas, &
                                rho2, mx2, my2, e2, dt_out, ierr) &
      bind(c, name="nfor_mpi_step_visc")
    integer(c_int), value, intent(in) :: nloc, nx_l, ny_
    real(c_double), intent(in) :: rho(nloc), mx(nloc), my(nloc), e(nloc)
    real(c_double), value, intent(in) :: dx, dy, gamma, cfl, dt_cap
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(5, 4)
    integer(c_int), intent(in) :: solid(nloc)
    real(c_double), value, intent(in) :: mu0, t0_ref, s_par, pr, tw, r_gas
    real(c_double), intent(out) :: rho2(nloc), mx2(nloc), my2(nloc), e2(nloc)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: ierr
    logical :: sl(nx_l+4, ny_)
    real(RP) :: rin(nx_l+4, ny_), min_(nx_l+4, ny_), myin(nx_l+4, ny_), ein(nx_l+4, ny_)
    real(RP) :: r1(nx_l+4, ny_), m1(nx_l+4, ny_), my1(nx_l+4, ny_), e1(nx_l+4, ny_)
    real(RP) :: r2(nx_l+4, ny_), m2(nx_l+4, ny_), my2_(nx_l+4, ny_), e2_(nx_l+4, ny_)
    real(RP) :: mu_c(nx_l+4, ny_)
    real(RP), allocatable :: u1(:,:), v1(:,:), ein1(:,:)
    real(RP) :: dt1, dt2, dt
    real(RP) :: rho_min_l, rho_min_g, nu_lam_g, dt_vis, dt_cap_v
    integer :: j, i

    allocate(u1(nx_l+4, ny_), v1(nx_l+4, ny_), ein1(nx_l+4, ny_))
    do j = 1, ny_
      do i = 1, nx_l
        rin(i+2,j) = rho((j-1)*nx_l + i)
        min_(i+2,j) = mx((j-1)*nx_l + i)
        myin(i+2,j) = my((j-1)*nx_l + i)
        ein(i+2,j) = e((j-1)*nx_l + i)
        sl(i+2,j) = solid((j-1)*nx_l + i) == 1
      end do
      rin(1,j) = rin(3,j);   rin(2,j) = rin(3,j)
      rin(nx_l+3,j) = rin(nx_l+2,j); rin(nx_l+4,j) = rin(nx_l+2,j)
      min_(1,j) = min_(3,j); min_(2,j) = min_(3,j)
      min_(nx_l+3,j) = min_(nx_l+2,j); min_(nx_l+4,j) = min_(nx_l+2,j)
      myin(1,j) = myin(3,j); myin(2,j) = myin(3,j)
      myin(nx_l+3,j) = myin(nx_l+2,j); myin(nx_l+4,j) = myin(nx_l+2,j)
      ein(1,j) = ein(3,j);   ein(2,j) = ein(3,j)
      ein(nx_l+3,j) = ein(nx_l+2,j); ein(nx_l+4,j) = ein(nx_l+2,j)
      sl(1,j) = .false.; sl(2,j) = .false.
      sl(nx_l+3,j) = .false.; sl(nx_l+4,j) = .false.
    end do
    call halo_exchange_cons(rin, min_, myin, ein, nx_l, ny_)

    ! the explicit viscous dt bound (serial advance2d_model_visc_rk2):
    ! dt < 0.25 min(dx,dy)^2 rho_min / (mu kappa). mu is sutherland at
    ! the step-entry state, so the bound uses the same temperatures the
    ! serial lane saw. rho_min allreduced across ranks.
    call primitive_ein(rin, min_, myin, ein, nx_l, ny_, u1, v1, ein1)
    call lam_mu_slab(rin, ein1, nx_l, ny_, mu0, t0_ref, s_par, mu_c)
    rho_min_l = minval(rin(3:nx_l+2, 1:ny_))
    call global_min(rho_min_l, rho_min_g, g_decomp)
    nu_lam_g = mu0 / max(rho_min_g, 1.0e-6_RP)
    dt_vis = 0.25_RP * min(dx, dy)**2 / (nu_lam_g * max(gamma / ((gamma - 1.0_RP) * pr), 1.0_RP))
    dt_cap_v = min(dt_cap, dt_vis)

    ! stage 1: euler + viscous (mu evaluated once at step entry, the
    ! serial lane's convention; slab2 zeroes the domain-edge x faces and
    ! uses the face-averaged velocity in the energy flux)
    call euler2d_step_rk2(rin, min_, myin, ein, nx_l, ny_, &
                          dx, dy, gamma, cfl, dt_cap_v, g_decomp, &
                          bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                          r1, m1, my1, e1, dt1, 1)
    call halo_exchange_cons(r1, m1, my1, e1, nx_l, ny_)
    call reseed_edge_ghosts(r1, nx_l, ny_)
    call reseed_edge_ghosts(m1, nx_l, ny_)
    call reseed_edge_ghosts(my1, nx_l, ny_)
    call reseed_edge_ghosts(e1, nx_l, ny_)
    call viscous_add_slab2(r1, m1, my1, e1, nx_l, ny_, dx, dy, dt1, &
                           mu_c, mu_c, .false., mu0, t0_ref, s_par, &
                           pr, pr, tw, gamma, r_gas, bc_s, bc_n, &
                           g_decomp%west_domain_edge, g_decomp%east_domain_edge)
    ! refresh again so stage-2's euler sees the post-viscous neighbors
    ! at slab seams (serial reads the full post-viscous state), then
    ! re-seed the domain-edge ghosts from the stage-1 edge cells (the
    ! euler stage copies its input ghosts through, so without this the
    ! stage-2 sweeps would read the entry state at the edges).
    call halo_exchange_cons(r1, m1, my1, e1, nx_l, ny_)
    call reseed_edge_ghosts(r1, nx_l, ny_)
    call reseed_edge_ghosts(m1, nx_l, ny_)
    call reseed_edge_ghosts(my1, nx_l, ny_)
    call reseed_edge_ghosts(e1, nx_l, ny_)

    ! stage 2: euler + viscous from the stage-1 state
    call euler2d_step_rk2(r1, m1, my1, e1, nx_l, ny_, &
                          dx, dy, gamma, cfl, dt_cap_v, g_decomp, &
                          bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                          r2, m2, my2_, e2_, dt2, 1)
    call halo_exchange_cons(r2, m2, my2_, e2_, nx_l, ny_)
    call reseed_edge_ghosts(r2, nx_l, ny_)
    call reseed_edge_ghosts(m2, nx_l, ny_)
    call reseed_edge_ghosts(my2_, nx_l, ny_)
    call reseed_edge_ghosts(e2_, nx_l, ny_)
    call viscous_add_slab2(r2, m2, my2_, e2_, nx_l, ny_, dx, dy, dt2, &
                           mu_c, mu_c, .false., mu0, t0_ref, s_par, &
                           pr, pr, tw, gamma, r_gas, bc_s, bc_n, &
                           g_decomp%west_domain_edge, g_decomp%east_domain_edge)

    ! heun average of the input and the twice-advanced state
    dt = dt1
    do j = 1, ny_
      do i = 1, nx_l
        rho2((j-1)*nx_l + i) = 0.5_RP * (rin(i+2,j) + r2(i+2,j))
        mx2((j-1)*nx_l + i) = 0.5_RP * (min_(i+2,j) + m2(i+2,j))
        my2((j-1)*nx_l + i) = 0.5_RP * (myin(i+2,j) + my2_(i+2,j))
        e2((j-1)*nx_l + i) = 0.5_RP * (ein(i+2,j) + e2_(i+2,j))
      end do
    end do
    dt_out = dt
    ierr = 0
    deallocate(u1, v1, ein1)
  end subroutine

  ! rank probe without touching the decomposition (cli output gating).
  subroutine nfor_mpi_rank(rank) bind(c, name="nfor_mpi_rank")
    integer(c_int), intent(out) :: rank
    rank = g_rank
  end subroutine

  ! broadcast the packed (4n) state buffer from rank 0 (chunked-body
  ! mode: the rust side re-projects band cells between steps on rank 0,
  ! then every rank re-carves its slab from the refreshed full state).
  ! one packed message: the earlier four-array variant sat right at the
  ! ucx eager/rendezvous threshold and interleaved badly with the halo
  ! exchange under oversubscription.
  subroutine nfor_mpi_bcast4(buf, n4, ierr) &
      bind(c, name="nfor_mpi_bcast4")
    integer(c_int), value, intent(in) :: n4
    real(c_double), intent(inout) :: buf(n4)
    integer(c_int), intent(out) :: ierr
    integer :: e2
    call mpi_bcast(buf, n4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, e2)
    ierr = e2
  end subroutine

  ! gather the full-domain state to rank 0 (rust then writes output).
  subroutine nfor_mpi_gather(rho, mx, my, e, nloc, nx, ny_, out, ierr) &
      bind(c, name="nfor_mpi_gather")
    integer(c_int), value, intent(in) :: nloc, nx, ny_
    real(c_double), intent(in) :: rho(nloc), mx(nloc), my(nloc), e(nloc)
    real(c_double), intent(out) :: out(*)
    integer(c_int), intent(out) :: ierr
    real(RP) :: send(4, nloc), recv(4, nx*ny_)
    integer :: counts(g_decomp%nranks), disps(g_decomp%nranks)
    integer :: r, base, rem, off, e2, k, k_src, k_dst, nxl_r, i0_r, i_l, j
    ! per-rank cell counts (mirror make_decomp)
    base = nx / g_decomp%nranks
    rem = mod(nx, g_decomp%nranks)
    do r = 1, g_decomp%nranks
      if (r-1 < rem) then
        counts(r) = (base+1) * ny_
      else
        counts(r) = base * ny_
      end if
    end do
    disps(1) = 0
    do r = 2, g_decomp%nranks
      disps(r) = disps(r-1) + 4*counts(r-1)
    end do
    do k = 1, nloc
      send(1,k) = rho(k); send(2,k) = mx(k); send(3,k) = my(k); send(4,k) = e(k)
    end do
    call mpi_gatherv(send, 4*nloc, MPI_DOUBLE_PRECISION, recv, 4*counts, disps, &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, e2)
    if (g_decomp%rank == 0) then
      ! recv arrives rank-major; scatter into the global row-major layout the
      ! rust side unpacks: global cell (j, i0_r + i_l - 1) = slot
      ! disps(r) + (j-1)*nx_local(r) + i_l.
      do r = 1, g_decomp%nranks
        nxl_r = counts(r) / ny_
        i0_r = (r-1)*base + min(r-1, rem) + 1
        do j = 1, ny_
          do i_l = 1, nxl_r
            k_src = disps(r)/4 + (j-1)*nxl_r + i_l
            k_dst = (j-1)*nx + (i0_r - 1) + i_l
            out(4*(k_dst-1) + 1) = recv(1, k_src)
            out(4*(k_dst-1) + 2) = recv(2, k_src)
            out(4*(k_dst-1) + 3) = recv(3, k_src)
            out(4*(k_dst-1) + 4) = recv(4, k_src)
          end do
        end do
      end do
    end if
    ierr = 0
  end subroutine

  ! 4-field halo exchange for the conserved slab (rho, mx, my, e)
  subroutine halo_exchange_cons(rho, mx, my, e, nx_l, ny_)
    real(RP), intent(inout) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                               my(nx_l+4, ny_), e(nx_l+4, ny_)
    integer, intent(in) :: nx_l, ny_
    real(RP) :: buf(4, nx_l+4, ny_)
    if (g_decomp%nranks == 1) return
    buf(1,:,:) = rho
    buf(2,:,:) = mx
    buf(3,:,:) = my
    buf(4,:,:) = e
    call halo_exchange_x(buf, 4, g_decomp)
    rho = buf(1,:,:)
    mx = buf(2,:,:)
    my = buf(3,:,:)
    e = buf(4,:,:)
  end subroutine

  ! halo exchange for the padded nu_tilde slab (1 variable).
  subroutine halo_exchange_nt(nt, nx_l, ny_)
    integer, intent(in) :: nx_l, ny_
    real(RP), intent(inout) :: nt(nx_l+4, ny_)
    real(RP) :: buf(1, nx_l+4, ny_)
    if (g_decomp%nranks == 1) return
    buf(1,:,:) = nt
    call halo_exchange_x(buf, 1, g_decomp)
    nt = buf(1,:,:)
  end subroutine

  ! one coupled heun step of the sa rans lane on this rank's slab,
  ! mirroring advance2d_sa_model_rk2 op for op. stage k:
  ! euler (heun half) -> halo -> eddy mu_t from the stage-entry
  ! nu_tilde at the post-euler rho -> viscous_add_slab2 with mu_cells
  ! evaluated (bisection route) at the stage-entry state -> halo ->
  ! sa transport with dt_k. stage 2 re-evaluates everything at the
  ! stage-1 output. heun average closes state and nu_tilde together.
  ! the sa diffusivity dt cap arrives from rust as dt_cap (global
  ! each step from the full field; fortran only enforces it inside
  ! the euler stage-1 half, which is where the march picks dt).
  subroutine nfor_mpi_step_sa(rho, mx, my, e, nt, d_slab, nloc, nx_l, ny_, &
                              dx, dy, gamma, cfl, dt_cap, &
                              bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                              mu0, t0_ref, s_par, pr, pr_t, tw, r_gas, &
                              rho2, mx2, my2, e2, nt2, dt_out, ierr) &
      bind(c, name="nfor_mpi_step_sa")
    integer(c_int), value, intent(in) :: nloc, nx_l, ny_
    real(c_double), intent(in) :: rho(nloc), mx(nloc), my(nloc), e(nloc), nt(nloc)
    real(c_double), intent(in) :: d_slab(nloc)
    real(c_double), value, intent(in) :: dx, dy, gamma, cfl, dt_cap
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(5, 4)
    integer(c_int), intent(in) :: solid(nloc)
    real(c_double), value, intent(in) :: mu0, t0_ref, s_par, pr, pr_t, tw, r_gas
    real(c_double), intent(out) :: rho2(nloc), mx2(nloc), my2(nloc), e2(nloc), nt2(nloc)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: ierr
    logical :: sl(nx_l+4, ny_)
    real(RP) :: rin(nx_l+4, ny_), min_(nx_l+4, ny_), myin(nx_l+4, ny_), ein(nx_l+4, ny_)
    real(RP) :: ntin(nx_l+4, ny_)
    real(RP) :: r1(nx_l+4, ny_), m1(nx_l+4, ny_), my1(nx_l+4, ny_), e1(nx_l+4, ny_)
    real(RP) :: r2(nx_l+4, ny_), m2(nx_l+4, ny_), my2_(nx_l+4, ny_), e2_(nx_l+4, ny_)
    real(RP) :: nt1(nx_l+4, ny_), nt2_(nx_l+4, ny_)
    real(RP) :: mu_c(nx_l+4, ny_), mu_t1(nx_l+4, ny_), mu_t2(nx_l+4, ny_)
    real(RP), allocatable :: u1(:,:), v1(:,:), ein1(:,:)
    real(RP), allocatable :: us(:,:), vs(:,:), es(:,:)
    real(RP) :: dslab(nx_l, ny_)
    real(RP) :: dt1, dt2, dt
    real(RP) :: rho_min_l, rho_min_g, nt_max_l, nt_max_g, dt_sa, dt_cap_sa
    integer :: j, i

    allocate(u1(nx_l+4, ny_), v1(nx_l+4, ny_), ein1(nx_l+4, ny_))
    allocate(us(nx_l+4, ny_), vs(nx_l+4, ny_), es(nx_l+4, ny_))

    ! slab carve: conserved state, nu_tilde, wall distance
    do j = 1, ny_
      do i = 1, nx_l
        dslab(i,j) = d_slab((j-1)*nx_l + i)
        rin(i+2,j) = rho((j-1)*nx_l + i)
        min_(i+2,j) = mx((j-1)*nx_l + i)
        myin(i+2,j) = my((j-1)*nx_l + i)
        ein(i+2,j) = e((j-1)*nx_l + i)
        ntin(i+2,j) = nt((j-1)*nx_l + i)
        sl(i+2,j) = solid((j-1)*nx_l + i) == 1
      end do
      rin(1,j) = rin(3,j);   rin(2,j) = rin(3,j)
      rin(nx_l+3,j) = rin(nx_l+2,j); rin(nx_l+4,j) = rin(nx_l+2,j)
      min_(1,j) = min_(3,j); min_(2,j) = min_(3,j)
      min_(nx_l+3,j) = min_(nx_l+2,j); min_(nx_l+4,j) = min_(nx_l+2,j)
      myin(1,j) = myin(3,j); myin(2,j) = myin(3,j)
      myin(nx_l+3,j) = myin(nx_l+2,j); myin(nx_l+4,j) = myin(nx_l+2,j)
      ein(1,j) = ein(3,j);   ein(2,j) = ein(3,j)
      ein(nx_l+3,j) = ein(nx_l+2,j); ein(nx_l+4,j) = ein(nx_l+2,j)
      ntin(1,j) = ntin(3,j); ntin(2,j) = ntin(3,j)
      ntin(nx_l+3,j) = ntin(nx_l+2,j); ntin(nx_l+4,j) = ntin(nx_l+2,j)
      sl(1,j) = .false.; sl(2,j) = .false.
      sl(nx_l+3,j) = .false.; sl(nx_l+4,j) = .false.
    end do
    call halo_exchange_cons(rin, min_, myin, ein, nx_l, ny_)
    call halo_exchange_nt(ntin, nx_l, ny_)

    ! ---- the sa diffusivity dt cap ----
    ! rust: rho_min over the whole field (max 1e-6), nt_max over the
    ! whole field floored at nt_inf; dt_sa = 0.25 min(dx,dy)^2 sigma /
    ! (mu/rho_min + nt_max), allreduce max/min over ranks.
    rho_min_l = minval(rin(3:nx_l+2, 1:ny_))
    rho_min_l = max(rho_min_l, 1.0e-6_RP)
    nt_max_l = maxval(ntin(3:nx_l+2, 1:ny_))
    call global_min(rho_min_l, rho_min_g, g_decomp)
    call global_max(nt_max_l, nt_max_g, g_decomp)
    dt_sa = 0.25_RP * min(dx, dy)**2 * SIGMA_SA / (mu0 / rho_min_g + nt_max_g)
    dt_cap_sa = min(dt_cap, dt_sa)

    ! ---- stage 1 ----
    ! mu_cells at the step-entry state (bisection route)
    call primitive_ein(rin, min_, myin, ein, nx_l, ny_, u1, v1, ein1)
    call sa_mu_slab(rin, ein1, nx_l, ny_, mu0, t0_ref, s_par, mu_c)
    call euler2d_step_rk2(rin, min_, myin, ein, nx_l, ny_, &
                          dx, dy, gamma, cfl, dt_cap_sa, g_decomp, &
                          bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                          r1, m1, my1, e1, dt1, 1)

    call halo_exchange_cons(r1, m1, my1, e1, nx_l, ny_)
    call reseed_edge_ghosts(r1, nx_l, ny_)
    call reseed_edge_ghosts(m1, nx_l, ny_)
    call reseed_edge_ghosts(my1, nx_l, ny_)
    call reseed_edge_ghosts(e1, nx_l, ny_)
    call eddy_mu_t_slab(ntin, r1, mu0, nx_l, ny_, mu_t1)
    call viscous_add_slab2(r1, m1, my1, e1, nx_l, ny_, dx, dy, dt1, &
                           mu_c, mu_t1, .true., mu0, t0_ref, s_par, &
                           pr, pr_t, tw, gamma, r_gas, bc_s, bc_n, &
                           g_decomp%west_domain_edge, g_decomp%east_domain_edge)
    call halo_exchange_cons(r1, m1, my1, e1, nx_l, ny_)
    call reseed_edge_ghosts(r1, nx_l, ny_)
    call reseed_edge_ghosts(m1, nx_l, ny_)
    call reseed_edge_ghosts(my1, nx_l, ny_)
    call reseed_edge_ghosts(e1, nx_l, ny_)
    ! transport on the post-viscous state (ntin -> nt1)
    call primitive_ein(r1, m1, my1, e1, nx_l, ny_, us, vs, es)
    nt1 = ntin
    call sa_transport_slab(nt1, r1, us, vs, dslab, nx_l, ny_, dx, dy, dt1, &
                           mu0, bc_s, bc_n)
    call halo_exchange_nt(nt1, nx_l, ny_)
    call reseed_edge_ghosts(nt1, nx_l, ny_)

    ! ---- stage 2 (re-evaluate at the stage-1 state) ----
    call reseed_edge_ghosts(r1, nx_l, ny_)
    call reseed_edge_ghosts(m1, nx_l, ny_)
    call reseed_edge_ghosts(my1, nx_l, ny_)
    call reseed_edge_ghosts(e1, nx_l, ny_)
    call primitive_ein(r1, m1, my1, e1, nx_l, ny_, u1, v1, ein1)
    call sa_mu_slab(r1, ein1, nx_l, ny_, mu0, t0_ref, s_par, mu_c)
    call euler2d_step_rk2(r1, m1, my1, e1, nx_l, ny_, &
                          dx, dy, gamma, cfl, dt_cap_sa, g_decomp, &
                          bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                          r2, m2, my2_, e2_, dt2, 1)
    call halo_exchange_cons(r2, m2, my2_, e2_, nx_l, ny_)
    call reseed_edge_ghosts(r2, nx_l, ny_)
    call reseed_edge_ghosts(m2, nx_l, ny_)
    call reseed_edge_ghosts(my2_, nx_l, ny_)
    call reseed_edge_ghosts(e2_, nx_l, ny_)
    ! stage-2 eddy viscosity: rust's mu_t2 = mu_t_of(&t1.nu_tilde, &s1.rho)
    ! runs AFTER the stage-2 euler mutated s1, so the rho field is the
    ! post-euler r2 and the nu_tilde field is the stage-1 transported nt1.
    call eddy_mu_t_slab(nt1, r2, mu0, nx_l, ny_, mu_t2)
    call viscous_add_slab2(r2, m2, my2_, e2_, nx_l, ny_, dx, dy, dt2, &
                           mu_c, mu_t2, .true., mu0, t0_ref, s_par, &
                           pr, pr_t, tw, gamma, r_gas, bc_s, bc_n, &
                           g_decomp%west_domain_edge, g_decomp%east_domain_edge)
    call halo_exchange_cons(r2, m2, my2_, e2_, nx_l, ny_)
    call reseed_edge_ghosts(r2, nx_l, ny_)
    call reseed_edge_ghosts(m2, nx_l, ny_)
    call reseed_edge_ghosts(my2_, nx_l, ny_)
    call reseed_edge_ghosts(e2_, nx_l, ny_)
    call primitive_ein(r2, m2, my2_, e2_, nx_l, ny_, us, vs, es)
    nt2_ = nt1
    call sa_transport_slab(nt2_, r2, us, vs, dslab, nx_l, ny_, dx, dy, dt2, &
                           mu0, bc_s, bc_n)
    call reseed_edge_ghosts(nt2_, nx_l, ny_)
    call halo_exchange_nt(nt2_, nx_l, ny_)

    ! ---- heun average of the entry state and the twice-advanced state ----
    dt = dt1
    do j = 1, ny_
      do i = 1, nx_l
        rho2((j-1)*nx_l + i) = 0.5_RP * (rin(i+2,j) + r2(i+2,j))
        mx2((j-1)*nx_l + i) = 0.5_RP * (min_(i+2,j) + m2(i+2,j))
        my2((j-1)*nx_l + i) = 0.5_RP * (myin(i+2,j) + my2_(i+2,j))
        e2((j-1)*nx_l + i) = 0.5_RP * (ein(i+2,j) + e2_(i+2,j))
        nt2((j-1)*nx_l + i) = 0.5_RP * (ntin(i+2,j) + nt2_(i+2,j))
      end do
    end do
    dt_out = dt
    ierr = 0
    deallocate(u1, v1, ein1, us, vs, es)
  end subroutine

  ! gather the full-domain nu_tilde field to rank 0 (one variable).
  subroutine nfor_mpi_gather_nt(nt, nloc, nx, ny_, out, ierr) &
      bind(c, name="nfor_mpi_gather_nt")
    integer(c_int), value, intent(in) :: nloc, nx, ny_
    real(c_double), intent(in) :: nt(nloc)
    real(c_double), intent(out) :: out(*)
    integer(c_int), intent(out) :: ierr
    real(RP) :: send(nloc), recv(nx*ny_)
    integer :: counts(g_decomp%nranks), disps(g_decomp%nranks)
    integer :: r, base, rem, off, e2, k, k_src, k_dst, nxl_r, i0_r, i_l, j
    base = nx / g_decomp%nranks
    rem = mod(nx, g_decomp%nranks)
    do r = 1, g_decomp%nranks
      if (r-1 < rem) then
        counts(r) = (base+1) * ny_
      else
        counts(r) = base * ny_
      end if
    end do
    disps(1) = 0
    do r = 2, g_decomp%nranks
      disps(r) = disps(r-1) + counts(r-1)
    end do
    do k = 1, nloc
      send(k) = nt(k)
    end do
    call mpi_gatherv(send, nloc, MPI_DOUBLE_PRECISION, recv, counts, disps, &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, e2)
    if (g_decomp%rank == 0) then
      do r = 1, g_decomp%nranks
        nxl_r = counts(r) / ny_
        i0_r = (r-1)*base + min(r-1, rem) + 1
        do j = 1, ny_
          do i_l = 1, nxl_r
            k_src = disps(r) + (j-1)*nxl_r + i_l
            k_dst = (j-1)*nx + (i0_r - 1) + i_l
            out(k_dst) = recv(k_src)
          end do
        end do
      end do
    end if
    ierr = 0
  end subroutine

end module nuforbridge

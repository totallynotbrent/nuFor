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
      nbk12(i) = 1   ! blocks with grabau terms; rust marks has_grabau
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
    real(RP) :: rout(nx_l+4, ny_), mout(nx_l+4, ny_), myout(nx_l+4, ny_), eout(nx_l+4, ny_)
    real(RP) :: dt
    integer :: j, i
    ! unpack rust flat (nx_l x ny, row-major) into the padded slab: interiors
    ! 3..nx_l+2, ghosts seeded as copies of the edge interiors (halo exchange
    ! overwrites at np>1)
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
    call euler2d_step_rk2(rin, min_, myin, ein, nx_l, ny_, &
                          dx, dy, gamma, cfl, dt_cap, g_decomp, &
                          bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                          rout, mout, myout, eout, dt, 1)
    ! pack back to flat (interior only)
    do j = 1, ny_
      do i = 1, nx_l
        rho2((j-1)*nx_l + i) = rout(i+2,j)
        mx2((j-1)*nx_l + i) = mout(i+2,j)
        my2((j-1)*nx_l + i) = myout(i+2,j)
        e2((j-1)*nx_l + i) = eout(i+2,j)
      end do
    end do
    dt_out = dt
    ierr = 0
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

end module nuforbridge

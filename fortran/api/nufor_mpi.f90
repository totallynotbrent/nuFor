! nufor_mpi.f90 — decomposition, halo exchange, and global reductions for the
! mpi-parallel 2d march. one module, no save-state surprises: every routine
! takes its context explicitly so ranks stay reentrant.
module nufor_mpi
  use mpi
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: DP = real64

  ! decomposition info for one rank (1d slabs along x for v1).
  type :: decomp_t
    integer :: nranks        ! total ranks
    integer :: rank          ! this rank 0-based
    integer :: i0, i1        ! cell x-range [i0, i1] inclusive, global indexing
    integer :: nx_local      ! i1 - i0 + 1
    integer :: ny            ! full y extent (no y decomposition in v1)
    integer :: west_rank, east_rank  ! neighbor ranks, MPI_PROC_NULL at edges
    logical :: west_domain_edge, east_domain_edge
  end type

contains

  ! build the 1d slab decomposition: contiguous x ranges, sizes as even as
  ! the remainder allows. ranks of the west/east neighbors come from rank
  ! order; domain-edge ranks get MPI_PROC_NULL and an edge flag.
  function make_decomp(nx, ny, nranks, rank) result(d)
    integer, intent(in) :: nx, ny, nranks, rank
    type(decomp_t) :: d
    integer :: base, rem, k
    if (nx < nranks) then
      ! more ranks than columns: early ranks get one column, the rest idle
      ! (legal; the caller should avoid it for performance).
      if (rank < nx) then
        d%i0 = rank + 1
        d%i1 = rank + 1
      else
        d%i0 = 1
        d%i1 = 0   ! empty range
      end if
    else
      base = nx / nranks
      rem = mod(nx, nranks)
      if (rank < rem) then
        d%i0 = rank * (base + 1) + 1
        d%i1 = d%i0 + base
      else
        d%i0 = rem * (base + 1) + (rank - rem) * base + 1
        d%i1 = d%i0 + base - 1
      end if
    end if
    d%nranks = nranks
    d%rank = rank
    d%nx_local = d%i1 - d%i0 + 1
    d%ny = ny
    d%west_rank = rank - 1
    d%east_rank = rank + 1
    if (d%west_rank < 0) d%west_rank = MPI_PROC_NULL
    if (d%east_rank > nranks - 1) d%east_rank = MPI_PROC_NULL
    d%west_domain_edge = (rank == 0)
    d%east_domain_edge = (rank == nranks - 1)
  end function

  ! exchange two x-face halo layers: each rank sends its two edge interior
  ! columns to each neighbor and receives the neighbor's two edge columns as
  ! ghosts. two layers are the minimum for the muscl limiter at subdomain
  ! faces. f is (nvar, nx_local+4, ny); interior columns are 3..nx_local+2;
  ! ghost columns 1,2 (west) and nx_local+3,4 (east).
  subroutine halo_exchange_x(f, nvar, d)
    real(DP), intent(inout) :: f(:,:,:)
    integer, intent(in) :: nvar
    type(decomp_t), intent(in) :: d
    integer, parameter :: TAG_W = 101, TAG_E = 102
    real(DP) :: send_w(nvar, 2*d%ny), send_e(nvar, 2*d%ny)
    real(DP) :: recv_w(nvar, 2*d%ny), recv_e(nvar, 2*d%ny)
    integer :: reqs(4), j, v, ierr

    if (d%nranks == 1) return   ! no neighbors: ghost columns unused at np=1

    do j = 1, d%ny
      do v = 1, nvar
        send_w(v, 2*j-1) = f(v, 3, j)               ! first interior column
        send_w(v, 2*j)   = f(v, 4, j)               ! second interior column
        send_e(v, 2*j-1) = f(v, d%nx_local+1, j)   ! second-to-last interior
        send_e(v, 2*j)   = f(v, d%nx_local+2, j)   ! last interior column
      end do
    end do

    reqs = MPI_REQUEST_NULL
    if (d%west_rank /= MPI_PROC_NULL) &
      call mpi_irecv(recv_w, nvar*2*d%ny, MPI_DOUBLE_PRECISION, d%west_rank, &
                     TAG_E, MPI_COMM_WORLD, reqs(1), ierr)
    if (d%east_rank /= MPI_PROC_NULL) &
      call mpi_irecv(recv_e, nvar*2*d%ny, MPI_DOUBLE_PRECISION, d%east_rank, &
                     TAG_W, MPI_COMM_WORLD, reqs(2), ierr)
    if (d%west_rank /= MPI_PROC_NULL) &
      call mpi_isend(send_w, nvar*2*d%ny, MPI_DOUBLE_PRECISION, d%west_rank, &
                     TAG_W, MPI_COMM_WORLD, reqs(3), ierr)
    if (d%east_rank /= MPI_PROC_NULL) &
      call mpi_isend(send_e, nvar*2*d%ny, MPI_DOUBLE_PRECISION, d%east_rank, &
                     TAG_E, MPI_COMM_WORLD, reqs(4), ierr)
    call mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)

    ! write received ghosts: west ghost cols 1,2 hold the west neighbor's
    ! last two interior columns (its col order preserved); east ghosts
    ! nx_local+3,4 hold the east neighbor's first two.
    if (d%west_rank /= MPI_PROC_NULL) then
      do j = 1, d%ny
        do v = 1, nvar
          f(v, 1, j) = recv_w(v, 2*j-1)
          f(v, 2, j) = recv_w(v, 2*j)
        end do
      end do
    end if
    if (d%east_rank /= MPI_PROC_NULL) then
      do j = 1, d%ny
        do v = 1, nvar
          f(v, d%nx_local+3, j) = recv_e(v, 2*j-1)
          f(v, d%nx_local+4, j) = recv_e(v, 2*j)
        end do
      end do
    end if
  end subroutine

  ! global max over all ranks (cfl smax, dt).
  subroutine global_max(local, global, d)
    real(DP), intent(in) :: local
    real(DP), intent(out) :: global
    type(decomp_t), intent(in) :: d
    integer :: ierr
    call mpi_allreduce(local, global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
  end subroutine

  subroutine global_min(local, global, d)
    real(DP), intent(in) :: local
    real(DP), intent(out) :: global
    type(decomp_t), intent(in) :: d
    integer :: ierr
    call mpi_allreduce(local, global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                       MPI_COMM_WORLD, ierr)
  end subroutine

  ! is this rank the output rank (rank 0)?
  pure function is_rank0(d) result(r)
    type(decomp_t), intent(in) :: d
    logical :: r
    r = d%rank == 0
  end function


  ! ---------------------------------------------------------------
  ! seam column agreement: publish the two columns that close one end
  ! of this slab to the neighbour on that side and compare them against
  ! the ghost columns the halo delivered. the two views of a shared
  ! global column must agree exactly, and the largest disagreement is
  ! reported rather than assumed.
  !
  ! the west ghosts 1,2 hold the west neighbour's last two interior
  ! columns; the east ghosts nx_local+3,nx_local+4 hold the east
  ! neighbour's first two. f is (nvar, nx_local+4, ny) with the interior
  ! columns at 3..nx_local+2, matching halo_exchange_x.
  subroutine seam_column_check(f, nvar, d, max_diff)
    real(DP), intent(in) :: f(:,:,:)
    integer, intent(in) :: nvar
    type(decomp_t), intent(in) :: d
    real(DP), intent(out) :: max_diff
    integer, parameter :: TAG_LOW = 431, TAG_HIGH = 432
    real(DP) :: snd_lo(2*nvar*d%ny), snd_hi(2*nvar*d%ny)
    real(DP) :: rcv_lo(2*nvar*d%ny), rcv_hi(2*nvar*d%ny)
    real(DP) :: loc_lo(2*nvar*d%ny), loc_hi(2*nvar*d%ny)
    real(DP) :: dd, gd
    integer :: reqs(4), j, v, ierr

    if (d%nranks == 1) then
      max_diff = 0.0_DP
      return
    end if

    do j = 1, d%ny
      do v = 1, nvar
        snd_lo(v+(j-1)*nvar) = f(v, 3, j)
        snd_lo(v+d%ny*nvar+(j-1)*nvar) = f(v, 4, j)
        snd_hi(v+(j-1)*nvar) = f(v, d%nx_local+1, j)
        snd_hi(v+d%ny*nvar+(j-1)*nvar) = f(v, d%nx_local+2, j)
        loc_lo(v+(j-1)*nvar) = f(v, 1, j)
        loc_lo(v+d%ny*nvar+(j-1)*nvar) = f(v, 2, j)
        loc_hi(v+(j-1)*nvar) = f(v, d%nx_local+3, j)
        loc_hi(v+d%ny*nvar+(j-1)*nvar) = f(v, d%nx_local+4, j)
      end do
    end do

    ! the low columns travel west and the high columns travel east; the
    ! receiver hears them from the opposite direction of travel.
    dd = 0.0_DP
    reqs = MPI_REQUEST_NULL
    if (d%east_rank /= MPI_PROC_NULL) &
      call mpi_irecv(rcv_lo, 2*nvar*d%ny, MPI_DOUBLE_PRECISION, d%east_rank, &
                     TAG_LOW, MPI_COMM_WORLD, reqs(1), ierr)
    if (d%west_rank /= MPI_PROC_NULL) &
      call mpi_irecv(rcv_hi, 2*nvar*d%ny, MPI_DOUBLE_PRECISION, d%west_rank, &
                     TAG_HIGH, MPI_COMM_WORLD, reqs(2), ierr)
    if (d%west_rank /= MPI_PROC_NULL) &
      call mpi_isend(snd_lo, 2*nvar*d%ny, MPI_DOUBLE_PRECISION, d%west_rank, &
                     TAG_LOW, MPI_COMM_WORLD, reqs(3), ierr)
    if (d%east_rank /= MPI_PROC_NULL) &
      call mpi_isend(snd_hi, 2*nvar*d%ny, MPI_DOUBLE_PRECISION, d%east_rank, &
                     TAG_HIGH, MPI_COMM_WORLD, reqs(4), ierr)
    if (d%east_rank /= MPI_PROC_NULL) then
      call mpi_wait(reqs(1), MPI_STATUS_IGNORE, ierr)
      dd = max(dd, max_abs_diff_arrays(rcv_lo, loc_hi, 2*nvar*d%ny))
    end if
    if (d%west_rank /= MPI_PROC_NULL) then
      call mpi_wait(reqs(2), MPI_STATUS_IGNORE, ierr)
      dd = max(dd, max_abs_diff_arrays(rcv_hi, loc_lo, 2*nvar*d%ny))
    end if
    call mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
    call global_max(dd, gd, d)
    max_diff = gd
  end subroutine seam_column_check

  ! largest absolute difference between two equal-length arrays.
  pure function max_abs_diff_arrays(a, b, n) result(m)
    real(DP) :: m
    integer, intent(in) :: n
    real(DP), intent(in) :: a(n), b(n)
    integer :: i
    m = 0.0_DP
    do i = 1, n
      m = max(m, abs(a(i)-b(i)))
    end do
  end function max_abs_diff_arrays

end module nufor_mpi

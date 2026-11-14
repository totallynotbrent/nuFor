! nuforcollect.f90 - the c-abi collectors: rank probe, chunk broadcast, and
! the gathers that reassemble the global row-major fields from the per-rank
! slabs (the conserved four, an arbitrary set of planes, and the turbulence
! plane).
module nuforcollect
  use iso_c_binding, only: c_int, c_double
  use iso_fortran_env, only: real64
  use mpi
  use nufor_mpi
  use nuforhalos, only: g_decomp, g_rank
  implicit none

  integer, parameter :: RP = real64

contains

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

  ! general nf-field gather: each of the `nf` planes (plane-major, plane*padded
  ! followed by cells) is scattered into the global row-major layout, field
  ! f occupying out((f-1)*nx*ny_ + k_dst). lets non-conserved planes (the
  ! species compositions) ride the same slab layout as the conserved four.
  subroutine nfor_mpi_gather_nf(fields, nf, nloc, nx, ny_, out, ierr) &
      bind(c, name="nfor_mpi_gather_nf")
    integer(c_int), value, intent(in) :: nf, nloc, nx, ny_
    real(c_double), intent(in) :: fields(*)
    real(c_double), intent(out) :: out(*)
    integer(c_int), intent(out) :: ierr
    real(RP), allocatable :: send(:), recv(:)
    integer :: counts(g_decomp%nranks), disps(g_decomp%nranks)
    integer :: r, base, rem, e2, f, k, l_src, l_recv, k_dst, nxl_r, i0_r, i_l, j
    allocate(send(nf*nloc), recv(nf*nx*ny_))
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
      disps(r) = disps(r-1) + nf*counts(r-1)
    end do
    ! flatten (k, f) column-major pairs into field-contiguous blocks: send
    ! holds field f's slab at [(f-1)*nloc + k], matching the rust caller's
    ! plane-major buffer.
    do f = 1, nf
      do k = 1, nloc
        send((f-1)*nloc + k) = fields((f-1)*nloc + k)
      end do
    end do
    call mpi_gatherv(send, nf*nloc, MPI_DOUBLE_PRECISION, recv, nf*counts, disps, &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, e2)
    if (g_decomp%rank == 0) then
      ! recv holds rank r's field f block at disps(r) + (f-1)*counts(r) + l;
      ! scatter into the global row-major, field-major out buffer.
      do r = 1, g_decomp%nranks
        nxl_r = counts(r) / ny_
        i0_r = (r-1)*base + min(r-1, rem) + 1
        do j = 1, ny_
          do i_l = 1, nxl_r
            l_src = (j-1)*nxl_r + i_l
            k_dst = (j-1)*nx + (i0_r - 1) + i_l
            do f = 1, nf
              l_recv = disps(r) + (f-1)*counts(r) + l_src
              out((f-1)*nx*ny_ + k_dst) = recv(l_recv)
            end do
          end do
        end do
      end do
    end if
    deallocate(send, recv)
    ierr = 0
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

end module nuforcollect

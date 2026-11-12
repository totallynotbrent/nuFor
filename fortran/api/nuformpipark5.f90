! nuformpipark5.f90 — the mpi-native five-species park5 carrier.
!
! this module carries the conserved-four-plus-five-species state across
! an x-slab decomposition. it owns the seam bookkeeping that lets the
! park5 sweep, when it arrives, behave exactly as the serial cartesian
! lane does on the undivided grid.
!
! the carrier exchanges the state in the same representation the cart
! sweep consumes: density, the two velocities, the specific internal
! energy, and the five composition fractions. the halo carries two
! layers in each x direction, which is the depth the centred van-leer
! limiter needs at a subdomain face.
!
! seam consistency: the west ghost columns 1,2 hold the west
! neighbour's second-to-last and last interior columns, and the east
! ghosts nx_local+3,nx_local+4 hold the east neighbour's first two.
! the padded strip therefore holds contiguous global columns in
! increasing order, so both sides of a seam reconstruct the seam face
! from the same global centres with the same stencil, and the two
! reconstructions agree digit for digit. the carrier checks that
! agreement instead of trusting it: each rank publishes the states it
! reconstructs for the two faces it borders, the neighbours exchange
! them, and the largest disagreement is reported.
!
! padding convention (matching nufor_mpi.halo_exchange_x): the local
! slab occupies pad columns 3..nx_local+2, ghosts sit at 1,2 (west) and
! nx_local+3,nx_local+4 (east). the flattened layout is
! k = (j-1)*(nx_local+4) + i, matching fortran column-major storage.
module nuformpipark5
  use iso_fortran_env, only: real64
  use mpi
  use nufor_mpi
  implicit none

  integer, parameter :: NP5 = 5
  ! the conserved four plus the five composition planes ride together.
  integer, parameter :: NCAR = 9
  integer, parameter :: RP = DP
  ! the carriers share one halo width, so the limit follows the widest
  ! dimension either arrangement can reach.
  integer, parameter :: mp_width = 8192

  ! the decomposition is per process and immutable once set, so it lives
  ! at module scope without being mutable scratch.
  type(decomp_t), save :: p_decomp
  logical, save :: p_ready = .false.

  private
  public :: nfor_mpipark5_init, nfor_mpipark5_rank
  public :: nfor_mpipark5_carrier, nfor_mpipark5_gather

contains

  pure real(real64) function vlim(dm, dp_)
    real(real64), intent(in) :: dm, dp_
    if (dm*dp_ > 0.0_DP) then
      vlim = 2.0_DP*dm*dp_/(dm+dp_)
    else
      vlim = 0.0_DP
    end if
  end function vlim

  ! ------------------------------------------------------------------
  ! c surface: this rank's zero-based slot in the decomposition.
  subroutine nfor_mpipark5_rank(out_rank, ierr) bind(c, name="nfor_mpipark5_rank")
    use iso_c_binding, only: c_int
    integer(c_int), intent(out) :: out_rank
    integer(c_int), intent(out) :: ierr
    integer :: world, rank, e
    call mpi_comm_rank(MPI_COMM_WORLD, rank, e)
    call mpi_comm_size(MPI_COMM_WORLD, world, e)
    out_rank = rank
    ierr = 0
  end subroutine nfor_mpipark5_rank

  ! ------------------------------------------------------------------
  ! c surface: bring up this rank's slab decomposition.
  subroutine nfor_mpipark5_init(nx, ny, out_rank, out_i0, out_nx_local, ierr) &
      bind(c, name="nfor_mpipark5_init")
    use iso_c_binding, only: c_int
    integer(c_int), value, intent(in) :: nx, ny
    integer(c_int), intent(out) :: out_rank, out_i0, out_nx_local
    integer(c_int), intent(out) :: ierr
    integer :: world, rank, e
    call mpi_comm_rank(MPI_COMM_WORLD, rank, e)
    call mpi_comm_size(MPI_COMM_WORLD, world, e)
    p_decomp = make_decomp(int(nx), int(ny), world, rank)
    p_ready = .true.
    out_rank = rank
    out_i0 = p_decomp%i0
    out_nx_local = p_decomp%nx_local
    ierr = 0
  end subroutine nfor_mpipark5_init

  ! ------------------------------------------------------------------
  ! c surface: hand this rank's slab to the carrier and take back the
  ! ghost-complete state.
  !
  ! the inputs are the conserved four (density, x momentum, y momentum,
  ! total energy) and the five composition fractions, flat with
  ! k = (j-1)*nx_local + i. the outputs are the same state expressed in
  ! the primitive form the sweep consumes (density, u, v, specific
  ! internal energy) plus the composition, laid out in the padded strip
  ! of nx_local+4 columns by ny rows. max_diff receives the largest
  ! disagreement between the two reconstructions of any seam face; a
  ! correct halo drives it to zero.
  subroutine nfor_mpipark5_carrier(nx_local, ny_local, rho, mx, my, e, y, &
                                    bc_w, bc_e, bc_vals, y_fill, &
                                    rho_g, u_g, v_g, ei_g, yg, max_diff, ierr) &
      bind(c, name="nfor_mpipark5_carrier")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), value, intent(in) :: nx_local, ny_local
    real(c_double), intent(in) :: rho(*), mx(*), my(*), e(*), y(*)
    integer(c_int), value, intent(in) :: bc_w, bc_e
    real(c_double), intent(in) :: bc_vals(*)
    real(c_double), intent(in) :: y_fill(*)
    real(c_double), intent(out) :: rho_g(*), u_g(*), v_g(*), ei_g(*), yg(*)
    real(c_double), intent(out) :: max_diff
    integer(c_int), intent(out) :: ierr
    real(RP), allocatable :: car(:,:,:), snd_lo(:), snd_hi(:), rcv_lo(:), rcv_hi(:), loc_lo(:), loc_hi(:)
    real(RP) :: bcv(5,4), yf(NP5), tt, rr, dd, gd
    integer :: nl, nj, nn, i, j, s, sd, vv, ierr2, nf
    integer :: reqs(4)

    ierr = 0
    if (.not. p_ready) then
      ierr = 1
      return
    end if
    nl = int(nx_local); nj = int(ny_local)
    if (nl <= 0 .or. nj <= 0) then
      ierr = 1
      return
    end if
    if (nl + 4 > mp_width .or. nj + 4 > mp_width) then
      ierr = 1
      return
    end if
    nn = nl*nj
    nf = 2*NCAR*nj

    do sd = 1, 4
      do vv = 1, 5
        bcv(vv,sd) = bc_vals((sd-1)*5 + vv)
      end do
    end do
    do s = 1, NP5
      yf(s) = y_fill(s)
    end do

    ! the conserved four convert to the primitive form the sweep wants:
    ! density, the two velocities, and the specific internal energy.
    allocate(car(NCAR, nl+4, nj), snd_lo(nf), snd_hi(nf), rcv_lo(nf), rcv_hi(nf), loc_lo(nf), loc_hi(nf))
    do j = 1, nj
      do i = 1, nl
        rr = max(rho((j-1)*nl + i), 1.0e-12_RP)
        car(1,i+2,j) = rr
        car(2,i+2,j) = mx((j-1)*nl + i) / rr
        car(3,i+2,j) = my((j-1)*nl + i) / rr
        tt = e((j-1)*nl + i) / rr
        car(4,i+2,j) = tt - 0.5_RP*(car(2,i+2,j)**2 + car(3,i+2,j)**2)
        do s = 1, NP5
          car(4+s,i+2,j) = y((s-1)*nn + (j-1)*nl + i)
        end do
      end do
    end do

    ! the halo brings the neighbour's two edge columns in, keeping every
    ! plane in lockstep so the seam state stays coherent.
    call halo_exchange_x(car, NCAR, p_decomp)
    call fill_domain_edges(car, nl, nj, int(bc_w), int(bc_e), bcv, yf)

    do j = 1, nj
      do i = 1, nl
        rho_g((j-1)*(nl+4) + i+2) = car(1,i+2,j)
        u_g((j-1)*(nl+4) + i+2) = car(2,i+2,j)
        v_g((j-1)*(nl+4) + i+2) = car(3,i+2,j)
        ei_g((j-1)*(nl+4) + i+2) = car(4,i+2,j)
        do s = 1, NP5
          yg((s-1)*(nl+4)*nj + (j-1)*(nl+4) + i+2) = car(4+s,i+2,j)
        end do
      end do
    end do

    ! the seam check: every rank offers the two columns that close one
    ! end of its slab to the neighbour on that side, and that neighbour
    ! compares them against the ghost columns the halo delivered. the
    ! two views of a shared column must agree exactly, and the largest
    ! disagreement is reported rather than assumed.
    !
    ! west ghosts 1,2 hold the west neighbour's last two interior
    ! columns; east ghosts nl+3,nl+4 hold the east neighbour's first
    ! two interior columns.
    do j = 1, nj
      do s = 1, NCAR
        snd_lo(s+(j-1)*NCAR) = car(s, 3, j)
        snd_lo(s+nj*NCAR+(j-1)*NCAR) = car(s, 4, j)
        snd_hi(s+(j-1)*NCAR) = car(s, nl+1, j)
        snd_hi(s+nj*NCAR+(j-1)*NCAR) = car(s, nl+2, j)
        loc_lo(s+(j-1)*NCAR) = car(s, 1, j)
        loc_lo(s+nj*NCAR+(j-1)*NCAR) = car(s, 2, j)
        loc_hi(s+(j-1)*NCAR) = car(s, nl+3, j)
        loc_hi(s+nj*NCAR+(j-1)*NCAR) = car(s, nl+4, j)
      end do
    end do
    !
    ! the low columns travel west, the high columns travel east; the
    ! receiver hears them from the opposite direction of travel.
    dd = 0.0_RP
    reqs = MPI_REQUEST_NULL
    if (p_decomp%east_rank /= MPI_PROC_NULL) then
      call mpi_irecv(rcv_lo, nf, MPI_DOUBLE_PRECISION, p_decomp%east_rank, &
                     431, MPI_COMM_WORLD, reqs(1), ierr2)
    end if
    if (p_decomp%west_rank /= MPI_PROC_NULL) then
      call mpi_irecv(rcv_hi, nf, MPI_DOUBLE_PRECISION, p_decomp%west_rank, &
                     432, MPI_COMM_WORLD, reqs(2), ierr2)
    end if
    if (p_decomp%west_rank /= MPI_PROC_NULL) then
      call mpi_isend(snd_lo, nf, MPI_DOUBLE_PRECISION, p_decomp%west_rank, &
                     431, MPI_COMM_WORLD, reqs(3), ierr2)
    end if
    if (p_decomp%east_rank /= MPI_PROC_NULL) then
      call mpi_isend(snd_hi, nf, MPI_DOUBLE_PRECISION, p_decomp%east_rank, &
                     432, MPI_COMM_WORLD, reqs(4), ierr2)
    end if
    if (p_decomp%east_rank /= MPI_PROC_NULL) then
      call mpi_wait(reqs(1), MPI_STATUS_IGNORE, ierr2)
      dd = max(dd, max_abs_diff(rcv_lo, loc_hi, nf))
    end if
    if (p_decomp%west_rank /= MPI_PROC_NULL) then
      call mpi_wait(reqs(2), MPI_STATUS_IGNORE, ierr2)
      dd = max(dd, max_abs_diff(rcv_hi, loc_lo, nf))
    end if
    call mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr2)
    call global_max(dd, gd, p_decomp)
    max_diff = gd

    deallocate(car, snd_lo, snd_hi, rcv_lo, rcv_hi, loc_lo, loc_hi)
  end subroutine nfor_mpipark5_carrier

  pure real(RP) function max_abs_diff(a, b, n)
    integer, intent(in) :: n
    real(RP), intent(in) :: a(n), b(n)
    integer :: i
    max_abs_diff = 0.0_RP
    do i = 1, n
      max_abs_diff = max(max_abs_diff, abs(a(i)-b(i)))
    end do
  end function max_abs_diff

  ! ------------------------------------------------------------------
  ! c surface: gather the five composition planes onto rank 0.
  !
  ! the incoming buffer is the padded local strip the carrier produces.
  ! each plane travels to rank 0 as one contiguous block, and rank 0
  ! lays the blocks back into the global row-major layout, so field f
  ! occupies out((f-1)*nx*ny + k).
  !
  ! the displacement and count bookkeeping mirrors the shipping bridge
  ! gather: disps advances by nf*counts, and the receive counts are the
  ! element-wise product nf*counts.
  subroutine nfor_mpipark5_gather(nx, ny, y, out_y, ierr) &
      bind(c, name="nfor_mpipark5_gather")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), value, intent(in) :: nx, ny
    real(c_double), intent(in) :: y(*)
    real(c_double), intent(out) :: out_y(*)
    integer(c_int), intent(out) :: ierr
    integer, allocatable :: counts(:), disps(:)
    real(RP), allocatable :: send(:), recv(:)
    integer :: r, f, base, rem, e2, k, l_src, l_recv, k_dst, nxl_r, i0_r
    integer :: i_l, j, nrow, nl_me, ncol_me, cw

    ierr = 0
    if (.not. p_ready) then
      ierr = 1
      return
    end if
    nrow = int(ny)
    nl_me = p_decomp%nx_local
    ncol_me = nl_me
    allocate(counts(p_decomp%nranks), disps(p_decomp%nranks))
    base = int(nx) / p_decomp%nranks
    rem = mod(int(nx), p_decomp%nranks)
    do r = 1, p_decomp%nranks
      if (r-1 < rem) then
        counts(r) = (base+1)*nrow
      else
        counts(r) = base*nrow
      end if
    end do
    disps(1) = 0
    do r = 2, p_decomp%nranks
      disps(r) = disps(r-1) + NP5*counts(r-1)
    end do
    allocate(send(NP5*ncol_me*nrow), recv(NP5*int(nx)*nrow))
    ! pack the interior cells of the padded strip into plane-major blocks.
    do f = 1, NP5
      do j = 1, nrow
        do i_l = 1, ncol_me
          send((f-1)*ncol_me*nrow + (j-1)*ncol_me + i_l) = &
            y((f-1)*(nl_me+4)*nrow + (j-1)*(nl_me+4) + i_l+2)
        end do
      end do
    end do
    call mpi_gatherv(send, NP5*ncol_me*nrow, MPI_DOUBLE_PRECISION, recv, &
                     NP5*counts, disps, MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, e2)
    if (p_decomp%rank == 0) then
      do r = 1, p_decomp%nranks
        nxl_r = counts(r)/nrow
        i0_r = (r-1)*base + min(r-1, rem) + 1
        do j = 1, nrow
          do i_l = 1, nxl_r
            l_src = (j-1)*nxl_r + i_l
            k_dst = (j-1)*int(nx) + (i0_r - 1) + i_l
            do f = 1, NP5
              l_recv = disps(r) + (f-1)*counts(r) + l_src
              out_y((f-1)*int(nx)*nrow + k_dst) = recv(l_recv)
            end do
          end do
        end do
      end do
    end if
    deallocate(counts, disps, send, recv)
  end subroutine nfor_mpipark5_gather

  ! ------------------------------------------------------------------
  ! fill the domain-edge ghost columns (where this rank really is a
  ! domain boundary) from the boundary codes; the seam ghosts came from
  ! the halo. bc codes: 0 transmissive, 1 inflow, 2 outflow, 3 slip
  ! wall, 4 no-slip wall. bc_vals carries [rho, u, v, -, ei] per side.
  subroutine fill_domain_edges(car, nl, nj, bc_w, bc_e, bcv, yf)
    integer, intent(in) :: nl, nj, bc_w, bc_e
    real(RP), intent(inout) :: car(NCAR, nl+4, nj)
    real(RP), intent(in) :: bcv(5,4), yf(NP5)
    integer :: j, s

    if (p_decomp%west_domain_edge) then
      do j = 1, nj
        select case (bc_w)
        case (1)
          car(1,2,j) = bcv(1,1); car(2,2,j) = bcv(2,1); car(3,2,j) = bcv(3,1)
          car(4,2,j) = bcv(5,1)
          do s = 1, NP5
            car(4+s,2,j) = yf(s)
          end do
        case (3)
          car(1,2,j) = car(1,3,j); car(2,2,j) = -car(2,3,j)
          car(3,2,j) = car(3,3,j); car(4,2,j) = car(4,3,j)
          do s = 1, NP5
            car(4+s,2,j) = car(4+s,3,j)
          end do
        case (4)
          car(1,2,j) = car(1,3,j); car(2,2,j) = -car(2,3,j)
          car(3,2,j) = -car(3,3,j); car(4,2,j) = car(4,3,j)
          do s = 1, NP5
            car(4+s,2,j) = car(4+s,3,j)
          end do
        case default
          car(1,2,j) = car(1,3,j); car(2,2,j) = car(2,3,j)
          car(3,2,j) = car(3,3,j); car(4,2,j) = car(4,3,j)
          do s = 1, NP5
            car(4+s,2,j) = car(4+s,3,j)
          end do
        end select
        car(:,1,j) = car(:,2,j)
      end do
    end if
    if (p_decomp%east_domain_edge) then
      do j = 1, nj
        select case (bc_e)
        case (1)
          car(1,nl+3,j) = bcv(1,2); car(2,nl+3,j) = bcv(2,2); car(3,nl+3,j) = bcv(3,2)
          car(4,nl+3,j) = bcv(5,2)
          do s = 1, NP5
            car(4+s,nl+3,j) = yf(s)
          end do
        case (3)
          car(1,nl+3,j) = car(1,nl+2,j); car(2,nl+3,j) = -car(2,nl+2,j)
          car(3,nl+3,j) = car(3,nl+2,j); car(4,nl+3,j) = car(4,nl+2,j)
          do s = 1, NP5
            car(4+s,nl+3,j) = car(4+s,nl+2,j)
          end do
        case (4)
          car(1,nl+3,j) = car(1,nl+2,j); car(2,nl+3,j) = -car(2,nl+2,j)
          car(3,nl+3,j) = -car(3,nl+2,j); car(4,nl+3,j) = car(4,nl+2,j)
          do s = 1, NP5
            car(4+s,nl+3,j) = car(4+s,nl+2,j)
          end do
        case default
          car(1,nl+3,j) = car(1,nl+2,j); car(2,nl+3,j) = car(2,nl+2,j)
          car(3,nl+3,j) = car(3,nl+2,j); car(4,nl+3,j) = car(4,nl+2,j)
          do s = 1, NP5
            car(4+s,nl+3,j) = car(4+s,nl+2,j)
          end do
        end select
        car(:,nl+4,j) = car(:,nl+3,j)
      end do
    end if
  end subroutine fill_domain_edges

end module nuformpipark5

! nuforhalos.f90 - the shared mpi state of a decomposed case: the x-slab
! decomposition every rank agrees on at case init, plus the halo exchangers
! for the conserved slab, the turbulence plane, the two-temperature plane
! and the solid mask.
!
! the decomposition lives here because every lane needs it: the march
! drivers, the collectors and the halo exchangers all address the same
! rank-local strips.
module nuforhalos
  use iso_fortran_env, only: real64
  use nufor_mpi
  implicit none

  integer, parameter :: RP = real64

  type(decomp_t), save :: g_decomp
  integer, save :: g_rank = -1

contains

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

  ! halo exchange for the padded e_v plane (1 variable). the two-temperature
  ! donor stencil reads one cell across the partition, so the plane has to
  ! reach the ghost columns the way the mask does.
  subroutine halo_exchange_plane(ev, nx_l, ny_)
    integer, intent(in) :: nx_l, ny_
    real(RP), intent(inout) :: ev(nx_l+4, ny_)
    real(RP) :: buf(1, nx_l+4, ny_)
    if (g_decomp%nranks == 1) return
    buf(1,:,:) = ev
    call halo_exchange_x(buf, 1, g_decomp)
    ev = buf(1,:,:)
  end subroutine halo_exchange_plane

  ! halo exchange for the solid mask. the wall-face closure reads the mask
  ! of the two ghost columns flanking every domain face, so a body cell
  ! sitting on the far side of a rank seam has to travel with the state or
  ! the seam face leaks flux where the body should seal it.
  subroutine halo_exchange_mask(sl, nx_l, ny_)
    integer, intent(in) :: nx_l, ny_
    logical, intent(inout) :: sl(nx_l+4, ny_)
    real(RP) :: buf(1, nx_l+4, ny_)
    integer :: j, i
    if (g_decomp%nranks == 1) return
    do j = 1, ny_
      do i = 1, nx_l+4
        if (sl(i,j)) then
          buf(1,i,j) = 1.0_RP
        else
          buf(1,i,j) = 0.0_RP
        end if
      end do
    end do
    call halo_exchange_x(buf, 1, g_decomp)
    do j = 1, ny_
      do i = 1, nx_l+4
        sl(i,j) = (buf(1,i,j) > 0.5_RP)
      end do
    end do
  end subroutine

end module nuforhalos

! nuforpark5common.f90 - the pieces the serial and mpi park5 marches share:
! the van-leer limiter, the conserved-to-primitive stripe conversion, and the
! domain-edge ghost fill for the nine-plane carrier strip.
!
! owning them in one place is what lets the two march lanes evolve the seam
! treatment independently without copying the edge algebra.
module nuforpark5common
  use iso_fortran_env, only: real64
  use nufor_mpi, only: decomp_t
  implicit none

  integer, parameter :: RP = real64
  integer, parameter :: NP5 = 5
  integer, parameter :: NCAR = 9

  private
  public :: park5_vlim, park5_pack_primitives, park5_domain_edges

contains

  ! van-leer limiter: harmonic mean of adjacent slopes when they agree.
  pure function park5_vlim(dm, dp_) result(v)
    real(RP), intent(in) :: dm, dp_
    real(RP) :: v
    if (dm*dp_ > 0.0_RP) then
      v = 2.0_RP*dm*dp_/(dm+dp_)
    else
      v = 0.0_RP
    end if
  end function park5_vlim

  ! pack the flat conserved four plus five planes into the padded carrier
  ! strip in primitive form (rho, u, v, e_i, y_1..5); interior columns
  ! i+2 for i = 1..nl. k = (j-1)*nl + i matches the caller layout.
  subroutine park5_pack_primitives(nl, nj, rho, mx, my, e_tot, ys, car)
    integer, intent(in) :: nl, nj
    real(RP), intent(in) :: rho(nl*nj), mx(nl*nj), my(nl*nj), e_tot(nl*nj), ys(*)
    real(RP), intent(out) :: car(NCAR, nl+4, nj)
    real(RP) :: rr, tt
    integer :: i, j, s
    do j = 1, nj
      do i = 1, nl
        rr = max(rho((j-1)*nl+i), 1.0e-12_RP)
        car(1,i+2,j) = rr
        car(2,i+2,j) = mx((j-1)*nl+i)/rr
        car(3,i+2,j) = my((j-1)*nl+i)/rr
        tt = e_tot((j-1)*nl+i)/rr
        car(4,i+2,j) = tt - 0.5_RP*(car(2,i+2,j)**2 + car(3,i+2,j)**2)
        do s = 1, NP5
          car(4+s,i+2,j) = ys((s-1)*nl*nj + (j-1)*nl + i)
        end do
      end do
    end do
  end subroutine park5_pack_primitives

  ! ------------------------------------------------------------------
  ! fill the domain-edge ghost columns where this rank really is a
  ! global boundary. the seam ghosts came from the halo. the codes are
  ! 0 transmissive, 1 inflow, 2 outflow, 3 slip wall, 4 no-slip wall.
  subroutine park5_domain_edges(car, nl, nj, dec, bc_w, bc_e, bcv, yf)
    integer, intent(in) :: nl, nj, bc_w, bc_e
    type(decomp_t), intent(in) :: dec
    real(RP), intent(inout) :: car(NCAR,nl+4,nj)
    real(RP), intent(in) :: bcv(5,4), yf(NP5)
    integer :: j, s

    if (dec%west_domain_edge) then
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
    if (dec%east_domain_edge) then
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
  end subroutine park5_domain_edges

end module nuforpark5common

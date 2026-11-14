! nuforboundary.f90 — outer-domain ghost-column fill for one conservative
! or primitive variable on the west/east domain edge of a padded strip.
module nuforboundary
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: DP = real64

contains


  ! apply boundary ghost cells for one padded strip at the domain edges.
  ! bc codes: 0=transmissive, 1=supersonic_inflow(rho,u,v,p), 2=supersonic_outflow,
  ! 3=slip_wall, 4=no_slip_wall. inflow carries (rho,u,v,p) in bc_vals.
  ! vars: 1=rho, 2=u, 3=v, 4=p (padded strips are primitive).
  subroutine bc_ghost(var, bc_code, bc_vals, n, padded, side)
    integer, intent(in) :: var, bc_code, n, side  ! side: 0=west, 1=east
    real(DP), intent(in) :: bc_vals(5)
    real(DP), intent(inout) :: padded(n+2)
    ! fills both ghost layers at the domain edge (the muscl limiter at the
    ! first interior face reads two cells beyond it)
    if (side == 0) then
      select case (bc_code)
      case (0)   ! transmissive: mirror the interior into both ghosts
        padded(1) = padded(3)
        padded(2) = padded(3)
      case (1)   ! supersonic inflow: fixed freestream primitive
        padded(1) = bc_vals(var)
        padded(2) = bc_vals(var)
      case (2)   ! supersonic outflow
        padded(1) = padded(3)
        padded(2) = padded(3)
      case (3)   ! slip wall: mirror normal velocity about the wall
        if (var == 2) then
          padded(1) = -padded(3)
          padded(2) = -padded(3)
        else
          padded(1) = padded(3)
          padded(2) = padded(3)
        end if
      case (4)   ! no-slip wall: both velocities reversed
        if (var == 2 .or. var == 3) then
          padded(1) = -padded(3)
          padded(2) = -padded(3)
        else
          padded(1) = padded(3)
          padded(2) = padded(3)
        end if
      end select
    else
      select case (bc_code)
      case (0)
        padded(n+1) = padded(n)
        padded(n+2) = padded(n)
      case (1)
        padded(n+1) = bc_vals(var)
        padded(n+2) = bc_vals(var)
      case (2)
        padded(n+1) = padded(n)
        padded(n+2) = padded(n)
      case (3)
        if (var == 2) then
          padded(n+1) = -padded(n)
          padded(n+2) = -padded(n)
        else
          padded(n+1) = padded(n)
          padded(n+2) = padded(n)
        end if
      case (4)
        if (var == 2 .or. var == 3) then
          padded(n+1) = -padded(n)
          padded(n+2) = -padded(n)
        else
          padded(n+1) = padded(n)
          padded(n+2) = padded(n)
        end if
      end select
    end if
  end subroutine

end module nuforboundary

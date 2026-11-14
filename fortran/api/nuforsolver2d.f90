! nuforsolver2d.f90 — the mpi-parallel 2d euler march driver: rk2 heun over
! 1d x-slabs with halo exchange, assembling per-row/per-column flux slices from
! nuforfluxes and applying them across the strip.
!
! discipline notes:
! - state layout matches rust exactly: flat row-major (j*nx + i), one array
!   per conserved variable, ghost columns bookended by the halo exchange.
! - the np=1 path runs this same code (one rank, no exchange) and is the
!   bit-identity anchor for the np>1 gates.
! - the tgas1 block coefficients and cea table are shipped once at init from
!   rust (which owns the parsed data), held by nufoeqair.
module nuforsolver2d
  use iso_fortran_env, only: real64
  use nufor_mpi, only: decomp_t, global_max
  use nuforfaces, only: MAXN, face_l, face_r, face_e_l, face_e_r, face_a_l, face_a_r, face_p_l, face_p_r, slope_y
  use nufoeqair, only: cea_pt_at, tgas1_pressure_at, tgas1_sound_at
  use nuforriemann, only: hllc_face_x
  use nuforfluxes, only: x_row_fluxes
  implicit none
  integer, parameter :: DP = real64

contains


  subroutine euler2d_step_rk2(rho_in, mx_in, my_in, e_in, nx_l, ny_, &
                              dx, dy, gamma, cfl, dt_cap, d, &
                              bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                              rho_out, mx_out, my_out, e_out, dt_out, steps_flag)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho_in(nx_l+4, ny_), mx_in(nx_l+4, ny_), &
                            my_in(nx_l+4, ny_), e_in(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, gamma, cfl, dt_cap
    type(decomp_t), intent(in) :: d
    integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(DP), intent(in) :: bc_vals(5, 4)
    logical, intent(in) :: solid(nx_l+4, ny_)
    real(DP), intent(out) :: rho_out(nx_l+4, ny_), mx_out(nx_l+4, ny_), &
                             my_out(nx_l+4, ny_), e_out(nx_l+4, ny_)
    real(DP), intent(out) :: dt_out
    integer, intent(in) :: steps_flag
    real(DP), allocatable :: u(:,:), v(:,:), p(:,:), a(:,:), ein(:,:)
    real(DP), allocatable :: fx(:,:,:), fy(:,:,:)
    real(DP) :: smax, s_local, dt, ge_ff
    real(DP) :: fm, fmx2, fmy2, fe, dm, dp_
    real(DP) :: pad_y(MAXN)
    integer :: i, j, f, var, ierr, nzf

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), p(nx_l+4, ny_), a(nx_l+4, ny_))
    allocate(ein(nx_l+4, ny_))
    allocate(fx(4, nx_l+2, ny_), fy(4, nx_l, ny_+1))

    ! primitive conversion over the whole padded slab; ein is the internal
    ! energy per unit mass (the eqair closure variable the rust march
    ! reconstructs on faces)
    do j = 1, ny_
      do i = 1, nx_l+4
        u(i,j) = mx_in(i,j) / max(rho_in(i,j), 1.0e-12_DP)
        v(i,j) = my_in(i,j) / max(rho_in(i,j), 1.0e-12_DP)
        ein(i,j) = e_in(i,j) / max(rho_in(i,j), 1.0e-12_DP) &
                 - 0.5_DP * (u(i,j)**2 + v(i,j)**2)
        call cea_pt_at(rho_in(i,j), ein(i,j), p(i,j))
        p(i,j) = max(p(i,j), 1.0e-12_DP)
        ge_ff = 1.0_DP + p(i,j) / (rho_in(i,j) * max(ein(i,j), 1.0e-3_DP))
        ge_ff = min(max(ge_ff, 1.05_DP), 1.67_DP)
        a(i,j) = sqrt(ge_ff * p(i,j) / max(rho_in(i,j), 1.0e-12_DP))
        a(i,j) = max(a(i,j), 50.0_DP)
      end do
    end do

    ! ---- x sweep per row: muscl over the padded strip, faces 2..nx_l+2 ----
    do j = 1, ny_
      call x_row_fluxes(j, nx_l, ny_, gamma, rho_in, ein, u, v, p, a, &
                        d, bc_w, bc_e, bc_vals, fx(1:4, 1:nx_l+2, j))
    end do

    ! ---- y sweep per column (y undecomposed in v1) ----
    do i = 3, nx_l+2
      do var = 1, 5
        do j = 1, ny_
          select case (var)
          case (1); pad_y(j+1) = rho_in(i,j)
          case (2); pad_y(j+1) = u(i,j)
          case (3); pad_y(j+1) = v(i,j)
          case (4); pad_y(j+1) = p(i,j)
          case (5); pad_y(j+1) = ein(i,j)
          end select
        end do
        ! domain-edge bcs (every rank holds the full y extent); the y strip
        ! carries ONE ghost slot per side (pad_y(1) south, pad_y(ny_+2)
        ! north), so fill them directly instead of the two-layer bc_ghost,
        ! which would clobber the first and last interior cells. var 5
        ! (ein) mirrors like a scalar in every bc case.
        select case (bc_s)
        case (1)   ! supersonic inflow: freestream primitive (ein preinverted)
          if (var == 5) then
            pad_y(1) = bc_vals(5,3)
          else
            pad_y(1) = bc_vals(var,3)
          end if
        case (3)   ! slip wall: reflect the y-normal velocity v
          if (var == 3) then
            pad_y(1) = -pad_y(2)
          else
            pad_y(1) = pad_y(2)
          end if
        case (4)   ! no-slip wall: both velocities reversed
          if (var == 2 .or. var == 3) then
            pad_y(1) = -pad_y(2)
          else
            pad_y(1) = pad_y(2)
          end if
        case default  ! transmissive / supersonic outflow: mirror the edge cell
          pad_y(1) = pad_y(2)
        end select
        select case (bc_n)
        case (1)
          if (var == 5) then
            pad_y(ny_+2) = bc_vals(5,4)
          else
            pad_y(ny_+2) = bc_vals(var,4)
          end if
        case (3)
          if (var == 3) then
            pad_y(ny_+2) = -pad_y(ny_+1)
          else
            pad_y(ny_+2) = pad_y(ny_+1)
          end if
        case (4)
          if (var == 2 .or. var == 3) then
            pad_y(ny_+2) = -pad_y(ny_+1)
          else
            pad_y(ny_+2) = pad_y(ny_+1)
          end if
        case default
          pad_y(ny_+2) = pad_y(ny_+1)
        end select
        ! van leer slopes: slope_y(f) is centered on cell f (slots f..f+2).
        ! face f separates cells (f-1) | f, mirroring the x-sweep's
        ! convention (face fx(i-1) is cell i's +x face): face 1 is the
        ! south domain-edge face (ghost | cell 1), face ny_+1 the north
        ! edge (cell ny_ | ghost).
        do f = 1, ny_
          dm = pad_y(f+1) - pad_y(f)
          dp_ = pad_y(f+2) - pad_y(f+1)
          if (dm * dp_ > 0.0_DP) then
            slope_y(f) = 2.0_DP * dm * dp_ / (dm + dp_)
          else
            slope_y(f) = 0.0_DP
          end if
        end do
        face_l(var, 1) = pad_y(1)
        face_r(var, 1) = pad_y(2) - 0.5_DP * slope_y(1)
        do f = 2, ny_+1
          face_l(var, f) = pad_y(f) + 0.5_DP * slope_y(f-1)
          face_r(var, f) = pad_y(f+1) - 0.5_DP * slope_y(f)
        end do
        face_r(var, ny_+1) = pad_y(ny_+2)
      end do
      do f = 1, ny_+1
        face_e_l(f) = face_l(5,f)
        face_e_r(f) = face_r(5,f)
        face_p_l(f) = tgas1_pressure_at(face_l(1,f), face_e_l(f))
        if (face_p_l(f) < 0.0_DP .or. face_p_l(f) /= face_p_l(f)) face_p_l(f) = 1.0e-12_DP
        face_p_r(f) = tgas1_pressure_at(face_r(1,f), face_e_r(f))
        if (face_p_r(f) < 0.0_DP .or. face_p_r(f) /= face_p_r(f)) face_p_r(f) = 1.0e-12_DP
        face_a_l(f) = tgas1_sound_at(face_l(1,f), face_e_l(f))
        face_a_r(f) = tgas1_sound_at(face_r(1,f), face_e_r(f))
        ! rust's f64::max(NaN, 1e-12) = 1e-12; floor NaN explicitly here
        if (face_a_l(f) /= face_a_l(f)) face_a_l(f) = 1.0e-12_DP
        if (face_a_r(f) /= face_a_r(f)) face_a_r(f) = 1.0e-12_DP
      end do
      do f = 1, ny_+1
        ! axis-1: pass normal momentum as arg2 (u-slot = v), tangential as arg3
        call hllc_face_x(face_l(1,f), face_l(3,f), face_l(2,f), face_p_l(f), &
                         face_e_l(f), face_a_l(f), &
                         face_r(1,f), face_r(3,f), face_r(2,f), face_p_r(f), &
                         face_e_r(f), face_a_r(f), gamma, fm, fmy2, fmx2, fe)
        fy(1, i-2, f) = fm; fy(2, i-2, f) = fmx2; fy(3, i-2, f) = fmy2; fy(4, i-2, f) = fe
      end do
    end do

    ! ---- global cfl: max(|u|,|v|) + a with the cea gamma_eff a field ----
    ! slip-wall enforcement: zero the flux across any face touching a
    ! solid cell, matching the rust march's wall closure (the riemann
    ! solve never runs across the wall). face fx(:,f,j) sits between
    ! padded cells f+1 and f+2; face fy(:,f,j) separates padded y-cells
    ! j and j+1 at x-slot f+2 (face 1 = ghost | cell 1).
    nzf = 0
    do j = 1, ny_
      do f = 1, nx_l+2
        if (solid(f+1,j) .or. solid(f+2,j)) then
          ! the mirrored slip-wall riemann, same as the rust wall
          ! closure: reflect the fluid state about the face and run the
          ! same hllc the interior uses. the solve returns the wall
          ! pressure the reflection builds, so the wall is impermeable
          ! but not force-free.
          ! the rust wall closure puts the mirror on the same side as the
          ! solid: left-solid runs (mirror, fluid), right-solid runs
          ! (fluid, mirror).
          if (solid(f+1,j)) then
            call hllc_face_x(rho_in(f+2,j), -u(f+2,j), -v(f+2,j), p(f+2,j), &
                             ein(f+2,j), a(f+2,j), &
                             rho_in(f+2,j), u(f+2,j), v(f+2,j), p(f+2,j), &
                             ein(f+2,j), a(f+2,j), &
                             gamma, fx(1,f,j), fx(2,f,j), fx(3,f,j), fx(4,f,j))
          else
            call hllc_face_x(rho_in(f+1,j), u(f+1,j), v(f+1,j), p(f+1,j), &
                             ein(f+1,j), a(f+1,j), &
                             rho_in(f+1,j), -u(f+1,j), -v(f+1,j), p(f+1,j), &
                             ein(f+1,j), a(f+1,j), &
                             gamma, fx(1,f,j), fx(2,f,j), fx(3,f,j), fx(4,f,j))
          end if
          nzf = nzf + 1
        end if
      end do
    end do
    do j = 1, ny_+1
      do f = 1, nx_l
        ! face j separates y-cells (j-1) | j (1-based); face 1 has a
        ! ghost below and face ny_+1 a ghost above, so only test the
        ! cell indices that exist (no out-of-bounds solid reads).
        if ((j > 1 .and. solid(f+2,j-1)) .or. (j <= ny_ .and. solid(f+2,j))) then
          ! mirrored slip-wall riemann on the y face: the u-slot takes
          ! the y-normal velocity v, the v-slot the tangential u, and
          ! the fluid-side normal velocity is reflected.
          if (j > 1 .and. solid(f+2,j-1)) then
            call hllc_face_x(rho_in(f+2,j), -v(f+2,j), -u(f+2,j), p(f+2,j), &
                             ein(f+2,j), a(f+2,j), &
                             rho_in(f+2,j), v(f+2,j), u(f+2,j), p(f+2,j), &
                             ein(f+2,j), a(f+2,j), &
                             gamma, fm, fmx2, fmy2, fe)
          else
            call hllc_face_x(rho_in(f+2,j-1), v(f+2,j-1), u(f+2,j-1), p(f+2,j-1), &
                             ein(f+2,j-1), a(f+2,j-1), &
                             rho_in(f+2,j-1), -v(f+2,j-1), -u(f+2,j-1), p(f+2,j-1), &
                             ein(f+2,j-1), a(f+2,j-1), &
                             gamma, fm, fmx2, fmy2, fe)
          end if
          ! hllc_face_x works in the face frame: fmx2 is the normal
          ! (y) momentum flux and fmy2 the tangential (x) momentum flux,
          ! so they land in swapped fy slots.
          fy(1,f,j) = fm
          fy(2,f,j) = fmy2
          fy(3,f,j) = fmx2
          fy(4,f,j) = fe
          nzf = nzf + 1
        end if
      end do
    end do
    s_local = 0.0_DP
    do j = 1, ny_
      do i = 3, nx_l+2
        s_local = max(s_local, max(abs(u(i,j)), abs(v(i,j))) + a(i,j))
      end do
    end do
    call global_max(s_local, smax, d)
    dt = min(cfl * min(dx, dy) / max(smax, 1.0e-12_DP), dt_cap)
    dt_out = dt

    ! ---- conservative update on interior cells (3..nx_l+2) ----
    do j = 1, ny_
      do i = 3, nx_l+2
        if (solid(i,j)) cycle
        rho_out(i,j) = rho_in(i,j) - dt * ( &
              (fx(1, i-1, j) - fx(1, i-2, j)) / dx &
            + (fy(1, i-2, j+1) - fy(1, i-2, j)) / dy )
        mx_out(i,j) = mx_in(i,j) - dt * ( &
              (fx(2, i-1, j) - fx(2, i-2, j)) / dx &
            + (fy(2, i-2, j+1) - fy(2, i-2, j)) / dy )
        my_out(i,j) = my_in(i,j) - dt * ( &
              (fx(3, i-1, j) - fx(3, i-2, j)) / dx &
            + (fy(3, i-2, j+1) - fy(3, i-2, j)) / dy )
        e_out(i,j) = e_in(i,j) - dt * ( &
              (fx(4, i-1, j) - fx(4, i-2, j)) / dx &
            + (fy(4, i-2, j+1) - fy(4, i-2, j)) / dy )
      end do
    end do
    ! solid cells + ghosts copy through
    do j = 1, ny_
      do i = 1, nx_l+4
        if (solid(i,j) .or. i < 3 .or. i > nx_l+2) then
          rho_out(i,j) = rho_in(i,j)
          mx_out(i,j) = mx_in(i,j)
          my_out(i,j) = my_in(i,j)
          e_out(i,j) = e_in(i,j)
        end if
      end do
    end do
    deallocate(u, v, p, a, ein, fx, fy)
  end subroutine


  ! primitives (u, v) and internal energy per mass from the padded
  ! conserved slab, mirroring cons_to_prim2d with the same floors.
  subroutine primitive_ein(rho, mx, my, e, nx_l, ny_, u, v, ein)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                            my(nx_l+4, ny_), e(nx_l+4, ny_)
    real(DP), intent(out) :: u(nx_l+4, ny_), v(nx_l+4, ny_), ein(nx_l+4, ny_)
    integer :: i, j
    do j = 1, ny_
      do i = 1, nx_l+4
        u(i,j) = mx(i,j) / max(rho(i,j), 1.0e-12_DP)
        v(i,j) = my(i,j) / max(rho(i,j), 1.0e-12_DP)
        ein(i,j) = e(i,j) / max(rho(i,j), 1.0e-12_DP) &
                 - 0.5_DP * (u(i,j)**2 + v(i,j)**2)
      end do
    end do
  end subroutine

end module nuforsolver2d

! nuforpg2d.f90 — perfect-gas 2d euler heun kernel on one rank's padded
! x-slab, mirroring the rust planar model march for ThermoModel::Perfect:
! pressure (not internal energy) is the muscl-reconstructed closure
! variable, cell p = (gamma-1)*(e - 0.5*(mx^2+my^2)/rho), face sound
! a = sqrt(gamma*p/rho) and the hllc floors live in nuforriemann.
module nuforpg2d
  use iso_fortran_env, only: real64
  use nufor_mpi, only: decomp_t, global_max
  use nuforhalos, only: halo_exchange_cons
  use nuforfaces, only: MAXN
  use nuforriemann, only: hllc_face_x
  use nuforboundary, only: bc_ghost
  implicit none
  integer, parameter :: DP = real64

contains

  ! primitives and cell pressure/sound over the padded slab, perfect gas.
  subroutine pg_primitives(rho, mx, my, e, nx_l, ny_, gamma, u, v, p, a)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                            my(nx_l+4, ny_), e(nx_l+4, ny_)
    real(DP), intent(in) :: gamma
    real(DP), intent(out) :: u(nx_l+4, ny_), v(nx_l+4, ny_), &
                             p(nx_l+4, ny_), a(nx_l+4, ny_)
    real(DP) :: rc
    integer :: i, j
    do j = 1, ny_
      do i = 1, nx_l+4
        rc = max(rho(i,j), 1.0e-12_DP)
        u(i,j) = mx(i,j) / rc
        v(i,j) = my(i,j) / rc
        p(i,j) = (gamma - 1.0_DP) * (e(i,j) - 0.5_DP * (mx(i,j)**2 + my(i,j)**2) / rc)
        a(i,j) = sqrt(gamma * max(p(i,j), 1.0e-12_DP) / rc)
      end do
    end do
  end subroutine pg_primitives

  ! one forward-euler stage on the padded slab: muscl (van leer) on
  ! rho/u/v/p, hllc per face, conservative update of the interior.
  subroutine pg2d_stage(rho_in, mx_in, my_in, e_in, nx_l, ny_, &
                        dx, dy, gamma, dt, d, &
                        bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                        rho_out, mx_out, my_out, e_out)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho_in(nx_l+4, ny_), mx_in(nx_l+4, ny_), &
                            my_in(nx_l+4, ny_), e_in(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, gamma, dt
    type(decomp_t), intent(in) :: d
    integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(DP), intent(in) :: bc_vals(5, 4)
    logical, intent(in) :: solid(nx_l+4, ny_)
    real(DP), intent(out) :: rho_out(nx_l+4, ny_), mx_out(nx_l+4, ny_), &
                             my_out(nx_l+4, ny_), e_out(nx_l+4, ny_)
    real(DP), allocatable :: u(:,:), v(:,:), p(:,:), a(:,:)
    real(DP), allocatable :: fx(:,:,:), fy(:,:,:)
    real(DP) :: fx_rho_l(nx_l+2,ny_), fx_rho_r(nx_l+2,ny_), fx_u_l(nx_l+2,ny_), fx_u_r(nx_l+2,ny_)
    real(DP) :: fx_v_l(nx_l+2,ny_), fx_v_r(nx_l+2,ny_), fx_p_l(nx_l+2,ny_), fx_p_r(nx_l+2,ny_)
    real(DP) :: fy_rho_l(nx_l,ny_+1), fy_rho_r(nx_l,ny_+1), fy_u_l(nx_l,ny_+1), fy_u_r(nx_l,ny_+1)
    real(DP) :: fy_v_l(nx_l,ny_+1), fy_v_r(nx_l,ny_+1), fy_p_l(nx_l,ny_+1), fy_p_r(nx_l,ny_+1)
    real(DP) :: pad(MAXN), sl(MAXN), fl(MAXN), fr(MAXN)
    real(DP) :: al, ar, fm, fmx, fmy, fe
    integer :: i, j, f, var

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), p(nx_l+4, ny_), a(nx_l+4, ny_))
    allocate(fx(4, nx_l+2, ny_), fy(4, nx_l, ny_+1))
    call pg_primitives(rho_in, mx_in, my_in, e_in, nx_l, ny_, gamma, u, v, p, a)

    ! x sweep per row: muscl on rho,u,v,p over the padded strip
    do j = 1, ny_
      do var = 1, 4
        do i = 1, nx_l+4
          select case (var)
          case (1); pad(i) = rho_in(i,j)
          case (2); pad(i) = u(i,j)
          case (3); pad(i) = v(i,j)
          case (4); pad(i) = p(i,j)
          end select
        end do
        if (d%west_domain_edge) then
          call bc_ghost(var, bc_w, bc_vals(1,1), nx_l+2, pad, 0)
        end if
        if (d%east_domain_edge) then
          call bc_ghost(var, bc_e, bc_vals(1,2), nx_l+2, pad, 1)
        end if
        do i = 1, nx_l+2
          sl(i) = vanleer(pad(i+1) - pad(i), pad(i+2) - pad(i+1))
        end do
        if (d%west_domain_edge) then
          fl(1) = pad(2)
        else
          fl(1) = pad(2) + 0.5_DP * sl(1)
        end if
        fr(1) = pad(3) - 0.5_DP * sl(2)
        do f = 2, nx_l
          fl(f) = pad(f+1) + 0.5_DP * sl(f)
          fr(f) = pad(f+2) - 0.5_DP * sl(f+1)
        end do
        fl(nx_l+1) = pad(nx_l+2) + 0.5_DP * sl(nx_l+1)
        if (d%east_domain_edge) then
          fr(nx_l+1) = pad(nx_l+3)
        else
          fr(nx_l+1) = pad(nx_l+3) - 0.5_DP * sl(nx_l+2)
        end if
        ! stash this var's face pair aside per face below
        do f = 1, nx_l+1
          select case (var)
          case (1); fx_rho_l(f,j) = fl(f); fx_rho_r(f,j) = fr(f)
          case (2); fx_u_l(f,j) = fl(f);   fx_u_r(f,j) = fr(f)
          case (3); fx_v_l(f,j) = fl(f);   fx_v_r(f,j) = fr(f)
          case (4); fx_p_l(f,j) = fl(f);   fx_p_r(f,j) = fr(f)
          end select
        end do
      end do
      do f = 1, nx_l+1
        if (solid(f+1,j) .or. solid(f+2,j)) then
          if (solid(f+1,j)) then
            call hllc_face_x(rho_in(f+2,j), -u(f+2,j), -v(f+2,j), p(f+2,j), &
                             ein_cell(e_in(f+2,j), mx_in(f+2,j), my_in(f+2,j), rho_in(f+2,j)), a(f+2,j), &
                             rho_in(f+2,j), u(f+2,j), v(f+2,j), p(f+2,j), &
                             ein_cell(e_in(f+2,j), mx_in(f+2,j), my_in(f+2,j), rho_in(f+2,j)), a(f+2,j), gamma, fx(1,f,j), fx(2,f,j), fx(3,f,j), fx(4,f,j))
          else
            call hllc_face_x(rho_in(f+1,j), u(f+1,j), v(f+1,j), p(f+1,j), &
                             ein_cell(e_in(f+1,j), mx_in(f+1,j), my_in(f+1,j), rho_in(f+1,j)), a(f+1,j), &
                             rho_in(f+1,j), -u(f+1,j), -v(f+1,j), p(f+1,j), &
                             ein_cell(e_in(f+1,j), mx_in(f+1,j), my_in(f+1,j), rho_in(f+1,j)), a(f+1,j), gamma, fx(1,f,j), fx(2,f,j), fx(3,f,j), fx(4,f,j))
          end if
        else
          al = sqrt(gamma * max(fx_p_l(f,j), 1.0e-12_DP) / max(fx_rho_l(f,j), 1.0e-12_DP))
          ar = sqrt(gamma * max(fx_p_r(f,j), 1.0e-12_DP) / max(fx_rho_r(f,j), 1.0e-12_DP))
          call hllc_face_x(fx_rho_l(f,j), fx_u_l(f,j), fx_v_l(f,j), fx_p_l(f,j), &
                           fx_p_l(f,j)/(fx_rho_l(f,j)*(gamma-1.0_DP)), al, &
                           fx_rho_r(f,j), fx_u_r(f,j), fx_v_r(f,j), fx_p_r(f,j), &
                           fx_p_r(f,j)/(fx_rho_r(f,j)*(gamma-1.0_DP)), ar, gamma, fm, fmx, fmy, fe)
          fx(1,f,j) = fm; fx(2,f,j) = fmx; fx(3,f,j) = fmy; fx(4,f,j) = fe
        end if
      end do
    end do

    ! y sweep per interior column: one ghost slot per side
    do i = 3, nx_l+2
      do var = 1, 4
        do j = 1, ny_
          select case (var)
          case (1); pad(j+1) = rho_in(i,j)
          case (2); pad(j+1) = u(i,j)
          case (3); pad(j+1) = v(i,j)
          case (4); pad(j+1) = p(i,j)
          end select
        end do
        select case (bc_s)
        case (1); pad(1) = bc_vals(var,3)
        case (3)
          if (var == 3) then; pad(1) = -pad(2); else; pad(1) = pad(2); end if
        case (4)
          if (var == 2 .or. var == 3) then; pad(1) = -pad(2); else; pad(1) = pad(2); end if
        case default; pad(1) = pad(2)
        end select
        select case (bc_n)
        case (1); pad(ny_+2) = bc_vals(var,4)
        case (3)
          if (var == 3) then; pad(ny_+2) = -pad(ny_+1); else; pad(ny_+2) = pad(ny_+1); end if
        case (4)
          if (var == 2 .or. var == 3) then; pad(ny_+2) = -pad(ny_+1); else; pad(ny_+2) = pad(ny_+1); end if
        case default; pad(ny_+2) = pad(ny_+1)
        end select
        do f = 1, ny_
          sl(f) = vanleer(pad(f+1) - pad(f), pad(f+2) - pad(f+1))
        end do
        fl(1) = pad(1)
        fr(1) = pad(2) - 0.5_DP * sl(1)
        do f = 2, ny_+1
          fl(f) = pad(f) + 0.5_DP * sl(f-1)
          fr(f) = pad(f+1) - 0.5_DP * sl(f)
        end do
        fr(ny_+1) = pad(ny_+2)
        do f = 1, ny_+1
          select case (var)
          case (1); fy_rho_l(i-2,f) = fl(f); fy_rho_r(i-2,f) = fr(f)
          case (2); fy_u_l(i-2,f) = fl(f);   fy_u_r(i-2,f) = fr(f)
          case (3); fy_v_l(i-2,f) = fl(f);   fy_v_r(i-2,f) = fr(f)
          case (4); fy_p_l(i-2,f) = fl(f);   fy_p_r(i-2,f) = fr(f)
          end select
        end do
      end do
      do f = 1, ny_+1
        if ((f > 1 .and. solid(i,f-1)) .or. (f <= ny_ .and. solid(i,f))) then
          if (f > 1 .and. solid(i,f-1)) then
            call hllc_face_x(rho_in(i,f), -v(i,f), -u(i,f), p(i,f), &
                             ein_cell(e_in(i,f), mx_in(i,f), my_in(i,f), rho_in(i,f)), a(i,f), &
                             rho_in(i,f), v(i,f), u(i,f), p(i,f), &
                             ein_cell(e_in(i,f), mx_in(i,f), my_in(i,f), rho_in(i,f)), a(i,f), gamma, fm, fmx, fmy, fe)
          else
            call hllc_face_x(rho_in(i,f-1), v(i,f-1), u(i,f-1), p(i,f-1), &
                             ein_cell(e_in(i,f-1), mx_in(i,f-1), my_in(i,f-1), rho_in(i,f-1)), a(i,f-1), &
                             rho_in(i,f-1), -v(i,f-1), -u(i,f-1), p(i,f-1), &
                             ein_cell(e_in(i,f-1), mx_in(i,f-1), my_in(i,f-1), rho_in(i,f-1)), a(i,f-1), gamma, fm, fmx, fmy, fe)
          end if
          fy(1,i-2,f) = fm
          fy(2,i-2,f) = fmy
          fy(3,i-2,f) = fmx
          fy(4,i-2,f) = fe
        else
          al = sqrt(gamma * max(fy_p_l(i-2,f), 1.0e-12_DP) / max(fy_rho_l(i-2,f), 1.0e-12_DP))
          ar = sqrt(gamma * max(fy_p_r(i-2,f), 1.0e-12_DP) / max(fy_rho_r(i-2,f), 1.0e-12_DP))
          call hllc_face_x(fy_rho_l(i-2,f), fy_v_l(i-2,f), fy_u_l(i-2,f), fy_p_l(i-2,f), &
                           fy_p_l(i-2,f)/(fy_rho_l(i-2,f)*(gamma-1.0_DP)), al, &
                           fy_rho_r(i-2,f), fy_v_r(i-2,f), fy_u_r(i-2,f), fy_p_r(i-2,f), &
                           fy_p_r(i-2,f)/(fy_rho_r(i-2,f)*(gamma-1.0_DP)), ar, gamma, fm, fmx, fmy, fe)
          fy(1,i-2,f) = fm
          fy(2,i-2,f) = fmy
          fy(3,i-2,f) = fmx
          fy(4,i-2,f) = fe
        end if
      end do
    end do

    ! conservative update on interior cells
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
    deallocate(u, v, p, a, fx, fy)
  end subroutine pg2d_stage


  ! internal energy per mass of a cell, given its padded conserveds
  pure function ein_cell(ev, mxv, myv, rv) result(eint)
    real(DP), intent(in) :: ev, mxv, myv, rv
    real(DP) :: eint, uu, vv
    uu = mxv / rv
    vv = myv / rv
    eint = ev / rv - 0.5_DP * (uu*uu + vv*vv)
  end function ein_cell

  ! the van leer limiter, matching rust axi::van_leer
  pure function vanleer(dm, dp_) result(s)
    real(DP), intent(in) :: dm, dp_
    real(DP) :: s
    if (dm * dp_ <= 0.0_DP) then
      s = 0.0_DP
    else
      s = 2.0_DP * dm * dp_ / (dm + dp_)
    end if
  end function vanleer

  ! two-stage heun on the padded slab: stage1 -> halo -> stage2 -> average,
  ! mirroring rust axi::advance2d_model_rk2 and the eqair bridge.
  subroutine pg2d_step_rk2(rho_in, mx_in, my_in, e_in, nx_l, ny_, &
                           dx, dy, gamma, cfl, dt_cap, d, &
                           bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                           rho_out, mx_out, my_out, e_out, dt_out)
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
    real(DP), allocatable :: u(:,:), v(:,:), p(:,:), a(:,:)
    real(DP) :: r1(nx_l+4, ny_), m1(nx_l+4, ny_), my1(nx_l+4, ny_), e1(nx_l+4, ny_)
    real(DP) :: r2(nx_l+4, ny_), m2(nx_l+4, ny_), my2(nx_l+4, ny_), e2(nx_l+4, ny_)
    real(DP) :: s_local, smax, dt
    integer :: i, j

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), p(nx_l+4, ny_), a(nx_l+4, ny_))
    call pg_primitives(rho_in, mx_in, my_in, e_in, nx_l, ny_, gamma, u, v, p, a)
    s_local = 0.0_DP
    do j = 1, ny_
      do i = 3, nx_l+2
        s_local = max(s_local, max(abs(u(i,j)), abs(v(i,j))) + a(i,j))
      end do
    end do
    call global_max(s_local, smax, d)
    dt = min(cfl * min(dx, dy) / max(smax, 1.0e-12_DP), dt_cap)
    dt_out = dt

    call pg2d_stage(rho_in, mx_in, my_in, e_in, nx_l, ny_, &
                    dx, dy, gamma, dt, d, bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                    r1, m1, my1, e1)
    call halo_exchange_cons(r1, m1, my1, e1, nx_l, ny_)
    ! the rust lane recomputes the timestep inside every advance2d_model
    ! call, so the second stage runs on its own dt taken from the
    ! stage-one state.
    call pg_primitives(r1, m1, my1, e1, nx_l, ny_, gamma, u, v, p, a)
    s_local = 0.0_DP
    do j = 1, ny_
      do i = 3, nx_l+2
        s_local = max(s_local, max(abs(u(i,j)), abs(v(i,j))) + a(i,j))
      end do
    end do
    call global_max(s_local, smax, d)
    dt = min(cfl * min(dx, dy) / max(smax, 1.0e-12_DP), dt_cap)
    call pg2d_stage(r1, m1, my1, e1, nx_l, ny_, &
                    dx, dy, gamma, dt, d, bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                    r2, m2, my2, e2)
    do j = 1, ny_
      do i = 1, nx_l+4
        rho_out(i,j) = 0.5_DP * rho_in(i,j) + 0.5_DP * r2(i,j)
        mx_out(i,j)  = 0.5_DP * mx_in(i,j)  + 0.5_DP * m2(i,j)
        my_out(i,j)  = 0.5_DP * my_in(i,j)  + 0.5_DP * my2(i,j)
        e_out(i,j)   = 0.5_DP * e_in(i,j)   + 0.5_DP * e2(i,j)
      end do
    end do
    deallocate(u, v, p, a)
  end subroutine pg2d_step_rk2

end module nuforpg2d

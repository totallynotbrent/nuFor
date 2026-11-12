! nuforcartpark5.f90 — the cartesian park5 species-carrying 2d march.
!
! this is the serial/cartesian counterpart of the mpi eqair lane in
! nuforsolver2d.f90, specialised for the five-species mixture: every
! face derives its pressure and sound speed from the park5 mixture
! closure at the muscl-reconstructed face composition (mirroring
! axi.rs face_pair for ThermoModel::Park5), the species are advected
! by the same hllc face mass fluxes the flow advanced under, and the
! operator-split chemistry runs after the transport update.
!
! padding convention: cell i (1-based) lives at pad index i+2, ghosts
! at 2 and nx+3, and one sacrificial slab at 1 and nx+4 so centred
! van-leer slopes never read out of bounds. face f (1..nx+1) sits
! between centres f+1 and f+2, i.e. cell f-1 and cell f.
module nuforcartpark5
  use iso_fortran_env, only: real64
  use nuforpark5, only: closure_triple, transport_planar, chem_substep_grid
  use nuforsolver2d, only: hllc_face_x_pub => hllc_face_x
  implicit none

  integer, parameter :: DP = real64
  integer, parameter :: NP5 = 5
  ! the slope scratch serves both sweeps, so it takes the larger extent.
  integer, parameter :: np_slope = 8192
  ! the same for the ghost-padding scratch.
  integer, parameter :: np_pad = 8200

  private
  public :: cart_park5_step, nfor_cart_park5_step

contains

  pure real(DP) function vlim(dm, dp_)
    real(DP), intent(in) :: dm, dp_
    if (dm*dp_ > 0.0_DP) then
      vlim = 2.0_DP*dm*dp_/(dm+dp_)
    else
      vlim = 0.0_DP
    end if
  end function vlim

  ! one forward-euler park5 step on a plain nx*ny cartesian grid.
  ! arrays are flat with k = (j-1)*nx + i (i 0-based externally).
  ! bcs: 0 transmissive, 1 supersonic inflow (bc_vals), 2 outflow,
  ! 3 slip wall, 4 no-slip wall. y_fill is the fill composition at
  ! inflow edges.
  subroutine cart_park5_step(nx, ny, rho, mx, my, e, y, dxs, dys, cfl, &
                             bc_w, bc_e, bc_s, bc_n, bc_vals, y_fill, &
                             rho2, mx2, my2, e2, y2, dt_out)
    integer, intent(in) :: nx, ny
    real(DP), intent(in) :: rho(nx*ny), mx(nx*ny), my(nx*ny), e(nx*ny)
    real(DP), intent(in) :: y(NP5, nx*ny)
    real(DP), intent(in) :: dxs(nx), dys(ny), cfl
    integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(DP), intent(in) :: bc_vals(5, 4)   ! (var, side), side 1=w 2=e 3=s 4=n
    real(DP), intent(in) :: y_fill(NP5)
    real(DP), intent(out) :: rho2(nx*ny), mx2(nx*ny), my2(nx*ny), e2(nx*ny)
    real(DP), intent(inout) :: y2(NP5, nx*ny)
    real(DP), intent(out) :: dt_out

    real(DP) :: u(nx*ny), v(nx*ny), ei(nx*ny), pp(nx*ny), aa(nx*ny)
    real(DP) :: fx4(4, nx+1, ny), fy4(4, ny+1, nx)
    real(DP) :: fxm(nx+1, ny), fym(ny+1, nx)
    ! the pads live on the heap: gfortran promotes large automatic
    ! arrays to static storage, and a shared static would let
    ! concurrent marches overwrite each other's faces.
    real(DP), allocatable :: thr(:), tu(:), tv(:), te(:), ty(:,:)
    real(DP), allocatable :: sxr(:), sxu(:), sxv(:), sxe(:), sxy(:,:)
    real(DP) :: yl(NP5), yr(NP5), ycc(NP5), sum_y
    real(DP) :: rl, rul, vul, el, al, pl, rr, rur, vur, er_, ar_, pr_
    real(DP) :: fm, fmx, fmy, fe, smax, dt, dum
    real(DP) :: rhe(nx*ny), dee(nx*ny), yc(NP5, nx*ny)
    integer :: i, j, f, s, k, c

    ! the composition advances incrementally, so it starts as a copy.
    y2 = y
    allocate(thr(np_pad), tu(np_pad), tv(np_pad), te(np_pad), ty(NP5, np_pad))
    allocate(sxr(np_slope), sxu(np_slope), sxv(np_slope), sxe(np_slope), sxy(NP5, np_slope))

    ! primitives + per-cell mixture closure
    smax = 0.0_DP
    do k = 1, nx*ny
      u(k) = mx(k) / max(rho(k), 1.0e-12_DP)
      v(k) = my(k) / max(rho(k), 1.0e-12_DP)
      ei(k) = e(k) / max(rho(k), 1.0e-12_DP) - 0.5_DP*(u(k)**2 + v(k)**2)
      do s = 1, NP5
        ycc(s) = y(s, k)
      end do
      call closure_triple(rho(k), ei(k), ycc, pp(k), dum, aa(k))
      pp(k) = max(pp(k), 1.0e-12_DP)
      aa(k) = max(aa(k), 50.0_DP)
      smax = max(smax, max(abs(u(k)), abs(v(k))) + aa(k))
    end do
    dt = cfl * min(minval(dxs), minval(dys)) / max(smax, 1.0e-12_DP)
    dt_out = dt

    ! ===================== x sweep (rows) =====================
    do j = 1, ny
      do i = 1, nx
        k = (j-1)*nx + i
        c = i + 2
        thr(c) = rho(k); tu(c) = u(k); tv(c) = v(k); te(c) = ei(k)
        do s = 1, NP5
          ty(s,c) = y(s,k)
        end do
      end do
      ! west ghosts at c=2 (real) and c=1 (sacrificial copy)
      select case (bc_w)
      case (1)
        thr(2) = bc_vals(1,1); tu(2) = bc_vals(2,1); tv(2) = bc_vals(3,1)
        te(2) = bc_vals(5,1)
        ty(:,2) = y_fill
      case (3)
        thr(2) = thr(3); tu(2) = -tu(3); tv(2) = tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case (4)
        thr(2) = thr(3); tu(2) = -tu(3); tv(2) = -tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case default
        thr(2) = thr(3); tu(2) = tu(3); tv(2) = tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      end select
      thr(1) = thr(2); tu(1) = tu(2); tv(1) = tv(2); te(1) = te(2); ty(:,1) = ty(:,2)
      ! east ghosts at c=nx+3 and c=nx+4
      select case (bc_e)
      case (1)
        thr(nx+3) = bc_vals(1,2); tu(nx+3) = bc_vals(2,2); tv(nx+3) = bc_vals(3,2)
        te(nx+3) = bc_vals(5,2)
        ty(:,nx+3) = y_fill
      case (3)
        thr(nx+3) = thr(nx+2); tu(nx+3) = -tu(nx+2); tv(nx+3) = tv(nx+2)
        te(nx+3) = te(nx+2); ty(:,nx+3) = ty(:,nx+2)
      case (4)
        thr(nx+3) = thr(nx+2); tu(nx+3) = -tu(nx+2); tv(nx+3) = -tv(nx+2)
        te(nx+3) = te(nx+2); ty(:,nx+3) = ty(:,nx+2)
      case default
        thr(nx+3) = thr(nx+2); tu(nx+3) = tu(nx+2); tv(nx+3) = tv(nx+2)
        te(nx+3) = te(nx+2); ty(:,nx+3) = ty(:,nx+2)
      end select
      thr(nx+4) = thr(nx+3); tu(nx+4) = tu(nx+3); tv(nx+4) = tv(nx+3)
      te(nx+4) = te(nx+3); ty(:,nx+4) = ty(:,nx+3)
      ! centred van-leer slopes at real centres c=2..nx+3
      do c = 2, nx+3
        sxr(c-1) = vlim(thr(c)-thr(c-1), thr(c+1)-thr(c))
        sxu(c-1) = vlim(tu(c)-tu(c-1), tu(c+1)-tu(c))
        sxv(c-1) = vlim(tv(c)-tv(c-1), tv(c+1)-tv(c))
        sxe(c-1) = vlim(te(c)-te(c-1), te(c+1)-te(c))
        do s = 1, NP5
          sxy(s,c-1) = vlim(ty(s,c)-ty(s,c-1), ty(s,c+1)-ty(s,c))
        end do
      end do
      ! face f = 1..nx+1 between centres cL=f+1 and cR=f+2
      do f = 1, nx+1
        rl = thr(f+1) + 0.5_DP*sxr(f);  rr = thr(f+2) - 0.5_DP*sxr(f+1)
        rul = tu(f+1) + 0.5_DP*sxu(f);  rur = tu(f+2) - 0.5_DP*sxu(f+1)
        vul = tv(f+1) + 0.5_DP*sxv(f);  vur = tv(f+2) - 0.5_DP*sxv(f+1)
        el = te(f+1) + 0.5_DP*sxe(f);   er_ = te(f+2) - 0.5_DP*sxe(f+1)
        do s = 1, NP5
          yl(s) = ty(s,f+1) + 0.5_DP*sxy(s,f)
          yr(s) = ty(s,f+2) - 0.5_DP*sxy(s,f+1)
        end do
        call closure_triple(max(rl,1.0e-12_DP), el, yl, pl, dum, al)
        call closure_triple(max(rr,1.0e-12_DP), er_, yr, pr_, dum, ar_)
        pl = max(pl, 1.0e-12_DP); pr_ = max(pr_, 1.0e-12_DP)
        al = max(al, 1.0e-12_DP); ar_ = max(ar_, 1.0e-12_DP)
        call hllc_face_x_pub(rl, rul, vul, pl, el, al, rr, rur, vur, pr_, er_, ar_, &
                             1.4_DP, fm, fmx, fmy, fe)
        fx4(1,f,j) = fm; fx4(2,f,j) = fmx; fx4(3,f,j) = fmy; fx4(4,f,j) = fe
        fxm(f,j) = fm
      end do
    end do

    ! ===================== y sweep (columns) =====================
    do i = 1, nx
      do j = 1, ny
        k = (j-1)*nx + i
        c = j + 2
        thr(c) = rho(k); tu(c) = v(k); tv(c) = u(k); te(c) = ei(k)
        do s = 1, NP5
          ty(s,c) = y(s,k)
        end do
      end do
      select case (bc_s)
      case (1)
        thr(2) = bc_vals(1,3); tu(2) = bc_vals(3,3); tv(2) = bc_vals(2,3)
        te(2) = bc_vals(5,3)
        ty(:,2) = y_fill
      case (3)
        thr(2) = thr(3); tu(2) = -tu(3); tv(2) = tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case (4)
        thr(2) = thr(3); tu(2) = -tu(3); tv(2) = -tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case default
        thr(2) = thr(3); tu(2) = tu(3); tv(2) = tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      end select
      thr(1) = thr(2); tu(1) = tu(2); tv(1) = tv(2); te(1) = te(2); ty(:,1) = ty(:,2)
      select case (bc_n)
      case (1)
        thr(ny+3) = bc_vals(1,4); tu(ny+3) = bc_vals(3,4); tv(ny+3) = bc_vals(2,4)
        te(ny+3) = bc_vals(5,4)
        ty(:,ny+3) = y_fill
      case (3)
        thr(ny+3) = thr(ny+2); tu(ny+3) = -tu(ny+2); tv(ny+3) = tv(ny+2)
        te(ny+3) = te(ny+2); ty(:,ny+3) = ty(:,ny+2)
      case (4)
        thr(ny+3) = thr(ny+2); tu(ny+3) = -tu(ny+2); tv(ny+3) = -tv(ny+2)
        te(ny+3) = te(ny+2); ty(:,ny+3) = ty(:,ny+2)
      case default
        thr(ny+3) = thr(ny+2); tu(ny+3) = tu(ny+2); tv(ny+3) = tv(ny+2)
        te(ny+3) = te(ny+2); ty(:,ny+3) = ty(:,ny+2)
      end select
      thr(ny+4) = thr(ny+3); tu(ny+4) = tu(ny+3); tv(ny+4) = tv(ny+3)
      te(ny+4) = te(ny+3); ty(:,ny+4) = ty(:,ny+3)
      do c = 2, ny+3
        sxr(c-1) = vlim(thr(c)-thr(c-1), thr(c+1)-thr(c))
        sxu(c-1) = vlim(tu(c)-tu(c-1), tu(c+1)-tu(c))
        sxv(c-1) = vlim(tv(c)-tv(c-1), tv(c+1)-tv(c))
        sxe(c-1) = vlim(te(c)-te(c-1), te(c+1)-te(c))
        do s = 1, NP5
          sxy(s,c-1) = vlim(ty(s,c)-ty(s,c-1), ty(s,c+1)-ty(s,c))
        end do
      end do
      do f = 1, ny+1
        rl = thr(f+1) + 0.5_DP*sxr(f);  rr = thr(f+2) - 0.5_DP*sxr(f+1)
        rul = tu(f+1) + 0.5_DP*sxu(f);  rur = tu(f+2) - 0.5_DP*sxu(f+1)   ! normal = v
        vul = tv(f+1) + 0.5_DP*sxv(f);  vur = tv(f+2) - 0.5_DP*sxv(f+1)   ! tangent = u
        el = te(f+1) + 0.5_DP*sxe(f);   er_ = te(f+2) - 0.5_DP*sxe(f+1)
        do s = 1, NP5
          yl(s) = ty(s,f+1) + 0.5_DP*sxy(s,f)
          yr(s) = ty(s,f+2) - 0.5_DP*sxy(s,f+1)
        end do
        call closure_triple(max(rl,1.0e-12_DP), el, yl, pl, dum, al)
        call closure_triple(max(rr,1.0e-12_DP), er_, yr, pr_, dum, ar_)
        pl = max(pl, 1.0e-12_DP); pr_ = max(pr_, 1.0e-12_DP)
        al = max(al, 1.0e-12_DP); ar_ = max(ar_, 1.0e-12_DP)
        ! normal momentum into arg2, tangential into arg3, outputs likewise
        call hllc_face_x_pub(rl, rul, vul, pl, el, al, rr, rur, vur, pr_, er_, ar_, &
                             1.4_DP, fm, fmy, fmx, fe)
        fy4(1,f,i) = fm; fy4(2,f,i) = fmx; fy4(3,f,i) = fmy; fy4(4,f,i) = fe
        fym(f,i) = fm
      end do
    end do

    ! ===================== conservative update =====================
    do j = 1, ny
      do i = 1, nx
        k = (j-1)*nx + i
        rho2(k) = rho(k) - dt*((fx4(1,i+1,j)-fx4(1,i,j))/dxs(i) &
                             + (fy4(1,j+1,i)-fy4(1,j,i))/dys(j))
        mx2(k) = mx(k) - dt*((fx4(2,i+1,j)-fx4(2,i,j))/dxs(i) &
                           + (fy4(2,j+1,i)-fy4(2,j,i))/dys(j))
        my2(k) = my(k) - dt*((fx4(3,i+1,j)-fx4(3,i,j))/dxs(i) &
                           + (fy4(3,j+1,i)-fy4(3,j,i))/dys(j))
        e2(k) = e(k) - dt*((fx4(4,i+1,j)-fx4(4,i,j))/dxs(i) &
                         + (fy4(4,j+1,i)-fy4(4,j,i))/dys(j))
      end do
    end do

    ! ===================== species transport + chemistry =====================
    call transport_planar(nx, ny, rho2, dt, dxs, dys, fxm, fym, y2, dum)
    do k = 1, nx*ny
      sum_y = 0.0_DP
      do s = 1, NP5
        sum_y = sum_y + y2(s,k)
      end do
      if (sum_y > 1.0e-30_DP) then
        do s = 1, NP5
          y2(s,k) = y2(s,k) / sum_y
        end do
      end if
      rhe(k) = max(rho2(k), 1.0e-12_DP)
      dee(k) = e2(k)/rhe(k) - 0.5_DP*((mx2(k)/rhe(k))**2 + (my2(k)/rhe(k))**2)
      do s = 1, NP5
        yc(s,k) = y2(s,k)
      end do
    end do
    call chem_substep_grid(nx*ny, rhe, dee, dt, yc, dum, dum)
    do k = 1, nx*ny
      do s = 1, NP5
        y2(s,k) = max(yc(s,k), 0.0_DP)
      end do
    end do
    deallocate(thr, tu, tv, te, ty, sxr, sxu, sxv, sxe, sxy)
  end subroutine cart_park5_step

  ! ------------------------------------------------------------------
  ! c surface: the marshalling for the rust caller. the buffers cross
  ! the boundary in the flat layout the rust side lays out.
  subroutine nfor_cart_park5_step(nx, ny, rho, mx, my, e, y, dxs, dys, cfl, &
                                  bc_w, bc_e, bc_s, bc_n, bc_vals, y_fill, &
                                  rho2, mx2, my2, e2, y2, dt_out, err) &
      bind(c, name="nfor_cart_park5_step")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), value, intent(in) :: nx, ny
    real(c_double), intent(in) :: rho(*), mx(*), my(*), e(*), y(*)
    real(c_double), intent(in) :: dxs(*), dys(*)
    ! the cfl scalar crosses by value, so it must be declared that way.
    real(c_double), value, intent(in) :: cfl
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(*)
    real(c_double), intent(in) :: y_fill(*)
    real(c_double), intent(out) :: rho2(*), mx2(*), my2(*), e2(*), y2(*)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: err
    real(c_double) :: bcv(5,4), yf(5)
    ! the core wants explicit-shape arrays, so the crossing buffers are
    ! copied into locals with a definite shape first.
    real(c_double), allocatable :: rh(:), mh(:), nh(:), eh(:), ro(:), mo(:), no(:), qo(:)
    real(c_double), allocatable :: dxn(:), dyn(:), ya(:,:), y2a(:,:)
    integer :: s, v, sd, nn, ni, nj

    if (nx <= 0 .or. ny <= 0) then
      err = 1
      return
    end if
    ! the scratch pads are fixed size, so an oversized grid is refused
    ! rather than scribbling past them.
    if (int(nx) > np_slope - 4 .or. int(ny) > np_pad - 4) then
      err = 1
      return
    end if
    ni = int(nx); nj = int(ny); nn = ni*nj
    allocate(rh(nn), mh(nn), nh(nn), eh(nn), ro(nn), mo(nn), no(nn), qo(nn))
    allocate(dxn(ni), dyn(nj), ya(5,nn), y2a(5,nn))
    rh = rho(1:nn); mh = mx(1:nn); nh = my(1:nn); eh = e(1:nn)
    dxn = dxs(1:ni); dyn = dys(1:nj)
    do sd = 1, 4
      do v = 1, 5
        bcv(v,sd) = bc_vals((sd-1)*5 + v)
      end do
    end do
    do s = 1, 5
      yf(s) = y_fill(s)
    end do
    do s = 1, 5
      ya(s,1:nn) = y((s-1)*nn + 1:s*nn)
    end do
    call cart_park5_step(ni, nj, rh, mh, nh, eh, ya, &
                         dxn, dyn, cfl, int(bc_w), int(bc_e), int(bc_s), int(bc_n), &
                         bcv, yf, ro, mo, no, qo, y2a, dt_out)
    rho2(1:nn) = ro; mx2(1:nn) = mo; my2(1:nn) = no; e2(1:nn) = qo
    do s = 1, 5
      y2((s-1)*nn + 1:s*nn) = y2a(s,1:nn)
    end do
    deallocate(rh, mh, nh, eh, ro, mo, no, qo, dxn, dyn, ya, y2a)
    err = 0
  end subroutine nfor_cart_park5_step

end module nuforcartpark5

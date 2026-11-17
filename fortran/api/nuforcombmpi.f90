! nuforcombmpi.f90 -- the mpi-native li2004 combustion march.
!
! this module advances the nine-species reactive mixture on an x-slab
! decomposition, transcribing the serial axisymmetric recipe the CLI's
! combustion route calls (advance2d_combustion_axi_rk2): heun average of
! two single-pass advances on the conserved vectors, radial y divergence
! on the face radii, the p/rc radial source, donor-cell species transport,
! per-cell renorm, and the operator-split chemistry on the post-Heun
! state.
!
! the layout mapping follows nuforpark5mpi: local cell c stands for pad
! column c+2, two-deep ghosts carry the neighbour columns, and the halo
! exchange plus the domain-edge prescription give both ranks of a seam
! the same global columns, hence the same face states and fluxes.
module nuforcombmpi
  use iso_fortran_env, only: real64
  use mpi
  use nufor_mpi
  use nuforcomb, only: comb_triple, chem_substep_grid_li, internal_energy_at, sound_speed_at
  use nuforriemann, only: hllc_face_x_pub => hllc_face_x
  implicit none

  integer, parameter :: RP = real64
  integer, parameter :: NC9 = 9

  private
  public :: nfor_mpicomb_march

contains

  pure function vlim(dm, dp_) result(v)
    real(RP), intent(in) :: dm, dp_
    real(RP) :: v
    if (dm*dp_ > 0.0_RP) then
      v = 2.0_RP*dm*dp_/(dm+dp_)
    else
      v = 0.0_RP
    end if
  end function vlim

  ! ------------------------------------------------------------------
  ! c surface: advance this rank's slab one heun combustion step.
  ! all carriers are conserved vectors, flat k = (j-1)*nx_local + i;
  ! faces_y carries the ny+1 face radii; ys has one plane per species.
  subroutine nfor_mpicomb_march(nx, ny, rho, mx, my, e_tot, ys, &
                                 dxs, dys, faces_y, cfl, react, &
                                 bc_w, bc_e, bc_s, bc_n, bc_vals, y_fill, &
                                 rho_o, mx_o, my_o, e_o, ys_o, dt_out, err) &
      bind(c, name="nfor_mpicomb_march")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), value, intent(in) :: nx, ny
    real(c_double), intent(in) :: rho(*), mx(*), my(*), e_tot(*), ys(*)
    real(c_double), intent(in) :: dxs(*), dys(*), faces_y(*)
    real(c_double), value, intent(in) :: cfl
    integer(c_int), value, intent(in) :: react
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(*), y_fill(*)
    real(c_double), intent(out) :: rho_o(*), mx_o(*), my_o(*), e_o(*), ys_o(*)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: err

    real(RP) :: bcv(5,4), yf(NC9)
    real(RP), allocatable :: c0(:,:,:), ca(:,:,:)          ! padded carriers, 13 wide
    real(RP), allocatable :: rho1(:), mx1(:), my1(:), e1(:), y1(:,:)
    real(RP), allocatable :: efloor(:), rh(:), eh(:), ycw(:,:)
    real(RP), allocatable :: scr(:,:,:)
    real(RP), allocatable :: fx4(:,:,:), fy4(:,:,:), fxm(:,:), fym(:,:)
    real(RP), allocatable :: sxr(:), sxu(:), sxv(:), sxe(:), sxy(:,:)
    real(RP), allocatable :: thr(:), tu(:), tv(:), te(:), ty(:,:)
    real(RP), allocatable :: yx(:,:,:), yy(:,:,:)
    real(RP) :: pl, al, smax, smax_g, dt, dt2
    real(RP) :: dmin_x, dmin_x_g, dmin_y, ysum, e_int_k, worst_ov
    type(decomp_t) :: dec
    integer :: world, rank, eerr, i, j, s, k, nl, nj, nn, ng, sd, vv, step, nbmax, nf ,ixof, iyof
    ixof = 0; iyof = 1
    block
      character(len=8) :: ev
      call get_environment_variable('NF_OX', ev); if (len_trim(ev)>0) read(ev,*) ixof
      call get_environment_variable('NF_OY', ev); if (len_trim(ev)>0) read(ev,*) iyof
    end block

    err = 0
    call mpi_comm_rank(MPI_COMM_WORLD, rank, eerr)
    call mpi_comm_size(MPI_COMM_WORLD, world, eerr)
    ng = int(nx); nj = int(ny)
    if (ng <= 0 .or. nj <= 0) then
      err = 1
      return
    end if
    dec = make_decomp(ng, nj, world, rank)
    nl = dec%nx_local
    if (nl <= 0) then
      dt_out = 0.0_c_double
      return
    end if
    nn = nl*nj
    nbmax = max(nl, nj) + 4
    nf = nl + 1

    do sd = 1, 4
      do vv = 1, 5
        bcv(vv,sd) = bc_vals((sd-1)*5 + vv)
      end do
    end do
    do s = 1, NC9
      yf(s) = y_fill(s)
    end do

    allocate(c0(13,nl+4,nj), ca(13,nl+4,nj), scr(13,nl+4,nj))
    allocate(fx4(4,0:nf,nj), fy4(4,0:nj+1,nl), fxm(0:nf,nj), fym(0:nj+1,nl))
    allocate(sxr(nbmax), sxu(nbmax), sxv(nbmax), sxe(nbmax), sxy(NC9,nbmax))
    allocate(thr(nbmax), tu(nbmax), tv(nbmax), te(nbmax), ty(NC9,nbmax))
    allocate(yx(NC9,0:nf,nj), yy(NC9,0:nj+1,nl))
    allocate(rho1(nn), mx1(nn), my1(nn), e1(nn), y1(NC9,nn))
    allocate(efloor(nn), rh(nn), eh(nn), ycw(NC9,nn))

    ! pack entry state into a padded primitive carrier
    do j = 1, nj
      do i = 1, nl
        k = (j-1)*nl + i
        c0(1,i+2,j) = max(rho(k), 1.0e-12_RP)
        c0(2,i+2,j) = mx(k)/c0(1,i+2,j)
        c0(3,i+2,j) = my(k)/c0(1,i+2,j)
        c0(4,i+2,j) = e_tot(k)/c0(1,i+2,j) - 0.5_RP*(c0(2,i+2,j)**2 + c0(3,i+2,j)**2)
        do s = 1, NC9
          c0(4+s,i+2,j) = ys((s-1)*nn + k)
        end do
        efloor(k) = internal_energy_at(200.0_RP, c0(5:13,i+2,j))
      end do
    end do
    ca = c0

    ! ---- global timestep probe on the entry state ----
    smax = 0.0_RP
    dmin_x = huge(1.0_RP)
    do i = 1, nl
      dmin_x = min(dmin_x, dxs(i))
    end do
    dmin_y = huge(1.0_RP)
    do j = 1, nj
      dmin_y = min(dmin_y, dys(j))
    end do
    do j = 1, nj
      do i = 1, nl
        k = (j-1)*nl + i
        al = sound_speed_at(max(rho(k),1.0e-12_RP), c0(4,i+2,j), c0(5:13,i+2,j))
        smax = max(smax, max(abs(c0(2,i+2,j)), abs(c0(3,i+2,j))) + al)
      end do
    end do
    call global_max(smax, smax_g, dec)
    call global_min(dmin_x, dmin_x_g, dec)
    ! heun pass 1 dt probes the entry state
    dt = cfl*min(dmin_x_g, dmin_y)/max(smax_g, 1.0e-12_RP)

    ! ============ heun on conserved vectors ============
    ! the serial rk2 advances the pass-1 result a second time (s1 += dt*l(s1))
    ! and then mixes 0.5*u0 + 0.5*s1, so pass 2 must start from the pass-1
    ! state, not from the heun average.
    do step = 1, 2
      call single_pass(dec, nl, nj, ca, dxs, dys, faces_y, dt, bcv, yf, &
                       int(bc_w), int(bc_e), int(bc_s), int(bc_n), efloor, &
                       rho1, mx1, my1, e1, y1, scr, fx4, fy4, fxm, fym, &
                       sxr, sxu, sxv, sxe, sxy, thr, tu, tv, te, ty, yx, yy, &
                       nbmax, nf, nn, err)
      if (err /= 0) return
      if (step == 1) then
        ! repack the pass-1 state as pass 2's entry
        do j = 1, nj
          do i = 1, nl
            k = (j-1)*nl + i
            ca(1,i+2,j) = max(rho1(k), 1.0e-12_RP)
            ca(2,i+2,j) = mx1(k)/ca(1,i+2,j)
            ca(3,i+2,j) = my1(k)/ca(1,i+2,j)
            ca(4,i+2,j) = e1(k)/ca(1,i+2,j) - 0.5_RP*(ca(2,i+2,j)**2 + ca(3,i+2,j)**2)
            do s = 1, NC9
              ca(4+s,i+2,j) = y1(s,k)
            end do
          end do
        end do
        ! pass 2 dt probes the pass-1 state, like the serial lane
        smax = 0.0_RP
        do j = 1, nj
          do i = 1, nl
            k = (j-1)*nl + i
            al = sound_speed_at(max(ca(1,i+2,j),1.0e-12_RP), ca(4,i+2,j), ca(5:13,i+2,j))
            smax = max(smax, max(abs(ca(2,i+2,j)), abs(ca(3,i+2,j))) + al)
          end do
        end do
        call global_max(smax, smax_g, dec)
        dt = cfl*min(dmin_x_g, dmin_y)/max(smax_g, 1.0e-12_RP)
      else
        do k = 1, nn
          rho_o(k) = 0.5_RP*rho(k) + 0.5_RP*rho1(k)
          mx_o(k)  = 0.5_RP*mx(k)  + 0.5_RP*mx1(k)
          my_o(k)  = 0.5_RP*my(k)  + 0.5_RP*my1(k)
          e_o(k)   = 0.5_RP*e_tot(k) + 0.5_RP*e1(k)
          do s = 1, NC9
            ys_o((s-1)*nn + k) = 0.5_RP*ys((s-1)*nn + k) + 0.5_RP*y1(s,k)
          end do
        end do
      end if
    end do

    ! ---- renormalise compositions on the post-heun state ----
    do k = 1, nn
      ysum = 0.0_RP
      do s = 1, NC9
        ysum = ysum + ys_o((s-1)*nn + k)
      end do
      if (ysum > 1.0e-30_RP) then
        do s = 1, NC9
          ys_o((s-1)*nn + k) = ys_o((s-1)*nn + k)/ysum
        end do
      end if
    end do

    ! ---- chemistry on the post-heun state, dt from the same probe ----
    if (react /= 0) then
      smax = 0.0_RP
      do k = 1, nn
        rh(k) = max(rho_o(k), 1.0e-12_RP)
        eh(k) = e_o(k)/rh(k) - 0.5_RP*((mx_o(k)/rh(k))**2 + (my_o(k)/rh(k))**2)
        al = sound_speed_at(rh(k), eh(k), ys_o_col(nn, k, ys_o))
        smax = max(smax, max(abs(mx_o(k)/rh(k)), abs(my_o(k)/rh(k))) + al)
      end do
      call global_max(smax, smax_g, dec)
      dt2 = cfl*min(dmin_x_g, dmin_y)/max(smax_g, 1.0e-12_RP)
      dt_out = dt2
      do k = 1, nn
        do s = 1, NC9
          ycw(s,k) = ys_o((s-1)*nn + k)
        end do
      end do
      call chem_substep_grid_li(nn, rh, eh, dt2, ycw, worst_ov)
      do k = 1, nn
        do s = 1, NC9
          ys_o((s-1)*nn + k) = max(ycw(s,k), 0.0_RP)
        end do
      end do
    else
      dt_out = dt
    end if

    deallocate(c0, ca, scr, fx4, fy4, fxm, fym, sxr, sxu, sxv, sxe, sxy)
    deallocate(thr, tu, tv, te, ty, yx, yy, rho1, mx1, my1, e1, y1)
    deallocate(efloor, rh, eh, ycw)
  end subroutine nfor_mpicomb_march

  ! composition column k viewed out of the flat species planes
  function ys_o_col(nn, k, ys) result(col)
    integer, intent(in) :: nn, k
    real(RP), intent(in) :: ys(*)
    real(RP) :: col(NC9)
    integer :: s
    do s = 1, NC9
      col(s) = ys((s-1)*nn + k)
    end do
  end function ys_o_col

  ! ------------------------------------------------------------------
  ! one single-pass advance of the padded carrier, producing conserved
  ! outputs. mirrors advance2d_combustion_axi step by step.
  subroutine single_pass(dec, nl, nj, ca, dxs, dys, faces_y, dt, bcv, yf, &
                         bc_w, bc_e, bc_s, bc_n, efloor, &
                         rho1, mx1, my1, e1, y1, scr, fx4, fy4, fxm, fym, &
                         sxr, sxu, sxv, sxe, sxy, thr, tu, tv, te, ty, yx, yy, &
                         nbmax, nf, nn, err)
    type(decomp_t), intent(in) :: dec
    integer, intent(in) :: nl, nj, bc_w, bc_e, bc_s, bc_n, nbmax, nf, nn
    real(RP), intent(inout) :: ca(13,nl+4,nj)
    real(RP), intent(in) :: dxs(nl), dys(nj), faces_y(nj+1), dt, bcv(5,4), yf(NC9), efloor(nn)
    real(RP), intent(out) :: rho1(nn), mx1(nn), my1(nn), e1(nn), y1(NC9,nn)
    real(RP), intent(inout) :: scr(13,nl+4,nj)
    real(RP), intent(inout) :: fx4(4,0:nf,nj), fy4(4,0:nj+1,nl), fxm(0:nf,nj), fym(0:nj+1,nl)
    real(RP), intent(inout) :: sxr(nbmax), sxu(nbmax), sxv(nbmax), sxe(nbmax), sxy(NC9,nbmax)
    real(RP), intent(inout) :: thr(nbmax), tu(nbmax), tv(nbmax), te(nbmax), ty(NC9,nbmax)
    real(RP), intent(inout) :: yx(NC9,0:nf,nj), yy(NC9,0:nj+1,nl)
    integer, intent(out) :: err

    real(RP) :: rl, rul, vul, el, pl, al, rr, rur, vur, er_, pr_, ar_
        real(RP) :: fm, fqx, fqy, fe
    real(RP) :: yleft(NC9), yright(NC9), ysum, rho_y, flux_x, flux_y
    integer :: ng = 0
    real(RP) :: dxi, rad, dflux, new_rho_y, rm, rhs_rc, rc, dr
    real(RP) :: pk, eik, uk, vk, tk_dummy, rpv
    integer :: i, j, f, s, k, c, ixof, iyof

    err = 0
    scr = ca
    ixof = 0; iyof = 1
    block
      character(len=8) :: evx
      call get_environment_variable('NF_OX', evx); if (len_trim(evx)>0) read(evx,*) ixof
      call get_environment_variable('NF_OY', evx); if (len_trim(evx)>0) read(evx,*) iyof
    end block

    ! seed ghosts as edge copies, then halo + domain edges
    do j = 1, nj
      do i = 1, 2
        scr(:,i,j) = scr(:,3,j)
      end do
      do i = nl+3, nl+4
        scr(:,i,j) = scr(:,nl+2,j)
      end do
    end do
    call halo_exchange_x(scr, 13, dec)
    call comb_domain_edges(scr, nl, nj, dec, bc_w, bc_e, bcv, yf)

    fx4 = 0.0_RP; fy4 = 0.0_RP; fxm = 0.0_RP; fym = 0.0_RP; yx = 0.0_RP; yy = 0.0_RP
    ! ===================== x sweep =====================
    do j = 1, nj
      thr(1:nl+4) = scr(1,1:nl+4,j)
      tu(1:nl+4) = scr(2,1:nl+4,j)
      tv(1:nl+4) = scr(3,1:nl+4,j)
      te(1:nl+4) = scr(4,1:nl+4,j)
      ty(:,1:nl+4) = scr(5:13,1:nl+4,j)
      do c = 2, nl+3
        sxr(c-1) = vlim(thr(c)-thr(c-1), thr(c+1)-thr(c))
        sxu(c-1) = vlim(tu(c)-tu(c-1), tu(c+1)-tu(c))
        sxv(c-1) = vlim(tv(c)-tv(c-1), tv(c+1)-tv(c))
        sxe(c-1) = vlim(te(c)-te(c-1), te(c+1)-te(c))
        do s = 1, NC9
          sxy(s,c-1) = vlim(ty(s,c)-ty(s,c-1), ty(s,c+1)-ty(s,c))
        end do
      end do
      do f = 1, nf
        ! boundary faces match the serial reconstruction exactly: the west
        ! face keeps its left state unsloped, the east face its right one
        if (f == 1) then
          rl = thr(2);             rr = thr(3) - 0.5_RP*sxr(2)
          rul = tu(2);             rur = tu(3) - 0.5_RP*sxu(2)
          vul = tv(2);             vur = tv(3) - 0.5_RP*sxv(2)
          el = te(2);              er_ = te(3) - 0.5_RP*sxe(2)
          do s = 1, NC9
            yleft(s) = ty(s,2);    yright(s) = ty(s,3) - 0.5_RP*sxy(s,2)
          end do
        else if (f == nf) then
          rl = thr(f+1) + 0.5_RP*sxr(f); rr = thr(f+2)
          rul = tu(f+1) + 0.5_RP*sxu(f); rur = tu(f+2)
          vul = tv(f+1) + 0.5_RP*sxv(f); vur = tv(f+2)
          el = te(f+1) + 0.5_RP*sxe(f);  er_ = te(f+2)
          do s = 1, NC9
            yleft(s) = ty(s,f+1) + 0.5_RP*sxy(s,f); yright(s) = ty(s,f+2)
          end do
        else
          rl = thr(f+1) + 0.5_RP*sxr(f);  rr = thr(f+2) - 0.5_RP*sxr(f+1)
          rul = tu(f+1) + 0.5_RP*sxu(f);  rur = tu(f+2) - 0.5_RP*sxu(f+1)
          vul = tv(f+1) + 0.5_RP*sxv(f);  vur = tv(f+2) - 0.5_RP*sxv(f+1)
          el = te(f+1) + 0.5_RP*sxe(f);   er_ = te(f+2) - 0.5_RP*sxe(f+1)
          do s = 1, NC9
            yleft(s) = ty(s,f+1) + 0.5_RP*sxy(s,f)
            yright(s) = ty(s,f+2) - 0.5_RP*sxy(s,f+1)
          end do
        end if
        call comb_triple(max(rl,1.0e-12_RP), el, yleft, pl, tk_dummy, al)
        call comb_triple(max(rr,1.0e-12_RP), er_, yright, pr_, tk_dummy, ar_)
        pl = max(pl, 1.0e-12_RP); pr_ = max(pr_, 1.0e-12_RP)
        al = max(al, 1.0e-12_RP); ar_ = max(ar_, 1.0e-12_RP)
        block
          character(len=8) :: dfxv
          call get_environment_variable('NF_DEBUG_FX', dfxv)
          if (trim(dfxv) == '1' .and. j == 8 .and. (f == 11 .or. f == 12 .or. f == 13)) then
            write(*,'(A,I3,8ES16.8)') 'FP ', f, rl, rul, el, pl, rr, rur, er_, pr_
          end if
        end block
        call hllc_face_x_pub(rl, rul, vul, pl, el, al, rr, rur, vur, pr_, er_, ar_, &
                             1.4_RP, fm, fqx, fqy, fe)
        fx4(1,f,j) = fm; fx4(2,f,j) = fqx; fx4(3,f,j) = fqy; fx4(4,f,j) = fe
        fxm(f,j) = fm
      end do
    end do

    ! ===================== y sweep =====================
    do i = 1, nl
      do j = 1, nj
        c = j + 2
        thr(c) = scr(1,i+2,j); tu(c) = scr(2,i+2,j)
        tv(c) = scr(3,i+2,j); te(c) = scr(4,i+2,j)
        ty(:,c) = scr(5:13,i+2,j)
      end do
      select case (bc_s)
      case (1)
        thr(2) = bcv(1,3); tu(2) = bcv(2,3); tv(2) = bcv(3,3)
        te(2) = bcv(5,3); ty(:,2) = yf
      case (3)
        thr(2) = thr(3); tu(2) = tu(3); tv(2) = -tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case (4)
        thr(2) = thr(3); tu(2) = -tu(3); tv(2) = -tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case default
        thr(2) = thr(3); tu(2) = tu(3); tv(2) = tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      end select
      thr(1) = thr(2); tu(1) = tu(2); tv(1) = tv(2); te(1) = te(2); ty(:,1) = ty(:,2)
      select case (bc_n)
      case (1)
        thr(nj+3) = bcv(1,4); tu(nj+3) = bcv(2,4); tv(nj+3) = bcv(3,4)
        te(nj+3) = bcv(5,4); ty(:,nj+3) = yf
      case (3)
        thr(nj+3) = thr(nj+2); tu(nj+3) = tu(nj+2); tv(nj+3) = -tv(nj+2)
        te(nj+3) = te(nj+2); ty(:,nj+3) = ty(:,nj+2)
      case (4)
        thr(nj+3) = thr(nj+2); tu(nj+3) = -tu(nj+2); tv(nj+3) = -tv(nj+2)
        te(nj+3) = te(nj+2); ty(:,nj+3) = ty(:,nj+2)
      case default
        thr(nj+3) = thr(nj+2); tu(nj+3) = tu(nj+2); tv(nj+3) = tv(nj+2)
        te(nj+3) = te(nj+2); ty(:,nj+3) = ty(:,nj+2)
      end select
      thr(nj+4) = thr(nj+3); tu(nj+4) = tu(nj+3); tv(nj+4) = tv(nj+3)
      te(nj+4) = te(nj+3); ty(:,nj+4) = ty(:,nj+3)
      do c = 2, nj+3
        sxr(c-1) = vlim(thr(c)-thr(c-1), thr(c+1)-thr(c))
        sxu(c-1) = vlim(tu(c)-tu(c-1), tu(c+1)-tu(c))
        sxv(c-1) = vlim(tv(c)-tv(c-1), tv(c+1)-tv(c))
        sxe(c-1) = vlim(te(c)-te(c-1), te(c+1)-te(c))
        do s = 1, NC9
          sxy(s,c-1) = vlim(ty(s,c)-ty(s,c-1), ty(s,c+1)-ty(s,c))
        end do
      end do
      do f = 1, nj+1
        if (f == 1) then
          rl = thr(2);                   rr = thr(3) - 0.5_RP*sxr(1)
          rul = tv(2);                   rur = tv(3) - 0.5_RP*sxv(1)
          vul = tu(2);                   vur = tu(3) - 0.5_RP*sxu(1)
          el = te(2);                    er_ = te(3) - 0.5_RP*sxe(1)
          do s = 1, NC9
            yleft(s) = ty(s,2);          yright(s) = ty(s,3) - 0.5_RP*sxy(s,1)
          end do
        else if (f == nj+1) then
          rl = thr(f+1) + 0.5_RP*sxr(f); rr = thr(f+2)
          rul = tv(f+1) + 0.5_RP*sxv(f); rur = tv(f+2)
          vul = tu(f+1) + 0.5_RP*sxu(f); vur = tu(f+2)
          el = te(f+1) + 0.5_RP*sxe(f);  er_ = te(f+2)
          do s = 1, NC9
            yleft(s) = ty(s,f+1) + 0.5_RP*sxy(s,f)
            yright(s) = ty(s,f+2)
          end do
        else
          rl = thr(f+1) + 0.5_RP*sxr(f);  rr = thr(f+2) - 0.5_RP*sxr(f+1)
          rul = tv(f+1) + 0.5_RP*sxv(f);  rur = tv(f+2) - 0.5_RP*sxv(f+1)
          vul = tu(f+1) + 0.5_RP*sxu(f);  vur = tu(f+2) - 0.5_RP*sxu(f+1)
          el = te(f+1) + 0.5_RP*sxe(f);   er_ = te(f+2) - 0.5_RP*sxe(f+1)
          do s = 1, NC9
            yleft(s) = ty(s,f+1) + 0.5_RP*sxy(s,f)
            yright(s) = ty(s,f+2) - 0.5_RP*sxy(s,f+1)
          end do
        end if
        call comb_triple(max(rl,1.0e-12_RP), el, yleft, pl, tk_dummy, al)
        call comb_triple(max(rr,1.0e-12_RP), er_, yright, pr_, tk_dummy, ar_)
        pl = max(pl, 1.0e-12_RP); pr_ = max(pr_, 1.0e-12_RP)
        al = max(al, 1.0e-12_RP); ar_ = max(ar_, 1.0e-12_RP)
        ! y-normal faces ride the x solver with velocities permuted
        call hllc_face_x_pub(rl, rul, vul, pl, el, al, rr, rur, vur, pr_, er_, ar_, &
                             1.4_RP, fm, fqy, fqx, fe)
        fy4(1,f,i) = fm; fy4(2,f,i) = fqx; fy4(3,f,i) = fqy; fy4(4,f,i) = fe
        fym(f,i) = fm
      end do
    end do

    block
      character(len=8) :: dfxv
      call get_environment_variable('NF_DEBUG_FX', dfxv)
      if (trim(dfxv) == '1') then
      do f = 1, nf
        write(*,'(A,I3,4ES24.15)') 'FXD ', f, fx4(1,f,8), fx4(2,f,8), fx4(3,f,8), fx4(4,f,8)
      end do
      end if
    end block
    ! ===================== conservative update (annular) =====================
    do j = 1, nj
      rm = faces_y(j); rpv = faces_y(j+1)
      rc = 0.5_RP*(rm + rpv); dr = dys(j)
      do i = 1, nl
        k = (j-1)*nl + i
        dxi = dt/dxs(i)
        rad = (rpv*fy4(1,j+1,i) - rm*fy4(1,j,i))/(rc*dr)
        rho1(k) = ca(1,i+2,j) - (dxi*(fx4(1,i,j)-fx4(1,i-1,j)) + dt*rad)
        rad = (rpv*fy4(2,j+1,i) - rm*fy4(2,j,i))/(rc*dr)
        mx1(k) = ca(1,i+2,j)*ca(2,i+2,j) - (dxi*(fx4(2,i,j)-fx4(2,i-1,j)) + dt*rad)
        ! radial pressure source on the y momentum
        call comb_triple(ca(1,i+2,j), ca(4,i+2,j), ca(5:13,i+2,j), pk, tk_dummy, al)
        rad = (rpv*fy4(3,j+1,i) - rm*fy4(3,j,i))/(rc*dr)
        my1(k) = ca(1,i+2,j)*ca(3,i+2,j) - (dxi*(fx4(3,i,j)-fx4(3,i-1,j)) + dt*rad) + dt*pk/rc
        rad = (rpv*fy4(4,j+1,i) - rm*fy4(4,j,i))/(rc*dr)
        e1(k) = ca(1,i+2,j)*(ca(4,i+2,j) + 0.5_RP*(ca(2,i+2,j)**2 + ca(3,i+2,j)**2)) &
                - (dxi*(fx4(4,i,j)-fx4(4,i-1,j)) + dt*rad)
        block
          character(len=8) :: duv3
          call get_environment_variable('NF_DEBUG_E', duv3)
          if (trim(duv3) == '1' .and. j == 1 .and. i <= 3) then
            write(*,'(A,I3,5ES16.8)') 'EE ', i, ca(4,i+2,j), dxi*(fx4(4,i,j)-fx4(4,i-1,j)), dt*rad, fx4(4,i,j)
          end if
        end block
        ! positivity guard mirroring the serial lane
        eik = e1(k)/max(rho1(k),1.0e-12_RP) - 0.5_RP*((mx1(k)/max(rho1(k),1.0e-12_RP))**2 &
              + (my1(k)/max(rho1(k),1.0e-12_RP))**2)
        if (.not. (e1(k) > 0.0_RP .and. rho1(k) > 0.0_RP) .or. eik < efloor(k)) then
          ng = ng + 1
          if (ng <= 3) write(*,'(A,I4,4ES16.8)') 'GRIP ', k, e1(k), rho1(k), eik, efloor(k)
          rho1(k) = max(rho1(k), 1.0e-10_RP)
          mx1(k) = 0.0_RP
          my1(k) = 0.0_RP
          e1(k) = efloor(k)*rho1(k) + 1.0e-6_RP
        end if
      end do
    end do

    block
      character(len=8) :: duv2
      integer :: jj, ff
      call get_environment_variable('NF_DEBUG_U', duv2)
      if (trim(duv2) == '1') then
        do i = 1, nl
          k = i
          write(*,'(A,I3,5ES16.8)') 'UU ', i, rho1(k), mx1(k), my1(k), e1(k), dt
        if (i == 13) then
          do jj = 1, nj+1
            write(*,'(A,I3,4ES16.8)') 'FY ', jj, fy4(1,jj,13), fy4(2,jj,13), fy4(3,jj,13), fy4(4,jj,13)
          end do
          do ff = 12, 15
            write(*,'(A,I3,4ES16.8)') 'FX ', ff, fx4(1,ff,1), fx4(2,ff,1), fx4(3,ff,1), fx4(4,ff,1)
          end do
          write(*,'(A,4ES16.8)') 'GEO ', faces_y(1), faces_y(2), faces_y(nj+1), dys(1)
          write(*,'(A,4ES16.8)') 'FYe ', fy4(4,1,1), fy4(4,2,1), fy4(4,2,2), fy4(4,1,2)
          write(*,'(A,4ES16.8)') 'SCe ', scr(4,3,1), scr(4,3,2), scr(4,4,1), scr(4,4,2)
        end if
        end do
      end if
    end block
    ! ===================== species transport (donor-cell) =====================
    do j = 1, nj
      do f = 1, nf
        if (fxm(f,j) >= 0.0_RP) then
          do s = 1, NC9
            if (f > 1) then
              yx(s, f,j) = scr(4+s,f+1,j)
            else if (dec%west_domain_edge) then
              yx(s, f,j) = scr(4+s,3,j)
            else
              yx(s, f,j) = scr(4+s,3,j)
            end if
          end do
        else
          do s = 1, NC9
            if (f < nf) then
              yx(s, f,j) = scr(4+s,f+1,j)
            else if (dec%east_domain_edge) then
              yx(s, f,j) = scr(4+s,nl+2,j)
            else
              yx(s, f,j) = scr(4+s,nl+2,j)
            end if
          end do
        end if
      end do
    end do
    do i = 1, nl
      do f = 1, nj+1
        if (fym(f,i) >= 0.0_RP) then
          do s = 1, NC9
            if (f > 1) then
              yy(s, f,i) = scr(4+s,i+2,f-1)
            else
              yy(s, f,i) = scr(4+s,i+2,3)
            end if
          end do
        else
          do s = 1, NC9
            if (f < nj+1) then
              yy(s, f,i) = scr(4+s,i+2,f)
            else
              yy(s, f,i) = scr(4+s,i+2,nj+2)
            end if
          end do
        end if
      end do
    end do

    do s = 1, NC9
      do j = 1, nj
        rm = faces_y(j); rpv = faces_y(j+1)
        rc = 0.5_RP*(rm + rpv); dr = dys(j)
        do i = 1, nl
          k = (j-1)*nl + i
          flux_x = yx(s,i,j)*fxm(i,j) - yx(s,i-1,j)*fxm(i-1,j)
          rad = (rpv*yy(s,j+iyof,i)*fym(j+iyof,i) - rm*yy(s,j+iyof-1,i)*fym(j+iyof-1,i))/(rc*dr)
          dxi = dt/dxs(i)
          rho_y = rho1(k)*scr(4+s,i+2,j)
          dflux = -(dxi*flux_x + dt*rad)
          new_rho_y = rho_y + dflux
          if (new_rho_y > 0.0_RP) then
            y1(s,k) = new_rho_y/rho1(k)
          else
            y1(s,k) = 0.0_RP
          end if
        end do
      end do
    end do
  end subroutine single_pass

  ! ------------------------------------------------------------------
  ! domain-edge ghost columns for the x faces
  subroutine comb_domain_edges(car, nl, nj, dec, bc_w, bc_e, bcv, yf)
    integer, intent(in) :: nl, nj, bc_w, bc_e
    type(decomp_t), intent(in) :: dec
    real(RP), intent(inout) :: car(13,nl+4,nj)
    real(RP), intent(in) :: bcv(5,4), yf(NC9)
    integer :: j, s
    if (dec%west_domain_edge) then
      do j = 1, nj
        select case (bc_w)
        case (1)
          car(1,2,j) = bcv(1,1); car(2,2,j) = bcv(2,1); car(3,2,j) = bcv(3,1)
          car(4,2,j) = bcv(5,1)
          do s = 1, NC9
            car(4+s,2,j) = yf(s)
          end do
        case (3)
          car(1,2,j) = car(1,3,j); car(2,2,j) = -car(2,3,j)
          car(3,2,j) = car(3,3,j); car(4,2,j) = car(4,3,j)
          do s = 1, NC9
            car(4+s,2,j) = car(4+s,3,j)
          end do
        case (4)
          car(1,2,j) = car(1,3,j); car(2,2,j) = -car(2,3,j)
          car(3,2,j) = -car(3,3,j); car(4,2,j) = car(4,3,j)
          do s = 1, NC9
            car(4+s,2,j) = car(4+s,3,j)
          end do
        case default
          car(1,2,j) = car(1,3,j); car(2,2,j) = car(2,3,j)
          car(3,2,j) = car(3,3,j); car(4,2,j) = car(4,3,j)
          do s = 1, NC9
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
          do s = 1, NC9
            car(4+s,nl+3,j) = yf(s)
          end do
        case (3)
          car(1,nl+3,j) = car(1,nl+2,j); car(2,nl+3,j) = -car(2,nl+2,j)
          car(3,nl+3,j) = car(3,nl+2,j); car(4,nl+3,j) = car(4,nl+2,j)
          do s = 1, NC9
            car(4+s,nl+3,j) = car(4+s,nl+2,j)
          end do
        case (4)
          car(1,nl+3,j) = car(1,nl+2,j); car(2,nl+3,j) = -car(2,nl+2,j)
          car(3,nl+3,j) = -car(3,nl+2,j); car(4,nl+3,j) = car(4,nl+2,j)
          do s = 1, NC9
            car(4+s,nl+3,j) = car(4+s,nl+2,j)
          end do
        case default
          car(1,nl+3,j) = car(1,nl+2,j); car(2,nl+3,j) = car(2,nl+2,j)
          car(3,nl+3,j) = car(3,nl+2,j); car(4,nl+3,j) = car(4,nl+2,j)
          do s = 1, NC9
            car(4+s,nl+3,j) = car(4+s,nl+2,j)
          end do
        end select
        car(:,nl+4,j) = car(:,nl+3,j)
      end do
    end if
  end subroutine comb_domain_edges

  ! gather the padded local strip to rank 0 with nf arbitrary plane blocks;
  ! same packing/displacement scheme as nfor_mpipark5_gather but field-count
  ! generic, so the combustion lane's 13 planes ride the same decomposition.
  subroutine nfor_mpicomb_gather(nf, nx, ny, y, out_y, ierr) &
      bind(c, name="nfor_mpicomb_gather")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), value, intent(in) :: nf, nx, ny
    real(c_double), intent(in) :: y(*)
    real(c_double), intent(out) :: out_y(*)
    integer(c_int), intent(out) :: ierr
    integer, allocatable :: counts(:), disps(:)
    real(RP), allocatable :: send(:), recv(:)
    integer :: r, f, base, rem, e2, nrow, nl_me, nw, wr, ws2
    integer :: i_l, j, nxl_r, i0_r, l_src, l_recv, k_dst
    type(decomp_t) :: gd
    ierr = 0
    call mpi_comm_rank(MPI_COMM_WORLD, wr, e2)
    call mpi_comm_size(MPI_COMM_WORLD, ws2, e2)
    gd = make_decomp(int(nx), int(ny), ws2, wr)
    nrow = int(ny)
    nl_me = gd%nx_local
    nw = nl_me + 4
    base = int(nx) / gd%nranks
    rem = mod(int(nx), gd%nranks)
    allocate(counts(gd%nranks), disps(gd%nranks))
    do r = 1, gd%nranks
      if (r-1 < rem) then
        counts(r) = (base+1)*nrow
      else
        counts(r) = base*nrow
      end if
    end do
    disps(1) = 0
    do r = 2, gd%nranks
      disps(r) = disps(r-1) + int(nf)*counts(r-1)
    end do
    allocate(send(int(nf)*nl_me*nrow), recv(int(nf)*int(nx)*nrow))
    do f = 1, int(nf)
      do j = 1, nrow
        do i_l = 1, nl_me
          send((f-1)*nl_me*nrow + (j-1)*nl_me + i_l) = &
            y((f-1)*nw*nrow + (j-1)*nw + i_l+2)
        end do
      end do
    end do
    call mpi_gatherv(send, int(nf)*nl_me*nrow, MPI_DOUBLE_PRECISION, recv, &
                     int(nf)*counts, disps, MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, e2)
    if (wr == 0) then
      do r = 1, gd%nranks
        nxl_r = counts(r)/nrow
        i0_r = (r-1)*base + min(r-1, rem) + 1
        do j = 1, nrow
          do i_l = 1, nxl_r
            l_src = (j-1)*nxl_r + i_l
            k_dst = (j-1)*int(nx) + (i0_r - 1) + i_l
            do f = 1, int(nf)
              l_recv = disps(r) + (f-1)*counts(r) + l_src
              out_y((f-1)*int(nx)*nrow + k_dst) = recv(l_recv)
            end do
          end do
        end do
      end do
    end if
    deallocate(counts, disps, send, recv)
  end subroutine nfor_mpicomb_gather
end module
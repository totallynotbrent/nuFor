! nuforpark5mpi.f90 — the mpi-native five-species park5 march.
!
! this module advances the five-species mixture on an x-slab decomposition.
! it transcribes the serial cartesian recipe in nuforcartpark5.f90 onto the
! slab layout, so the production mpi lane owns the same physics the serial
! lane owns instead of falling back to a rank-zero serial march.
!
! the layout mapping is the load-bearing part. the serial lane puts global
! cell g at pad column g+2, with one ghost column on each side. the slab
! puts its local cell c at pad column c+2 and holds global cells i0-2 and
! i0-1 in the west ghosts and i1+1 and i1+2 in the east ghosts. so every
! slab pad column stands for one global pad column, the shift being the
! rank's i0, and the two-depth halo supplies exactly the columns the
! serial lane reads at a seam.
!
! seam consistency: the face states at a seam are reconstructed from the
! same global centres on both sides of it, because both ranks hold the same
! global columns in their ghost cells and evaluate the same expressions in
! the same order on the same doubles. the shared face therefore carries the
! same mass flux and the same upwind donor on both ranks, so the species
! divergence cancels across the interface.
!
! at a global domain edge the rank applies the boundary prescription to its
! ghost column and copies that column into the sacrificial one, exactly as
! the serial lane does, which keeps the np=1 slab identical to the serial
! lane.
!
! the flow faces come from the fortran hllc kernel and the mixture closure
! from the fortran park5 closure, so this march shares its numerics with
! the serial park5 lane rather than approximating it.
module nuforpark5mpi
  use iso_fortran_env, only: real64
  use mpi
  use nufor_mpi
  use nuforpark5, only: closure_triple, chem_substep_grid
  use nuforprops, only: park5_viscous_update
  use nuforriemann, only: hllc_face_x_pub => hllc_face_x
  use nuforpark5common, only: park5_vlim, park5_domain_edges
  implicit none

  integer, parameter :: RP = real64
  ! the primitive four travel with the five composition planes.
  integer, parameter :: NCAR = 9
  integer, parameter :: NP5 = 5

  private
  public :: nfor_mpipark5_march

contains

  ! ------------------------------------------------------------------
  ! c surface: advance this rank's slab one park5 step.
  !
  ! the inputs are the conserved four (density, x momentum, y momentum,
  ! total energy) and the five composition fractions, flat with
  ! k = (j-1)*nx_local + i. the outputs are the same state after the step.
  subroutine nfor_mpipark5_march(nx, ny, rho, mx, my, e_tot, ys, &
                                 dxs, dys, faces_y, geom, visc, cfl, &
                                 bc_w, bc_e, bc_s, bc_n, bc_vals, y_fill, &
                                 rho_o, mx_o, my_o, e_o, ys_o, dt_out, err) &
      bind(c, name="nfor_mpipark5_march")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), value, intent(in) :: nx, ny
    real(c_double), intent(in) :: rho(*), mx(*), my(*), e_tot(*), ys(*)
    real(c_double), intent(in) :: dxs(*), dys(*), faces_y(*)
    integer(c_int), value, intent(in) :: geom, visc
    ! the cfl scalar crosses by value, so it must be declared that way.
    real(c_double), value, intent(in) :: cfl
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(*), y_fill(*)
    real(c_double), intent(out) :: rho_o(*), mx_o(*), my_o(*), e_o(*), ys_o(*)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: err

    real(RP) :: bcv(5,4), yf(NP5)
    real(RP), allocatable :: car(:,:,:)
    real(RP), allocatable :: fx4(:,:,:), fy4(:,:,:), fxm(:,:), fym(:,:)
    real(RP), allocatable :: sxr(:), sxu(:), sxv(:), sxe(:), sxy(:,:)
    real(RP), allocatable :: thr(:), tu(:), tv(:), te(:), ty(:,:)
    real(RP), allocatable :: yx(:,:,:), yy(:,:,:)
    real(RP), allocatable :: rhe(:), dee(:), ycw(:,:), pf(:)
    real(RP) :: rl, rul, vul, el, pl, al
    real(RP) :: rr, rur, vur, er_, pr_, ar_
    real(RP) :: fm, fqx, fqy, fe, dummy, smax, smax_g, dt
    real(RP) :: dmin_x, dmin_x_g, dmin_y, dt_den
    real(RP) :: rm, rpv, rc, pk
    real(RP) :: yleft(NP5), yright(NP5), ysum, rho_y, worst_dy, worst_ov
    real(RP) :: flux_x, flux_y, dxi, dyi, dflux, new_rho_y
    real(RP) :: dwx, dex, pri
    real(RP), allocatable :: cwl(:,:), cel(:,:)
    logical :: annular
    type(decomp_t) :: dec
    integer :: world, rank, eerr, ierr, i, j, f, s, k, c, nl, nj, nn, ng, nb, nf, nbmax
    integer :: gi, sd, vv

    err = 0
    call mpi_comm_rank(MPI_COMM_WORLD, rank, eerr)
    call mpi_comm_size(MPI_COMM_WORLD, world, eerr)
    ng = int(nx); nj = int(ny)
    if (ng <= 0 .or. nj <= 0) then
      err = 1
      return
    end if
    ! the shared viscous kernel is planar-metric only; refusing annular +
    ! viscous here keeps the r-weighted lane from taking a wrong-metric pass.
    if (geom /= 0 .and. visc /= 0) then
      err = 3
      dt_out = 0.0_c_double
      return
    end if
    dec = make_decomp(ng, nj, world, rank)
    nl = dec%nx_local
    if (nl <= 0) then
      ! a rank with no columns has nothing to advance, and it must not
      ! desynchronise the collective calls the other ranks are making.
      dt_out = 0.0_c_double
      return
    end if
    nn = nl*nj
    nb = nl + 4
    nf = nl + 1
    ! the slope scratch serves both sweeps, so it takes the larger extent.
    nbmax = max(nl, nj) + 4

    do sd = 1, 4
      do vv = 1, 5
        bcv(vv,sd) = bc_vals((sd-1)*5 + vv)
      end do
    end do
    do s = 1, NP5
      yf(s) = y_fill(s)
    end do

    allocate(car(NCAR,nb,nj))
    allocate(fx4(4,nf+1,nj), fy4(4,nj+1,nl))
    allocate(fxm(nf+1,nj), fym(nj+1,nl))
    allocate(sxr(nbmax), sxu(nbmax), sxv(nbmax), sxe(nbmax), sxy(NP5,nbmax))
    allocate(thr(nbmax), tu(nbmax), tv(nbmax), te(nbmax), ty(NP5,nbmax))
    allocate(yx(NP5,nf+1,nj), yy(NP5,nj+1,nl))
    allocate(rhe(nn), dee(nn), ycw(NP5,nn), pf(nn))

    annular = (geom /= 0)
    if (annular) then
      ! the southernmost face of an annulus must sit off the axis for
      ! rc to be invertible.
      if (faces_y(1) <= 0.0_RP) then
        err = 2
        return
      end if
    end if

    ! ---------- primitives and the per-cell closure ----------
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
        car(1,i+2,j) = max(rho(k), 1.0e-12_RP)
        car(2,i+2,j) = mx(k) / car(1,i+2,j)
        car(3,i+2,j) = my(k) / car(1,i+2,j)
        car(4,i+2,j) = e_tot(k) / car(1,i+2,j) &
                       - 0.5_RP*(car(2,i+2,j)**2 + car(3,i+2,j)**2)
        do s = 1, NP5
          car(4+s,i+2,j) = ys((s-1)*nn + k)
        end do
        call closure_triple(rho(k), car(4,i+2,j), car(5:NCAR,i+2,j), &
                            pl, dummy, al)
        pf(k) = max(pl, 1.0e-12_RP)
        ! the serial lane floors the mixture sound speed for the cfl
        ! bound, which is what holds the timestep down at low speed.
        al = max(al, 50.0_RP)
        smax = max(smax, max(abs(car(2,i+2,j)), abs(car(3,i+2,j))) + al)
      end do
    end do
    call global_max(smax, smax_g, dec)
    call global_min(dmin_x, dmin_x_g, dec)
    dt_den = max(smax_g, 1.0e-12_RP)
    dt = cfl * min(dmin_x_g, dmin_y) / dt_den
    dt_out = dt

    ! ---------- the x halo carries the primitives and the species -----
    ! the ghost columns start as copies of the slab edge; the exchange
    ! overwrites them with the real neighbours, and the boundary
    ! prescription replaces them at a global domain edge.
    do j = 1, nj
      do i = 1, 2
        car(:,i,j) = car(:,3,j)
      end do
      do i = nl+3, nl+4
        car(:,i,j) = car(:,nl+2,j)
      end do
    end do
    call halo_exchange_x(car, NCAR, dec)
    call park5_domain_edges(car, nl, nj, dec, int(bc_w), int(bc_e), bcv, yf)

    ! ===================== x sweep (rows) =====================
    ! face f (1..nl+1) sits between pad columns f+1 and f+2, so local
    ! cell c spans the faces c and c+1, as in the serial lane.
    do j = 1, nj
      ! the strip carries the padded row, ghost columns included, so the
      ! slope stencil reads the same columns the serial lane reads.
      thr(1:nl+4) = car(1,1:nl+4,j)
      tu(1:nl+4) = car(2,1:nl+4,j)
      tv(1:nl+4) = car(3,1:nl+4,j)
      te(1:nl+4) = car(4,1:nl+4,j)
      ty(:,1:nl+4) = car(5:NCAR,1:nl+4,j)
      do c = 2, nl+3
        sxr(c-1) = park5_vlim(thr(c)-thr(c-1), thr(c+1)-thr(c))
        sxu(c-1) = park5_vlim(tu(c)-tu(c-1), tu(c+1)-tu(c))
        sxv(c-1) = park5_vlim(tv(c)-tv(c-1), tv(c+1)-tv(c))
        sxe(c-1) = park5_vlim(te(c)-te(c-1), te(c+1)-te(c))
        do s = 1, NP5
          sxy(s,c-1) = park5_vlim(ty(s,c)-ty(s,c-1), ty(s,c+1)-ty(s,c))
        end do
      end do
      do f = 1, nf
        rl = thr(f+1) + 0.5_RP*sxr(f);  rr = thr(f+2) - 0.5_RP*sxr(f+1)
        rul = tu(f+1) + 0.5_RP*sxu(f);  rur = tu(f+2) - 0.5_RP*sxu(f+1)
        vul = tv(f+1) + 0.5_RP*sxv(f);  vur = tv(f+2) - 0.5_RP*sxv(f+1)
        el = te(f+1) + 0.5_RP*sxe(f);   er_ = te(f+2) - 0.5_RP*sxe(f+1)
        do s = 1, NP5
          yleft(s) = ty(s,f+1) + 0.5_RP*sxy(s,f)
          yright(s) = ty(s,f+2) - 0.5_RP*sxy(s,f+1)
        end do
        call closure_triple(max(rl,1.0e-12_RP), el, yleft, pl, dummy, al)
        call closure_triple(max(rr,1.0e-12_RP), er_, yright, pr_, dummy, ar_)
        pl = max(pl, 1.0e-12_RP); pr_ = max(pr_, 1.0e-12_RP)
        al = max(al, 1.0e-12_RP); ar_ = max(ar_, 1.0e-12_RP)
        call hllc_face_x_pub(rl, rul, vul, pl, el, al, rr, rur, vur, pr_, er_, ar_, &
                             1.4_RP, fm, fqx, fqy, fe)
        fx4(1,f,j) = fm; fx4(2,f,j) = fqx; fx4(3,f,j) = fqy; fx4(4,f,j) = fe
        fxm(f,j) = fm
      end do
    end do

    ! ===================== y sweep (columns) =====================
    ! the rows are not decomposed, so every column owns the full vertical
    ! extent and the boundary prescription applies at the global edges.
    do i = 1, nl
      do j = 1, nj
        k = (j-1)*nl + i
        c = j + 2
        thr(c) = car(1,i+2,j); tu(c) = car(2,i+2,j)
        tv(c) = car(3,i+2,j); te(c) = car(4,i+2,j)
        ty(:,c) = car(5:NCAR,i+2,j)
      end do
      select case (int(bc_s))
      case (1)
        ! the strip holds u in tu and v in tv and the face normal is the
        ! y direction, so the normal ghost goes into tv.
        thr(2) = bcv(1,3); tu(2) = bcv(2,3); tv(2) = bcv(3,3)
        te(2) = bcv(5,3)
        ty(:,2) = yf
      case (3)
        thr(2) = thr(3); tu(2) = tu(3); tv(2) = -tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case (4)
        thr(2) = thr(3); tu(2) = -tu(3); tv(2) = -tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      case default
        thr(2) = thr(3); tu(2) = tu(3); tv(2) = tv(3); te(2) = te(3); ty(:,2) = ty(:,3)
      end select
      thr(1) = thr(2); tu(1) = tu(2); tv(1) = tv(2); te(1) = te(2); ty(:,1) = ty(:,2)
      select case (int(bc_n))
      case (1)
        thr(nj+3) = bcv(1,4); tu(nj+3) = bcv(2,4); tv(nj+3) = bcv(3,4)
        te(nj+3) = bcv(5,4)
        ty(:,nj+3) = yf
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
        sxr(c-1) = park5_vlim(thr(c)-thr(c-1), thr(c+1)-thr(c))
        sxu(c-1) = park5_vlim(tu(c)-tu(c-1), tu(c+1)-tu(c))
        sxv(c-1) = park5_vlim(tv(c)-tv(c-1), tv(c+1)-tv(c))
        sxe(c-1) = park5_vlim(te(c)-te(c-1), te(c+1)-te(c))
        do s = 1, NP5
          sxy(s,c-1) = park5_vlim(ty(s,c)-ty(s,c-1), ty(s,c+1)-ty(s,c))
        end do
      end do
      do f = 1, nj+1
        rl = thr(f+1) + 0.5_RP*sxr(f);  rr = thr(f+2) - 0.5_RP*sxr(f+1)
        ! the serial lane swaps the velocities for the y faces: the
        ! normal component goes into the solver's normal slot.
        rul = tv(f+1) + 0.5_RP*sxv(f);  rur = tv(f+2) - 0.5_RP*sxv(f+1)
        vul = tu(f+1) + 0.5_RP*sxu(f);  vur = tu(f+2) - 0.5_RP*sxu(f+1)
        el = te(f+1) + 0.5_RP*sxe(f);   er_ = te(f+2) - 0.5_RP*sxe(f+1)
        do s = 1, NP5
          yleft(s) = ty(s,f+1) + 0.5_RP*sxy(s,f)
          yright(s) = ty(s,f+2) - 0.5_RP*sxy(s,f+1)
        end do
        call closure_triple(max(rl,1.0e-12_RP), el, yleft, pl, dummy, al)
        call closure_triple(max(rr,1.0e-12_RP), er_, yright, pr_, dummy, ar_)
        pl = max(pl, 1.0e-12_RP); pr_ = max(pr_, 1.0e-12_RP)
        al = max(al, 1.0e-12_RP); ar_ = max(ar_, 1.0e-12_RP)
        call hllc_face_x_pub(rl, rul, vul, pl, el, al, rr, rur, vur, pr_, er_, ar_, &
                             1.4_RP, fm, fqy, fqx, fe)
        fy4(1,f,i) = fm; fy4(2,f,i) = fqx; fy4(3,f,i) = fqy; fy4(4,f,i) = fe
        fym(f,i) = fm
      end do
    end do

    ! ===================== conservative update =====================
    do j = 1, nj
      if (annular) then
        rm = faces_y(j); rpv = faces_y(j+1)
        rc = 0.5_RP*(rm + rpv)
      end if
      do i = 1, nl
        k = (j-1)*nl + i
        gi = dec%i0 + i - 1
        if (annular) then
          ! annular: r-weighted radial divergence plus the p/rc source on
          ! radial momentum, grouped exactly like the serial axi oracle.
          rho_o(k) = rho(k) - dt*((fx4(1,i+1,j)-fx4(1,i,j))/dxs(i) &
                     + (rpv*fy4(1,j+1,i)-rm*fy4(1,j,i))/(rc*dys(j)))
          mx_o(k) = mx(k) - dt*((fx4(2,i+1,j)-fx4(2,i,j))/dxs(i) &
                   + (rpv*fy4(2,j+1,i)-rm*fy4(2,j,i))/(rc*dys(j)))
          my_o(k) = my(k) - dt*((fx4(3,i+1,j)-fx4(3,i,j))/dxs(i) &
                   + (rpv*fy4(3,j+1,i)-rm*fy4(3,j,i))/(rc*dys(j))) &
                   + dt*pf(k)/rc
          e_o(k) = e_tot(k) - dt*((fx4(4,i+1,j)-fx4(4,i,j))/dxs(i) &
                    + (rpv*fy4(4,j+1,i)-rm*fy4(4,j,i))/(rc*dys(j)))
        else
          rho_o(k) = rho(k) - dt*((fx4(1,i+1,j)-fx4(1,i,j))/dxs(i) &
                                  + (fy4(1,j+1,i)-fy4(1,j,i))/dys(j))
          mx_o(k) = mx(k) - dt*((fx4(2,i+1,j)-fx4(2,i,j))/dxs(i) &
                                + (fy4(2,j+1,i)-fy4(2,j,i))/dys(j))
          my_o(k) = my(k) - dt*((fx4(3,i+1,j)-fx4(3,i,j))/dxs(i) &
                                + (fy4(3,j+1,i)-fy4(3,j,i))/dys(j))
          e_o(k) = e_tot(k) - dt*((fx4(4,i+1,j)-fx4(4,i,j))/dxs(i) &
                                  + (fy4(4,j+1,i)-fy4(4,j,i))/dys(j))
        end if
      end do
    end do

    ! ===================== species transport =====================
    ! the transport is the serial lane's donor-cell upwind rule, laid on
    ! the slab. the face donors come from the padded strip, so the west
    ! ghost supplies the donor for the first seam face and the east ghost
    ! supplies it for the last one.
    !
    ! face f carries the donor from the left cell f and the right cell
    ! f+1, exactly as the serial lane reads them, with the serial lane's
    ! boundary convention at the two global ends: the outermost faces fall
    ! back to the interior cell.
    do j = 1, nj
      do f = 1, nf
        if (fxm(f,j) >= 0.0_RP) then
          do s = 1, NP5
            if (f > 1) then
              yx(s,f,j) = car(4+s,f+1,j)
            else if (dec%west_domain_edge) then
              yx(s,f,j) = car(4+s,3,j)
            else
              yx(s,f,j) = car(4+s,2,j)
            end if
          end do
        else
          do s = 1, NP5
            if (f < nf) then
              yx(s,f,j) = car(4+s,f+2,j)
            else if (dec%east_domain_edge) then
              yx(s,f,j) = car(4+s,nl+2,j)
            else
              yx(s,f,j) = car(4+s,nl+3,j)
            end if
          end do
        end if
      end do
    end do
    do i = 1, nl
      do f = 1, nj+1
        if (fym(f,i) >= 0.0_RP) then
          do s = 1, NP5
            if (f > 1) then
              yy(s,f,i) = car(4+s,i+2,f-1)
            else
              yy(s,f,i) = car(4+s,i+2,1)
            end if
          end do
        else
          do s = 1, NP5
            if (f < nj+1) then
              yy(s,f,i) = car(4+s,i+2,f)
            else
              yy(s,f,i) = car(4+s,i+2,nj)
            end if
          end do
        end if
      end do
    end do

    do s = 1, NP5
      do j = 1, nj
        if (annular) then
          rm = faces_y(j); rpv = faces_y(j+1)
          rc = 0.5_RP*(rm + rpv)
        end if
        do i = 1, nl
          k = (j-1)*nl + i
          ! the serial lane's arithmetic: the mass flux difference across
          ! the cell's two x faces and its two y faces, weighted by the
          ! upwind donors and scaled by the updated density.
          flux_x = yx(s,i+1,j)*fxm(i+1,j) - yx(s,i,j)*fxm(i,j)
          dxi = dt / dxs(i)
          if (annular) then
            flux_y = (rpv*yy(s,j+1,i)*fym(j+1,i) - rm*yy(s,j,i)*fym(j,i))/(rc*dys(j))
            dyi = dt
          else
            flux_y = yy(s,j+1,i)*fym(j+1,i) - yy(s,j,i)*fym(j,i)
            dyi = dt / dys(j)
          end if
          rho_y = rho_o(k) * car(4+s,i+2,j)
          dflux = -(dxi*flux_x + dyi*flux_y)
          new_rho_y = rho_y + dflux
          if (new_rho_y > 0.0_RP) then
            ys_o((s-1)*nn + k) = new_rho_y / rho_o(k)
          else
            ys_o((s-1)*nn + k) = 0.0_RP
          end if
        end do
      end do
    end do

    ! opt-in laminar pass, the same shared kernel the serial oracle calls.
    ! the kernel's x-face stencil reaches one cell past the slab, so the
    ! post-transport primitives are repacked into the padded carrier, the
    ! seam halo is refreshed, and the neighbour edge widths ride along in
    ! a scalar exchange; global domain edges pass zero flux as before.
    if (visc /= 0) then
      do j = 1, nj
        do i = 1, nl
          k = (j-1)*nl + i
          pri = max(rho_o(k), 1.0e-12_RP)
          car(1,i+2,j) = pri
          car(2,i+2,j) = mx_o(k) / pri
          car(3,i+2,j) = my_o(k) / pri
          car(4,i+2,j) = e_o(k) / pri - 0.5_RP*(car(2,i+2,j)**2 + car(3,i+2,j)**2)
          do s = 1, NP5
            car(4+s,i+2,j) = ys_o((s-1)*nn + k)
          end do
        end do
      end do
      call halo_exchange_x(car, NCAR, dec)
      allocate(cwl(NCAR,nj), cel(NCAR,nj))
      do j = 1, nj
        cwl(:,j) = car(:,2,j)
        cel(:,j) = car(:,nl+3,j)
      end do
      ! west neighbour's last cell width and east neighbour's first.
      dwx = 0.0_RP; dex = 0.0_RP
      call mpi_sendrecv(dxs(nl), 1, MPI_DOUBLE_PRECISION, dec%east_rank, 103, &
                        dwx, 1, MPI_DOUBLE_PRECISION, dec%west_rank, 103, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_sendrecv(dxs(1), 1, MPI_DOUBLE_PRECISION, dec%west_rank, 104, &
                        dex, 1, MPI_DOUBLE_PRECISION, dec%east_rank, 104, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call park5_viscous_update(nl, nj, rho_o, mx_o, my_o, e_o, ys_o, dxs, dys, dt, &
                                dec%west_domain_edge, dec%east_domain_edge, &
                                dwx, dex, cwl, cel)
      deallocate(cwl, cel)
    end if

    ! ===================== normalisation and chemistry =====================
    do k = 1, nn
      ysum = 0.0_RP
      do s = 1, NP5
        ysum = ysum + ys_o((s-1)*nn + k)
      end do
      if (ysum > 1.0e-30_RP) then
        do s = 1, NP5
          ys_o((s-1)*nn + k) = ys_o((s-1)*nn + k) / ysum
        end do
      end if
      rhe(k) = max(rho_o(k), 1.0e-12_RP)
      dee(k) = e_o(k)/rhe(k) &
               - 0.5_RP*((mx_o(k)/rhe(k))**2 + (my_o(k)/rhe(k))**2)
      do s = 1, NP5
        ycw(s,k) = ys_o((s-1)*nn + k)
      end do
    end do
    call chem_substep_grid(nn, rhe, dee, dt, ycw, worst_dy, worst_ov)
    do k = 1, nn
      do s = 1, NP5
        ys_o((s-1)*nn + k) = max(ycw(s,k), 0.0_RP)
      end do
    end do

    deallocate(car, fx4, fy4, fxm, fym)
    deallocate(sxr, sxu, sxv, sxe, sxy)
    deallocate(thr, tu, tv, te, ty)
    deallocate(yx, yy, rhe, dee, ycw, pf)
    err = 0
  end subroutine nfor_mpipark5_march

  ! ------------------------------------------------------------------
  ! fill the domain-edge ghost columns where this rank really is a
  ! global boundary. the seam ghosts came from the halo. the codes are
  ! 0 transmissive, 1 inflow, 2 outflow, 3 slip wall, 4 no-slip wall.
  subroutine apply_domain_edges(car, nl, nj, dec, bc_w, bc_e, bcv, yf)
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
  end subroutine apply_domain_edges

end module nuforpark5mpi

! nuforvisc.f90 — viscous slab contributions: molecular lambda/mu per cell
! (direct-read and sagutherland-style) and the staggered-grid flux-divergence
! adds for the conservative update, with no-slip wall shear and conduction.
module nuforvisc
  use iso_fortran_env, only: real64
  use nufoeqair, only: cea_t_at
  implicit none
  integer, parameter :: DP = real64

contains


  ! the viscous half-step on this rank's padded slab, mirroring
  ! viscous2d::add_viscous_cells_mu on a uniform grid: closure
  ! temperature from the cea table, sutherland mu per cell, central
  ! cell gradients (y mirrored at the domain rows, x from the halo
  ! which carries the pad2d copy semantics), face-averaged gradients,
  ! no-slip wall shear + isothermal-wall conduction, divergence add.
  ! call right after a halo refresh of the post-euler state.
  subroutine viscous_add_slab(rho, mx, my, e, nx_l, ny_, dx, dy, dt, &
                              mu0, t0_ref, s_par, pr, tw, gamma, r_gas, &
                              bc_s, bc_n, solid)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(inout) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                               my(nx_l+4, ny_), e(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, dt, mu0, t0_ref, s_par, pr, tw, gamma, r_gas
    integer, intent(in) :: bc_s, bc_n
    logical, intent(in) :: solid(nx_l+4, ny_)
    real(DP), allocatable :: u(:,:), v(:,:), tt(:,:), mu(:,:)
    real(DP), allocatable :: fxv(:,:,:), fyv(:,:,:), gu(:,:,:)
    real(DP) :: ein, kap, dtdn, q_wall
    real(DP) :: gux, guy, gvx, gvy, gtx, gty
    real(DP) :: muf, divw, txx, txy, tyy, q
    integer :: i, j, f, k, a, b

    if (mu0 <= 0.0_DP) return
    ! kappa = mu-independent cp/pr with real units, mirroring the rust
    ! closure lane: k = gamma r_gas / ((gamma-1) pr).
    kap = gamma * r_gas / ((gamma - 1.0_DP) * pr)

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), tt(nx_l+4, ny_), mu(nx_l+4, ny_))
    allocate(fxv(4, 0:nx_l+2, ny_), fyv(3, 0:nx_l, 0:ny_+1))
    allocate(gu(6, nx_l+4, ny_))

    ! primitives + closure temperature + sutherland mu per padded cell
    do j = 1, ny_
      do i = 1, nx_l+4
        u(i,j) = mx(i,j) / max(rho(i,j), 1.0e-12_DP)
        v(i,j) = my(i,j) / max(rho(i,j), 1.0e-12_DP)
        ein = e(i,j) / max(rho(i,j), 1.0e-12_DP) - 0.5_DP * (u(i,j)**2 + v(i,j)**2)
        call cea_t_at(rho(i,j), ein, tt(i,j))
        tt(i,j) = max(tt(i,j), 1.0_DP)
        mu(i,j) = mu0 * (tt(i,j) / t0_ref)**1.5_DP * (t0_ref + s_par) / (tt(i,j) + s_par)
      end do
    end do

    ! cell-centred gradients with the wall override on no-slip sides,
    ! stored per padded cell: gu(6, i, j) = ux, uy, vx, vy, tx, ty.
    ! y rows 1 and ny_ mirror their edge row (pad2d); no-slip wall
    ! cells then override uy/vy with the quadratic-consistent stencil.
    do j = 1, ny_
      do i = 2, nx_l+3
        gu(1,i,j) = (u(i+1,j) - u(i-1,j)) / (2.0_DP * dx)
        gu(3,i,j) = (v(i+1,j) - v(i-1,j)) / (2.0_DP * dx)
        gu(5,i,j) = (tt(i+1,j) - tt(i-1,j)) / (2.0_DP * dx)
        if (j == 1) then
          ! mirrored ghost row (pad2d): (u2 - u1)/(2 dy)
          gu(2,i,j) = (u(i,2) - u(i,1)) / (2.0_DP * dy)
          gu(4,i,j) = (v(i,2) - v(i,1)) / (2.0_DP * dy)
          gu(6,i,j) = (tt(i,2) - tt(i,1)) / (2.0_DP * dy)
        else if (j == ny_) then
          gu(2,i,j) = (u(i,ny_) - u(i,ny_-1)) / (2.0_DP * dy)
          gu(4,i,j) = (v(i,ny_) - v(i,ny_-1)) / (2.0_DP * dy)
          gu(6,i,j) = (tt(i,ny_) - tt(i,ny_-1)) / (2.0_DP * dy)
        else
          gu(2,i,j) = (u(i,j+1) - u(i,j-1)) / (2.0_DP * dy)
          gu(4,i,j) = (v(i,j+1) - v(i,j-1)) / (2.0_DP * dy)
          gu(6,i,j) = (tt(i,j+1) - tt(i,j-1)) / (2.0_DP * dy)
        end if
      end do
      ! no-slip wall-cell override: uy = flip*(3 u0 + u1)/(3 dy)
      if (bc_s == 4) then
        do i = 2, nx_l+3
          gu(2,i,1) = (3.0_DP * u(i,1) + u(i,2)) / (3.0_DP * dy)
          gu(4,i,1) = (3.0_DP * v(i,1) + v(i,2)) / (3.0_DP * dy)
        end do
      end if
      if (bc_n == 4) then
        do i = 2, nx_l+3
          gu(2,i,ny_) = -(3.0_DP * u(i,ny_) + u(i,ny_-1)) / (3.0_DP * dy)
          gu(4,i,ny_) = -(3.0_DP * v(i,ny_) + v(i,ny_-1)) / (3.0_DP * dy)
        end do
      end if
    end do

    ! x-face fluxes: face k between padded cells k+1 | k+2, k=1..nx_l+1
    ! (k=0 and k=nx_l+2 are the ghost-edge faces: zero, matching the
    ! rust pass whose domain-edge faces carry zero flux).
    fxv = 0.0_DP
    do j = 1, ny_
      do k = 1, nx_l+1
        a = k+1; b = k+2
        muf = 0.5_DP * (mu(a,j) + mu(b,j))
        gux = 0.5_DP * (gu(1,a,j) + gu(1,b,j))
        guy = 0.5_DP * (gu(2,a,j) + gu(2,b,j))
        gvx = 0.5_DP * (gu(3,a,j) + gu(3,b,j))
        gvy = 0.5_DP * (gu(4,a,j) + gu(4,b,j))
        gtx = 0.5_DP * (gu(5,a,j) + gu(5,b,j))
        divw = gux + gvy
        txx = muf * (gux * 2.0_DP - divw * 2.0_DP / 3.0_DP)
        txy = muf * (guy + gvx)
        q = muf * kap * gtx
        fxv(1,k,j) = txx
        fxv(2,k,j) = txy
        fxv(4,k,j) = u(a,j) * txx + v(a,j) * txy + q
      end do
    end do

    ! y-face fluxes: face f below row f (between rows f-1 | f),
    ! f=1..ny_+1 (f=1 ghost edge, f=ny_+1 north edge); interior
    ! f=2..ny_. x slots 1..nx_l map to fyv(:, i-2, :) for cell i.
    fyv = 0.0_DP
    do i = 3, nx_l+2
      do f = 2, ny_
        muf = 0.5_DP * (mu(i,f-1) + mu(i,f))
        gux = 0.5_DP * (gu(1,i,f-1) + gu(1,i,f))
        guy = 0.5_DP * (gu(2,i,f-1) + gu(2,i,f))
        gvx = 0.5_DP * (gu(3,i,f-1) + gu(3,i,f))
        gvy = 0.5_DP * (gu(4,i,f-1) + gu(4,i,f))
        gty = 0.5_DP * (gu(6,i,f-1) + gu(6,i,f))
        divw = gux + gvy
        txy = muf * (guy + gvx)
        tyy = muf * (gvy * 2.0_DP - divw * 2.0_DP / 3.0_DP)
        q = muf * kap * gty
        fyv(1,i-2,f) = txy
        fyv(2,i-2,f) = tyy
        fyv(3,i-2,f) = u(i,f-1) * txy + v(i,f-1) * tyy + q
      end do
      ! south no-slip wall face f=1: one-sided shear (9 u0 - u1)/(3 dy)
      if (bc_s == 4) then
        muf = mu(i,1)
        txy = muf * (9.0_DP * u(i,1) - u(i,2)) / (3.0_DP * dy)
        tyy = muf * (9.0_DP * v(i,1) - v(i,2)) / (3.0_DP * dy)
        q_wall = 0.0_DP
        if (tw > 0.0_DP) then
          dtdn = (-3.0_DP * tw + 4.0_DP * tt(i,1) - tt(i,2)) / (2.0_DP * dy)
          q_wall = muf * kap * dtdn
        end if
        fyv(1,i-2,1) = txy
        fyv(2,i-2,1) = tyy
        fyv(3,i-2,1) = u(i,1) * txy + v(i,1) * tyy + q_wall
      end if
      ! north no-slip wall face f=ny_+1: outward normal -y negates.
      if (bc_n == 4) then
        muf = mu(i,ny_)
        txy = muf * (9.0_DP * u(i,ny_) - u(i,ny_-1)) / (3.0_DP * dy)
        tyy = muf * (9.0_DP * v(i,ny_) - v(i,ny_-1)) / (3.0_DP * dy)
        q_wall = 0.0_DP
        if (tw > 0.0_DP) then
          dtdn = (3.0_DP * tw - 4.0_DP * tt(i,ny_) + tt(i,ny_-1)) / (2.0_DP * dy)
          q_wall = muf * kap * dtdn
        end if
        fyv(1,i-2,ny_+1) = -txy
        fyv(2,i-2,ny_+1) = -tyy
        fyv(3,i-2,ny_+1) = -(u(i,ny_) * txy + v(i,ny_) * tyy + q_wall)
      end if
    end do

    ! divergence into interior cells: cell (i,j) adds
    ! dt/dx * (fx(:,i-1,j) - fx(:,i-2,j)) + dt/dy * (fy(:,i-2,j+1) - fy(:,i-2,j))
    do j = 1, ny_
      do i = 3, nx_l+2
        if (solid(i,j)) cycle
        mx(i,j) = mx(i,j) + dt / dx * (fxv(1,i-1,j) - fxv(1,i-2,j)) &
                            + dt / dy * (fyv(1,i-2,j+1) - fyv(1,i-2,j))
        my(i,j) = my(i,j) + dt / dx * (fxv(2,i-1,j) - fxv(2,i-2,j)) &
                            + dt / dy * (fyv(2,i-2,j+1) - fyv(2,i-2,j))
        e(i,j) = e(i,j) + dt / dx * (fxv(4,i-1,j) - fxv(4,i-2,j)) &
                       + dt / dy * (fyv(3,i-2,j+1) - fyv(3,i-2,j))
      end do
    end do

    deallocate(u, v, tt, mu, fxv, fyv, gu)
  end subroutine


  ! per-cell molecular viscosity, direct-read route (fused_state_mu):
  ! t from the cea table at (rho, e), then sutherland.
  subroutine lam_mu_slab(rho, ein, nx_l, ny_, mu0, t0_ref, s_par, mu_out)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho(nx_l+4, ny_), ein(nx_l+4, ny_)
    real(DP), intent(in) :: mu0, t0_ref, s_par
    real(DP), intent(out) :: mu_out(nx_l+4, ny_)
    real(DP) :: tt
    integer :: i, j
    do j = 1, ny_
      do i = 1, nx_l+4
        call cea_t_at(rho(i,j), ein(i,j), tt)
        tt = max(tt, 1.0_DP)
        mu_out(i,j) = mu0 * (tt / t0_ref)**1.5_DP * (t0_ref + s_par) / (tt + s_par)
      end do
    end do
  end subroutine


  ! the viscous half-step with an optional eddy-viscosity coupling,
  ! mirroring viscous2d::add_viscous_cells_mu on the uniform branch:
  ! closure temperature from the cea table, central cell gradients on
  ! pad2d/pad2d_vel ghosts (no-slip y flip; the second row's centered
  ! stencil reads the flipped ghost exactly as the rust padded arrays
  ! do), wall-cell uy/vy override, face-averaged face velocity in the
  ! energy flux, zero flux on the x domain-edge faces, one-sided
  ! no-slip wall shear with the eddy mu folded in, molecular kappa in
  ! every heat-flux term. mu_cells arrives from the caller at the
  ! state the serial lane evaluated it at (step entry), mu_t from the
  ! caller's eddy pass (empty flag = laminar).
  subroutine viscous_add_slab2(rho, mx, my, e, nx_l, ny_, dx, dy, dt, &
                               mu_cells, mu_t, has_turb, mu0, t0_ref, s_par, &
                               pr, pr_t, tw, gamma, r_gas, bc_s, bc_n, &
                               west_dom, east_dom)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(inout) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                               my(nx_l+4, ny_), e(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, dt
    real(DP), intent(in) :: mu_cells(nx_l+4, ny_)
    real(DP), intent(in) :: mu_t(nx_l+4, ny_)
    logical, intent(in) :: has_turb
    real(DP), intent(in) :: mu0, t0_ref, s_par, pr, pr_t, tw, gamma, r_gas
    integer, intent(in) :: bc_s, bc_n
    logical, intent(in) :: west_dom, east_dom
    real(DP), allocatable :: u(:,:), v(:,:), tt(:,:)
    real(DP), allocatable :: u2(:,:), v2(:,:), t2(:,:)
    real(DP), allocatable :: fxv(:,:,:), fyv(:,:,:), gu(:,:,:)
    real(DP) :: ein, kap, kap_t, dtdn, q_wall
    real(DP) :: gux, guy, gvx, gvy, gtx, gty
    real(DP) :: muf, mu_tf, divw, txx, txy, tyy, q, uf, vf
    integer :: i, j, f, k

    if (mu0 <= 0.0_DP) return
    kap = gamma * r_gas / ((gamma - 1.0_DP) * pr)
    kap_t = gamma / ((gamma - 1.0_DP) * pr_t)

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), tt(nx_l+4, ny_))
    allocate(u2(nx_l+4, ny_+2), v2(nx_l+4, ny_+2), t2(nx_l+4, ny_+2))
    allocate(fxv(4, 0:nx_l+2, ny_), fyv(3, 0:nx_l, 0:ny_+1))
    allocate(gu(6, nx_l+4, ny_))

    ! primitives + closure temperature per padded cell
    do j = 1, ny_
      do i = 1, nx_l+4
        u(i,j) = mx(i,j) / max(rho(i,j), 1.0e-12_DP)
        v(i,j) = my(i,j) / max(rho(i,j), 1.0e-12_DP)
        ein = e(i,j) / max(rho(i,j), 1.0e-12_DP) &
            - 0.5_DP * (u(i,j)**2 + v(i,j)**2)
        call cea_t_at(rho(i,j), ein, tt(i,j))
        tt(i,j) = max(tt(i,j), 1.0_DP)
      end do
    end do

    ! padded rows: pad2d copies for tt; pad2d_vel for u, v (no-slip
    ! flip on the y ghosts). the x columns carry the state's own ghost
    ! columns (copies at domain edges, true neighbors at seams).
    do i = 1, nx_l+4
      u2(i,1) = u(i,1); v2(i,1) = v(i,1); t2(i,1) = tt(i,1)
      u2(i,ny_+2) = u(i,ny_); v2(i,ny_+2) = v(i,ny_); t2(i,ny_+2) = tt(i,ny_)
      if (bc_s == 4) then
        u2(i,1) = -u(i,1); v2(i,1) = -v(i,1)
      end if
      if (bc_n == 4) then
        u2(i,ny_+2) = -u(i,ny_); v2(i,ny_+2) = -v(i,ny_)
      end if
      do j = 1, ny_
        u2(i,j+1) = u(i,j); v2(i,j+1) = v(i,j); t2(i,j+1) = tt(i,j)
      end do
    end do

    ! cell-centred gradients over the padded range the face loops read
    ! (cols 2..nx_l+3: the seam faces at np>1 average neighbor gradients,
    ! and the domain-edge faces never enter the divergence but the
    ! arrays must not carry garbage anywhere they are read)
    do j = 1, ny_
      do i = 2, nx_l+3
        gu(1,i,j) = (u2(i+1,j+1) - u2(i-1,j+1)) / (2.0_DP * dx)
        gu(3,i,j) = (v2(i+1,j+1) - v2(i-1,j+1)) / (2.0_DP * dx)
        gu(5,i,j) = (t2(i+1,j+1) - t2(i-1,j+1)) / (2.0_DP * dx)
        gu(2,i,j) = (u2(i,j+2) - u2(i,j)) / (2.0_DP * dy)
        gu(4,i,j) = (v2(i,j+2) - v2(i,j)) / (2.0_DP * dy)
        gu(6,i,j) = (t2(i,j+2) - t2(i,j)) / (2.0_DP * dy)
      end do
    end do
    ! no-slip wall-cell uy/vy override (uniform wall_cell_grad). the
    ! ghost columns carry the neighbor rank's wall cells, whose grads
    ! rust also overrides, so the seam face averages two overridden
    ! values exactly like the serial lane.
    if (bc_s == 4) then
      do i = 1, nx_l+4
        gu(2,i,1) = (3.0_DP * u(i,1) + u(i,2)) / (3.0_DP * dy)
        gu(4,i,1) = (3.0_DP * v(i,1) + v(i,2)) / (3.0_DP * dy)
      end do
    end if
    if (bc_n == 4) then
      do i = 1, nx_l+4
        gu(2,i,ny_) = -(3.0_DP * u(i,ny_) + u(i,ny_-1)) / (3.0_DP * dy)
        gu(4,i,ny_) = -(3.0_DP * v(i,ny_) + v(i,ny_-1)) / (3.0_DP * dy)
      end do
    end if

    ! x-face fluxes: face k sits between padded cols k+1 | k+2. rust
    ! faces f=1..nx-1 map to interior-interior pairs; at np>1 the slab
    ! seam faces (k=1 west, k=nx_l+1 east) are true interior faces of
    ! the global grid and carry flux (the halo holds the neighbor); the
    ! x domain-edge faces stay zero exactly like rust, whose viscous
    ! x-face loop runs 1..nx-1 and never touches them.
    fxv = 0.0_DP
    do j = 1, ny_
      do k = 1, nx_l+1
        if (west_dom .and. k == 1) cycle
        if (east_dom .and. k == nx_l+1) cycle
        i = k + 1
        muf = 0.5_DP * (mu_cells(i,j) + mu_cells(i+1,j))
        if (has_turb) then
          mu_tf = 0.5_DP * (mu_t(i,j) + mu_t(i+1,j))
        else
          mu_tf = 0.0_DP
        end if
        gux = 0.5_DP * (gu(1,i,j) + gu(1,i+1,j))
        guy = 0.5_DP * (gu(2,i,j) + gu(2,i+1,j))
        gvx = 0.5_DP * (gu(3,i,j) + gu(3,i+1,j))
        gvy = 0.5_DP * (gu(4,i,j) + gu(4,i+1,j))
        gtx = 0.5_DP * (gu(5,i,j) + gu(5,i+1,j))
        divw = gux + gvy
        txx = (muf + mu_tf) * (gux * 2.0_DP - divw * 2.0_DP / 3.0_DP)
        txy = (muf + mu_tf) * (guy + gvx)
        q = (muf * kap + mu_tf * kap_t) * gtx
        uf = 0.5_DP * (u(i,j) + u(i+1,j))
        vf = 0.5_DP * (v(i,j) + v(i+1,j))
        fxv(1,k,j) = txx
        fxv(2,k,j) = txy
        fxv(4,k,j) = uf * txx + vf * txy + q
      end do
    end do

    ! y-face fluxes: interior f=2..ny; south face f=1 and north face
    ! f=ny_+1 carry the one-sided no-slip wall treatment, else zero.
    fyv = 0.0_DP
    do i = 3, nx_l+2
      do f = 2, ny_
        muf = 0.5_DP * (mu_cells(i,f-1) + mu_cells(i,f))
        if (has_turb) then
          mu_tf = 0.5_DP * (mu_t(i,f-1) + mu_t(i,f))
        else
          mu_tf = 0.0_DP
        end if
        gux = 0.5_DP * (gu(1,i,f-1) + gu(1,i,f))
        guy = 0.5_DP * (gu(2,i,f-1) + gu(2,i,f))
        gvx = 0.5_DP * (gu(3,i,f-1) + gu(3,i,f))
        gvy = 0.5_DP * (gu(4,i,f-1) + gu(4,i,f))
        gty = 0.5_DP * (gu(6,i,f-1) + gu(6,i,f))
        divw = gux + gvy
        txy = (muf + mu_tf) * (guy + gvx)
        tyy = (muf + mu_tf) * (gvy * 2.0_DP - divw * 2.0_DP / 3.0_DP)
        q = (muf * kap + mu_tf * kap_t) * gty
        uf = 0.5_DP * (u(i,f-1) + u(i,f))
        vf = 0.5_DP * (v(i,f-1) + v(i,f))
        fyv(1,i-2,f) = txy
        fyv(2,i-2,f) = tyy
        fyv(3,i-2,f) = uf * txy + vf * tyy + q
      end do
      if (bc_s == 4) then
        if (has_turb) then
          muf = mu_cells(i,1) + mu_t(i,1)
        else
          muf = mu_cells(i,1)
        end if
        txy = muf * (9.0_DP * u(i,1) - u(i,2)) / (3.0_DP * dy)
        if (i == 3) then
        end if
        if (i == 49 .or. i == 50) then
        end if
        tyy = muf * (9.0_DP * v(i,1) - v(i,2)) / (3.0_DP * dy)
        q_wall = 0.0_DP
        if (tw > 0.0_DP) then
          dtdn = (-3.0_DP * tw + 4.0_DP * tt(i,1) - tt(i,2)) / (2.0_DP * dy)
          q_wall = muf * kap * dtdn
        end if
        fyv(1,i-2,1) = txy
        fyv(2,i-2,1) = tyy
        fyv(3,i-2,1) = u(i,1) * txy + v(i,1) * tyy + q_wall
      end if
      if (bc_n == 4) then
        if (has_turb) then
          muf = mu_cells(i,ny_) + mu_t(i,ny_)
        else
          muf = mu_cells(i,ny_)
        end if
        txy = muf * (9.0_DP * u(i,ny_) - u(i,ny_-1)) / (3.0_DP * dy)
        tyy = muf * (9.0_DP * v(i,ny_) - v(i,ny_-1)) / (3.0_DP * dy)
        q_wall = 0.0_DP
        if (tw > 0.0_DP) then
          dtdn = (3.0_DP * tw - 4.0_DP * tt(i,ny_) + tt(i,ny_-1)) / (2.0_DP * dy)
          q_wall = muf * kap * dtdn
        end if
        fyv(1,i-2,ny_+1) = -txy
        fyv(2,i-2,ny_+1) = -tyy
        fyv(3,i-2,ny_+1) = -(u(i,ny_) * txy + v(i,ny_) * tyy + q_wall)
      end if
    end do

    ! divergence into interior cells
    do j = 1, ny_
      do i = 3, nx_l+2
        mx(i,j) = mx(i,j) + dt / dx * (fxv(1,i-1,j) - fxv(1,i-2,j)) &
                            + dt / dy * (fyv(1,i-2,j+1) - fyv(1,i-2,j))
        my(i,j) = my(i,j) + dt / dx * (fxv(2,i-1,j) - fxv(2,i-2,j)) &
                            + dt / dy * (fyv(2,i-2,j+1) - fyv(2,i-2,j))
        e(i,j) = e(i,j) + dt / dx * (fxv(4,i-1,j) - fxv(4,i-2,j)) &
                       + dt / dy * (fyv(3,i-2,j+1) - fyv(3,i-2,j))
      end do
    end do

    deallocate(u, v, tt, u2, v2, t2, fxv, fyv, gu)
  end subroutine

end module nuforvisc

! nufoorturb.f90 — spalart-allmaras turbulence (sa.rs 1994 formulation): the
! closure constants, the per-cell mu/sigma quantities and the nu-tilde
! transport sweep. shipped equation-of-state scalars come from nufoeqair.
module nufoorturb
  use iso_fortran_env, only: real64
  use nufoeqair, only: cea_pt_at, cea_t_at
  implicit none
  integer, parameter :: DP = real64

  real(DP), parameter :: CB1_SA   = 0.1355_DP
  real(DP), parameter :: CB2_SA   = 0.622_DP
  real(DP), parameter :: SIGMA_SA = 2.0_DP / 3.0_DP
  real(DP), parameter :: KAPPA_SA = 0.41_DP
  real(DP), parameter :: CW1_SA   = CB1_SA / (KAPPA_SA*KAPPA_SA) + (1.0_DP + CB2_SA) / SIGMA_SA
  real(DP), parameter :: CW2_SA   = 0.3_DP
  real(DP), parameter :: CW3_SA   = 2.0_DP
  real(DP), parameter :: CV1_SA   = 7.1_DP
  real(DP), parameter :: CS_SA    = 0.3_DP
  real(DP), parameter :: RMAX_SA  = 10.0_DP

contains


  ! sa source closures, one cell: mirrors sa.rs stilde/r/g/fw/source in the
  ! same operation order (stilde floor, r clip, g, fw powf).
  subroutine sa_source_closure(nt_c, s, dd, nu, s_tilde, rr, prod, destr, src)
    real(DP), intent(in) :: nt_c, s, dd, nu
    real(DP), intent(out) :: s_tilde, rr, prod, destr, src
    real(DP) :: d2, chi_v, c3, fv1v, fv2v, st_kd, r6, g6, gg, c36
    d2 = max(dd*dd, 1.0e-24_DP)
    chi_v = nt_c / nu
    c3 = chi_v * chi_v * chi_v
    fv1v = c3 / (c3 + CV1_SA*CV1_SA*CV1_SA)
    fv2v = 1.0_DP - chi_v / (1.0_DP + chi_v * fv1v)
    s_tilde = max(s + nt_c * fv2v / (KAPPA_SA * KAPPA_SA * d2), CS_SA * s)
    rr = min(nt_c / (max(s_tilde, 1.0e-10_DP) * KAPPA_SA * KAPPA_SA * d2), RMAX_SA)
    r6 = rr**6
    gg = rr + CW2_SA * (r6 - rr)
    g6 = gg**6
    c36 = CW3_SA**6
    prod = CB1_SA * s_tilde * nt_c
    destr = CW1_SA * gg * ((1.0_DP + c36) / (g6 + c36))**(1.0_DP/6.0_DP) * nt_c * nt_c / d2
    src = prod - destr
  end subroutine


  ! per-cell molecular viscosity, bisection route (the sa lane's
  ! thermo::pressure + thermo::temperature): p from the table at the
  ! cell's (rho, e), then 24 log-e bisection iterations recover the
  ! e whose p matches, then t from that e, then sutherland.
  subroutine sa_mu_slab(rho, ein, nx_l, ny_, mu0, t0_ref, s_par, mu_out)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho(nx_l+4, ny_), ein(nx_l+4, ny_)
    real(DP), intent(in) :: mu0, t0_ref, s_par
    real(DP), intent(out) :: mu_out(nx_l+4, ny_)
    real(DP) :: pp, tt, target_p, lo, hi, mid, emid, p_mid, best_e, r
    integer :: i, j, it
    do j = 1, ny_
      do i = 1, nx_l+4
        r = max(rho(i,j), 1.0e-12_DP)
        call cea_pt_at(r, ein(i,j), pp)
        target_p = max(pp, 1.0e-3_DP)
        lo = 0.005_DP
        hi = 9.5_DP
        best_e = 1.0_DP
        do it = 1, 24
          mid = 0.5_DP * (lo + hi)
          emid = 10.0_DP**mid
          call cea_pt_at(r, emid, p_mid)
          if (p_mid < target_p) then
            lo = mid
          else
            hi = mid
          end if
          best_e = emid
        end do
        call cea_t_at(r, best_e, tt)
        tt = max(tt, 1.0_DP)
        mu_out(i,j) = mu0 * (tt / t0_ref)**1.5_DP * (t0_ref + s_par) / (tt + s_par)
      end do
    end do
  end subroutine


  ! eddy viscosity per padded cell: mu_t = rho nu_tilde fv1(chi),
  ! chi = nu_tilde / (mu/rho) — sa.rs eddy_viscosity op order.
  subroutine eddy_mu_t_slab(nt, rho, mu_lam, nx_l, ny_, mu_t_out)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: nt(nx_l+4, ny_), rho(nx_l+4, ny_), mu_lam
    real(DP), intent(out) :: mu_t_out(nx_l+4, ny_)
    real(DP) :: chi_v, c3
    integer :: i, j
    do j = 1, ny_
      do i = 1, nx_l+4
        chi_v = nt(i,j) / (mu_lam / max(rho(i,j), 1.0e-12_DP))
        c3 = chi_v * chi_v * chi_v
        mu_t_out(i,j) = rho(i,j) * nt(i,j) * c3 / (c3 + CV1_SA*CV1_SA*CV1_SA)
      end do
    end do
  end subroutine


  ! the sa transport step on this rank's slab, mirroring
  ! turb2d::advance_turb on the uniform branch op for op. nt arrives
  ! with its x ghosts already rule-filled and halo-exchanged by the
  ! caller; the y ghosts follow pad_scalar here (solid walls zero,
  ! open sides copy). the mean-flow (rho,u,v) slab is the padded
  ! post-viscous state; y ghost rows for u/v follow pad2d_vel
  ! (no-slip flip, open copy). dd is the per-cell wall distance on
  ! the local interior. advection is conservative in rho*nu_tilde
  ! (divided by rho), the domain-edge faces are single-sided on the
  ! mean-flow velocity exactly as the rust f==0 / f==nx branches.
  subroutine sa_transport_slab(nt, rho, u, v, dd, nx_l, ny_, dx, dy, dt, &
                               mu_lam, bc_s, bc_n)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(inout) :: nt(nx_l+4, ny_)
    real(DP), intent(in) :: rho(nx_l+4, ny_), u(nx_l+4, ny_), v(nx_l+4, ny_)
    real(DP), intent(in) :: dd(nx_l, ny_)
    real(DP), intent(in) :: dx, dy, dt, mu_lam
    integer, intent(in) :: bc_s, bc_n
    real(DP), allocatable :: pnt(:,:), up(:,:), vp(:,:)
    real(DP) :: nt_c, nu_c, src, rhs, s_vort
    real(DP) :: adv_x_w, adv_x_e, adv_y_s, adv_y_n
    real(DP) :: dif_x_w, dif_x_e, dif_y_s, dif_y_n
    real(DP) :: rho_f, v_f, u_f, nu_f, nt_f, nt_up
    real(DP) :: left, right, below, above
    real(DP) :: gx, gy, cb2_term
    real(DP) :: s_tilde, rr, prod, destr
    integer :: i, j, c, jp
    logical :: wall_s, wall_n

    if (dt <= 0.0_DP) return
    wall_s = (bc_s == 3 .or. bc_s == 4)
    wall_n = (bc_n == 3 .or. bc_n == 4)

    allocate(pnt(nx_l+4, ny_+2)); pnt = 0.0_DP
    allocate(up(nx_l+4, ny_+2));  up = 0.0_DP
    allocate(vp(nx_l+4, ny_+2));  vp = 0.0_DP
    do j = 1, ny_
      do i = 1, nx_l+4
        pnt(i,j+1) = nt(i,j)
        up(i,j+1) = u(i,j)
        vp(i,j+1) = v(i,j)
      end do
    end do
    ! y ghost rows: pad_scalar for nt (wall zero / copy), pad2d_vel for
    ! u, v (no-slip flip / copy), pad2d for rho (always copy — rho is
    ! read through the state slab directly at rows 1..ny_, and the face
    ! averages at the edge faces are single-sided there, so no padded
    ! rho is needed).
    do i = 3, nx_l+2
      if (wall_s) then
        pnt(i,1) = 0.0_DP
      else
        pnt(i,1) = pnt(i,2)
      end if
      if (wall_n) then
        pnt(i,ny_+2) = 0.0_DP
      else
        pnt(i,ny_+2) = pnt(i,ny_+1)
      end if
      ! rust vorticity() pads u, v with plain copies (pad, not pad2d_vel),
      ! so the centered vorticity stencil at the wall row sees the
      ! mirrored velocity, not the no-slip flip.
      up(i,1) = up(i,2); vp(i,1) = vp(i,2)
      up(i,ny_+2) = up(i,ny_+1); vp(i,ny_+2) = vp(i,ny_+1)
    end do

    do j = 1, ny_
      do i = 3, nx_l+2
        c = i
        jp = j + 1
        nt_c = nt(c,j)
        nu_c = mu_lam / max(rho(c,j), 1.0e-12_DP)

        ! west x face (between cols c-1 | c)
        left = pnt(c-1,jp); right = pnt(c,jp)
        u_f = 0.5_DP * (u(c-1,j) + u(c,j))
        rho_f = 0.5_DP * (rho(c-1,j) + rho(c,j))
        nt_up = merge(left, right, u_f >= 0.0_DP)
        adv_x_w = rho_f * u_f * nt_up
        nt_f = 0.5_DP * (left + right)
        nu_f = mu_lam / max(rho_f, 1.0e-12_DP)
        dif_x_w = (nu_f + nt_f) / SIGMA_SA * (right - left) / dx

        ! east x face (between cols c | c+1)
        left = pnt(c,jp); right = pnt(c+1,jp)
        u_f = 0.5_DP * (u(c,j) + u(c+1,j))
        rho_f = 0.5_DP * (rho(c,j) + rho(c+1,j))
        nt_up = merge(left, right, u_f >= 0.0_DP)
        adv_x_e = rho_f * u_f * nt_up
        nt_f = 0.5_DP * (left + right)
        nu_f = mu_lam / max(rho_f, 1.0e-12_DP)
        dif_x_e = (nu_f + nt_f) / SIGMA_SA * (right - left) / dx

        ! south y face (between rows j-1 | j); the domain-edge face is
        ! single-sided on the mean flow (rust f==0 branch), so the
        ! padded ghost rows never enter it.
        if (j == 1) then
          v_f = v(c,1); rho_f = rho(c,1)
        else
          v_f = 0.5_DP * (vp(c,jp-1) + vp(c,jp))
          rho_f = 0.5_DP * (rho(c,j-1) + rho(c,j))
        end if
        below = pnt(c,jp-1); above = pnt(c,jp)
        nt_up = merge(below, above, v_f >= 0.0_DP)
        adv_y_s = rho_f * v_f * nt_up
        nt_f = 0.5_DP * (below + above)
        nu_f = mu_lam / max(rho_f, 1.0e-12_DP)
        dif_y_s = (nu_f + nt_f) / SIGMA_SA * (above - below) / dy

        ! north y face (between rows j | j+1)
        if (j == ny_) then
          v_f = v(c,ny_); rho_f = rho(c,ny_)
        else
          v_f = 0.5_DP * (vp(c,jp) + vp(c,jp+1))
          rho_f = 0.5_DP * (rho(c,j) + rho(c,j+1))
        end if
        below = pnt(c,jp); above = pnt(c,jp+1)
        nt_up = merge(below, above, v_f >= 0.0_DP)
        adv_y_n = rho_f * v_f * nt_up
        nt_f = 0.5_DP * (below + above)
        nu_f = mu_lam / max(rho_f, 1.0e-12_DP)
        dif_y_n = (nu_f + nt_f) / SIGMA_SA * (above - below) / dy

        ! conservative flux differences: advection out minus in,
        ! diffusion in minus out (the laplacian sign is plus).
        rhs = -((adv_x_e - adv_x_w) / dx + (adv_y_n - adv_y_s) / dy) &
              / max(rho(c,j), 1.0e-12_DP) &
              + (dif_x_e - dif_x_w) / dx + (dif_y_n - dif_y_s) / dy

        ! c_b2/sigma cross term on centered grads of the padded field
        gx = (pnt(c+1,jp) - pnt(c-1,jp)) / (2.0_DP * dx)
        gy = (pnt(c,jp+1) - pnt(c,jp-1)) / (2.0_DP * dy)
        cb2_term = CB2_SA / SIGMA_SA * (gx * gx + gy * gy)

        ! vorticity s = |dv/dx - du/dy| from centered velocity grads
        s_vort = abs( (vp(c+1,jp) - vp(c-1,jp)) / (2.0_DP * dx) &
                    - (up(c,jp+1) - up(c,jp-1)) / (2.0_DP * dy) )

        call sa_source_closure(nt_c, s_vort, dd(c-2,j), nu_c, &
                               s_tilde, rr, prod, destr, src)
        rhs = rhs + cb2_term + src
        nt(c,j) = nt_c + dt * rhs
      end do
    end do

    ! positivity: the transported field is physically non-negative
    do j = 1, ny_
      do i = 3, nx_l+2
        nt(i,j) = max(nt(i,j), 0.0_DP)
      end do
    end do
    deallocate(pnt, up, vp)
  end subroutine

end module nufoorturb

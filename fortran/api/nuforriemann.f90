! nuforriemann.f90 — the hllc approximate riemann solver for one face,
! taking left/right primitive-plus-closure states and returning the flux.
module nuforriemann
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: DP = real64

contains


  ! hllc flux across one face, x-normal: a faithful port of rust
  ! hllc2d::hllc_flux (ql/qr pressure sharpening, toro star states via
  ! rankine-hugoniot, contact-wave selection). primitives in, flux out.
  subroutine hllc_face_x(rl, ul, vl, pl, el_t, al, rr, ur, vr, pr, er_t, ar, gam, &
                         fmass, fmx, fmy, fe)
    real(DP), intent(in) :: rl, ul, vl, pl, el_t, al, rr, ur, vr, pr, er_t, ar, gam
    real(DP), intent(out) :: fmass, fmx, fmy, fe
    real(DP) :: rl_, pl_, al_, rr_, pr_, ar_
    real(DP) :: gl, gr, rbar, abar, pstar, ql, qr, sl, sr, denom, sstar
    real(DP) :: etl, etr, elt, ert
    real(DP) :: fl(4), fr(4), usl(4), usr(4), fsl(4), fsr(4), f(4)
    rl_ = max(rl, 1.0e-12_DP); pl_ = max(pl, 1.0e-12_DP); al_ = max(al, 1.0e-12_DP)
    rr_ = max(rr, 1.0e-12_DP); pr_ = max(pr, 1.0e-12_DP); ar_ = max(ar, 1.0e-12_DP)
    if (pl_ > 1.0e-9_DP .and. rl_ > 1.0e-9_DP) then
      gl = al_ * al_ * rl_ / pl_
    else
      gl = gam
    end if
    if (pr_ > 1.0e-9_DP .and. rr_ > 1.0e-9_DP) then
      gr = ar_ * ar_ * rr_ / pr_
    else
      gr = gam
    end if
    rbar = 0.5_DP * (rl_ + rr_)
    abar = 0.5_DP * (al_ + ar_)
    pstar = max(0.5_DP * (pl_ + pr_) - 0.5_DP * (ur - ul) * rbar * abar, 0.0_DP)
    if (pstar <= pl_) then
      ql = 1.0_DP
    else
      ql = sqrt(1.0_DP + (gl + 1.0_DP) / (2.0_DP * gl) * (pstar / pl_ - 1.0_DP))
    end if
    if (pstar <= pr_) then
      qr = 1.0_DP
    else
      qr = sqrt(1.0_DP + (gr + 1.0_DP) / (2.0_DP * gr) * (pstar / pr_ - 1.0_DP))
    end if
    sl = min(ul - al_ * ql, ur - ar_ * qr)
    sr = max(ul + al_ * ql, ur + ar_ * qr)
    denom = rl_ * (sl - ul) - rr_ * (sr - ur)
    sstar = (pr_ - pl_ + rl_ * ul * (sl - ul) - rr_ * ur * (sr - ur)) / denom
    ! total energies (per volume) and physical fluxes
    elt = rl_ * (el_t + 0.5_DP * (ul**2 + vl**2))
    ert = rr_ * (er_t + 0.5_DP * (ur**2 + vr**2))
    fl = [rl_ * ul, rl_ * ul**2 + pl_, rl_ * ul * vl, ul * (elt + pl_)]
    fr = [rr_ * ur, rr_ * ur**2 + pr_, rr_ * ur * vr, ur * (ert + pr_)]
    ! star states and star fluxes: f* = f + s (u* - u)
    usl = star_state(rl_, ul, vl, pl_, elt, sl, sstar)
    usr = star_state(rr_, ur, vr, pr_, ert, sr, sstar)
    fsl = fl + sl * (usl - [rl_, rl_ * ul, rl_ * vl, elt])
    fsr = fr + sr * (usr - [rr_, rr_ * ur, rr_ * vr, ert])
    if (sl > 0.0_DP) then
      f = fl
    else if (sr < 0.0_DP) then
      f = fr
    else if (sstar > 0.0_DP) then
      f = fsl
    else
      f = fsr
    end if
    fmass = f(1); fmx = f(2); fmy = f(3); fe = f(4)
  end subroutine


  ! toro star conserved state for one side
  pure function star_state(r, un, ut, p, e_tot, s_side, sstar) result(u)
    real(DP), intent(in) :: r, un, ut, p, e_tot, s_side, sstar
    real(DP) :: u(4), fac
    fac = r * (s_side - un) / (s_side - sstar)
    u = [fac, fac * sstar, fac * ut, &
         fac * (e_tot / r + (sstar - un) * (sstar + p / (r * (s_side - un))))]
  end function

end module nuforriemann

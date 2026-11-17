! nuforcomb.f90 -- the li2004 9-species mixture closure and per-cell
! chemistry integrator for the combustion lane.
!
! ports crates/nufor-core/src/combustion_closure.rs (temperature_at,
! pressure_at, sound_speed_at, internal_energy_at, mix_gas_constant) and
! the li arm of crates/nufor-core/src/chem_kernel.rs +
! crates/nufor-core/src/combustion.rs omega/k_f/reaction_rate. same fit
! ordering, same operation order, so a fortran cell agrees with a rust
! cell to the last rounding.
module nuforcomb
  use iso_fortran_env, only: real64
  use nuforcombconst, only: RP, NS, NR, RU, P0, NASA_LOW, NASA_HIGH, SMASS, RATE_KF, RATE_N, RATE_E, NUT, R14A, R14N, R14E, R19A, R19N, R19E, FAL9_A0, FAL9_N0, FAL9_E0, FAL9_AINF, FAL9_NINF, FAL9_EINF, FAL9_FC, FAL15_A0, FAL15_N0, FAL15_E0, FAL15_AINF, FAL15_NINF, FAL15_EINF, FAL15_FC
  implicit none
  real(RP), parameter :: T_SWITCH = 1000.0_RP, T_FLOOR = 200.0_RP

contains

  ! cp/R for species s (1-based) at t, floored at 200 K as in rust.
  pure function cp_r(s, t)
    integer, intent(in) :: s
    real(RP), intent(in) :: t
    real(RP) :: cp_r
    real(RP) :: tt, a(7)
    tt = max(t, T_FLOOR)
    if (tt < T_SWITCH) then
      a = NASA_LOW(:, s)
    else
      a = NASA_HIGH(:, s)
    end if
    cp_r = a(1) + a(2)*tt + a(3)*tt**2 + a(4)*tt**3 + a(5)*tt**4
  end function cp_r

  ! h/(R T) for species s at t, formation included.
  pure function h_rt(s, t)
    integer, intent(in) :: s
    real(RP), intent(in) :: t
    real(RP) :: h_rt
    real(RP) :: tt, a(7)
    tt = max(t, T_FLOOR)
    if (tt < T_SWITCH) then
      a = NASA_LOW(:, s)
    else
      a = NASA_HIGH(:, s)
    end if
    h_rt = a(1) + a(2)*tt/2.0_RP + a(3)*tt**2/3.0_RP + a(4)*tt**3/4.0_RP &
           + a(5)*tt**4/5.0_RP + a(6)/tt
  end function h_rt

  ! s/R for species s at t, 1 atm standard state.
  pure function s_r(s, t)
    integer, intent(in) :: s
    real(RP), intent(in) :: t
    real(RP) :: s_r
    real(RP) :: tt, a(7)
    tt = max(t, T_FLOOR)
    if (tt < T_SWITCH) then
      a = NASA_LOW(:, s)
    else
      a = NASA_HIGH(:, s)
    end if
    s_r = a(1)*log(tt) + a(2)*tt + a(3)*tt**2/2.0_RP + a(4)*tt**3/3.0_RP &
          + a(5)*tt**4/4.0_RP + a(7)
  end function s_r

  pure function g_rt(s, t)
    integer, intent(in) :: s
    real(RP), intent(in) :: t
    real(RP) :: g_rt
    g_rt = h_rt(s, t) - s_r(s, t)
  end function g_rt

  ! static enthalpy per mass for species s at t, formation included.
  pure function h_mass(s, t)
    integer, intent(in) :: s
    real(RP), intent(in) :: t
    real(RP) :: h_mass
    h_mass = RU * t * h_rt(s, t) / SMASS(s)
  end function h_mass

  ! cp per mass, j/(kg k).
  pure function cp_mass(s, t)
    integer, intent(in) :: s
    real(RP), intent(in) :: t
    real(RP) :: cp_mass
    cp_mass = RU * cp_r(s, t) / SMASS(s)
  end function cp_mass

  pure function mix_gas_constant(y)
    real(RP), intent(in) :: y(9)
    real(RP) :: mix_gas_constant
    integer :: s
    mix_gas_constant = 0.0_RP
    do s = 1, 9
      mix_gas_constant = mix_gas_constant + y(s)*RU/SMASS(s)
    end do
  end function mix_gas_constant

  pure function mix_cp_mass(y, t)
    real(RP), intent(in) :: y(9), t
    real(RP) :: mix_cp_mass
    integer :: s
    mix_cp_mass = 0.0_RP
    do s = 1, 9
      mix_cp_mass = mix_cp_mass + y(s)*cp_mass(s, t)
    end do
  end function mix_cp_mass

  pure function mix_h_mass(y, t)
    real(RP), intent(in) :: y(9), t
    real(RP) :: mix_h_mass
    integer :: s
    mix_h_mass = 0.0_RP
    do s = 1, 9
      mix_h_mass = mix_h_mass + y(s)*h_mass(s, t)
    end do
  end function mix_h_mass

  ! internal energy per mass at (t, y); the clamp mirrors rust.
  pure function internal_energy_at(t, y)
    real(RP), intent(in) :: t, y(9)
    real(RP) :: internal_energy_at
    real(RP) :: tt
    tt = max(t, T_FLOOR)
    internal_energy_at = mix_h_mass(y, tt) - mix_gas_constant(y)*tt
  end function internal_energy_at

  ! temperature from (e_int, y): newton on the bracket [200,6000] with
  ! the analytic cv derivative, bisection fallback. mirrors rust exactly.
  function temperature_at(e_int, y) result(t)
    real(RP), intent(in) :: e_int, y(9)
    real(RP) :: t
    real(RP), parameter :: LO = 200.0_RP, HI = 6000.0_RP
    real(RP) :: f_lo, f_hi, tt, fv, dv, nxt, aa, bb, mid
    integer :: it
    f_lo = internal_energy_at(LO, y) - e_int
    f_hi = internal_energy_at(HI, y) - e_int
    if (f_lo > 0.0_RP) then
      t = LO
      return
    end if
    if (f_hi < 0.0_RP) then
      t = HI
      return
    end if
    tt = 0.5_RP*(LO + HI)
    do it = 1, 12
      fv = internal_energy_at(tt, y) - e_int
      dv = mix_cp_mass(y, tt) - mix_gas_constant(y)
      if (dv <= 1.0_RP) exit
      nxt = tt - fv/dv
      if (nxt < LO) nxt = LO
      if (nxt > HI) nxt = HI
      if (abs(nxt - tt) < 1e-6_RP) then
        t = nxt
        return
      end if
      tt = nxt
    end do
    if (abs(internal_energy_at(tt, y) - e_int) < 1e-3_RP) then
      t = tt
      return
    end if
    aa = LO; bb = HI
    do it = 1, 24
      mid = 0.5_RP*(aa + bb)
      if (internal_energy_at(mid, y) - e_int < 0.0_RP) then
        aa = mid
      else
        bb = mid
      end if
    end do
    t = 0.5_RP*(aa + bb)
  end function temperature_at

  function pressure_at(rho, e_int, y)
    real(RP), intent(in) :: rho, e_int, y(9)
    real(RP) :: pressure_at
    pressure_at = rho*mix_gas_constant(y)*temperature_at(e_int, y)
  end function pressure_at

  ! frozen-composition sound speed, floored like the rust closure.
  function sound_speed_at(rho, e_int, y)
    real(RP), intent(in) :: rho, e_int, y(9)
    real(RP) :: sound_speed_at
    real(RP) :: cp, cv, gam, pp
    cp = mix_cp_mass(y, temperature_at(e_int, y))
    cv = max(cp - mix_gas_constant(y), 1.0_RP)
    gam = max(cp/cv, 1.01_RP)
    pp = max(pressure_at(rho, e_int, y), 1e-9_RP)
    sound_speed_at = sqrt(gam*pp/rho)
  end function sound_speed_at

  ! one closure tuple: pressure, temperature, sound speed for a cell.
  subroutine comb_triple(rho, e_int, y, pp, tt, aa)
    real(RP), intent(in) :: rho, e_int, y(9)
    real(RP), intent(out) :: pp, tt, aa
    tt = temperature_at(e_int, y)
    pp = rho*mix_gas_constant(y)*tt
    aa = sound_speed_at(rho, e_int, y)
  end subroutine comb_triple

  ! arrhenius k = a t^n exp(-e/(ru t)); negative e is legal.
  pure function arrh(a, n, e, t)
    real(RP), intent(in) :: a, n, e, t
    real(RP) :: arrh
    arrh = a*t**n*exp(-e/(RU*t))
  end function arrh

  ! forward rate constant for non-falloff reaction r (0-based).
  pure function k_forward(r, t)
    integer, intent(in) :: r
    real(RP), intent(in) :: t
    real(RP) :: k_forward
    integer :: rf
    rf = r + 1
    if (rf == 9 .or. rf == 15) then
      k_forward = 0.0_RP
      return
    end if
    k_forward = arrh(RATE_KF(rf), RATE_N(rf), RATE_E(rf), t)
    if (rf == 14) k_forward = k_forward + arrh(R14A, R14N, R14E, t)
    if (rf == 19) k_forward = k_forward + arrh(R19A, R19N, R19E, t)
  end function k_forward

  ! ln Kc for reaction r (0-based) from the nasa fits, p0 standard.
  pure function ln_kc(r, t)
    integer, intent(in) :: r
    real(RP), intent(in) :: t
    real(RP) :: ln_kc
    real(RP) :: dnu, lp
    integer :: s, rf
    rf = r + 1
    dnu = 0.0_RP
    lp = 0.0_RP
    do s = 1, 9
      lp = lp - NUT(rf, s)*g_rt(s, t)
      dnu = dnu + NUT(rf, s)
    end do
    ln_kc = lp + dnu*log(P0/(RU*t))
  end function ln_kc

  ! troe falloff constant for rx (0-based 8 -> pair 9, 14 -> pair 15).
  pure function k_falloff(r, t, cm)
    integer, intent(in) :: r
    real(RP), intent(in) :: t, cm
    real(RP) :: k_falloff
    real(RP) :: k0, kinf, pr, a0, n0, e0, ai, ni, ei, fc
    if (r == 8) then
      a0 = FAL9_A0; n0 = FAL9_N0; e0 = FAL9_E0
      ai = FAL9_AINF; ni = FAL9_NINF; ei = FAL9_EINF; fc = FAL9_FC
    else
      a0 = FAL15_A0; n0 = FAL15_N0; e0 = FAL15_E0
      ai = FAL15_AINF; ni = FAL15_NINF; ei = FAL15_EINF; fc = FAL15_FC
    end if
    k0 = arrh(a0, n0, e0, t)
    kinf = arrh(ai, ni, ei, t)
    if (kinf == 0.0_RP) then
      k_falloff = 0.0_RP
      return
    end if
    pr = k0*cm/kinf
    k_falloff = kinf*(pr/(1.0_RP + pr))*fc
  end function k_falloff

  ! third-body efficiency, li footnotes a/b/c/f (r,s 1-based).
  pure function tb_eff(r, s)
    integer, intent(in) :: r, s
    real(RP) :: tb_eff
    if ((r == 5 .and. s == 3) .or. (r == 6 .and. s == 3) .or. (r == 7 .and. s == 3)) then
      tb_eff = 12.0_RP
    else if ((r == 5 .and. s == 1) .or. (r == 6 .and. s == 1) .or. (r == 7 .and. s == 1)) then
      tb_eff = 2.5_RP
    else if (r == 8 .and. s == 3) then
      tb_eff = 12.0_RP
    else if (r == 9 .and. s == 3) then
      tb_eff = 11.0_RP
    else if (r == 9 .and. s == 1) then
      tb_eff = 2.0_RP
    else if (r == 9 .and. s == 2) then
      tb_eff = 0.78_RP
    else if (r == 15 .and. s == 3) then
      tb_eff = 12.0_RP
    else if (r == 15 .and. s == 1) then
      tb_eff = 2.5_RP
    else
      tb_eff = 1.0_RP
    end if
  end function tb_eff

  pure function conc_m(r, nd) ! r 1-based, nd mol/m3
    integer, intent(in) :: r
    real(RP), intent(in) :: nd(9)
    real(RP) :: conc_m
    integer :: s
    conc_m = 0.0_RP
    do s = 1, 9
      conc_m = conc_m + tb_eff(r, s)*nd(s)
    end do
  end function conc_m

  ! net rate of reaction r (0-based), mol/m3/s, forward minus reverse.
  function reaction_rate(r, nd, t) result(q)
    integer, intent(in) :: r
    real(RP), intent(in) :: nd(9), t
    real(RP) :: q
    real(RP) :: kc, kf, fwd, rev, cm
    integer :: s, rf
    rf = r + 1
    kc = exp(ln_kc(r, t))
    if (rf == 9 .or. rf == 15) then
      kf = k_falloff(r, t, conc_m(rf, nd))
    else
      kf = k_forward(r, t)
    end if
    if (kf == 0.0_RP) then
      q = 0.0_RP
      return
    end if
    fwd = kf
    rev = kf/kc
    do s = 1, 9
      if (NUT(rf, s) < 0.0_RP) then
        fwd = fwd*nd(s)**(-NUT(rf, s))
      else if (NUT(rf, s) > 0.0_RP) then
        rev = rev*nd(s)**NUT(rf, s)
      end if
    end do
    if (rf >= 5 .and. rf <= 8) then
      cm = conc_m(rf, nd)
      fwd = fwd*cm
      rev = rev*cm
    end if
    q = fwd - rev
  end function reaction_rate

  ! per-species mass production dy/dt (1/s) for one cell.
  subroutine omega_li(rho, y, t, dydt)
    real(RP), intent(in) :: rho, y(9), t
    real(RP), intent(out) :: dydt(9)
    real(RP) :: nd(9), q
    integer :: r, s
    nd = rho*y/SMASS
    dydt = 0.0_RP
    do r = 0, NR - 1
      q = reaction_rate(r, nd, t)
      if (q == 0.0_RP) cycle
      do s = 1, 9
        if (NUT(r+1, s) /= 0.0_RP) then
          dydt(s) = dydt(s) + NUT(r+1, s)*q*SMASS(s)/rho
        end if
      end do
    end do
  end subroutine omega_li

  ! the li arm of the per-cell chemistry step: adaptive explicit stepping
  ! at fixed total energy with rate scaling for positivity, mirroring
  ! chem_kernel::li_arm. returns the worst oversub measure.
  subroutine chem_cell_li(e_int, rho, y, dt, oversub)
    real(RP), intent(in) :: e_int, rho, dt
    real(RP), intent(inout) :: y(9)
    real(RP), intent(out) :: oversub
    real(RP) :: a(9), dydt(9), t, remaining, rate_scale, tau, hh, scale, allowed, dd
    integer :: s, steps
    a = y
    oversub = 0.0_RP
    t = temperature_at(e_int, a)
    remaining = dt
    steps = 0
    do while (remaining > 1e-14_RP .and. steps < 512)
      call omega_li(rho, a, t, dydt)
      rate_scale = 0.0_RP
      do s = 1, 9
        rate_scale = max(rate_scale, abs(dydt(s)))
      end do
      if (rate_scale > 1e-30_RP) then
        tau = 1.0_RP/rate_scale
      else
        tau = remaining
      end if
      hh = min(remaining, tau*0.2_RP)
      scale = 1.0_RP
      do s = 1, 9
        if (dydt(s) < 0.0_RP .and. a(s) > 0.0_RP) then
          allowed = a(s)/(-dydt(s)*hh)
          if (allowed < scale) scale = allowed
        end if
      end do
      do s = 1, 9
        dd = dydt(s)*hh*scale
        a(s) = a(s) + dd
      end do
      t = temperature_at(e_int, a)
      remaining = remaining - hh
      steps = steps + 1
      if (steps == 1) oversub = max(oversub, dt/(tau*0.2_RP))
    end do
    y = a
  end subroutine chem_cell_li

  ! chemistry over a local slab of cells; y is (species, cell).
  subroutine chem_substep_grid_li(nn, rho, e_int, dt, y, worst_ov)
    integer, intent(in) :: nn
    real(RP), intent(in) :: rho(nn), e_int(nn), dt
    real(RP), intent(inout) :: y(9, nn)
    real(RP), intent(out) :: worst_ov
    real(RP) :: ov
    integer :: k
    worst_ov = 0.0_RP
    do k = 1, nn
      call chem_cell_li(e_int(k), rho(k), y(:, k), dt, ov)
      if (ov > worst_ov) worst_ov = ov
    end do
  end subroutine chem_substep_grid_li

end module nuforcomb

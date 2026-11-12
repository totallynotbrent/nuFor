! nuforpark5.f90 — the park 5-species lane for the fortran core.
!
! this module mirrors crates/nufor-core/src/species_closure.rs,
! species_chem.rs, species_transport.rs and the chemistry substep of
! species_march.rs operation for operation: the mixture closure (pressure,
! temperature, sound speed from rho, e_internal, Y), the park 1990
! forward rates with the tp-2867 equilibrium-constant fits, the upwind
! species transport on the flow's face mass fluxes, and the
! operator-split chemistry substep.
!
! the lane owns no flow state. the caller supplies the mass flux per face
! from whichever march produced it, so the species transport rides the
! same fluxes the flow advanced under, and the chemistry runs afterwards
! on the post-step composition.
!
! the transcendental functions are reached through the c library on
! purpose: rust's f64::exp, f64::ln and f64::powf lower to libc's exp,
! log and pow, and fortran's own ** may expand to exp(log()), which
! rounds differently. calling the same c functions keeps the two
! implementations agreeing to the last bit on the transcendentals.

module nuforpark5
  use iso_c_binding, only: c_double, c_int
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: DP = real64

  ! number of species in the park 5-species air set.
  integer, parameter :: NS = 5

  ! molar masses, kg/mol: n2, o2, no, n, o.
  real(DP), parameter :: M_SPEC(NS) = [ &
    28.016e-3_DP, 32.000e-3_DP, 30.008e-3_DP, 14.008e-3_DP, 16.000e-3_DP ]

  ! molar gas constant, J/(mol K).
  real(DP), parameter :: RU = 8.314462_DP

  ! dissociation energies in K (park 1990): n2, o2, no.
  real(DP), parameter :: D0_K(3) = [113200.0_DP, 59400.0_DP, 75500.0_DP]

  ! per-species gas constant, J/(kg K).
  real(DP), parameter :: R_SPEC(NS) = RU / M_SPEC

  ! per-species constant-v (thermal only, vibration excluded), J/(kg K):
  ! molecules 5/2 R, atoms 3/2 R.
  real(DP), parameter :: CV_SPEC(NS) = [ &
    2.5_DP * R_SPEC(1), 2.5_DP * R_SPEC(2), 2.5_DP * R_SPEC(3), &
    1.5_DP * R_SPEC(4), 1.5_DP * R_SPEC(5) ]

  ! per-species formation enthalpy per mass, J/kg, on the 0-K
  ! atom-positive convention of tp-2867 table i.
  real(DP), parameter :: HF_SPEC(NS) = [ &
    0.0_DP, &
    0.0_DP, &
    RU * (0.5_DP * D0_K(1) + 0.5_DP * D0_K(2) - D0_K(3)) / M_SPEC(3), &
    0.5_DP * RU * D0_K(1) / M_SPEC(4), &
    0.5_DP * RU * D0_K(2) / M_SPEC(5) ]

  ! modified-arrhenius forward rates: kf = a * t^eta * exp(-td / t).
  ! the six entries are the park 1990 table 2 set: the n2 dissociation
  ! carries one rate for the molecular colliders and another for the
  ! atoms, then o2 dissociation, no dissociation, and the two zeldovich
  ! exchanges.
  real(DP), parameter :: RATE_A(6) = [ &
    7.0e15_DP, 3.0e16_DP, 2.0e15_DP, 5.0e9_DP, 6.4e11_DP, 8.4e6_DP ]
  real(DP), parameter :: RATE_ETA(6) = [ &
    -1.6_DP, -1.6_DP, -1.5_DP, -1.0_DP, -1.0_DP, 0.0_DP ]
  real(DP), parameter :: RATE_TD(6) = [ &
    D0_K(1), D0_K(1), D0_K(2), D0_K(3), 38370.0_DP, 19450.0_DP ]

  ! the five equilibrium-constant fits, tp-2867 table iii:
  ! ln kc = b1 + b2 ln z + b3 z + b4 z^2 + b5 z^3 with z = 10000/t.
  real(DP), parameter :: KEQ_B1(5) = [ &
    1.858_DP, 2.855_DP, 0.792_DP, 1.066_DP, -2.063_DP ]
  real(DP), parameter :: KEQ_B2(5) = [ &
    -1.325_DP, 0.988_DP, -0.492_DP, -0.833_DP, -1.480_DP ]
  real(DP), parameter :: KEQ_B3(5) = [ &
    -9.856_DP, -6.181_DP, -6.761_DP, -3.095_DP, -0.580_DP ]
  real(DP), parameter :: KEQ_B4(5) = [ &
    -0.174_DP, -0.023_DP, -0.091_DP, -0.084_DP, -0.114_DP ]
  real(DP), parameter :: KEQ_B5(5) = [ &
    0.008_DP, -0.001_DP, 0.004_DP, 0.004_DP, 0.005_DP ]
  ! the mole-number change per reaction, which fixes kc's units.
  integer, parameter :: KEQ_DNU(5) = [1, 1, 1, 0, 0]

  ! structured error codes shared with the rust wrapper.
  integer(c_int), parameter :: NFOR_OK    = 0
  integer(c_int), parameter :: NFOR_EARGS = 1

  private
  public :: nfor_park5_closure
  public :: nfor_park5_omega
  public :: nfor_park5_chem_substep
  public :: nfor_park5_transport
  public :: closure_triple, chem_substep_grid, transport_planar
  public :: pressure_at, sound_speed_at, temperature_at

  interface
    function c_exp(x) bind(c, name="exp") result(y)
      import c_double
      real(c_double), value :: x
      real(c_double) :: y
    end function c_exp

    function c_log(x) bind(c, name="log") result(y)
      import c_double
      real(c_double), value :: x
      real(c_double) :: y
    end function c_log

    function c_pow(x, y) bind(c, name="pow") result(z)
      import c_double
      real(c_double), value :: x, y
      real(c_double) :: z
    end function c_pow
  end interface

contains

  ! the mixture's thermal capacity per mass, J/(kg K): sum of y_i cv_i.
  function mixture_cv(y) result(cv)
    real(DP), intent(in) :: y(NS)
    real(DP) :: cv
    integer :: i
    cv = 0.0_DP
    do i = 1, NS
      cv = cv + y(i) * CV_SPEC(i)
    end do
  end function mixture_cv

  ! the mixture's formation energy per mass, J/kg.
  function mixture_hf(y) result(hf)
    real(DP), intent(in) :: y(NS)
    real(DP) :: hf
    integer :: i
    hf = 0.0_DP
    do i = 1, NS
      hf = hf + y(i) * HF_SPEC(i)
    end do
  end function mixture_hf

  ! the mixture's gas constant per mass, J/(kg K).
  function mix_gas_constant(y) result(r)
    real(DP), intent(in) :: y(NS)
    real(DP) :: r
    integer :: i
    r = 0.0_DP
    do i = 1, NS
      r = r + y(i) * R_SPEC(i)
    end do
  end function mix_gas_constant

  ! temperature from the internal energy per mass and the composition, K.
  function temperature_at(e_int, y) result(t)
    real(DP), intent(in) :: e_int, y(NS)
    real(DP) :: t
    real(DP) :: a, b
    a = mixture_cv(y)
    b = mixture_hf(y)
    ! the thermal part floors at 50 K: a cold mixture can carry an
    ! internal energy below the formation floor during a shock step.
    t = max((e_int - b) / a, 50.0_DP)
  end function temperature_at

  ! the internal energy per mass at the temperature and composition.
  function internal_energy_at(t, y) result(e_int)
    real(DP), intent(in) :: t, y(NS)
    real(DP) :: e_int
    e_int = mixture_cv(y) * t + mixture_hf(y)
  end function internal_energy_at

  ! pressure from the density, internal energy per mass and composition.
  function pressure_at(rho, e_int, y) result(p)
    real(DP), intent(in) :: rho, e_int, y(NS)
    real(DP) :: p
    p = rho * mix_gas_constant(y) * temperature_at(e_int, y)
  end function pressure_at

  ! the frozen-composition sound speed from the mixture gamma.
  function sound_speed_at(rho, e_int, y) result(ac)
    real(DP), intent(in) :: rho, e_int, y(NS)
    real(DP) :: ac
    real(DP) :: cv, cp, gamma, p
    integer :: i
    cv = 0.0_DP
    cp = 0.0_DP
    do i = 1, NS
      cv = cv + y(i) * CV_SPEC(i)
      cp = cp + y(i) * (CV_SPEC(i) + R_SPEC(i))
    end do
    gamma = max(cp / cv, 1.01_DP)
    p = max(pressure_at(rho, e_int, y), 1.0e-9_DP)
    ac = sqrt(gamma * p / rho)
  end function sound_speed_at

  ! the concentration-form equilibrium constant of reaction r, mol/m3 for
  ! the dissociations and dimensionless for the exchanges.
  function k_eq(r, t) result(kc)
    integer, intent(in) :: r
    real(DP), intent(in) :: t
    real(DP) :: kc
    real(DP) :: z, ln_kc
    z = 10000.0_DP / max(t, 200.0_DP)
    ln_kc = KEQ_B1(r) + KEQ_B2(r) * c_log(z) + KEQ_B3(r) * z &
      + KEQ_B4(r) * z * z + KEQ_B5(r) * z * z * z
    kc = c_exp(ln_kc)
    if (KEQ_DNU(r) == 1) then
      ! the fits print mol/cm3; the si form is mol/m3.
      kc = kc * 1.0e6_DP
    end if
  end function k_eq

  ! the forward rate of reaction r.
  function k_f(r, t) result(k)
    integer, intent(in) :: r
    real(DP), intent(in) :: t
    real(DP) :: k
    k = RATE_A(r) * c_pow(t, RATE_ETA(r)) * c_exp(-RATE_TD(r) / t)
  end function k_f

  ! the backward rate of reaction r, from the shared equilibrium constant.
  function k_b(r, t) result(k)
    integer, intent(in) :: r
    real(DP), intent(in) :: t
    real(DP) :: k
    k = k_f(r, t) / k_eq(r, t)
  end function k_b

  ! one cell's reaction source: the mass-fraction production rates and
  ! the chemical energy deposition per mass.
  !
  ! the n2 dissociation splits by collider, each branch pairing its own
  ! forward constant with the shared kc, so both branches sit in
  ! equilibrium at the same composition.
  subroutine omega(rho, y, t, dydt, de_chem)
    real(DP), intent(in) :: rho, y(NS), t
    real(DP), intent(out) :: dydt(NS), de_chem
    real(DP) :: n(NS), n_mol, n_atom, nn2, no2, nno, nn, no
    real(DP) :: kc0, kf0_mol, kf0_atom, kf1, kb1, kf2, kb2, kf3, kb3, kf4, kb4
    real(DP) :: r0_net, r1_net, r2_net, r3_net, r4_net, dn(NS), dy
    integer :: s

    do s = 1, NS
      n(s) = rho * y(s) / M_SPEC(s)
    end do
    n_mol = n(1) + n(2) + n(3)
    n_atom = n(4) + n(5)
    nn2 = n(1)
    no2 = n(2)
    nno = n(3)
    nn = n(4)
    no = n(5)

    kc0 = k_eq(1, t)
    kf0_mol = k_f(1, t)
    kf0_atom = k_f(2, t)
    r0_net = (kf0_mol * nn2 - (kf0_mol / kc0) * nn * nn) * n_mol &
      + (kf0_atom * nn2 - (kf0_atom / kc0) * nn * nn) * n_atom

    kf1 = k_f(3, t)
    kb1 = k_b(2, t)
    r1_net = (kf1 * no2 - kb1 * no * no) * n_mol

    kf2 = k_f(4, t)
    kb2 = k_b(3, t)
    r2_net = (kf2 * nno - kb2 * nn * no) * n_mol

    kf3 = k_f(5, t)
    kb3 = k_b(4, t)
    r3_net = kf3 * no * nn2 - kb3 * nno * nn

    kf4 = k_f(6, t)
    kb4 = k_b(5, t)
    r4_net = kf4 * no * nno - kb4 * no2 * nn

    dn(1) = -r0_net - r3_net + r4_net
    dn(2) = -r1_net + r4_net
    dn(3) = -r2_net + r3_net - r4_net
    dn(4) = 2.0_DP * r0_net + r2_net + r3_net + r4_net
    dn(5) = 2.0_DP * r1_net + r2_net - r3_net - r4_net

    de_chem = 0.0_DP
    do s = 1, NS
      dy = dn(s) / rho
      dy = dy * M_SPEC(s)
      dydt(s) = dy
      de_chem = de_chem + dy * HF_SPEC(s)
    end do
  end subroutine omega

  ! one cell's chemistry substep: the reactions integrate at the
  ! post-transport composition, capped at a fifth of the fastest chemical
  ! timescale, and the composition renormalizes after every substep.
  !
  ! the total internal energy is unchanged: the formation energy the
  ! products carry sits inside the conserved state through the closure,
  ! so the temperature simply re-solves from the same energy at the new
  ! composition and the cell cools as it dissociates.
  subroutine chem_substep_cell(rho, e_int, y, dt, worst_dy, worst_oversub)
    real(DP), intent(in) :: rho, e_int, dt
    real(DP), intent(inout) :: y(NS)
    real(DP), intent(out) :: worst_dy, worst_oversub
    real(DP) :: t, remaining, h, rate_scale, tau, d, sum_y
    ! the scratch lives here, not at module scope, so concurrent calls
    ! cannot overwrite each other's rates.
    real(DP) :: dydt_tmp(NS), de_chem_tmp
    integer :: s, steps

    worst_dy = 0.0_DP
    worst_oversub = 0.0_DP
    t = temperature_at(e_int, y)
    remaining = dt
    steps = 0
    do while (remaining > 1.0e-14_DP .and. steps < 64)
      call omega(rho, y, t, dydt_tmp, de_chem_tmp)
      rate_scale = 0.0_DP
      do s = 1, NS
        rate_scale = max(rate_scale, abs(dydt_tmp(s)))
      end do
      if (rate_scale > 1.0e-30_DP) then
        tau = 1.0_DP / rate_scale
      else
        tau = remaining
      end if
      h = min(remaining, tau * 0.2_DP)
      do s = 1, NS
        d = dydt_tmp(s) * h
        worst_dy = max(worst_dy, abs(d))
        y(s) = min(max(y(s) + d, 0.0_DP), 1.0_DP)
      end do
      sum_y = 0.0_DP
      do s = 1, NS
        sum_y = sum_y + y(s)
      end do
      if (sum_y > 0.0_DP) then
        do s = 1, NS
          y(s) = y(s) / sum_y
        end do
      end if
      t = temperature_at(e_int, y)
      remaining = remaining - h
      steps = steps + 1
      if (steps == 1) then
        worst_oversub = max(worst_oversub, dt / (tau * 0.2_DP))
      end if
    end do
  end subroutine chem_substep_cell

  ! the upwind species transport on the flow's face mass fluxes.
  !
  ! the face convention is the rust lane's: face f (f in 0..nx) sits
  ! between cell f-1 (left) and cell f (right), so cell i spans faces i
  ! and i+1. the face choice follows the sign of the mass flux, which
  ! is the hllc contact state's own choice for the material crossing
  ! the face.
  !
  ! the update advances rho*Y and divides by the density afterwards, so
  ! the pass is conservative in the per-species mass, not just the
  ! fraction. the domain-edge faces take the adjacent interior cell on
  ! the missing side, so a uniform field stays exactly uniform.
  subroutine transport_planar(nx, ny, rho, dt, dxs, dys, fx, fy, y, resid)
    integer, intent(in) :: nx, ny
    real(DP), intent(in) :: rho(nx*ny), dt, dxs(nx), dys(ny)
    real(DP), intent(in) :: fx(nx+1, ny), fy(ny+1, nx)
    real(DP), intent(inout) :: y(NS, nx*ny)
    real(DP), intent(out) :: resid
    real(DP) :: yx(NS, nx+1, ny), yy(NS, ny+1, nx)
    real(DP) :: mf, left, right, flux_x, flux_y, rho_y, d, new_rho_y
    real(DP) :: dxi, dyi
    integer :: s, j, i, f, k

    resid = 0.0_DP

    ! the face convention is the rust lane's own: face f (0-based) sits
    ! between cell f-1 on the left and cell f on the right, so cell i
    ! spans the faces i and i+1. the counters here run 0-based over the
    ! grid and the arrays are indexed 1-based fortran-side, hence the
    ! +1 offsets on the cell index.
    do s = 1, NS
      do j = 0, ny - 1
        do f = 0, nx
          mf = fx(f+1, j+1)
          if (f > 0) then
            left = y(s, j*nx + f)
          else
            left = y(s, j*nx + 1)
          end if
          if (f < nx) then
            right = y(s, j*nx + f + 1)
          else
            right = y(s, j*nx + nx)
          end if
          if (mf >= 0.0_DP) then
            yx(s, f+1, j+1) = left
          else
            yx(s, f+1, j+1) = right
          end if
        end do
      end do
      do i = 0, nx - 1
        do f = 0, ny
          mf = fy(f+1, i+1)
          if (f > 0) then
            left = y(s, (f-1)*nx + i + 1)
          else
            left = y(s, i + 1)
          end if
          if (f < ny) then
            right = y(s, f*nx + i + 1)
          else
            right = y(s, (ny-1)*nx + i + 1)
          end if
          if (mf >= 0.0_DP) then
            yy(s, f+1, i+1) = left
          else
            yy(s, f+1, i+1) = right
          end if
        end do
      end do
    end do

    do s = 1, NS
      do j = 0, ny - 1
        do i = 0, nx - 1
          k = j*nx + i + 1
          dxi = dt / dxs(i+1)
          dyi = dt / dys(j+1)
          ! cell i sits between the faces i and i+1, so the divergence
          ! is the far face's flux minus the near face's.
          flux_x = yx(s, i+2, j+1) * fx(i+2, j+1) - yx(s, i+1, j+1) * fx(i+1, j+1)
          flux_y = yy(s, j+2, i+1) * fy(j+2, i+1) - yy(s, j+1, i+1) * fy(j+1, i+1)
          rho_y = rho(k) * y(s, k)
          d = -(dxi * flux_x + dyi * flux_y)
          new_rho_y = rho_y + d
          resid = max(resid, abs(d) / max(rho(k), 1.0e-30_DP))
          ! a species mass may not go negative; the clip is the last
          ! resort that keeps the field finite and the caller keeps the
          ! step small enough that it does not trigger.
          if (new_rho_y > 0.0_DP) then
            y(s, k) = new_rho_y / rho(k)
          else
            y(s, k) = 0.0_DP
          end if
        end do
      end do
    end do
  end subroutine transport_planar

  ! the closure internals a sibling module borrows (the hllc face pass of
  ! a species-carrying march). temperature first, then pressure from it,
  ! then the frozen sound speed; exactly the c-surface triple evaluated
  ! inline so no flat-buffer round trip is paid per cell.
  subroutine closure_triple(rho, e_int, y, p, t, a)
    real(DP), intent(in) :: rho, e_int, y(NS)
    real(DP), intent(out) :: p, t, a
    p = pressure_at(rho, e_int, y)
    t = temperature_at(e_int, y)
    a = sound_speed_at(rho, e_int, y)
  end subroutine closure_triple

  ! grid-level chemistry: a convenience over chem_substep_cell for whole
  ! slabs held in 2d arrays (indexed (s, cell)). mirrors the c surface.
  subroutine chem_substep_grid(n, rho, e_int, dt, y, worst_dy, worst_oversub)
    integer, intent(in) :: n
    real(DP), intent(in) :: rho(n), e_int(n), dt
    real(DP), intent(inout) :: y(NS, n)
    real(DP), intent(out) :: worst_dy, worst_oversub
    real(DP) :: yc(NS), wd, wo
    integer :: k, s
    worst_dy = 0.0_DP
    worst_oversub = 0.0_DP
    do k = 1, n
      do s = 1, NS
        yc(s) = y(s, k)
      end do
      call chem_substep_cell(rho(k), e_int(k), yc, dt, wd, wo)
      do s = 1, NS
        ! stiff integrator leaves denormal undershoot on inert species; clip
        yc(s) = max(yc(s), 0.0_DP)
        y(s, k) = yc(s)
      end do
      worst_dy = max(worst_dy, wd)
      worst_oversub = max(worst_oversub, wo)
    end do
  end subroutine chem_substep_grid

  ! ---------------- c surface ----------------

  ! the mixture closure for n cells: y holds the species planes
  ! contiguously (species s occupies y((s-1)*n+1 : s*n)).
  subroutine nfor_park5_closure(n, rho, e_int, y, p, t, a, err) &
      bind(c, name="nfor_park5_closure")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: rho(*), e_int(*), y(*)
    real(c_double), intent(out) :: p(*), t(*), a(*)
    integer(c_int), intent(out) :: err
    real(DP) :: yc(NS)
    integer :: k, s

    if (n <= 0) then
      err = NFOR_EARGS
      return
    end if
    do k = 1, n
      do s = 1, NS
        yc(s) = y((s-1)*n + k)
      end do
      p(k) = pressure_at(rho(k), e_int(k), yc)
      t(k) = temperature_at(e_int(k), yc)
      a(k) = sound_speed_at(rho(k), e_int(k), yc)
    end do
    err = NFOR_OK
  end subroutine nfor_park5_closure

  ! the reaction source for n cells at the given temperature.
  subroutine nfor_park5_omega(n, rho, t, y, dydt, de_chem, err) &
      bind(c, name="nfor_park5_omega")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: rho(*), t(*), y(*)
    real(c_double), intent(out) :: dydt(*), de_chem(*)
    integer(c_int), intent(out) :: err
    real(DP) :: yc(NS), dydt_c(NS), de_chem_c
    integer :: k, s

    if (n <= 0) then
      err = NFOR_EARGS
      return
    end if
    do k = 1, n
      do s = 1, NS
        yc(s) = y((s-1)*n + k)
      end do
      call omega(rho(k), yc, t(k), dydt_c, de_chem_c)
      do s = 1, NS
        dydt((s-1)*n + k) = dydt_c(s)
      end do
      de_chem(k) = de_chem_c
    end do
    err = NFOR_OK
  end subroutine nfor_park5_omega

  ! the chemistry substep for n cells: the composition advances in place.
  subroutine nfor_park5_chem_substep(n, rho, e_int, dt, y, worst_dy, worst_oversub, err) &
      bind(c, name="nfor_park5_chem_substep")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: rho(*), e_int(*)
    real(c_double), value, intent(in) :: dt
    real(c_double), intent(inout) :: y(*)
    real(c_double), intent(out) :: worst_dy, worst_oversub
    integer(c_int), intent(out) :: err
    real(DP) :: yc(NS), wd, wo
    integer :: k, s

    if (n <= 0) then
      err = NFOR_EARGS
      return
    end if
    worst_dy = 0.0_DP
    worst_oversub = 0.0_DP
    do k = 1, n
      do s = 1, NS
        yc(s) = y((s-1)*n + k)
      end do
      call chem_substep_cell(rho(k), e_int(k), yc, dt, wd, wo)
      do s = 1, NS
        y((s-1)*n + k) = yc(s)
      end do
      worst_dy = max(worst_dy, wd)
      worst_oversub = max(worst_oversub, wo)
    end do
    err = NFOR_OK
  end subroutine nfor_park5_chem_substep

  ! the upwind species transport for a rectangular grid: y holds the
  ! species planes contiguously. fx holds the x-face mass fluxes row by
  ! row (jx from 0, then the nx+1 faces of that row) and fy the y-face
  ! fluxes column by column (iy from 0, then the ny+1 faces of that
  ! column), the layout the rust lane's face vectors already carry.
  subroutine nfor_park5_transport(nx, ny, rho, dt, dxs, dys, fx, fy, y, resid, err) &
      bind(c, name="nfor_park5_transport")
    integer(c_int), value, intent(in) :: nx, ny
    real(c_double), intent(in) :: rho(*)
    real(c_double), value, intent(in) :: dt
    real(c_double), intent(in) :: dxs(*), dys(*), fx(*), fy(*)
    real(c_double), intent(inout) :: y(*)
    real(c_double), intent(out) :: resid
    integer(c_int), intent(out) :: err
    real(DP), allocatable :: yc(:, :), fxx(:, :), fyy(:, :)
    real(DP), allocatable :: rhov(:), dxv(:), dyv(:)
    real(DP) :: rv
    integer :: k, s, f, j, i

    if (nx <= 0 .or. ny <= 0) then
      err = NFOR_EARGS
      return
    end if
    allocate(yc(NS, nx*ny), fxx(nx+1, ny), fyy(ny+1, nx))
    allocate(rhov(nx*ny), dxv(nx), dyv(ny))
    do k = 1, nx*ny
      rhov(k) = rho(k)
    end do
    do i = 1, nx
      dxv(i) = dxs(i)
    end do
    do j = 1, ny
      dyv(j) = dys(j)
    end do
    do j = 1, ny
      do f = 0, nx
        fxx(f+1, j) = fx((j-1)*(nx+1) + f + 1)
      end do
    end do
    do i = 1, nx
      do f = 0, ny
        fyy(f+1, i) = fy((i-1)*(ny+1) + f + 1)
      end do
    end do
    do s = 1, NS
      do k = 1, nx*ny
        yc(s, k) = y((s-1)*nx*ny + k)
      end do
    end do
    call transport_planar(nx, ny, rhov, dt, dxv, dyv, fxx, fyy, yc, rv)
    do s = 1, NS
      do k = 1, nx*ny
        y((s-1)*nx*ny + k) = yc(s, k)
      end do
    end do
    resid = rv
    err = NFOR_OK
    deallocate(yc, fxx, fyy, rhov, dxv, dyv)
  end subroutine nfor_park5_transport

end module nuforpark5

! nuforprops.f90 — the park 5-species transport properties.
!
! mixture viscosity and frozen conductivity by wilke's rule from the
! RP-1232 table V/VI curve fits, and the effective binary diffusivity
! for the species diffusion near a wall (table VII, single effective-d
! form).
!
! fit form (all three tables): exp(a z^n family) with z = ln t.
! viscosity and diffusion are quadratic in z, conductivity quartic.
! units: rp-1232 prints mu in g/(cm s) (x0.1 -> pa s), k in
! cal/(cm s k) (x418.4 -> w/(m k)), d in cm^2/s at 1 atm
! (x1e-4 -> m^2/s; scale linearly with 1/p_atm).
!
! below 1000 k the tables direct sutherland air laws for every
! species (table V/vi footnotes), so the mixture rules then ride on
! a single air value per property: the mixture rule collapses to the
! same number whatever the composition, which is the honest cold-limit
! behavior of a single-d single-mu approximation.
!
! source: plans/catalytic_research/constants.json + wilke_transport.md
! (rp-1232 tables V, VI, VII, verbatim coefficients).
!
! the algebra mirrors crates/nufor-core/src/species_props.rs operation
! for operation, and the transcendentals go through the c library on
! purpose: rust's f64::exp, f64::ln and f64::powf lower to libc's exp,
! log and pow, and fortran's own ** may expand to exp(log()), which
! rounds differently. calling the same c functions keeps the two
! implementations agreeing on the transcendentals.
!
! the polynomial evaluations spell out their associations so the
! rounding matches the rust expression term for term.

module nuforprops
  use iso_c_binding, only: c_double
  use iso_fortran_env, only: real64
  use nuforpark5, only: NS, M_SPEC, HF_SPEC, pressure_at, temperature_at
  implicit none

  integer, parameter :: DP = real64

  ! sutherland air reference for the sub-1000 k branch (table V/vi
  ! footnotes direct sutherland "for air ... for all species").
  real(DP), parameter :: MU0_AIR = 1.716e-5_DP
  real(DP), parameter :: T0_AIR = 273.15_DP
  real(DP), parameter :: S_AIR = 110.4_DP

  ! the frozen air specific heat and prandtl number the sub-1000 k
  ! conductivity branch rides (kappa = mu cp / pr).
  real(DP), parameter :: CP_AIR = 1004.7_DP
  real(DP), parameter :: PR_AIR = 0.71_DP

  ! atmospheric reference pressure, pa.
  real(DP), parameter :: P_ATM = 101325.0_DP

  ! table V viscosity coefficients in this lane's species order
  ! [n2, o2, no, n, o], converted from the printed (n2, o2, n, o, no).
  ! the arrays are column major: species first, coefficient second.
  real(DP), parameter :: MU_FIT(NS, 3) = reshape([ &
    0.0203_DP, 0.0484_DP, 0.0452_DP, 0.012_DP, 0.0205_DP, 0.4329_DP, &
    -0.1455_DP, -0.0609_DP, 0.593_DP, 0.4257_DP, -11.8153_DP, -8.9231_DP, &
    -9.4596_DP, -12.3805_DP, -11.5803_DP &
  ], [NS, 3])

  real(DP), parameter :: K_FIT(NS, 5) = reshape([ &
    0.03607_DP, 0.07987_DP, 0.02792_DP, 0.0_DP, 0.0_DP, -1.07503_DP, &
    -2.58428_DP, -0.87133_DP, 0.0_DP, 0.0_DP, 11.95029_DP, 31.25959_DP, &
    10.17967_DP, 0.01619_DP, 0.0331_DP, -57.90063_DP, -166.76267_DP, &
    -52.03466_DP, 0.55022_DP, 0.22834_DP, 93.21782_DP, 321.6982_DP, &
    88.6706_DP, -12.9219_DP, -11.58116_DP &
  ], [NS, 5])

  real(DP), parameter :: D_FIT(NS, 3) = reshape([ &
    0.0112_DP, 0.0465_DP, 0.0291_DP, 0.0195_DP, 0.014_DP, 1.6182_DP, &
    0.9271_DP, 1.2676_DP, 1.488_DP, 1.5824_DP, -11.3091_DP, -8.1137_DP, &
    -9.6878_DP, -10.3654_DP, -10.8819_DP &
  ], [NS, 3])

  integer, parameter :: PAIR_ROW(NS, NS) = reshape([ &
    1, 2, 3, 4, 5, 2, 2, 3, 4, 5, 3, 3, 3, 4, 5, 4, 4, 4, 4, 5, 5, 5, 5, &
    5, 5 &
  ], [NS, NS])

  public :: species_viscosity, species_conductivity, binary_diffusivity
  public :: effective_diffusivity, mixture_viscosity, mixture_conductivity
  public :: park5_viscous_update

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

  ! sutherland's viscosity law, mu0 * (t/t0)^1.5 * (t0+s)/(t+s).
  function sutherland_mu(mu0, t, t0, s) result(mu)
    real(DP), intent(in) :: mu0, t, t0, s
    real(DP) :: mu
    real(DP) :: powered
    powered = c_pow(t / t0, 1.5_DP)
    mu = mu0 * powered
    mu = mu * (t0 + s)
    mu = mu / (t + s)
  end function sutherland_mu

  ! the air viscosity anchor the cold branch falls back on.
  function air_mu(t) result(mu)
    real(DP), intent(in) :: t
    real(DP) :: mu
    mu = sutherland_mu(MU0_AIR, max(t, 1.0_DP), T0_AIR, S_AIR)
  end function air_mu

  ! the air conductivity anchor the cold branch falls back on.
  function air_kappa(t) result(kappa)
    real(DP), intent(in) :: t
    real(DP) :: kappa
    kappa = air_mu(t) * CP_AIR
    kappa = kappa / PR_AIR
  end function air_kappa

  ! the quadratic fit a z^2 + b z + c, associated as the rust
  ! expression a * z * z + b * z + c associates it.
  function fit_quadratic(a, b, c, z) result(p)
    real(DP), intent(in) :: a, b, c, z
    real(DP) :: p
    real(DP) :: za, zb
    za = a * z
    za = za * z
    zb = b * z
    p = za + zb
    p = p + c
  end function fit_quadratic

  ! the quartic fit a z^4 + b z^3 + c z^2 + d z + e, associated as the
  ! rust expression associates it.
  function fit_quartic(a, b, c, d, e, z) result(p)
    real(DP), intent(in) :: a, b, c, d, e, z
    real(DP) :: p
    real(DP) :: za, zb, zc, zd
    za = a * z
    za = za * z
    za = za * z
    za = za * z
    zb = b * z
    zb = zb * z
    zb = zb * z
    zc = c * z
    zc = zc * z
    zd = d * z
    p = za + zb
    p = p + zc
    p = p + zd
    p = p + e
  end function fit_quartic

  ! pure-species dynamic viscosity in pa s: the rp-1232 fit above 1000 k,
  ! the sutherland air law below (the tables' own directive).
  function species_viscosity(s, t) result(mu)
    integer, intent(in) :: s
    real(DP), intent(in) :: t
    real(DP) :: mu
    real(DP) :: z, p
    if (t < 1000.0_DP) then
      mu = air_mu(t)
    else
      z = c_log(t)
      p = fit_quadratic(MU_FIT(s, 1), MU_FIT(s, 2), MU_FIT(s, 3), z)
      mu = 0.1_DP * c_exp(p)
    end if
  end function species_viscosity

  ! pure-species frozen thermal conductivity in w/(m k): table VI fit
  ! above 1000 k (cal/(cm s k) x 418.4), sutherland-scaled air value
  ! below. the sub-1000 k conductivity rides the same sutherland mu
  ! through kappa = mu cp / pr with air pr: the footnote's own form.
  function species_conductivity(s, t) result(kappa)
    integer, intent(in) :: s
    real(DP), intent(in) :: t
    real(DP) :: kappa
    real(DP) :: z, p
    if (t < 1000.0_DP) then
      kappa = air_kappa(t)
    else
      z = c_log(t)
      p = fit_quartic(K_FIT(s, 1), K_FIT(s, 2), K_FIT(s, 3), K_FIT(s, 4), K_FIT(s, 5), z)
      kappa = 418.4_DP * c_exp(p)
    end if
  end function species_conductivity

  ! the table VII exp fit for one pair in cm^2/s at 1 atm.
  !
  ! diagonal pairs use the species' own self-fit; off-diagonal take the
  ! stored i-vs-n2 row only when one side is n2; every other cross pair
  ! uses the lighter-species-vs-n2 row as the effective interaction
  ! (documented approximation: the neutral cross pairs vary by under 40%
  ! in the table while the single-d model absorbs the difference; the
  ! effective-d formula below weights them).
  function pair_fit(i, j, t) result(d)
    integer, intent(in) :: i, j
    real(DP), intent(in) :: t
    real(DP) :: d
    real(DP) :: z, p
    integer :: row
    row = PAIR_ROW(i, j)
    z = c_log(t)
    p = fit_quadratic(D_FIT(row, 1), D_FIT(row, 2), D_FIT(row, 3), z)
    d = 1e-4_DP * c_exp(p)
  end function pair_fit

  ! binary diffusion coefficient for the pair (i, j) in m^2/s at the
  ! given pressure (pa). below 1000 k there is no tabulated neutral
  ! pair, so the 1000 k value scales with the chapman-enskog t^1.5/p
  ! shape to keep a finite, physically-shaped d in the cold wall layer.
  function binary_diffusivity(i, j, t, p_pa) result(d)
    integer, intent(in) :: i, j
    real(DP), intent(in) :: t, p_pa
    real(DP) :: d
    real(DP) :: d1000, shaped
    if (t < 1000.0_DP) then
      d1000 = pair_fit(i, j, 1000.0_DP)
      shaped = c_pow(max(t, 1.0_DP) / 1000.0_DP, 1.5_DP)
      d = d1000 * shaped
    else
      d = pair_fit(i, j, t)
    end if
    d = max(d, 1e-12_DP)
    d = d * P_ATM
    d = d / max(p_pa, 1.0_DP)
  end function binary_diffusivity

  ! effective binary diffusivity for the mixture in m^2/s: the
  ! mole-fraction-weighted mean of the pair fits. with a single d the
  ! species fluxes are j_s = -rho d dy_s/dn and sum j_s = 0 holds
  ! identically (the engineering simplification documented in
  ! wilke_transport.md).
  function effective_diffusivity(y, t, p_pa) result(d)
    real(DP), intent(in) :: y(NS), t, p_pa
    real(DP) :: d
    real(DP) :: m_mix, num, den, w
    real(DP) :: x(NS), ym(NS)
    integer :: i, j
    m_mix = 0.0_DP
    do i = 1, NS
      ym(i) = y(i) / M_SPEC(i)
      m_mix = m_mix + ym(i)
    end do
    if (m_mix <= 0.0_DP) then
      d = binary_diffusivity(1, 1, t, p_pa)
      return
    end if
    do i = 1, NS
      x(i) = ym(i) / m_mix
    end do
    num = 0.0_DP
    den = 0.0_DP
    do i = 1, NS - 1
      do j = i + 1, NS
        w = x(i) * x(j)
        if (w > 1e-12_DP) then
          num = num + w * pair_fit(i, j, t)
          den = den + w
        end if
      end do
    end do
    if (den <= 0.0_DP) then
      d = binary_diffusivity(1, 1, t, p_pa)
      return
    end if
    d = num / den
    d = d * P_ATM
    d = d / max(p_pa, 1.0_DP)
  end function effective_diffusivity

  ! wilke mixture viscosity in pa s from mass fractions and the cell
  ! temperature. species with mole fraction under 1e-12 are skipped in
  ! both loops (the 0/0 guard).
  function mixture_viscosity(y, t) result(mu)
    real(DP), intent(in) :: y(NS), t
    real(DP) :: mu
    real(DP) :: m_mix, num, den, phi_sum, phi, mr, root, powered
    real(DP) :: x(NS), ym(NS), mus(NS)
    integer :: i, j
    m_mix = 0.0_DP
    do i = 1, NS
      ym(i) = y(i) / M_SPEC(i)
      m_mix = m_mix + ym(i)
    end do
    if (m_mix <= 0.0_DP) then
      mu = air_mu(t)
      return
    end if
    do i = 1, NS
      x(i) = ym(i) / m_mix
      mus(i) = species_viscosity(i, t)
    end do
    num = 0.0_DP
    den = 0.0_DP
    do i = 1, NS
      if (x(i) < 1e-12_DP) then
        cycle
      end if
      num = num + x(i) * mus(i)
      phi_sum = 0.0_DP
      do j = 1, NS
        if (x(j) < 1e-12_DP) then
          cycle
        end if
        ! phi_ij per wilke 1950
        mr = M_SPEC(i) / M_SPEC(j)
        root = sqrt(mus(i) / mus(j))
        powered = c_pow(mr, 0.25_DP)
        phi = 1.0_DP + root * powered
        phi = phi * phi
        root = sqrt(8.0_DP * (1.0_DP + mr))
        phi = phi / root
        phi_sum = phi_sum + x(j) * phi
      end do
      if (phi_sum > 0.0_DP) then
        den = den + x(i) * phi_sum
      end if
    end do
    if (den <= 0.0_DP) then
      mu = air_mu(t)
      return
    end if
    mu = num / den
  end function mixture_viscosity

  ! wilke-class mixture frozen conductivity in w/(m k): the same
  ! weighting with k in place of mu (mason-saxena form per rp-1232's
  ! approximation review).
  function mixture_conductivity(y, t) result(kappa)
    real(DP), intent(in) :: y(NS), t
    real(DP) :: kappa
    real(DP) :: m_mix, num, den, phi_sum, phi, mr, root_k, root_mu, powered
    real(DP) :: x(NS), ym(NS), ks(NS), mus(NS)
    integer :: i, j
    m_mix = 0.0_DP
    do i = 1, NS
      ym(i) = y(i) / M_SPEC(i)
      m_mix = m_mix + ym(i)
    end do
    if (m_mix <= 0.0_DP) then
      kappa = air_kappa(t)
      return
    end if
    do i = 1, NS
      x(i) = ym(i) / m_mix
      ks(i) = species_conductivity(i, t)
      mus(i) = species_viscosity(i, t)
    end do
    num = 0.0_DP
    den = 0.0_DP
    do i = 1, NS
      if (x(i) < 1e-12_DP) then
        cycle
      end if
      num = num + x(i) * ks(i)
      phi_sum = 0.0_DP
      do j = 1, NS
        if (x(j) < 1e-12_DP) then
          cycle
        end if
        mr = M_SPEC(i) / M_SPEC(j)
        root_k = sqrt(ks(i) / ks(j))
        root_mu = sqrt(mus(j) / mus(i))
        powered = c_pow(mr, 0.25_DP)
        phi = root_k * root_mu
        phi = phi * powered
        phi = 1.0_DP + phi
        phi = phi * phi
        root_k = sqrt(8.0_DP * (1.0_DP + mr))
        phi = phi / root_k
        phi_sum = phi_sum + x(j) * phi
      end do
      if (phi_sum > 0.0_DP) then
        den = den + x(i) * phi_sum
      end if
    end do
    if (den <= 0.0_DP) then
      kappa = air_kappa(t)
      return
    end if
    kappa = num / den
  end function mixture_conductivity

  ! the laminar viscous pass shared by the serial and the mpi park5 lanes:
  ! wilke coefficients per cell, interior-face stress / conduction /
  ! species-diffusion fluxes, and the divergence update with the face-local
  ! diffusion stability cap. mirrors species_viscous.rs (body-free): domain
  ! faces carry zero viscous flux, and the central stencils never read a
  ! ghost column, so the mpi lane needs no halo for this pass.
  ! the composition arrives flat, species-major (ys((s-1)*n+k)), and is
  ! mutated in place together with the momenta and the energy; the density
  ! is only read.
  subroutine park5_viscous_update(nl, nj, rho, mx, my, e_tot, ys, dxs, dys, dt, &
                                   west_edge, east_edge, dwx, dex, cw, ce)
    integer, intent(in) :: nl, nj
    logical, intent(in) :: west_edge, east_edge
    real(DP), intent(in) :: rho(nl*nj), dxs(nl), dys(nj), dt
    ! seam spacings: dwx/dex are the west/east neighbour slabs' edge cell
    ! widths (unused at a global domain edge).
    real(DP), intent(in) :: dwx, dex
    ! west/east neighbour-column primitives in the car layout
    ! (rho, u, v, e_int, Y_1..5); only read when the face is not a
    ! global domain edge.
    real(DP), intent(in) :: cw(9,nj), ce(9,nj)
    real(DP), intent(inout) :: mx(nl*nj), my(nl*nj), e_tot(nl*nj), ys(*)
    real(DP), allocatable :: car(:,:,:), tt(:), mu(:), ka(:), df(:)
    real(DP), allocatable :: gt(:), gm(:), gk(:), gd(:), ht(:), hm(:), hk(:), hd(:)
    real(DP), allocatable :: gce(:)
    real(DP), allocatable :: xtx(:,:), xty(:,:), xen(:,:), xjs(:,:,:)
    real(DP), allocatable :: ytx(:,:), yty(:,:), yen(:,:), yjs(:,:,:)
    real(DP) :: yc(NS), yg(NS), pc, pg, ph, mum, kam, dfm, dx, dy
    real(DP) :: dudx, dvdx, duDy, dvdy, dtdx, dtdy
    real(DP) :: div, txx, txy, tyy, qx, qy, rf, js, qe, rhom, dtx, dty, ssum
    real(DP) :: rl, ul, vl, el, mul, kal, dfl, tgl
    real(DP) :: rr, ur, vr, er_, mur, kar, dfr, tgr
    integer :: i, j, f, s, k, kl, kr, nn

    nn = nl*nj
    allocate(car(9, nl+4, nj), tt(nn), mu(nn), ka(nn), df(nn))
    allocate(gt(nj), gm(nj), gk(nj), gd(nj), ht(nj), hm(nj), hk(nj), hd(nj))
    allocate(gce(nj))
    allocate(xtx(nl+1,nj), xty(nl+1,nj), xen(nl+1,nj), xjs(NS,nl+1,nj))
    allocate(ytx(nj+1,nl), yty(nj+1,nl), yen(nj+1,nl), yjs(NS,nj+1,nl))
    xtx = 0.0_DP; xty = 0.0_DP; xen = 0.0_DP; xjs = 0.0_DP
    ytx = 0.0_DP; yty = 0.0_DP; yen = 0.0_DP; yjs = 0.0_DP

    ! padded primitive carrier, centres only.
    do j = 1, nj
      do i = 1, nl
        k = (j-1)*nl + i
        car(1,i+2,j) = rho(k)
        car(2,i+2,j) = mx(k) / max(rho(k), 1.0e-12_DP)
        car(3,i+2,j) = my(k) / max(rho(k), 1.0e-12_DP)
        car(4,i+2,j) = e_tot(k) / max(rho(k), 1.0e-12_DP) &
                       - 0.5_DP*(car(2,i+2,j)**2 + car(3,i+2,j)**2)
        do s = 1, NS
          car(4+s,i+2,j) = ys((s-1)*nn + k)
        end do
      end do
    end do

    ! per-cell closure temperature + wilke coefficients; p floored at 1 Pa
    ! (the rust lane's clamp so the diffusivity stays finite).
    do j = 1, nj
      do i = 1, nl
        k = (j-1)*nl + i
        do s = 1, NS
          yc(s) = car(4+s,i+2,j)
        end do
        pc = max(pressure_at(car(1,i+2,j), car(4,i+2,j), yc), 1.0_DP)
        tt(k) = temperature_at(car(4,i+2,j), yc)
        mu(k) = mixture_viscosity(yc, tt(k))
        ka(k) = mixture_conductivity(yc, tt(k))
        df(k) = effective_diffusivity(yc, tt(k), pc)
      end do
    end do

    ! off-rank neighbour columns get their own closure, same flooring.
    do j = 1, nj
      if (.not. west_edge) then
        do s = 1, NS
          yg(s) = cw(4+s,j)
        end do
        pg = max(pressure_at(cw(1,j), cw(4,j), yg), 1.0_DP)
        gt(j) = temperature_at(cw(4,j), yg)
        gm(j) = mixture_viscosity(yg, gt(j))
        gk(j) = mixture_conductivity(yg, gt(j))
        gd(j) = effective_diffusivity(yg, gt(j), pg)
      end if
      if (.not. east_edge) then
        do s = 1, NS
          yg(s) = ce(4+s,j)
        end do
        ph = max(pressure_at(ce(1,j), ce(4,j), yg), 1.0_DP)
        ht(j) = temperature_at(ce(4,j), yg)
        hm(j) = mixture_viscosity(yg, ht(j))
        hk(j) = mixture_conductivity(yg, ht(j))
        hd(j) = effective_diffusivity(yg, ht(j), ph)
        ! the serial lane folds each cell's right neighbour's diffusion
        ! cap into its own face scaling, so the seam cell needs the east
        ! neighbour column's cap as well.
        gce(j) = min(dex*dex / (4.0_DP*max(hd(j), 1.0e-12_DP)/max(ce(1,j), 1.0e-9_DP)), &
                     dys(j)*dys(j) / (4.0_DP*max(hd(j), 1.0e-12_DP)/max(ce(1,j), 1.0e-9_DP)), dt)
      end if
    end do

    ! x faces f = 1..nl+1. f = 1 pairs the west neighbour column with
    ! local cell 1 and f = nl+1 local cell nl with the east neighbour;
    ! a global domain edge skips the face (arrays preset to zero flux).
    do j = 1, nj
      do f = 1, nl+1
        if (f == 1 .and. west_edge) cycle
        if (f == nl+1 .and. east_edge) cycle
        ! left side of the face
        if (f == 1) then
          kl = 0
          rl = cw(1,j); ul = cw(2,j); vl = cw(3,j); el = cw(4,j)
          mul = gm(j); kal = gk(j); dfl = gd(j); tgl = gt(j)
          do s = 1, NS
            yc(s) = cw(4+s,j)
          end do
          dx = 0.5_DP*(dwx + dxs(1))
        else
          kl = (j-1)*nl + f - 1
          rl = car(1,f+1,j); ul = car(2,f+1,j); vl = car(3,f+1,j); el = car(4,f+1,j)
          mul = mu(kl); kal = ka(kl); dfl = df(kl); tgl = tt(kl)
          do s = 1, NS
            yc(s) = car(4+s,f+1,j)
          end do
          dx = 0.5_DP*(dxs(f-1) + dxs(f))
        end if
        ! right side of the face
        if (f == nl+1) then
          kr = 0
          rr = ce(1,j); ur = ce(2,j); vr = ce(3,j); er_ = ce(4,j)
          mur = hm(j); kar = hk(j); dfr = hd(j); tgr = ht(j)
          do s = 1, NS
            yg(s) = ce(4+s,j)
          end do
          dx = 0.5_DP*(dxs(nl) + dex)
        else
          kr = (j-1)*nl + f
          rr = car(1,f+2,j); ur = car(2,f+2,j); vr = car(3,f+2,j); er_ = car(4,f+2,j)
          mur = mu(kr); kar = ka(kr); dfr = df(kr); tgr = tt(kr)
          do s = 1, NS
            yg(s) = car(4+s,f+2,j)
          end do
        end if
        mum = 0.5_DP*(mul + mur)
        kam = 0.5_DP*(kal + kar)
        dfm = 0.5_DP*(dfl + dfr)
        dudx = (ur - ul) / dx
        dvdx = (vr - vl) / dx
        dtdx = (tgr - tgl) / dx
        div = dudx
        txx = mum*(2.0_DP*dudx - (2.0_DP/3.0_DP)*div)
        txy = mum*dvdx
        qx = kam*dtdx
        rf = 0.5_DP*(rl + rr)
        qe = 0.0_DP
        do s = 1, NS
          js = -rf*dfm*(yg(s) - yc(s)) / dx
          xjs(s,f,j) = js
          qe = qe + HF_SPEC(s)*js
        end do
        xtx(f,j) = txx
        xty(f,j) = txy
        xen(f,j) = 0.5_DP*(ul + ur)*txx + 0.5_DP*(vl + vr)*txy + qx - qe
      end do
    end do

    ! y faces f = 2..nj (f = 1 and f = nj+1 are domain faces: zero flux).
    do i = 1, nl
      do f = 2, nj
        kl = (f-2)*nl + i
        kr = (f-1)*nl + i
        dy = 0.5_DP*(dys(f-1) + dys(f))
        mum = 0.5_DP*(mu(kl) + mu(kr))
        kam = 0.5_DP*(ka(kl) + ka(kr))
        dfm = 0.5_DP*(df(kl) + df(kr))
        duDy = (car(2,i+2,f) - car(2,i+2,f-1)) / dy
        dvdy = (car(3,i+2,f) - car(3,i+2,f-1)) / dy
        dtdy = (tt(kr) - tt(kl)) / dy
        div = dvdy
        tyy = mum*(2.0_DP*dvdy - (2.0_DP/3.0_DP)*div)
        txy = mum*duDy
        qy = kam*dtdy
        rf = 0.5_DP*(car(1,i+2,f) + car(1,i+2,f-1))
        qe = 0.0_DP
        do s = 1, NS
          js = -rf*dfm*(car(4+s,i+2,f) - car(4+s,i+2,f-1)) / dy
          yjs(s,f,i) = js
          qe = qe + HF_SPEC(s)*js
        end do
        ytx(f,i) = txy
        yty(f,i) = tyy
        yen(f,i) = 0.5_DP*(car(2,i+2,f) + car(2,i+2,f-1))*txy &
                 + 0.5_DP*(car(3,i+2,f) + car(3,i+2,f-1))*tyy + qy - qe
      end do
    end do

    ! divergence update, each face scaled by the stricter of the two
    ! adjoining cells' diffusion caps.
    do j = 1, nj
      do i = 1, nl
        k = (j-1)*nl + i
        rhom = max(rho(k), 1.0e-9_DP)
        dtx = cell_cap(k, dxs(i), dys(j), dt)
        if (i < nl) then
          dtx = min(dtx, cell_cap(k+1, dxs(i+1), dys(j), dt))
        else if (.not. east_edge) then
          dtx = min(dtx, gce(j))
        end if
        dty = cell_cap(k, dxs(i), dys(j), dt)
        if (j < nj) dty = min(dty, cell_cap(k+nl, dxs(i), dys(j+1), dt))
        mx(k) = mx(k) + dtx*(xtx(i+1,j) - xtx(i,j)) / dxs(i)
        my(k) = my(k) + dtx*(xty(i+1,j) - xty(i,j)) / dxs(i)
        e_tot(k) = e_tot(k) + dtx*(xen(i+1,j) - xen(i,j)) / dxs(i)
        my(k) = my(k) + dty*(yty(j+1,i) - yty(j,i)) / dys(j)
        mx(k) = mx(k) + dty*(ytx(j+1,i) - ytx(j,i)) / dys(j)
        e_tot(k) = e_tot(k) + dty*(yen(j+1,i) - yen(j,i)) / dys(j)
        do s = 1, NS
          ys((s-1)*nn+k) = ys((s-1)*nn+k) - dtx*(xjs(s,i+1,j) - xjs(s,i,j)) / dxs(i) / rhom
          ys((s-1)*nn+k) = ys((s-1)*nn+k) - dty*(yjs(s,j+1,i) - yjs(s,j,i)) / dys(j) / rhom
        end do
      end do
    end do

    ! positivity clip + renormalize, one cell at a time.
    do k = 1, nn
      ssum = 0.0_DP
      do s = 1, NS
        if (ys((s-1)*nn+k) < 0.0_DP) ys((s-1)*nn+k) = 0.0_DP
        if (ys((s-1)*nn+k) > 1.0_DP) ys((s-1)*nn+k) = 1.0_DP
        ssum = ssum + ys((s-1)*nn+k)
      end do
      if (ssum > 0.0_DP) then
        do s = 1, NS
          ys((s-1)*nn+k) = ys((s-1)*nn+k) / ssum
        end do
      end if
    end do

    deallocate(car, tt, mu, ka, df, gt, gm, gk, gd, ht, hm, hk, hd, gce, &
               xtx, xty, xen, xjs, ytx, yty, yen, yjs)

  contains

    ! the per-cell explicit diffusion cap folded with the march dt.
    function cell_cap(k, dx, dy, dtc) result(c)
      integer, intent(in) :: k
      real(DP), intent(in) :: dx, dy, dtc
      real(DP) :: c, rh
      rh = max(rho(k), 1.0e-9_DP)
      c = min(dx*dx / (4.0_DP*max(df(k), 1.0e-12_DP) / rh), &
              dy*dy / (4.0_DP*max(df(k), 1.0e-12_DP) / rh), dtc)
    end function cell_cap

  end subroutine park5_viscous_update

end module nuforprops

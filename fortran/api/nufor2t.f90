! nufor2t.f90 — the two-temperature lane in the mpi core. mirrors
! crates/nufor-core/src/park2t.rs operation for operation: the vibrational
! energy plane e_v is advected with donor cells on the cell-centred
! velocities, then relaxed toward the local translational temperature with
! a backward-euler landau-teller solve.
!
! the substep draws its own pressure and temperature from a constant gamma
! and the nitrogen gas constant, so it does not touch the closed-form
! closures the euler march uses. it also never writes the conserved state,
! so carrying the lane alongside the march cannot disturb the euler gates.
!
! the slab layout matches the rest of the bridge: ghost columns 1,2 (west)
! and nx_l+3,nx_l+4 (east), interior columns 3..nx_l+2. the donor stencil
! spans one cell, so one ghost column on each side is what the seam needs.
! the cell-centred velocities are owned by the conserving march; the plane
! crosses the rank seam through the same halo machinery the solid mask
! uses.
module nufor2t
  use iso_c_binding, only: c_double
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: DP = real64

  ! n2 molar mass, kg/mol.
  real(DP), parameter :: M_N2 = 28.016e-3_DP
  ! n2 characteristic vibrational temperature, K.
  real(DP), parameter :: THETA_N2 = 3393.0_DP
  ! specific gas constant of n2, J/(kg K).
  real(DP), parameter :: R_N2 = 8.314462_DP / M_N2
  ! avogadro constant scaled as the rust lane scales it.
  real(DP), parameter :: AVOGADRO = 6.022e23_DP

  private
  public :: park2t_substep_slab

  ! pi to double precision, spelled as the literal the rust reference
  ! uses rather than acos(-1), so the two agree bit for bit.
  real(DP), parameter :: PI_CONST = 3.141592653589793116_DP

  interface
    ! cbrt and pow are not fortran intrinsics, and rust's f64::cbrt and
    ! f64::powf lower to the c library's cbrt and pow. calling the same c
    ! functions keeps the two implementations bit-identical. fortran's own
    ! fractional ** may expand to exp(log()), which rounds differently.
    function c_cbrt(x) bind(c, name="cbrt") result(y)
      import c_double
      real(c_double), value :: x
      real(c_double) :: y
    end function c_cbrt

    function c_pow(x, y) bind(c, name="pow") result(z)
      import c_double
      real(c_double), value :: x, y
      real(c_double) :: z
    end function c_pow
  end interface

contains

  ! harmonic-oscillator vibrational energy per mass at vibrational
  ! temperature t_v, J/kg. mirrors park2t.rs e_vib_of_t_v.
  function e_vib_of_t_v(t_v) result(e_v)
    real(DP), intent(in) :: t_v
    real(DP) :: e_v
    real(DP) :: x, ex
    x = THETA_N2 / max(t_v, 1.0_DP)
    if (x > 40.0_DP) then
      e_v = 0.0_DP
      return
    end if
    ex = exp(x)
    e_v = R_N2 * THETA_N2 / max(ex - 1.0_DP, 1e-12_DP)
  end function e_vib_of_t_v

  ! millikan-white vibrational relaxation time for n2-n2 collisions with
  ! park's collision-limited cutoff. mirrors park2t.rs tau_vt_n2.
  function tau_vt_n2(t, p_pa, n) result(tau)
    real(DP), intent(in) :: t, p_pa, n
    real(DP) :: tau
    real(DP) :: mu_gmol, inv_t, a_s, b_s, expo, tau_mw, r_pair, c_bar, tau_cut
    ! reduced mass of the colliding pair n2-n2, g/mol.
    mu_gmol = M_N2 * 1e3_DP / 2.0_DP
    inv_t = 1.0_DP / max(t, 1.0_DP)
    a_s = 1.16e-3_DP * sqrt(mu_gmol) * c_pow(THETA_N2, 4.0_DP / 3.0_DP)
    b_s = 0.015_DP * c_pow(mu_gmol, 0.25_DP)
    expo = a_s * (c_cbrt(inv_t) - b_s) - 18.42_DP
    expo = min(max(expo, -50.0_DP), 50.0_DP)
    tau_mw = exp(expo) / max(p_pa, 1.0_DP)
    ! park cutoff: sigma_v = 1e-21 * t^2 m^2, mean thermal speed of the
    ! reduced pair.
    r_pair = 8.314462_DP / (mu_gmol * 1e-3_DP)
    c_bar = sqrt(8.0_DP * r_pair * t / PI_CONST)
    tau_cut = 1.0_DP / (1e-21_DP * t * t * max(n, 1.0_DP) * c_bar)
    tau = tau_mw + tau_cut
  end function tau_vt_n2

  ! one advect+relax substep on the vibrational plane of one rank's slab.
  !
  ! rho, mx, my, e carry the CONSERVED state after the march: the caller
  ! passes the post-euler state. ev carries the plane entering the step,
  ! with its ghost columns already holding the neighbouring ranks' planes,
  ! because the donor stencil reads one cell across the seam.
  !
  ! the x faces at a GLOBAL domain edge carry no neighbour, while a rank
  ! seam is interior: the halo supplies the neighbour's plane there. the
  ! caller passes the decomposition's edge flags so the twin can tell the
  ! two apart, which is the difference between a sealed boundary and a
  ! leak across the partition.
  subroutine park2t_substep_slab(rho, mx, my, e, ev, solid, nx_l, ny_, &
                                  dx, dy, gamma, dt, west_edge, east_edge, ev_out)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                            my(nx_l+4, ny_), e(nx_l+4, ny_)
    real(DP), intent(in) :: ev(nx_l+4, ny_)
    logical, intent(in) :: solid(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, gamma, dt
    logical, intent(in) :: west_edge, east_edge
    real(DP), intent(out) :: ev_out(nx_l+4, ny_)

    real(DP) :: ev_new(nx_l+4, ny_)
    real(DP) :: p(nx_l+4, ny_), t(nx_l+4, ny_)
    real(DP) :: u, v, rho_raw, rho_c, et, ein_full, ein_tr, d, f
    real(DP) :: num_den, tau, e_eq, w, ev_relax
    logical :: has_west, has_east
    integer :: i, j

    ev_new = ev

    ! translational pressure and temperature from the post-euler state:
    ! the internal pool holds e_v, so the split removes it before the gas
    ! constant sees it.
    do j = 1, ny_
      do i = 3, nx_l+2
        rho_raw = rho(i,j)
        u = mx(i,j) / max(rho_raw, 1.0e-12_DP)
        v = my(i,j) / max(rho_raw, 1.0e-12_DP)
        et = e(i,j) / max(rho_raw, 1.0e-12_DP)
        ein_full = et - 0.5_DP * (u*u + v*v)
        ein_tr = max(ein_full - ev(i,j), 1.0_DP)
        p(i,j) = (gamma - 1.0_DP) * rho_raw * ein_tr
        t(i,j) = p(i,j) / (max(rho_raw, 1.0e-12_DP) * R_N2)
      end do
    end do

    ! donor-cell advection of e_v on the cell-centred velocities. solid
    ! faces carry no mass flux, so the upwind donor there is the same cell:
    ! no transfer.
    do j = 1, ny_
      do i = 3, nx_l+2
        if (solid(i,j)) cycle
        rho_c = max(rho(i,j), 1.0e-12_DP)
        u = mx(i,j) / rho_c
        v = my(i,j) / rho_c
        d = 0.0_DP
        ! the west face has a neighbour everywhere except at the global
        ! west edge: on a rank seam the ghost column holds it.
        has_west = .not. (west_edge .and. i == 3)
        if (has_west) then
          if (.not. solid(i-1,j)) then
            if (u >= 0.0_DP) then
              f = rho_c * u * ev(i-1,j)
            else
              f = rho_c * u * ev(i,j)
            end if
            d = d - f / dx
          end if
        end if
        ! the east face likewise has a neighbour everywhere except at the
        ! global east edge.
        has_east = .not. (east_edge .and. i == nx_l+2)
        if (has_east) then
          if (.not. solid(i+1,j)) then
            if (u >= 0.0_DP) then
              f = rho_c * u * ev(i,j)
            else
              f = rho_c * u * ev(i+1,j)
            end if
            d = d + f / dx
          end if
        end if
        ! the y extent is not partitioned, so the southern and northern
        ! edges are always global edges.
        if (j > 1) then
          if (.not. solid(i,j-1)) then
            if (v >= 0.0_DP) then
              f = rho_c * v * ev(i,j-1)
            else
              f = rho_c * v * ev(i,j)
            end if
            d = d - f / dy
          end if
        end if
        if (j < ny_) then
          if (.not. solid(i,j+1)) then
            if (v >= 0.0_DP) then
              f = rho_c * v * ev(i,j)
            else
              f = rho_c * v * ev(i,j+1)
            end if
            d = d + f / dy
          end if
        end if
        ev_new(i,j) = max(ev(i,j) + dt * d, 0.0_DP)
      end do
    end do

    ! backward-euler relaxation with the translational draw. the
    ! vibrational pool lives inside the conserved total, so the draw moves
    ! energy between pools and the total is conserved by construction.
    do j = 1, ny_
      do i = 3, nx_l+2
        if (solid(i,j)) cycle
        rho_c = max(rho(i,j), 1.0e-12_DP)
        num_den = rho_c / M_N2 * AVOGADRO
        tau = tau_vt_n2(t(i,j), p(i,j), num_den)
        e_eq = e_vib_of_t_v(t(i,j))
        w = dt / max(tau, 1.0e-9_DP)
        ev_relax = (ev_new(i,j) + w * e_eq) / (1.0_DP + w)
        ev_new(i,j) = ev_relax
      end do
    end do

    ev_out = ev_new
  end subroutine park2t_substep_slab

end module nufor2t

! nuforpark5abi.f90 — the c surface for the park 5-species lane.
!
! the physics kernels live in nuforpark5; this module is the abi wrapper
! the rust core calls. it converts between the flat species-major buffers
! the rust side passes and the per-cell vectors the kernels want, and it
! owns the error codes.
!
! the split keeps nuforpark5.f90 a physics module and this one a boundary
! module, so neither grows past the size ceiling as the lane gains exports.

module nuforpark5abi
  use iso_c_binding, only: c_double, c_int
  use iso_fortran_env, only: real64
  use nuforpark5
  use nuforprops, only: mixture_viscosity, mixture_conductivity, &
                       effective_diffusivity
  implicit none

  integer, parameter :: DP = real64

  ! structured error codes shared with the rust wrapper.
  integer(c_int), parameter :: NFOR_OK    = 0
  integer(c_int), parameter :: NFOR_EARGS = 1

  public :: nfor_park5_closure
  public :: nfor_park5_omega
  public :: nfor_park5_chem_substep
  public :: nfor_park5_transport
  public :: nfor_park5_formation_enthalpy
  public :: nfor_park5_seed_energy
  public :: nfor_park5_internal_energy
  public :: nfor_park5_equilibria
  public :: nfor_park5_props

contains

  ! ---------------- c surface ----------------

  ! the batch transport-property kernel for n cells: the wilke mixture
  ! viscosity, the frozen mixture conductivity and the effective binary
  ! diffusivity, one value per cell.
  !
  ! t and p arrive per cell because the diffusivity scales with pressure
  ! and the fits with temperature; y holds the species planes
  ! contiguously (species s occupies y((s-1)*n+1 : s*n)).
  subroutine nfor_park5_props(n, y, t, p, mu, kappa, d, err) &
      bind(c, name="nfor_park5_props")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: y(*), t(*), p(*)
    real(c_double), intent(out) :: mu(*), kappa(*), d(*)
    integer(c_int), intent(out) :: err
    real(DP) :: yc(NS)
    integer :: k, s
    if (n <= 0) then
      err = NFOR_EARGS
      return
    end if
    do k = 1, int(n)
      do s = 1, NS
        yc(s) = y((s - 1) * int(n) + k)
      end do
      mu(k) = mixture_viscosity(yc, t(k))
      kappa(k) = mixture_conductivity(yc, t(k))
      d(k) = effective_diffusivity(yc, t(k), p(k))
    end do
    err = NFOR_OK
  end subroutine nfor_park5_props

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

  ! the per-species formation enthalpy per mass, J/kg. the rust side
  ! reads this table instead of keeping its own copy, so the two
  ! implementations cannot drift apart on the constants.
  subroutine nfor_park5_formation_enthalpy(hf, err) &
      bind(c, name="nfor_park5_formation_enthalpy")
    real(c_double), intent(out) :: hf(NS)
    integer(c_int), intent(out) :: err
    integer :: s
    do s = 1, NS
      hf(s) = HF_SPEC(s)
    end do
    err = NFOR_OK
  end subroutine nfor_park5_formation_enthalpy

  ! the mixture gas constant and the internal energy at a temperature,
  ! for the initial-condition seeding path. the rust side reads these
  ! from the fortran core rather than keeping its own formulas.
  subroutine nfor_park5_seed_energy(n, rho, p, y, t_out, e_out, err) &
      bind(c, name="nfor_park5_seed_energy")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: rho(*), p(*), y(*)
    real(c_double), intent(out) :: t_out(*), e_out(*)
    integer(c_int), intent(out) :: err
    real(DP) :: yc(NS), t
    integer :: k, s

    if (n <= 0) then
      err = NFOR_EARGS
      return
    end if
    do k = 1, n
      do s = 1, NS
        yc(s) = y((s-1)*n + k)
      end do
      t = p(k) / max(rho(k) * mix_gas_constant(yc), 1.0e-30_DP)
      t_out(k) = t
      e_out(k) = internal_energy_at(t, yc)
    end do
    err = NFOR_OK
  end subroutine nfor_park5_seed_energy

  ! the internal energy per mass at a temperature and composition, for
  ! the boundary conditions that state a temperature instead of a
  ! pressure. the rust side reads this from the fortran core.
  subroutine nfor_park5_internal_energy(n, t, y, e_out, err) &
      bind(c, name="nfor_park5_internal_energy")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: t(*), y(*)
    real(c_double), intent(out) :: e_out(*)
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
      e_out(k) = internal_energy_at(t(k), yc)
    end do
    err = NFOR_OK
  end subroutine nfor_park5_internal_energy

  ! the equilibrium constants of the five reactions at a temperature. the
  ! rust side reads these from the fortran core rather than keeping its
  ! own fits.
  subroutine nfor_park5_equilibria(n, t, kc, err) &
      bind(c, name="nfor_park5_equilibria")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: t(*)
    real(c_double), intent(out) :: kc(5, *)
    integer(c_int), intent(out) :: err
    integer :: k, r

    if (n <= 0) then
      err = NFOR_EARGS
      return
    end if
    do k = 1, n
      do r = 1, 5
        kc(r, k) = k_eq(r, t(k))
      end do
    end do
    err = NFOR_OK
  end subroutine nfor_park5_equilibria

end module nuforpark5abi

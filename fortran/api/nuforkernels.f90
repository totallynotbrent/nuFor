! c-interoperable kernel surface of the fortran numerical layer (spec 39).

module nuforkernels
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: nfor_version
  public :: nfor_saxpy
  public :: nfor_grid1d_init
  public :: nfor_prim_to_cons
  public :: nfor_cons_to_prim
  public :: nfor_eos_pressure
  public :: nfor_eos_sound_speed
  public :: nfor_eos_mach
  public :: nfor_eos_temperature
  public :: nfor_hll_flux
  public :: nfor_cfl_dt

  ! structured error codes shared with the rust wrapper.
  integer(c_int), parameter :: NFOR_OK    = 0  ! success
  integer(c_int), parameter :: NFOR_EARGS = 1  ! invalid argument (count <= 0)
  integer(c_int), parameter :: NFOR_EDATA = 2  ! numerical failure in kernel
contains

  ! returns the kernel library version as a nul-terminated string.
  subroutine nfor_version(ver, ver_len, err) bind(c, name="nfor_version")
    character(kind=c_char, len=1), intent(out) :: ver(*)
    integer(c_int), value, intent(in)          :: ver_len
    integer(c_int), intent(out)                :: err
    character(kind=c_char, len=*), parameter   :: v = "0.1.0"
    integer :: n, i
    if (ver_len <= 0) then
       err = NFOR_EARGS
       return
    end if
    n = min(len(v), ver_len - 1)
    do i = 1, n
       ver(i) = v(i:i)
    end do
    ver(n + 1) = c_null_char
    err = NFOR_OK
  end subroutine nfor_version

  ! y = alpha*x + y elementwise over count contiguous doubles.
  subroutine nfor_saxpy(count, alpha, x, y, err) bind(c, name="nfor_saxpy")
    integer(c_int), value, intent(in)  :: count
    real(c_double), value, intent(in)  :: alpha
    real(c_double), intent(in)         :: x(*)
    real(c_double), intent(inout)      :: y(*)
    integer(c_int), intent(out)        :: err
    integer :: i
    if (count <= 0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, count
       y(i) = alpha * x(i) + y(i)
    end do
    err = NFOR_OK
  end subroutine nfor_saxpy

  ! uniform 1D control-volume geometry: n cells over [xmin, xmax], faces close the domain.
  subroutine nfor_grid1d_init(n, xmin, xmax, centers, faces, dx, err) bind(c, name="nfor_grid1d_init")
    integer(c_int), value, intent(in) :: n
    real(c_double), value, intent(in) :: xmin
    real(c_double), value, intent(in) :: xmax
    real(c_double), intent(out)       :: centers(*)
    real(c_double), intent(out)       :: faces(*)
    real(c_double), intent(out)       :: dx
    integer(c_int), intent(out)       :: err
    real(c_double) :: h
    integer :: i
    if (n < 2 .or. xmax <= xmin) then
       err = NFOR_EARGS
       return
    end if
    h = (xmax - xmin) / real(n, c_double)
    do i = 1, n
       centers(i) = xmin + (real(i, c_double) - 0.5d0) * h
    end do
    do i = 1, n
       faces(i) = xmin + real(i - 1, c_double) * h
    end do
    faces(n + 1) = xmax
    dx = h
    err = NFOR_OK
  end subroutine nfor_grid1d_init

  ! primitives (rho,u,e_t) to conserved (rho,m,E); density must be positive.
  subroutine nfor_prim_to_cons(n, rho, u, et, m, E, err) bind(c, name="nfor_prim_to_cons")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho(*)
    real(c_double), intent(in)        :: u(*)
    real(c_double), intent(in)        :: et(*)
    real(c_double), intent(out)       :: m(*)
    real(c_double), intent(out)       :: E(*)
    integer(c_int), intent(out)       :: err
    integer :: i
    if (n <= 0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       if (rho(i) <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       m(i) = rho(i) * u(i)
       E(i) = rho(i) * et(i)
    end do
    err = NFOR_OK
  end subroutine nfor_prim_to_cons

  ! conserved (rho,m,E) to primitives (rho,u,e_t), the inverse of the above.
  subroutine nfor_cons_to_prim(n, rho, m, E, u, et, err) bind(c, name="nfor_cons_to_prim")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho(*)
    real(c_double), intent(in)        :: m(*)
    real(c_double), intent(in)        :: E(*)
    real(c_double), intent(out)       :: u(*)
    real(c_double), intent(out)       :: et(*)
    integer(c_int), intent(out)       :: err
    integer :: i
    if (n <= 0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       if (rho(i) <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       u(i) = m(i) / rho(i)
       et(i) = E(i) / rho(i)
    end do
    err = NFOR_OK
  end subroutine nfor_cons_to_prim

  ! ideal-gas pressure p = (gamma-1)*rho*e_int; non-positive density/internal energy invalid.
  subroutine nfor_eos_pressure(gamma, n, rho, et, u, p, err) bind(c, name="nfor_eos_pressure")
    real(c_double), value, intent(in) :: gamma
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho(*)
    real(c_double), intent(in)        :: et(*)
    real(c_double), intent(in)        :: u(*)
    real(c_double), intent(out)       :: p(*)
    integer(c_int), intent(out)       :: err
    real(c_double) :: e_int
    integer :: i
    if (n <= 0 .or. .not. ieee_is_finite(gamma) .or. gamma <= 1.0d0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       e_int = et(i) - 0.5d0 * u(i) * u(i)
       if (.not. ieee_is_finite(rho(i)) .or. .not. ieee_is_finite(e_int) .or. &
           rho(i) <= 0.0d0 .or. e_int <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       p(i) = (gamma - 1.0d0) * rho(i) * e_int
    end do
    err = NFOR_OK
  end subroutine nfor_eos_pressure

  ! ideal-gas sound speed a = sqrt(gamma*p/rho); needs positive density and pressure.
  subroutine nfor_eos_sound_speed(gamma, n, rho, p, a, err) bind(c, name="nfor_eos_sound_speed")
    real(c_double), value, intent(in) :: gamma
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho(*)
    real(c_double), intent(in)        :: p(*)
    real(c_double), intent(out)       :: a(*)
    integer(c_int), intent(out)       :: err
    integer :: i
    if (n <= 0 .or. .not. ieee_is_finite(gamma) .or. gamma <= 1.0d0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       if (.not. ieee_is_finite(rho(i)) .or. .not. ieee_is_finite(p(i)) .or. &
           rho(i) <= 0.0d0 .or. p(i) <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       a(i) = sqrt(gamma * p(i) / rho(i))
    end do
    err = NFOR_OK
  end subroutine nfor_eos_sound_speed

  ! mach number M = |u|/a from flow velocity and sound speed.
  subroutine nfor_eos_mach(n, u, a, mach, err) bind(c, name="nfor_eos_mach")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: u(*)
    real(c_double), intent(in)        :: a(*)
    real(c_double), intent(out)       :: mach(*)
    integer(c_int), intent(out)       :: err
    integer :: i
    if (n <= 0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       if (.not. ieee_is_finite(u(i)) .or. .not. ieee_is_finite(a(i)) .or. a(i) <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       mach(i) = abs(u(i)) / a(i)
    end do
    err = NFOR_OK
  end subroutine nfor_eos_mach

  ! ideal-gas temperature from p = rho*R*T; all arguments must be positive.
  subroutine nfor_eos_temperature(r, n, rho, p, t, err) bind(c, name="nfor_eos_temperature")
    real(c_double), value, intent(in) :: r
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho(*)
    real(c_double), intent(in)        :: p(*)
    real(c_double), intent(out)       :: t(*)
    integer(c_int), intent(out)       :: err
    integer :: i
    if (n <= 0 .or. .not. ieee_is_finite(r) .or. r <= 0.0d0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       if (.not. ieee_is_finite(rho(i)) .or. .not. ieee_is_finite(p(i)) .or. &
           rho(i) <= 0.0d0 .or. p(i) <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       t(i) = p(i) / (rho(i) * r)
    end do
    err = NFOR_OK
  end subroutine nfor_eos_temperature

  ! HLL flux for n 1D euler faces using davis wave-speed estimates s_l, s_r.
  subroutine nfor_hll_flux(gamma, n, rho_l, m_l, e_l, rho_r, m_r, e_r, f_rho, f_m, f_e, err) bind(c, name="nfor_hll_flux")
    real(c_double), value, intent(in) :: gamma
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho_l(*)
    real(c_double), intent(in)        :: m_l(*)
    real(c_double), intent(in)        :: e_l(*)
    real(c_double), intent(in)        :: rho_r(*)
    real(c_double), intent(in)        :: m_r(*)
    real(c_double), intent(in)        :: e_r(*)
    real(c_double), intent(out)       :: f_rho(*)
    real(c_double), intent(out)       :: f_m(*)
    real(c_double), intent(out)       :: f_e(*)
    integer(c_int), intent(out)       :: err
    real(c_double) :: u_l, p_l, a_l, u_r, p_r, a_r, s_l, s_r
    real(c_double) :: fl_rho, fl_m, fl_e, fr_rho, fr_m, fr_e
    real(c_double) :: e_int, denom
    integer :: i
    if (n <= 0 .or. .not. ieee_is_finite(gamma) .or. gamma <= 1.0d0) then
       err = NFOR_EARGS
       return
    end if
    do i = 1, n
       ! recover primitives and sound speeds on both states.
       e_int = e_l(i) / rho_l(i) - 0.5d0 * (m_l(i) / rho_l(i))**2
       if (.not. ieee_is_finite(rho_l(i)) .or. .not. ieee_is_finite(m_l(i)) .or. &
           .not. ieee_is_finite(e_l(i)) .or. rho_l(i) <= 0.0d0 .or. e_int <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       e_int = e_r(i) / rho_r(i) - 0.5d0 * (m_r(i) / rho_r(i))**2
       if (.not. ieee_is_finite(rho_r(i)) .or. .not. ieee_is_finite(m_r(i)) .or. &
           .not. ieee_is_finite(e_r(i)) .or. rho_r(i) <= 0.0d0 .or. e_int <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       u_l = m_l(i) / rho_l(i)
       p_l = (gamma - 1.0d0) * (e_l(i) - 0.5d0 * rho_l(i) * u_l * u_l)
       a_l = sqrt(gamma * p_l / rho_l(i))
       u_r = m_r(i) / rho_r(i)
       p_r = (gamma - 1.0d0) * (e_r(i) - 0.5d0 * rho_r(i) * u_r * u_r)
       a_r = sqrt(gamma * p_r / rho_r(i))
       s_l = min(u_l - a_l, u_r - a_r)
       s_r = max(u_l + a_l, u_r + a_r)
       fl_rho = rho_l(i) * u_l
       fl_m = rho_l(i) * u_l * u_l + p_l
       fl_e = u_l * (e_l(i) + p_l)
       fr_rho = rho_r(i) * u_r
       fr_m = rho_r(i) * u_r * u_r + p_r
       fr_e = u_r * (e_r(i) + p_r)
       if (s_l >= 0.0d0) then
          f_rho(i) = fl_rho
          f_m(i) = fl_m
          f_e(i) = fl_e
       else if (s_r <= 0.0d0) then
          f_rho(i) = fr_rho
          f_m(i) = fr_m
          f_e(i) = fr_e
       else
          denom = 1.0d0 / (s_r - s_l)
          f_rho(i) = denom * (s_r * fl_rho - s_l * fr_rho + s_l * s_r * (rho_r(i) - rho_l(i)))
          f_m(i) = denom * (s_r * fl_m - s_l * fr_m + s_l * s_r * (m_r(i) - m_l(i)))
          f_e(i) = denom * (s_r * fl_e - s_l * fr_e + s_l * s_r * (e_r(i) - e_l(i)))
       end if
    end do
    err = NFOR_OK
  end subroutine nfor_hll_flux

  ! explicit time step from the CFL condition dt = cfl*dx/s_max; conserved states in, s_max and dt out.
  subroutine nfor_cfl_dt(gamma, cfl, dx, n, rho, m, e, s_max, dt, err) bind(c, name="nfor_cfl_dt")
    real(c_double), value, intent(in) :: gamma
    real(c_double), value, intent(in) :: cfl
    real(c_double), value, intent(in) :: dx
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in)        :: rho(*)
    real(c_double), intent(in)        :: m(*)
    real(c_double), intent(in)        :: e(*)
    real(c_double), intent(out)       :: s_max
    real(c_double), intent(out)       :: dt
    integer(c_int), intent(out)       :: err
    real(c_double) :: u, e_int, p, a, speed
    integer :: i
    if (n <= 0 .or. .not. ieee_is_finite(gamma) .or. gamma <= 1.0d0 .or. &
        .not. ieee_is_finite(cfl) .or. cfl <= 0.0d0 .or. cfl > 1.0d0 .or. &
        .not. ieee_is_finite(dx) .or. dx <= 0.0d0) then
       err = NFOR_EARGS
       return
    end if
    s_max = 0.0d0
    do i = 1, n
       if (.not. ieee_is_finite(rho(i)) .or. .not. ieee_is_finite(m(i)) .or. &
           .not. ieee_is_finite(e(i)) .or. rho(i) <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       u = m(i) / rho(i)
       e_int = e(i) / rho(i) - 0.5d0 * u * u
       if (e_int <= 0.0d0) then
          err = NFOR_EDATA
          return
       end if
       p = (gamma - 1.0d0) * rho(i) * e_int
       a = sqrt(gamma * p / rho(i))
       speed = abs(u) + a
       if (speed > s_max) s_max = speed
    end do
    dt = cfl * dx / s_max
    err = NFOR_OK
  end subroutine nfor_cfl_dt

end module nuforkernels
! C-interoperable kernel surface of the Fortran numerical layer (spec 39).
! Every exported routine is bind(C) with contiguous arrays and explicit
! lengths; no derived types cross the boundary; errors are structured codes.

module nuforkernels
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: nfor_version
  public :: nfor_saxpy
  public :: nfor_grid1d_init
  public :: nfor_prim_to_cons
  public :: nfor_cons_to_prim

  ! Structured error codes shared with the Rust wrapper (src/lib.rs codes).
  integer(c_int), parameter :: NFOR_OK    = 0  ! success
  integer(c_int), parameter :: NFOR_EARGS = 1  ! invalid argument (count <= 0)
  integer(c_int), parameter :: NFOR_EDATA = 2  ! numerical failure in kernel
contains

  ! Returns the kernel library version as a NUL-terminated C string.
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

  ! y = alpha * x + y elementwise over `count` contiguous real(c_double).
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

  ! Uniform 1D control-volume geometry (spec 18): `n` cells over [xmin, xmax].
  ! centers(i) = xmin + (i - 1/2)*dx; face j sits at xmin + (j - 1)*dx with
  ! faces(1) = xmin and faces(n + 1) = xmax closing the domain exactly.
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

  ! Primitives (rho, u, e_t) -> conserved (rho, m, E): m = rho*u, E = rho*e_t.
  ! e_t is the total specific energy; the ideal-gas internal/kinetic split and
  ! pressure recovery belong to the EOS step. Density must be positive.
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

  ! Conserved (rho, m, E) -> primitives (rho, u, e_t): u = m/rho, e_t = E/rho.
  ! Inverse of nfor_prim_to_cons; density must be positive.
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

end module nuforkernels
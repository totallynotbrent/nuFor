! C-interoperable kernel surface of the Fortran numerical layer (spec 39).
! Every exported routine is bind(C) with contiguous arrays and explicit
! lengths; no derived types cross the boundary; errors are structured codes.

module nuforkernels
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: nfor_version
  public :: nfor_saxpy

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

end module nuforkernels
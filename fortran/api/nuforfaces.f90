! nuforfaces.f90 — the reusable per-strip face workbank: muscl face states,
! per-face closure properties and limited slopes, indexed [face]. callers
! fill it strip by strip (one row or one column at a time), so no rank-local
! heap allocation is needed inside the sweep.
module nuforfaces
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: DP = real64

  integer, parameter :: MAXN = 4096
  real(DP) :: face_l(5, MAXN), face_r(5, MAXN)
  real(DP) :: face_e_l(MAXN), face_e_r(MAXN), face_a_l(MAXN), face_a_r(MAXN)
  real(DP) :: face_p_l(MAXN), face_p_r(MAXN)
  real(DP) :: slope(MAXN), slope_y(MAXN)
end module nuforfaces

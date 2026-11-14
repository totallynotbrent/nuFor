! nufoeqair.f90 — the eqair thermochemical closure: bilinear cea (rho, e)->(t, p)
! table plus the 12-block tgas1 face fits shipped from rust at init, and the
! point evaluations (pressure, temperature, sound speed) used on every face.
module nufoeqair
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: DP = real64

  ! cea bilinear table (t, p per node), shipped from rust at init.
  integer :: cea_nr = 0, cea_ne = 0
  real(DP), allocatable :: cea_lr(:), cea_le(:)   ! log10 rho / e axes
  real(DP), allocatable :: cea_t(:,:), cea_p(:,:) ! node values (nr, ne)

  ! tgas1 face-fit coefficients for the wake, in the eqair.rs block layout:
  ! g1..g4 polynomial parts plus the grabau sigmoid set, 12 blocks.
  integer, parameter :: NB = 12
  integer :: nblocks = NB
  real(DP) :: blk_g1(NB,2), blk_g2(NB,3), blk_g3(NB,3), blk_g4(NB,3)
  real(DP) :: blk_g5(NB,2), blk_g6(NB,3), blk_g7(NB,3), blk_g8(NB,3)
  real(DP) :: blk_s(NB,4)
  logical :: blk_minus(NB)
  integer :: nbk(NB)          ! grabau terms per block
  real(DP) :: rho0, e0, ln10  ! normalization constants
  real(DP) :: z_breaks(3,5)   ! region z boundaries (rust Z_BREAKS)
  real(DP) :: cold_gamm(3)    ! per-region cold-limit gamma
  logical :: tgas1_ready = .false.

contains


  ! init: receive the table and tgas1 coefficients from rust through the
  ! c ABI (nufor_mpi_init, see nuforbridge below). stored module-side.
  subroutine set_cea_table(nr, ne, lr, le, t, p)
    integer, intent(in) :: nr, ne
    real(DP), intent(in) :: lr(nr), le(ne), t(nr,ne), p(nr,ne)
    cea_nr = nr; cea_ne = ne
    if (allocated(cea_lr)) deallocate(cea_lr, cea_le, cea_t, cea_p)
    allocate(cea_lr(nr), cea_le(ne), cea_t(nr,ne), cea_p(nr,ne))
    cea_lr = lr; cea_le = le; cea_t = t; cea_p = p
  end subroutine


  subroutine set_tgas1(g1,g2,g3,g4,g5,g6,g7,g8,s,minus,nbk_,rho0_,e0_)
    real(DP), intent(in) :: g1(NB,2), g2(NB,3), g3(NB,3), g4(NB,3)
    real(DP), intent(in) :: g5(NB,2), g6(NB,3), g7(NB,3), g8(NB,3), s(NB,4)
    logical, intent(in) :: minus(NB)
    integer, intent(in) :: nbk_(NB)
    real(DP), intent(in) :: rho0_, e0_
    blk_g1 = g1; blk_g2 = g2; blk_g3 = g3; blk_g4 = g4
    blk_g5 = g5; blk_g6 = g6; blk_g7 = g7; blk_g8 = g8
    blk_s = s; blk_minus = minus; nbk = nbk_
    rho0 = rho0_; e0 = e0_; ln10 = log(10.0_DP)
  end subroutine


  ! bilinear cea table lookup for all cells: t and p from (rho, e).
  subroutine cea_pt_all(rho, e, n, tout, pout)
    integer, intent(in) :: n
    real(DP), intent(in) :: rho(n), e(n)
    real(DP), intent(out) :: tout(n), pout(n)
    integer :: k, i, j
    real(DP) :: lr, le, fi, fj, t, p
    do k = 1, n
      lr = log10(max(rho(k), 1.0e-12_DP))
      lr = min(max(lr, cea_lr(1)), cea_lr(cea_nr))
      le = log10(max(e(k), 1.0e-3_DP))
      le = min(max(le, cea_le(1)), cea_le(cea_ne))
      i = 1
      do while (i + 1 < cea_nr .and. cea_lr(i+1) < lr)
        i = i + 1
      end do
      j = 1
      do while (j + 1 < cea_ne .and. cea_le(j+1) < le)
        j = j + 1
      end do
      if (cea_lr(i+1) > cea_lr(i)) then
        fi = (lr - cea_lr(i)) / (cea_lr(i+1) - cea_lr(i))
      else
        fi = 0.0_DP
      end if
      if (cea_le(j+1) > cea_le(j)) then
        fj = (le - cea_le(j)) / (cea_le(j+1) - cea_le(j))
      else
        fj = 0.0_DP
      end if
      fi = min(max(fi, 0.0_DP), 1.0_DP)
      fj = min(max(fj, 0.0_DP), 1.0_DP)
      t = cea_t(i,j)   * (1-fi) * (1-fj) &
        + cea_t(i+1,j) * fi     * (1-fj) &
        + cea_t(i,j+1) * (1-fi) * fj     &
        + cea_t(i+1,j+1) * fi   * fj
      p = cea_p(i,j)   * (1-fi) * (1-fj) &
        + cea_p(i+1,j) * fi     * (1-fj) &
        + cea_p(i,j+1) * (1-fi) * fj     &
        + cea_p(i+1,j+1) * fi   * fj
      tout(k) = t
      pout(k) = p
    end do
  end subroutine


  ! single-point cea lookup (p only)
  subroutine cea_pt_at(rho, e, p)
    real(DP), intent(in) :: rho, e
    real(DP), intent(out) :: p
    real(DP) :: lr, le, fi, fj
    integer :: i, j
    lr = log10(max(rho, 1.0e-12_DP))
    lr = min(max(lr, cea_lr(1)), cea_lr(cea_nr))
    le = log10(max(e, 1.0e-3_DP))
    le = min(max(le, cea_le(1)), cea_le(cea_ne))
    i = 1
    do while (i + 1 < cea_nr .and. cea_lr(i+1) < lr)
      i = i + 1
    end do
    j = 1
    do while (j + 1 < cea_ne .and. cea_le(j+1) < le)
      j = j + 1
    end do
    if (cea_lr(i+1) > cea_lr(i)) then
      fi = (lr - cea_lr(i)) / (cea_lr(i+1) - cea_lr(i))
    else
      fi = 0.0_DP
    end if
    if (cea_le(j+1) > cea_le(j)) then
      fj = (le - cea_le(j)) / (cea_le(j+1) - cea_le(j))
    else
      fj = 0.0_DP
    end if
    fi = min(max(fi, 0.0_DP), 1.0_DP)
    fj = min(max(fj, 0.0_DP), 1.0_DP)
    p = cea_p(i,j) * (1-fi) * (1-fj) + cea_p(i+1,j) * fi * (1-fj) &
      + cea_p(i,j+1) * (1-fi) * fj + cea_p(i+1,j+1) * fi * fj
  end subroutine


  ! tgas1 gamm + partials at (y, z) = (log10 rho/rho0, log10 e/e0):
  ! a faithful port of rust eqair::gamm_and_partials (12 blocks, grabau
  ! sigmoids, cold-limit branches).
  subroutine gamm_and_partials_tgas1(y, z, gamm, gr, ge)
    real(DP), intent(in) :: y, z
    real(DP), intent(out) :: gamm, gr, ge
    integer :: region, blk
    real(DP) :: g1, g2, g3, g4
    real(DP) :: g5, g6, g7, g8, num, numr, nume, s, sr, se, e_s, deno
    if (y <= -4.5_DP) then
      region = 1
      if (z <= z_breaks(1,1)) then
        blk = -1
      else if (z <= z_breaks(1,2)) then
        blk = 1
      else if (z <= z_breaks(1,3)) then
        blk = 2
      else if (z <= z_breaks(1,4)) then
        blk = 3
      else if (z <= z_breaks(1,5)) then
        blk = 4
      else
        blk = 5
      end if
    else if (y <= -0.5_DP) then
      region = 2
      if (z <= z_breaks(2,1)) then
        blk = -1
      else if (z <= z_breaks(2,2)) then
        blk = 6
      else if (z <= z_breaks(2,3)) then
        blk = 7
      else if (z <= z_breaks(2,4)) then
        blk = 8
      else
        blk = 9
      end if
    else
      region = 3
      if (z <= z_breaks(3,1)) then
        blk = -1
      else if (z <= z_breaks(3,2)) then
        blk = 10
      else if (z <= z_breaks(3,3)) then
        blk = 11
      else
        blk = 12
      end if
    end if
    if (blk < 0) then
      gamm = cold_gamm(region)
      gr = 0.0_DP
      ge = 0.0_DP
      return
    end if
    g1 = blk_g1(blk,1) + blk_g1(blk,2) * y
    g2 = (blk_g2(blk,1) + blk_g2(blk,2) * y) * z
    g3 = (blk_g3(blk,1) + blk_g3(blk,2) * z + blk_g3(blk,3) * y) * y * y
    g4 = (blk_g4(blk,1) + blk_g4(blk,2) * y + blk_g4(blk,3) * z) * z * z
    gamm = g1 + g2 + g3 + g4
    gr = blk_g1(blk,2) + blk_g2(blk,2) * z + blk_g3(blk,3) * y * y &
         + 2.0_DP * y * (blk_g3(blk,1) + blk_g3(blk,2) * z) &
         + blk_g4(blk,2) * z * z
    ge = blk_g2(blk,1) + blk_g2(blk,2) * y + blk_g3(blk,2) * y * y &
         + 2.0_DP * z * (blk_g4(blk,1) + blk_g4(blk,2) * y) &
         + 3.0_DP * blk_g4(blk,3) * z * z
    if (nbk(blk) > 0) then
      g5 = blk_g5(blk,1) + blk_g5(blk,2) * y
      g6 = (blk_g6(blk,1) + blk_g6(blk,2) * y) * z
      g7 = (blk_g7(blk,1) + blk_g7(blk,2) * z + blk_g7(blk,3) * y) * y * y
      g8 = (blk_g8(blk,1) + blk_g8(blk,2) * y + blk_g8(blk,3) * z) * z * z
      num = g5 + g6 + g7 + g8
      numr = blk_g5(blk,2) + blk_g6(blk,2) * z + blk_g7(blk,3) * y * y &
             + 2.0_DP * y * (blk_g7(blk,1) + blk_g7(blk,2) * z) &
             + blk_g8(blk,2) * z * z
      nume = blk_g6(blk,1) + blk_g6(blk,2) * y + blk_g7(blk,2) * y * y &
             + 2.0_DP * z * (blk_g8(blk,1) + blk_g8(blk,2) * y) &
             + 3.0_DP * blk_g8(blk,3) * z * z
      s = blk_s(blk,1) + blk_s(blk,2) * y + blk_s(blk,3) * z + blk_s(blk,4) * y * z
      sr = blk_s(blk,2) + blk_s(blk,4) * z
      se = blk_s(blk,3) + blk_s(blk,4) * y
      if (blk_minus(blk)) then
        e_s = exp(s)
        deno = 1.0_DP - e_s
        gamm = gamm + num / deno
        gr = gr + (numr + num * e_s * sr / deno) / deno
        ge = ge + (nume + num * e_s * se / deno) / deno
      else
        s = min(max(s, -30.0_DP), 30.0_DP)
        e_s = exp(s)
        deno = 1.0_DP + e_s
        gamm = gamm + num / deno
        gr = gr + (numr - num * e_s * sr / deno) / deno
        ge = ge + (nume - num * e_s * se / deno) / deno
      end if
    end if
  end subroutine


  ! tgas1 face pressure from (rho, e) — mirrors rust eqair_pressure_at
  function tgas1_pressure_at(rho_, e_) result(p)
    real(DP), intent(in) :: rho_, e_
    real(DP) :: p
    real(DP) :: y, z, gamm, gr, ge
    if (rho_ <= 0.0_DP .or. e_ <= 0.0_DP) then
      p = -1.0_DP
      return
    end if
    y = log10(rho_ / rho0)
    z = log10(e_ / e0)
    call gamm_and_partials_tgas1(y, z, gamm, gr, ge)
    p = (gamm - 1.0_DP) * e_ * rho_
  end function


  ! tgas1 face sound speed from (rho, e) — mirrors rust eqair_sound_at
  function tgas1_sound_at(rho_, e_) result(a)
    real(DP), intent(in) :: rho_, e_
    real(DP) :: a
    real(DP) :: y, z, gamm, gr, ge, asq
    ! mirror rust eqair_sound_at: invalid states return NaN and the
    ! call site floors it to 1e-12 exactly like rust's f64::max does.
    if (rho_ <= 0.0_DP .or. e_ <= 0.0_DP) then
      a = ieee_value(a, ieee_quiet_nan)
      return
    end if
    y = log10(rho_ / rho0)
    z = log10(e_ / e0)
    call gamm_and_partials_tgas1(y, z, gamm, gr, ge)
    asq = e_ * ((gamm - 1.0_DP) * (gamm + ge / ln10) + gr / ln10)
    if (asq > 0.0_DP) then
      a = sqrt(asq)
    else
      a = ieee_value(a, ieee_quiet_nan)
    end if
  end function


  ! ship the z-breaks + cold gamma + normalization (rust holds them)
  subroutine set_tgas1_meta(zb, cg, rho0_, e0_)
    real(DP), intent(in) :: zb(3,5), cg(3), rho0_, e0_
    z_breaks = zb
    cold_gamm = cg
    rho0 = rho0_
    e0 = e0_
    ln10 = log(10.0_DP)
    tgas1_ready = .true.
  end subroutine


  ! temperature at (rho, internal energy) from the cea table: the same
  ! bilinear read cea_pt_at does for pressure, on the temperature array.
  subroutine cea_t_at(rho, e, t)
    real(DP), intent(in) :: rho, e
    real(DP), intent(out) :: t
    real(DP) :: lr, le, fi, fj
    integer :: i, j
    lr = log10(max(rho, 1.0e-12_DP))
    lr = min(max(lr, cea_lr(1)), cea_lr(cea_nr))
    le = log10(max(e, 1.0e-3_DP))
    le = min(max(le, cea_le(1)), cea_le(cea_ne))
    i = 1
    do while (i + 1 < cea_nr .and. cea_lr(i+1) < lr)
      i = i + 1
    end do
    j = 1
    do while (j + 1 < cea_ne .and. cea_le(j+1) < le)
      j = j + 1
    end do
    if (cea_lr(i+1) > cea_lr(i)) then
      fi = (lr - cea_lr(i)) / (cea_lr(i+1) - cea_lr(i))
    else
      fi = 0.0_DP
    end if
    if (cea_le(j+1) > cea_le(j)) then
      fj = (le - cea_le(j)) / (cea_le(j+1) - cea_le(j))
    else
      fj = 0.0_DP
    end if
    fi = min(max(fi, 0.0_DP), 1.0_DP)
    fj = min(max(fj, 0.0_DP), 1.0_DP)
    t = cea_t(i,j) * (1-fi) * (1-fj) + cea_t(i+1,j) * fi * (1-fj) &
      + cea_t(i,j+1) * (1-fi) * fj + cea_t(i+1,j+1) * fi * fj
  end subroutine

end module nufoeqair

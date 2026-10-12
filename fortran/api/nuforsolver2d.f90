! nuforsolver2d.f90 — the mpi-parallel 2d euler march: rk2 heun over 1d x-slabs
! with halo exchange, muscl reconstruction, hllc fluxes, and the whole-array
! eqair closure kernels (bilinear cea table + tgas1 face evaluation).
!
! discipline notes:
! - state layout matches rust exactly: flat row-major (j*nx + i), one array
!   per conserved variable, ghost columns bookended by the halo exchange.
! - the np=1 path runs this same code (one rank, no exchange) and is the
!   bit-identity anchor for the np>1 gates.
! - the tgas1 block coefficients and cea table are shipped once at init from
!   rust (which owns the parsed data), not re-read from files here.
module nuforsolver2d
  use iso_fortran_env, only: real64
  use nufor_mpi, only: decomp_t, make_decomp, halo_exchange_x, global_max, is_rank0
  implicit none

  integer, parameter :: DP = real64

  ! cea bilinear table (t, p per node), shipped from rust at init.
  integer :: cea_nr = 0, cea_ne = 0
  real(DP), allocatable :: cea_lr(:), cea_le(:)   ! log10 rho / e axes
  real(DP), allocatable :: cea_t(:,:), cea_p(:,:) ! node values (nr, ne)

  ! tgas1 face-fit coefficients for the wake: shipped from rust at init.
  ! layout matches eqair.rs blocks: g1..g4 as (a,b,c) polynomial parts plus
  ! the grabau sigmoid set, 12 blocks, minus-flag for the sigmoid denominator.
  integer, parameter :: NB = 12
  integer :: nblocks = NB
  real(DP) :: blk_g1(NB,2), blk_g2(NB,3), blk_g3(NB,3), blk_g4(NB,3)
  real(DP) :: blk_g5(NB,2), blk_g6(NB,3), blk_g7(NB,3), blk_g8(NB,3)
  real(DP) :: blk_s(NB,4)
  logical :: blk_minus(NB)
  integer :: nbk(NB)          ! how many grabau terms each block carries
  real(DP) :: rho0, e0, ln10  ! normalization constants
  real(DP) :: z_breaks(3,5)   ! region z boundaries (rust Z_BREAKS)
  real(DP) :: cold_gamm(3)    ! per-region cold-limit gamma
  logical :: tgas1_ready = .false.

  ! per-row face work arrays (muscl states): (var, face) for the current row,
  ! plus per-face energy (conserved, from the closure-variable reconstruction
  ! used by the rust path) and sound speed.
  integer, parameter :: MAXN = 4096
  real(DP) :: face_l(5, MAXN), face_r(5, MAXN)
  real(DP) :: face_e_l(MAXN), face_e_r(MAXN), face_a_l(MAXN), face_a_r(MAXN)
  real(DP) :: face_p_l(MAXN), face_p_r(MAXN)
  real(DP) :: slope(MAXN), slope_y(MAXN)

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

  ! --- whole-array eqair kernels -------------------------------------------

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

  ! --- the march machinery --------------------------------------------------

  ! muscl strip face states for one row: van leer limiter, same as
  ! rust face_states (padded strip of n+2 -> n+1 left/right pairs).
  subroutine face_states_row(padded, n, left, right)
    integer, intent(in) :: n
    real(DP), intent(in) :: padded(n+2)
    real(DP), intent(out) :: left(n+1), right(n+1)
    real(DP) :: dm, dp_
    integer :: i
    do i = 1, n+1
      dm = padded(i+1) - padded(i)
      dp_ = padded(i+2) - padded(i+1)
      ! van Leer: (dm*dp_ > 0) ? 2*dm*dp_/(dm+dp_) : 0
      if (dm * dp_ > 0.0_DP) then
        left(i) = padded(i+1) + 0.5_DP * (2.0_DP * dm * dp_ / (dm + dp_))
      else
        left(i) = padded(i+1)
      end if
      right(i) = left(i)
    end do
  end subroutine

  ! hllc flux across one face, x-normal: a faithful port of rust
  ! hllc2d::hllc_flux (ql/qr pressure sharpening, toro star states via
  ! rankine-hugoniot, contact-wave selection). primitives in, flux out.
  subroutine hllc_face_x(rl, ul, vl, pl, el_t, al, rr, ur, vr, pr, er_t, ar, gam, &
                         fmass, fmx, fmy, fe)
    real(DP), intent(in) :: rl, ul, vl, pl, el_t, al, rr, ur, vr, pr, er_t, ar, gam
    real(DP), intent(out) :: fmass, fmx, fmy, fe
    real(DP) :: rl_, pl_, al_, rr_, pr_, ar_
    real(DP) :: gl, gr, rbar, abar, pstar, ql, qr, sl, sr, denom, sstar
    real(DP) :: etl, etr, elt, ert
    real(DP) :: fl(4), fr(4), usl(4), usr(4), fsl(4), fsr(4), f(4)
    rl_ = max(rl, 1.0e-12_DP); pl_ = max(pl, 1.0e-12_DP); al_ = max(al, 1.0e-12_DP)
    rr_ = max(rr, 1.0e-12_DP); pr_ = max(pr, 1.0e-12_DP); ar_ = max(ar, 1.0e-12_DP)
    if (pl_ > 1.0e-9_DP .and. rl_ > 1.0e-9_DP) then
      gl = al_ * al_ * rl_ / pl_
    else
      gl = gam
    end if
    if (pr_ > 1.0e-9_DP .and. rr_ > 1.0e-9_DP) then
      gr = ar_ * ar_ * rr_ / pr_
    else
      gr = gam
    end if
    rbar = 0.5_DP * (rl_ + rr_)
    abar = 0.5_DP * (al_ + ar_)
    pstar = max(0.5_DP * (pl_ + pr_) - 0.5_DP * (ur - ul) * rbar * abar, 0.0_DP)
    if (pstar <= pl_) then
      ql = 1.0_DP
    else
      ql = sqrt(1.0_DP + (gl + 1.0_DP) / (2.0_DP * gl) * (pstar / pl_ - 1.0_DP))
    end if
    if (pstar <= pr_) then
      qr = 1.0_DP
    else
      qr = sqrt(1.0_DP + (gr + 1.0_DP) / (2.0_DP * gr) * (pstar / pr_ - 1.0_DP))
    end if
    sl = min(ul - al_ * ql, ur - ar_ * qr)
    sr = max(ul + al_ * ql, ur + ar_ * qr)
    denom = rl_ * (sl - ul) - rr_ * (sr - ur)
    sstar = (pr_ - pl_ + rl_ * ul * (sl - ul) - rr_ * ur * (sr - ur)) / denom
    ! total energies (per volume) and physical fluxes
    elt = rl_ * (el_t + 0.5_DP * (ul**2 + vl**2))
    ert = rr_ * (er_t + 0.5_DP * (ur**2 + vr**2))
    fl = [rl_ * ul, rl_ * ul**2 + pl_, rl_ * ul * vl, ul * (elt + pl_)]
    fr = [rr_ * ur, rr_ * ur**2 + pr_, rr_ * ur * vr, ur * (ert + pr_)]
    ! star states and star fluxes: f* = f + s (u* - u)
    usl = star_state(rl_, ul, vl, pl_, elt, sl, sstar)
    usr = star_state(rr_, ur, vr, pr_, ert, sr, sstar)
    fsl = fl + sl * (usl - [rl_, rl_ * ul, rl_ * vl, elt])
    fsr = fr + sr * (usr - [rr_, rr_ * ur, rr_ * vr, ert])
    if (sl > 0.0_DP) then
      f = fl
    else if (sr < 0.0_DP) then
      f = fr
    else if (sstar > 0.0_DP) then
      f = fsl
    else
      f = fsr
    end if
    fmass = f(1); fmx = f(2); fmy = f(3); fe = f(4)
  end subroutine

  ! toro star conserved state for one side
  pure function star_state(r, un, ut, p, e_tot, s_side, sstar) result(u)
    real(DP), intent(in) :: r, un, ut, p, e_tot, s_side, sstar
    real(DP) :: u(4), fac
    fac = r * (s_side - un) / (s_side - sstar)
    u = [fac, fac * sstar, fac * ut, &
         fac * (e_tot / r + (sstar - un) * (sstar + p / (r * (s_side - un))))]
  end function

  ! --- the x-sweep with muscl + hllc + bc ------------------------------------

  ! apply boundary ghost cells for one padded strip at the domain edges.
  ! bc codes: 0=transmissive, 1=supersonic_inflow(rho,u,v,p), 2=supersonic_outflow,
  ! 3=slip_wall, 4=no_slip_wall. inflow carries (rho,u,v,p) in bc_vals.
  ! vars: 1=rho, 2=u, 3=v, 4=p (padded strips are primitive).
  subroutine bc_ghost(var, bc_code, bc_vals, n, padded, side)
    integer, intent(in) :: var, bc_code, n, side  ! side: 0=west, 1=east
    real(DP), intent(in) :: bc_vals(5)
    real(DP), intent(inout) :: padded(n+2)
    ! fills both ghost layers at the domain edge (the muscl limiter at the
    ! first interior face reads two cells beyond it)
    if (side == 0) then
      select case (bc_code)
      case (0)   ! transmissive: mirror the interior into both ghosts
        padded(1) = padded(3)
        padded(2) = padded(3)
      case (1)   ! supersonic inflow: fixed freestream primitive
        padded(1) = bc_vals(var)
        padded(2) = bc_vals(var)
      case (2)   ! supersonic outflow
        padded(1) = padded(3)
        padded(2) = padded(3)
      case (3)   ! slip wall: mirror normal velocity about the wall
        if (var == 2) then
          padded(1) = -padded(3)
          padded(2) = -padded(3)
        else
          padded(1) = padded(3)
          padded(2) = padded(3)
        end if
      case (4)   ! no-slip wall: both velocities reversed
        if (var == 2 .or. var == 3) then
          padded(1) = -padded(3)
          padded(2) = -padded(3)
        else
          padded(1) = padded(3)
          padded(2) = padded(3)
        end if
      end select
    else
      select case (bc_code)
      case (0)
        padded(n+1) = padded(n)
        padded(n+2) = padded(n)
      case (1)
        padded(n+1) = bc_vals(var)
        padded(n+2) = bc_vals(var)
      case (2)
        padded(n+1) = padded(n)
        padded(n+2) = padded(n)
      case (3)
        if (var == 2) then
          padded(n+1) = -padded(n)
          padded(n+2) = -padded(n)
        else
          padded(n+1) = padded(n)
          padded(n+2) = padded(n)
        end if
      case (4)
        if (var == 2 .or. var == 3) then
          padded(n+1) = -padded(n)
          padded(n+2) = -padded(n)
        else
          padded(n+1) = padded(n)
          padded(n+2) = padded(n)
        end if
      end select
    end if
  end subroutine

  ! --- the full euler rk2 step over this rank's slab -------------------------
  ! slab layout (nx_l+4, ny): ghost cols 1,2 (west) and nx_l+3,4 (east);
  ! interior cols 3..nx_l+2. at np=1 the ghosts are unused (bc fills the
  ! strip edges directly).

  subroutine euler2d_step_rk2(rho_in, mx_in, my_in, e_in, nx_l, ny_, &
                              dx, dy, gamma, cfl, dt_cap, d, &
                              bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                              rho_out, mx_out, my_out, e_out, dt_out, steps_flag)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(in) :: rho_in(nx_l+4, ny_), mx_in(nx_l+4, ny_), &
                            my_in(nx_l+4, ny_), e_in(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, gamma, cfl, dt_cap
    type(decomp_t), intent(in) :: d
    integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(DP), intent(in) :: bc_vals(5, 4)
    logical, intent(in) :: solid(nx_l+4, ny_)
    real(DP), intent(out) :: rho_out(nx_l+4, ny_), mx_out(nx_l+4, ny_), &
                             my_out(nx_l+4, ny_), e_out(nx_l+4, ny_)
    real(DP), intent(out) :: dt_out
    integer, intent(in) :: steps_flag
    real(DP), allocatable :: u(:,:), v(:,:), p(:,:), a(:,:), ein(:,:)
    real(DP), allocatable :: fx(:,:,:), fy(:,:,:)
    real(DP) :: smax, s_local, dt, ge_ff
    real(DP) :: fm, fmx2, fmy2, fe, dm, dp_
    real(DP) :: pad_y(MAXN)
    integer :: i, j, f, var, ierr, nzf

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), p(nx_l+4, ny_), a(nx_l+4, ny_))
    allocate(ein(nx_l+4, ny_))
    allocate(fx(4, nx_l+2, ny_), fy(4, nx_l, ny_+1))

    ! primitive conversion over the whole padded slab; ein is the internal
    ! energy per unit mass (the eqair closure variable the rust march
    ! reconstructs on faces)
    do j = 1, ny_
      do i = 1, nx_l+4
        u(i,j) = mx_in(i,j) / max(rho_in(i,j), 1.0e-12_DP)
        v(i,j) = my_in(i,j) / max(rho_in(i,j), 1.0e-12_DP)
        ein(i,j) = e_in(i,j) / max(rho_in(i,j), 1.0e-12_DP) &
                 - 0.5_DP * (u(i,j)**2 + v(i,j)**2)
        call cea_pt_at(rho_in(i,j), ein(i,j), p(i,j))
        p(i,j) = max(p(i,j), 1.0e-12_DP)
        ge_ff = 1.0_DP + p(i,j) / (rho_in(i,j) * max(ein(i,j), 1.0e-3_DP))
        ge_ff = min(max(ge_ff, 1.05_DP), 1.67_DP)
        a(i,j) = sqrt(ge_ff * p(i,j) / max(rho_in(i,j), 1.0e-12_DP))
        a(i,j) = max(a(i,j), 50.0_DP)
      end do
    end do

    ! ---- x sweep per row: muscl over the padded strip, faces 2..nx_l+2 ----
    do j = 1, ny_
      call x_row_fluxes(j, nx_l, ny_, gamma, rho_in, ein, u, v, p, a, &
                        d, bc_w, bc_e, bc_vals, fx(1:4, 1:nx_l+2, j))
    end do

    ! ---- y sweep per column (y undecomposed in v1) ----
    do i = 3, nx_l+2
      do var = 1, 5
        do j = 1, ny_
          select case (var)
          case (1); pad_y(j+1) = rho_in(i,j)
          case (2); pad_y(j+1) = u(i,j)
          case (3); pad_y(j+1) = v(i,j)
          case (4); pad_y(j+1) = p(i,j)
          case (5); pad_y(j+1) = ein(i,j)
          end select
        end do
        ! domain-edge bcs (every rank holds the full y extent); the y strip
        ! carries ONE ghost slot per side (pad_y(1) south, pad_y(ny_+2)
        ! north), so fill them directly instead of the two-layer bc_ghost,
        ! which would clobber the first and last interior cells. var 5
        ! (ein) mirrors like a scalar in every bc case.
        select case (bc_s)
        case (1)   ! supersonic inflow: freestream primitive (ein preinverted)
          if (var == 5) then
            pad_y(1) = bc_vals(5,3)
          else
            pad_y(1) = bc_vals(var,3)
          end if
        case (3)   ! slip wall: reflect the y-normal velocity v
          if (var == 3) then
            pad_y(1) = -pad_y(2)
          else
            pad_y(1) = pad_y(2)
          end if
        case (4)   ! no-slip wall: both velocities reversed
          if (var == 2 .or. var == 3) then
            pad_y(1) = -pad_y(2)
          else
            pad_y(1) = pad_y(2)
          end if
        case default  ! transmissive / supersonic outflow: mirror the edge cell
          pad_y(1) = pad_y(2)
        end select
        select case (bc_n)
        case (1)
          if (var == 5) then
            pad_y(ny_+2) = bc_vals(5,4)
          else
            pad_y(ny_+2) = bc_vals(var,4)
          end if
        case (3)
          if (var == 3) then
            pad_y(ny_+2) = -pad_y(ny_+1)
          else
            pad_y(ny_+2) = pad_y(ny_+1)
          end if
        case (4)
          if (var == 2 .or. var == 3) then
            pad_y(ny_+2) = -pad_y(ny_+1)
          else
            pad_y(ny_+2) = pad_y(ny_+1)
          end if
        case default
          pad_y(ny_+2) = pad_y(ny_+1)
        end select
        ! van leer slopes: slope_y(f) is centered on cell f (slots f..f+2).
        ! face f separates cells (f-1) | f, mirroring the x-sweep's
        ! convention (face fx(i-1) is cell i's +x face): face 1 is the
        ! south domain-edge face (ghost | cell 1), face ny_+1 the north
        ! edge (cell ny_ | ghost).
        do f = 1, ny_
          dm = pad_y(f+1) - pad_y(f)
          dp_ = pad_y(f+2) - pad_y(f+1)
          if (dm * dp_ > 0.0_DP) then
            slope_y(f) = 2.0_DP * dm * dp_ / (dm + dp_)
          else
            slope_y(f) = 0.0_DP
          end if
        end do
        face_l(var, 1) = pad_y(1)
        face_r(var, 1) = pad_y(2) - 0.5_DP * slope_y(1)
        do f = 2, ny_+1
          face_l(var, f) = pad_y(f) + 0.5_DP * slope_y(f-1)
          face_r(var, f) = pad_y(f+1) - 0.5_DP * slope_y(f)
        end do
        face_r(var, ny_+1) = pad_y(ny_+2)
      end do
      do f = 1, ny_+1
        face_e_l(f) = face_l(5,f)
        face_e_r(f) = face_r(5,f)
        face_p_l(f) = tgas1_pressure_at(face_l(1,f), face_e_l(f))
        if (face_p_l(f) < 0.0_DP) face_p_l(f) = 1.0e-12_DP
        face_p_r(f) = tgas1_pressure_at(face_r(1,f), face_e_r(f))
        if (face_p_r(f) < 0.0_DP) face_p_r(f) = 1.0e-12_DP
        face_a_l(f) = tgas1_sound_at(face_l(1,f), face_e_l(f))
        face_a_r(f) = tgas1_sound_at(face_r(1,f), face_e_r(f))
      end do
      do f = 1, ny_+1
        ! axis-1: pass normal momentum as arg2 (u-slot = v), tangential as arg3
        call hllc_face_x(face_l(1,f), face_l(3,f), face_l(2,f), face_p_l(f), &
                         face_e_l(f), face_a_l(f), &
                         face_r(1,f), face_r(3,f), face_r(2,f), face_p_r(f), &
                         face_e_r(f), face_a_r(f), gamma, fm, fmy2, fmx2, fe)
        fy(1, i-2, f) = fm; fy(2, i-2, f) = fmx2; fy(3, i-2, f) = fmy2; fy(4, i-2, f) = fe
      end do
    end do

    ! ---- global cfl: max(|u|,|v|) + a with the cea gamma_eff a field ----
    ! slip-wall enforcement: zero the flux across any face touching a
    ! solid cell, matching the rust march's wall closure (the riemann
    ! solve never runs across the wall). face fx(:,f,j) sits between
    ! padded cells f+1 and f+2; face fy(:,f,j) separates padded y-cells
    ! j and j+1 at x-slot f+2 (face 1 = ghost | cell 1).
    nzf = 0
    do j = 1, ny_
      do f = 1, nx_l+2
        if (solid(f+1,j) .or. solid(f+2,j)) then
          fx(1,f,j) = 0.0_DP
          fx(2,f,j) = 0.0_DP
          fx(3,f,j) = 0.0_DP
          fx(4,f,j) = 0.0_DP
          nzf = nzf + 1
        end if
      end do
    end do
    do j = 1, ny_+1
      do f = 1, nx_l
        ! face j separates y-cells (j-1) | j (1-based); face 1 has a
        ! ghost below and face ny_+1 a ghost above, so only test the
        ! cell indices that exist (no out-of-bounds solid reads).
        if ((j > 1 .and. solid(f+2,j-1)) .or. (j <= ny_ .and. solid(f+2,j))) then
          fy(1,f,j) = 0.0_DP
          fy(2,f,j) = 0.0_DP
          fy(3,f,j) = 0.0_DP
          fy(4,f,j) = 0.0_DP
          nzf = nzf + 1
        end if
      end do
    end do
    s_local = 0.0_DP
    do j = 1, ny_
      do i = 3, nx_l+2
        s_local = max(s_local, max(abs(u(i,j)), abs(v(i,j))) + a(i,j))
      end do
    end do
    call global_max(s_local, smax, d)
    dt = min(cfl * min(dx, dy) / max(smax, 1.0e-12_DP), dt_cap)
    dt_out = dt

    ! ---- conservative update on interior cells (3..nx_l+2) ----
    do j = 1, ny_
      do i = 3, nx_l+2
        if (solid(i,j)) cycle
        rho_out(i,j) = rho_in(i,j) - dt * ( &
              (fx(1, i-1, j) - fx(1, i-2, j)) / dx &
            + (fy(1, i-2, j+1) - fy(1, i-2, j)) / dy )
        mx_out(i,j) = mx_in(i,j) - dt * ( &
              (fx(2, i-1, j) - fx(2, i-2, j)) / dx &
            + (fy(2, i-2, j+1) - fy(2, i-2, j)) / dy )
        my_out(i,j) = my_in(i,j) - dt * ( &
              (fx(3, i-1, j) - fx(3, i-2, j)) / dx &
            + (fy(3, i-2, j+1) - fy(3, i-2, j)) / dy )
        e_out(i,j) = e_in(i,j) - dt * ( &
              (fx(4, i-1, j) - fx(4, i-2, j)) / dx &
            + (fy(4, i-2, j+1) - fy(4, i-2, j)) / dy )
      end do
    end do
    ! solid cells + ghosts copy through
    do j = 1, ny_
      do i = 1, nx_l+4
        if (solid(i,j) .or. i < 3 .or. i > nx_l+2) then
          rho_out(i,j) = rho_in(i,j)
          mx_out(i,j) = mx_in(i,j)
          my_out(i,j) = my_in(i,j)
          e_out(i,j) = e_in(i,j)
        end if
      end do
    end do

    deallocate(u, v, p, a, ein, fx, fy)
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

  ! the whole x-flux sweep for one row: muscl over the padded strip with the
  ! 2-ghost layout, tgas1 face p/a, hllc per face. outputs fx(4, nx_l+3)
  ! where fx(:,k,j) is the face between local cols k+1 and k+2 (i.e. global
  ! face k+1 of the slab interior).
  subroutine x_row_fluxes(j, nx_l, ny_, gamma, rho_in, ein, u, v, p, a, &
                          d, bc_w, bc_e, bc_vals, fx)
    integer, intent(in) :: j, nx_l, ny_
    real(DP), intent(in) :: gamma
    real(DP), intent(in) :: rho_in(nx_l+4, ny_), ein(nx_l+4, ny_)
    real(DP), intent(in) :: u(nx_l+4, ny_), v(nx_l+4, ny_), &
                            p(nx_l+4, ny_), a(nx_l+4, ny_)
    type(decomp_t), intent(in) :: d
    integer, intent(in) :: bc_w, bc_e
    real(DP), intent(in) :: bc_vals(5, 4)
    real(DP), intent(out) :: fx(4, nx_l+2)
    real(DP) :: pad(MAXN)
    real(DP) :: eip(MAXN)
    real(DP) :: dm, dp_, fm, fmx, fmy, fe
    integer :: i, f, var

    ! van leer slopes over the padded strip (d(f) for the cell at pad(f+1)),
    ! then rust face_states: fl(f)=pad(f+1)+0.5*d(f-1) from the left cell,
    ! fr(f)=pad(f+2)-0.5*d(f) from the right cell, boundary faces take the
    ! ghost/extrapolated endpoint directly.
    do var = 1, 4
      do i = 1, nx_l+4
        select case (var)
        case (1); pad(i) = rho_in(i,j)
        case (2); pad(i) = u(i,j)
        case (3); pad(i) = v(i,j)
        case (4); pad(i) = p(i,j)
        end select
      end do
      if (d%west_domain_edge) then
        call bc_ghost(var, bc_w, bc_vals(1,1), nx_l+2, pad, 0)
      end if
      if (d%east_domain_edge) then
        call bc_ghost(var, bc_e, bc_vals(1,2), nx_l+2, pad, 1)
      end if
      do i = 1, nx_l+2
        dm = pad(i+1) - pad(i)
        dp_ = pad(i+2) - pad(i+1)
        if (dm * dp_ > 0.0_DP) then
          slope(i) = 2.0_DP * dm * dp_ / (dm + dp_)
        else
          slope(i) = 0.0_DP
        end if
      end do
      ! faces 1..nx_l+1: face f separates cells pad(f+1) | pad(f+2). the
      ! DOMAIN boundary faces take the ghost state directly (rust fl[0]/
      ! fr[n]); interior slab faces — including the subdomain interface
      ! faces at np>1 — reconstruct with the van leer slopes on both sides.
      if (d%west_domain_edge) then
        face_l(var, 1) = pad(2)
      else
        face_l(var, 1) = pad(2) + 0.5_DP * slope(1)
      end if
      face_r(var, 1) = pad(3) - 0.5_DP * slope(2)
      do f = 2, nx_l
        face_l(var, f) = pad(f+1) + 0.5_DP * slope(f)
        face_r(var, f) = pad(f+2) - 0.5_DP * slope(f+1)
      end do
      face_l(var, nx_l+1) = pad(nx_l+2) + 0.5_DP * slope(nx_l+1)
      if (d%east_domain_edge) then
        face_r(var, nx_l+1) = pad(nx_l+3)
      else
        face_r(var, nx_l+1) = pad(nx_l+3) - 0.5_DP * slope(nx_l+2)
      end if
    end do

    ! muscl the closure variable ein with the same fl/fr scheme. inflow
    ! domain edges take the closure-inverted inflow ein (bc_vals(5, side))
    ! in both ghost layers, mirroring the prim-var ghost seeding.
    do i = 1, nx_l+4
      eip(i) = ein(i,j)
    end do
    if (d%west_domain_edge .and. bc_w == 1) then
      eip(1) = bc_vals(5,1)
      eip(2) = bc_vals(5,1)
    end if
    if (d%east_domain_edge .and. bc_e == 1) then
      eip(nx_l+3) = bc_vals(5,2)
      eip(nx_l+4) = bc_vals(5,2)
    end if
    do i = 1, nx_l+2
      dm = eip(i+1) - eip(i)
      dp_ = eip(i+2) - eip(i+1)
      if (dm * dp_ > 0.0_DP) then
        slope(i) = 2.0_DP * dm * dp_ / (dm + dp_)
      else
        slope(i) = 0.0_DP
      end if
    end do
    if (d%west_domain_edge) then
      face_e_l(1) = eip(2)
    else
      face_e_l(1) = eip(2) + 0.5_DP * slope(1)
    end if
    face_e_r(1) = eip(3) - 0.5_DP * slope(2)
    do f = 2, nx_l
      face_e_l(f) = eip(f+1) + 0.5_DP * slope(f)
      face_e_r(f) = eip(f+2) - 0.5_DP * slope(f+1)
    end do
    face_e_l(nx_l+1) = eip(nx_l+2) + 0.5_DP * slope(nx_l+1)
    if (d%east_domain_edge) then
      face_e_r(nx_l+1) = eip(nx_l+3)
    else
      face_e_r(nx_l+1) = eip(nx_l+3) - 0.5_DP * slope(nx_l+2)
    end if

    ! hllc per face: tgas1 face p/a from the reconstructed closure
    ! states, then the riemann flux into this row's fx slice.
    do f = 1, nx_l+1
      face_p_l(f) = tgas1_pressure_at(face_l(1,f), face_e_l(f))
      if (face_p_l(f) < 0.0_DP) face_p_l(f) = 1.0e-12_DP
      face_p_r(f) = tgas1_pressure_at(face_r(1,f), face_e_r(f))
      if (face_p_r(f) < 0.0_DP) face_p_r(f) = 1.0e-12_DP
      face_a_l(f) = tgas1_sound_at(face_l(1,f), face_e_l(f))
      face_a_r(f) = tgas1_sound_at(face_r(1,f), face_e_r(f))
      call hllc_face_x(face_l(1,f), face_l(2,f), face_l(3,f), face_p_l(f), &
                       face_e_l(f), face_a_l(f), &
                       face_r(1,f), face_r(2,f), face_r(3,f), face_p_r(f), &
                       face_e_r(f), face_a_r(f), gamma, fm, fmx, fmy, fe)
      fx(1,f) = fm; fx(2,f) = fmx; fx(3,f) = fmy; fx(4,f) = fe
    end do
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
    if (rho_ <= 0.0_DP .or. e_ <= 0.0_DP) then
      a = 300.0_DP
      return
    end if
    y = log10(rho_ / rho0)
    z = log10(e_ / e0)
    call gamm_and_partials_tgas1(y, z, gamm, gr, ge)
    asq = e_ * ((gamm - 1.0_DP) * (gamm + ge / ln10) + gr / ln10)
    if (asq > 0.0_DP) then
      a = sqrt(asq)
    else
      a = 300.0_DP
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

  ! the viscous half-step on this rank's padded slab, mirroring
  ! viscous2d::add_viscous_cells_mu on a uniform grid: closure
  ! temperature from the cea table, sutherland mu per cell, central
  ! cell gradients (y mirrored at the domain rows, x from the halo
  ! which carries the pad2d copy semantics), face-averaged gradients,
  ! no-slip wall shear + isothermal-wall conduction, divergence add.
  ! call right after a halo refresh of the post-euler state.
  subroutine viscous_add_slab(rho, mx, my, e, nx_l, ny_, dx, dy, dt, &
                              mu0, t0_ref, s_par, pr, tw, &
                              bc_s, bc_n, solid)
    integer, intent(in) :: nx_l, ny_
    real(DP), intent(inout) :: rho(nx_l+4, ny_), mx(nx_l+4, ny_), &
                               my(nx_l+4, ny_), e(nx_l+4, ny_)
    real(DP), intent(in) :: dx, dy, dt, mu0, t0_ref, s_par, pr, tw
    integer, intent(in) :: bc_s, bc_n
    logical, intent(in) :: solid(nx_l+4, ny_)
    real(DP), allocatable :: u(:,:), v(:,:), tt(:,:), mu(:,:)
    real(DP), allocatable :: fxv(:,:,:), fyv(:,:,:), gu(:,:,:)
    real(DP) :: ein, kap, dtdn, q_wall
    real(DP) :: gux, guy, gvx, gvy, gtx, gty
    real(DP) :: muf, divw, txx, txy, tyy, q
    integer :: i, j, f, k, a, b

    if (mu0 <= 0.0_DP) return
    ! kappa = gamma/((gamma-1) pr)/r0 with r0 = 78408.4/273.15, the
    ! closure-path form add_viscous_cells_mu takes under a thermo model.
    kap = 1.4_DP / (0.4_DP * pr) / (78408.4_DP / 273.15_DP)

    allocate(u(nx_l+4, ny_), v(nx_l+4, ny_), tt(nx_l+4, ny_), mu(nx_l+4, ny_))
    allocate(fxv(4, 0:nx_l+2, ny_), fyv(3, 0:nx_l-1, 0:ny_+1))
    allocate(gu(6, nx_l+4, ny_))

    ! primitives + closure temperature + sutherland mu per padded cell
    do j = 1, ny_
      do i = 1, nx_l+4
        u(i,j) = mx(i,j) / max(rho(i,j), 1.0e-12_DP)
        v(i,j) = my(i,j) / max(rho(i,j), 1.0e-12_DP)
        ein = e(i,j) / max(rho(i,j), 1.0e-12_DP) - 0.5_DP * (u(i,j)**2 + v(i,j)**2)
        call cea_t_at(rho(i,j), ein, tt(i,j))
        tt(i,j) = max(tt(i,j), 1.0_DP)
        mu(i,j) = mu0 * (tt(i,j) / t0_ref)**1.5_DP * (t0_ref + s_par) / (tt(i,j) + s_par)
      end do
    end do

    ! cell-centred gradients with the wall override on no-slip sides,
    ! stored per padded cell: gu(6, i, j) = ux, uy, vx, vy, tx, ty.
    ! y rows 1 and ny_ mirror their edge row (pad2d); no-slip wall
    ! cells then override uy/vy with the quadratic-consistent stencil.
    do j = 1, ny_
      do i = 2, nx_l+3
        gu(1,i,j) = (u(i+1,j) - u(i-1,j)) / (2.0_DP * dx)
        gu(3,i,j) = (v(i+1,j) - v(i-1,j)) / (2.0_DP * dx)
        gu(5,i,j) = (tt(i+1,j) - tt(i-1,j)) / (2.0_DP * dx)
        if (j == 1) then
          ! mirrored ghost row (pad2d): (u2 - u1)/(2 dy)
          gu(2,i,j) = (u(i,2) - u(i,1)) / (2.0_DP * dy)
          gu(4,i,j) = (v(i,2) - v(i,1)) / (2.0_DP * dy)
          gu(6,i,j) = (tt(i,2) - tt(i,1)) / (2.0_DP * dy)
        else if (j == ny_) then
          gu(2,i,j) = (u(i,ny_) - u(i,ny_-1)) / (2.0_DP * dy)
          gu(4,i,j) = (v(i,ny_) - v(i,ny_-1)) / (2.0_DP * dy)
          gu(6,i,j) = (tt(i,ny_) - tt(i,ny_-1)) / (2.0_DP * dy)
        else
          gu(2,i,j) = (u(i,j+1) - u(i,j-1)) / (2.0_DP * dy)
          gu(4,i,j) = (v(i,j+1) - v(i,j-1)) / (2.0_DP * dy)
          gu(6,i,j) = (tt(i,j+1) - tt(i,j-1)) / (2.0_DP * dy)
        end if
      end do
      ! no-slip wall-cell override: uy = flip*(3 u0 + u1)/(3 dy)
      if (bc_s == 4) then
        do i = 2, nx_l+3
          gu(2,i,1) = (3.0_DP * u(i,1) + u(i,2)) / (3.0_DP * dy)
          gu(4,i,1) = (3.0_DP * v(i,1) + v(i,2)) / (3.0_DP * dy)
        end do
      end if
      if (bc_n == 4) then
        do i = 2, nx_l+3
          gu(2,i,ny_) = -(3.0_DP * u(i,ny_) + u(i,ny_-1)) / (3.0_DP * dy)
          gu(4,i,ny_) = -(3.0_DP * v(i,ny_) + v(i,ny_-1)) / (3.0_DP * dy)
        end do
      end if
    end do

    ! x-face fluxes: face k between padded cells k+1 | k+2, k=1..nx_l+1
    ! (k=0 and k=nx_l+2 are the ghost-edge faces: zero, matching the
    ! rust pass whose domain-edge faces carry zero flux).
    fxv = 0.0_DP
    do j = 1, ny_
      do k = 1, nx_l+1
        a = k+1; b = k+2
        muf = 0.5_DP * (mu(a,j) + mu(b,j))
        gux = 0.5_DP * (gu(1,a,j) + gu(1,b,j))
        guy = 0.5_DP * (gu(2,a,j) + gu(2,b,j))
        gvx = 0.5_DP * (gu(3,a,j) + gu(3,b,j))
        gvy = 0.5_DP * (gu(4,a,j) + gu(4,b,j))
        gtx = 0.5_DP * (gu(5,a,j) + gu(5,b,j))
        divw = gux + gvy
        txx = muf * (gux * 2.0_DP - divw * 2.0_DP / 3.0_DP)
        txy = muf * (guy + gvx)
        q = muf * kap * gtx
        fxv(1,k,j) = txx
        fxv(2,k,j) = txy
        fxv(4,k,j) = u(a,j) * txx + v(a,j) * txy + q
      end do
    end do

    ! y-face fluxes: face f below row f (between rows f-1 | f),
    ! f=1..ny_+1 (f=1 ghost edge, f=ny_+1 north edge); interior
    ! f=2..ny_. x slots 1..nx_l map to fyv(:, i-2, :) for cell i.
    fyv = 0.0_DP
    do i = 3, nx_l+2
      do f = 2, ny_
        muf = 0.5_DP * (mu(i,f-1) + mu(i,f))
        gux = 0.5_DP * (gu(1,i,f-1) + gu(1,i,f))
        guy = 0.5_DP * (gu(2,i,f-1) + gu(2,i,f))
        gvx = 0.5_DP * (gu(3,i,f-1) + gu(3,i,f))
        gvy = 0.5_DP * (gu(4,i,f-1) + gu(4,i,f))
        gty = 0.5_DP * (gu(6,i,f-1) + gu(6,i,f))
        divw = gux + gvy
        txy = muf * (guy + gvx)
        tyy = muf * (gvy * 2.0_DP - divw * 2.0_DP / 3.0_DP)
        q = muf * kap * gty
        fyv(1,i-2,f) = txy
        fyv(2,i-2,f) = tyy
        fyv(3,i-2,f) = u(i,f-1) * txy + v(i,f-1) * tyy + q
      end do
      ! south no-slip wall face f=1: one-sided shear (9 u0 - u1)/(3 dy)
      if (bc_s == 4) then
        muf = mu(i,1)
        txy = muf * (9.0_DP * u(i,1) - u(i,2)) / (3.0_DP * dy)
        tyy = muf * (9.0_DP * v(i,1) - v(i,2)) / (3.0_DP * dy)
        q_wall = 0.0_DP
        if (tw > 0.0_DP) then
          dtdn = (-3.0_DP * tw + 4.0_DP * tt(i,1) - tt(i,2)) / (2.0_DP * dy)
          q_wall = muf * kap * dtdn
        end if
        fyv(1,i-2,1) = txy
        fyv(2,i-2,1) = tyy
        fyv(3,i-2,1) = u(i,1) * txy + v(i,1) * tyy + q_wall
      end if
      ! north no-slip wall face f=ny_+1: outward normal -y negates.
      if (bc_n == 4) then
        muf = mu(i,ny_)
        txy = muf * (9.0_DP * u(i,ny_) - u(i,ny_-1)) / (3.0_DP * dy)
        tyy = muf * (9.0_DP * v(i,ny_) - v(i,ny_-1)) / (3.0_DP * dy)
        q_wall = 0.0_DP
        if (tw > 0.0_DP) then
          dtdn = (3.0_DP * tw - 4.0_DP * tt(i,ny_) + tt(i,ny_-1)) / (2.0_DP * dy)
          q_wall = muf * kap * dtdn
        end if
        fyv(1,i-2,ny_+1) = -txy
        fyv(2,i-2,ny_+1) = -tyy
        fyv(3,i-2,ny_+1) = -(u(i,ny_) * txy + v(i,ny_) * tyy + q_wall)
      end if
    end do

    ! divergence into interior cells: cell (i,j) adds
    ! dt/dx * (fx(:,i-1,j) - fx(:,i-2,j)) + dt/dy * (fy(:,i-2,j+1) - fy(:,i-2,j))
    do j = 1, ny_
      do i = 3, nx_l+2
        if (solid(i,j)) cycle
        mx(i,j) = mx(i,j) + dt / dx * (fxv(1,i-1,j) - fxv(1,i-2,j)) &
                            + dt / dy * (fyv(1,i-2,j+1) - fyv(1,i-2,j))
        my(i,j) = my(i,j) + dt / dx * (fxv(2,i-1,j) - fxv(2,i-2,j)) &
                            + dt / dy * (fyv(2,i-2,j+1) - fyv(2,i-2,j))
        e(i,j) = e(i,j) + dt / dx * (fxv(4,i-1,j) - fxv(4,i-2,j)) &
                       + dt / dy * (fyv(3,i-2,j+1) - fyv(3,i-2,j))
      end do
    end do

    deallocate(u, v, tt, mu, fxv, fyv, gu)
  end subroutine

end module nuforsolver2d

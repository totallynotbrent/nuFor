! nuforfluxes.f90 — muscl strip face states and the assembled per-row x-flux
! slice: van-leer slopes, domain-ghost seeding, tgas1 face p/a, hllc per face.
! the strip reuses the shared face workbank in nuforfaces, filled per row.
module nuforfluxes
  use iso_fortran_env, only: real64
  use nufor_mpi, only: decomp_t
  use nuforfaces, only: DP => DP, MAXN, face_l, face_r, face_e_l, face_e_r, face_a_l, face_a_r, face_p_l, face_p_r, slope
  use nufoeqair, only: tgas1_pressure_at, tgas1_sound_at
  use nuforriemann, only: hllc_face_x
  use nuforboundary, only: bc_ghost
  implicit none

contains


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
      if (face_p_l(f) < 0.0_DP .or. face_p_l(f) /= face_p_l(f)) face_p_l(f) = 1.0e-12_DP
      face_p_r(f) = tgas1_pressure_at(face_r(1,f), face_e_r(f))
      if (face_p_r(f) < 0.0_DP .or. face_p_r(f) /= face_p_r(f)) face_p_r(f) = 1.0e-12_DP
      face_a_l(f) = tgas1_sound_at(face_l(1,f), face_e_l(f))
      face_a_r(f) = tgas1_sound_at(face_r(1,f), face_e_r(f))
      ! rust's f64::max(NaN, 1e-12) = 1e-12; floor NaN explicitly here
      if (face_a_l(f) /= face_a_l(f)) face_a_l(f) = 1.0e-12_DP
      if (face_a_r(f) /= face_a_r(f)) face_a_r(f) = 1.0e-12_DP
      call hllc_face_x(face_l(1,f), face_l(2,f), face_l(3,f), face_p_l(f), &
                       face_e_l(f), face_a_l(f), &
                       face_r(1,f), face_r(2,f), face_r(3,f), face_p_r(f), &
                       face_e_r(f), face_a_r(f), gamma, fm, fmx, fmy, fe)
      fx(1,f) = fm; fx(2,f) = fmx; fx(3,f) = fmy; fx(4,f) = fe
    end do
end subroutine

end module nuforfluxes

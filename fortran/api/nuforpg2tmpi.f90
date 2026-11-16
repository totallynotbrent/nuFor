! nuforpg2tmpi.f90 — c-bound mpi entry for the perfect-gas two-temperature
! (park2t) lane: the perfect-gas heun march over this rank's slab, then the
! landau-teller substep on the post-march state. mirrors nuforbridge's
! nfor_mpi_step_2t but the march kernel is the perfect-gas one.
module nuforpg2tmpi
  use iso_c_binding, only: c_int, c_double
  use iso_fortran_env, only: real64
  use nufor_mpi, only: decomp_t
  use nuforhalos, only: g_decomp, halo_exchange_cons, halo_exchange_mask, halo_exchange_plane
  use nuforpg2d, only: pg2d_step_rk2
  use nufor2t, only: park2t_substep_slab
  implicit none
  integer, parameter :: RP = real64

contains

  ! one park2t mpi step: perfect-gas heun march then the 2t substep on the
  ! marched state; ev enters/exits as the rank-local vibrational plane.
  subroutine nfor_mpi_step_pg2t(rho, mx, my, e, ev, nloc, nx_l, ny_, &
                                dx, dy, gamma, cfl, dt_cap, &
                                bc_w, bc_e, bc_s, bc_n, bc_vals, solid, &
                                rho2, mx2, my2, e2, ev2, dt_out, ierr) &
      bind(c, name="nfor_mpi_step_pg2t")
    integer(c_int), value, intent(in) :: nloc, nx_l, ny_
    real(c_double), intent(in) :: rho(nloc), mx(nloc), my(nloc), e(nloc), ev(nloc)
    real(c_double), value, intent(in) :: dx, dy, gamma, cfl, dt_cap
    integer(c_int), value, intent(in) :: bc_w, bc_e, bc_s, bc_n
    real(c_double), intent(in) :: bc_vals(5, 4)
    integer(c_int), intent(in) :: solid(nloc)   ! 0/1
    real(c_double), intent(out) :: rho2(nloc), mx2(nloc), my2(nloc), e2(nloc), ev2(nloc)
    real(c_double), intent(out) :: dt_out
    integer(c_int), intent(out) :: ierr
    logical :: sl(nx_l+4, ny_)
    real(RP) :: rf(nx_l+4, ny_), mf(nx_l+4, ny_), myf(nx_l+4, ny_), ef(nx_l+4, ny_)
    real(RP) :: r2(nx_l+4, ny_), m2(nx_l+4, ny_), my2f(nx_l+4, ny_), e2f(nx_l+4, ny_)
    real(RP) :: evin(nx_l+4, ny_), evout(nx_l+4, ny_)
    real(RP) :: dt
    integer :: j, i

    ierr = 0
    do j = 1, ny_
      do i = 1, nx_l
        rf(i+2,j) = rho((j-1)*nx_l + i)
        mf(i+2,j) = mx((j-1)*nx_l + i)
        myf(i+2,j) = my((j-1)*nx_l + i)
        ef(i+2,j) = e((j-1)*nx_l + i)
        evin(i+2,j) = ev((j-1)*nx_l + i)
        sl(i+2,j) = solid((j-1)*nx_l + i) == 1
      end do
      rf(1,j) = rf(3,j); rf(2,j) = rf(3,j)
      rf(nx_l+3,j) = rf(nx_l+2,j); rf(nx_l+4,j) = rf(nx_l+2,j)
      mf(1,j) = mf(3,j); mf(2,j) = mf(3,j)
      mf(nx_l+3,j) = mf(nx_l+2,j); mf(nx_l+4,j) = mf(nx_l+2,j)
      myf(1,j) = myf(3,j); myf(2,j) = myf(3,j)
      myf(nx_l+3,j) = myf(nx_l+2,j); myf(nx_l+4,j) = myf(nx_l+2,j)
      ef(1,j) = ef(3,j); ef(2,j) = ef(3,j)
      ef(nx_l+3,j) = ef(nx_l+2,j); ef(nx_l+4,j) = ef(nx_l+2,j)
      evin(1,j) = evin(3,j); evin(2,j) = evin(3,j)
      evin(nx_l+3,j) = evin(nx_l+2,j); evin(nx_l+4,j) = evin(nx_l+2,j)
      sl(1,j) = .false.; sl(2,j) = .false.
      sl(nx_l+3,j) = .false.; sl(nx_l+4,j) = .false.
    end do
    call halo_exchange_cons(rf, mf, myf, ef, nx_l, ny_)

    call pg2d_step_rk2(rf, mf, myf, ef, nx_l, ny_, &
                       dx, dy, gamma, cfl, dt_cap, g_decomp, &
                       bc_w, bc_e, bc_s, bc_n, bc_vals, sl, &
                       r2, m2, my2f, e2f, dt)
    dt_out = dt

    ! post-march padded state for the substep, with halos refreshed
    call halo_exchange_cons(r2, m2, my2f, e2f, nx_l, ny_)
    do j = 1, ny_
      do i = 1, nx_l
        rho2((j-1)*nx_l + i) = r2(i+2,j)
        mx2((j-1)*nx_l + i)  = m2(i+2,j)
        my2((j-1)*nx_l + i)  = my2f(i+2,j)
        e2((j-1)*nx_l + i)   = e2f(i+2,j)
      end do
    end do

    call halo_exchange_mask(sl, nx_l, ny_)
    call halo_exchange_plane(evin, nx_l, ny_)
    call park2t_substep_slab(r2, m2, my2f, e2f, evin, sl, nx_l, ny_, &
                             dx, dy, gamma, dt, &
                             g_decomp%west_domain_edge, g_decomp%east_domain_edge, &
                             evout)
    do j = 1, ny_
      do i = 1, nx_l
        ev2((j-1)*nx_l + i) = evout(i+2,j)
      end do
    end do
  end subroutine nfor_mpi_step_pg2t

end module nuforpg2tmpi

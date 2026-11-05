!*****************************************************************************************
!>
!  Unit tests for the Krylov (iterative) linear solvers: each one solves a
!  symmetric positive definite tridiagonal system with and without a
!  Jacobi preconditioner, and the result is checked against the exact
!  solution.

    module test_krylov

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUN_SUCCESS, SUN_PREC_NONE, SUN_PREC_LEFT, SUN_PREC_RIGHT, &
                                SUN_PREC_BOTH, SUN_CLASSICAL_GS, SUN_ERR_ARG_INCOMPATIBLE
    use fvode_linsol, only: sun_linear_solver, ls_problem
    use fvode_linsol_spgmr, only: sun_linsol_spgmr
    use fvode_linsol_spfgmr, only: sun_linsol_spfgmr
    use fvode_linsol_spbcgs, only: sun_linsol_spbcgs
    use fvode_linsol_sptfqmr, only: sun_linsol_sptfqmr
    use fvode_linsol_pcg, only: sun_linsol_pcg
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

    private

    public :: collect_krylov

    integer,parameter :: n = 20  !! system size
    real(wp),parameter :: tol = 1.0e-10_wp  !! linear solver tolerance

    type,extends(ls_problem) :: tridiag_problem
        !! A = tridiag(-1, d_i, -1) with d_i = 3 + i/n (SPD), and the Jacobi
        !! preconditioner P = diag(A)
        integer :: natimes = 0  !! number of A*v products
        integer :: npsolve = 0  !! number of preconditioner solves
    contains
        procedure :: ls_atimes => tri_atimes  !! z = A*v
        procedure :: ls_psetup => tri_psetup  !! nothing to set up
        procedure :: ls_psolve => tri_psolve  !! z = P^{-1} r
    end type tridiag_problem

    contains
!*****************************************************************************************

    subroutine collect_krylov(testsuite)
    !! collect all tests in this suite
    type(unittest_type),allocatable,intent(out) :: testsuite(:)  !! the unit tests of this suite
    testsuite = [ new_unittest("spgmr",         test_spgmr),   &
                  new_unittest("spgmr_classical_gs", test_spgmr_cgs), &
                  new_unittest("spfgmr",        test_spfgmr),  &
                  new_unittest("spbcgs",        test_spbcgs),  &
                  new_unittest("sptfqmr",       test_sptfqmr), &
                  new_unittest("pcg",           test_pcg),     &
                  new_unittest("right_prec_needs_zero_guess", test_right_guess) ]
    end subroutine collect_krylov

    pure real(wp) function diag(i)
    !! diagonal element i of A
    integer,intent(in) :: i  !! row index
    diag = 3.0_wp + real(i, wp) / n
    end function diag

    function tri_atimes(me, v, z) result(flag)
    !! z = A*v
    class(tridiag_problem),intent(inout) :: me  !! the problem
    real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
    real(wp),dimension(:),intent(out) :: z  !! the product A*v
    integer :: flag  !! 0
    integer :: i  !! row index
    do i = 1, n
        z(i) = diag(i) * v(i)
        if (i > 1) z(i) = z(i) - v(i-1)
        if (i < n) z(i) = z(i) - v(i+1)
    end do
    me%natimes = me%natimes + 1
    flag = 0
    end function tri_atimes

    function tri_psetup(me) result(flag)
    !! nothing to set up
    class(tridiag_problem),intent(inout) :: me  !! the problem
    integer :: flag  !! 0
    flag = 0
    end function tri_psetup

    function tri_psolve(me, r, z, tol, lr) result(flag)
    !! z = P^{-1} r with P = diag(A)
    class(tridiag_problem),intent(inout) :: me  !! the problem
    real(wp),dimension(:),intent(in) :: r  !! the right-hand side
    real(wp),dimension(:),intent(out) :: z  !! the solution
    real(wp),intent(in) :: tol  !! tolerance (not used)
    integer,intent(in) :: lr  !! left or right (not used)
    integer :: flag  !! 0
    integer :: i  !! row index
    do i = 1, n
        z(i) = r(i) / diag(i)
    end do
    me%npsolve = me%npsolve + 1
    flag = 0
    end function tri_psolve

!*****************************************************************************************
!>
!  Solve A x = b (b = A*xtrue) with the solver S (already configured)
!  starting from a zero guess, and check the result.

    subroutine solve_and_check(S, error, scaled)
    class(sun_linear_solver),intent(inout) :: S  !! the solver
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    logical,intent(in),optional :: scaled  !! use unit scaling vectors?
    type(tridiag_problem) :: prob  !! the problem
    real(wp),dimension(n) :: x, b, xtrue, ones  !! solution, right-hand side, exact solution, unit weights
    integer :: i, flag  !! index and return flag
    do i = 1, n
        xtrue(i) = sin(real(i, wp))
    end do
    flag = prob%ls_atimes(xtrue, b)
    S%has_atimes = .true.
    S%has_psetup = .true.
    S%has_psolve = .true.
    call check(error, S%initialize(), SUN_SUCCESS, 'initialize')
    if (allocated(error)) return
    if (present(scaled)) then
        if (scaled) then
            ones = 1.0_wp
            call check(error, S%set_scaling_vectors(ones, ones), SUN_SUCCESS, 'scaling vectors')
            if (allocated(error)) return
        end if
    end if
    call check(error, S%setup(prob=prob), SUN_SUCCESS, 'setup')
    if (allocated(error)) return
    x = 0.0_wp
    flag = S%set_zero_guess(.true.)
    flag = S%solve(x=x, b=b, tol=tol, prob=prob)
    call check(error, flag, SUN_SUCCESS, 'solve')
    if (allocated(error)) return
    call check(error, S%num_iters() > 0, 'iterations counted')
    if (allocated(error)) return
    call check(error, S%res_norm() <= tol, 'residual norm below tol')
    if (allocated(error)) return
    call check(error, maxval(abs(x - xtrue)) < 1.0e-8_wp, 'solution accurate')
    end subroutine solve_and_check

!*****************************************************************************************
    subroutine test_spgmr(error)
    !! GMRES: no, left, right and two-sided preconditioning, with scaling
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    integer :: pretype  !! preconditioning type
    do pretype = SUN_PREC_NONE, SUN_PREC_BOTH
        block
            type(sun_linsol_spgmr) :: S  !! the solver
            S = sun_linsol_spgmr(n, pretype, n)
            call solve_and_check(S, error, scaled=.true.)
            if (allocated(error)) return
        end block
    end do
    end subroutine test_spgmr

    subroutine test_spgmr_cgs(error)
    !! GMRES with classical Gram-Schmidt and restarts (small Krylov dimension)
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(sun_linsol_spgmr) :: S  !! the solver
    integer :: flag  !! return flag
    S = sun_linsol_spgmr(n, SUN_PREC_LEFT, 4)
    flag = S%set_gs_type(SUN_CLASSICAL_GS)
    flag = S%set_max_restarts(20)
    call solve_and_check(S, error)
    end subroutine test_spgmr_cgs

    subroutine test_spfgmr(error)
    !! flexible GMRES, without and with (right) preconditioning
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(sun_linsol_spfgmr) :: S  !! the solver
    S = sun_linsol_spfgmr(n, SUN_PREC_NONE, n)
    call solve_and_check(S, error, scaled=.true.)
    if (allocated(error)) return
    S = sun_linsol_spfgmr(n, SUN_PREC_LEFT, n)  ! means right preconditioning
    call check(error, S%pretype, SUN_PREC_RIGHT, 'SPFGMR uses right preconditioning')
    if (allocated(error)) return
    call solve_and_check(S, error)
    end subroutine test_spfgmr

    subroutine test_spbcgs(error)
    !! Bi-CGStab: no, left, right and two-sided preconditioning
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    integer :: pretype  !! preconditioning type
    do pretype = SUN_PREC_NONE, SUN_PREC_BOTH
        block
            type(sun_linsol_spbcgs) :: S  !! the solver
            S = sun_linsol_spbcgs(n, pretype, 2*n)
            call solve_and_check(S, error, scaled=.true.)
            if (allocated(error)) return
        end block
    end do
    end subroutine test_spbcgs

    subroutine test_sptfqmr(error)
    !! TFQMR: no, left, right and two-sided preconditioning
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    integer :: pretype  !! preconditioning type
    do pretype = SUN_PREC_NONE, SUN_PREC_BOTH
        block
            type(sun_linsol_sptfqmr) :: S  !! the solver
            S = sun_linsol_sptfqmr(n, pretype, 2*n)
            call solve_and_check(S, error, scaled=.true.)
            if (allocated(error)) return
        end block
    end do
    end subroutine test_sptfqmr

    subroutine test_pcg(error)
    !! conjugate gradients, without and with preconditioning
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(sun_linsol_pcg) :: S  !! the solver
    S = sun_linsol_pcg(n, SUN_PREC_NONE, n)
    call solve_and_check(S, error, scaled=.true.)
    if (allocated(error)) return
    S = sun_linsol_pcg(n, SUN_PREC_LEFT, n)
    call solve_and_check(S, error)
    end subroutine test_pcg

    subroutine test_right_guess(error)
    !! right preconditioning with a nonzero initial guess is rejected by
    !! SPBCGS and SPTFQMR
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(sun_linsol_spbcgs) :: S1  !! Bi-CGStab solver
    type(sun_linsol_sptfqmr) :: S2  !! TFQMR solver
    type(tridiag_problem) :: prob  !! the problem
    real(wp),dimension(n) :: x, b  !! solution and right-hand side
    b = 1.0_wp
    x = 1.0_wp
    S1 = sun_linsol_spbcgs(n, SUN_PREC_RIGHT, n)
    S1%has_atimes = .true.; S1%has_psolve = .true.
    call check(error, S1%solve(x=x, b=b, tol=tol, prob=prob), SUN_ERR_ARG_INCOMPATIBLE, 'SPBCGS')
    if (allocated(error)) return
    S2 = sun_linsol_sptfqmr(n, SUN_PREC_RIGHT, n)
    S2%has_atimes = .true.; S2%has_psolve = .true.
    call check(error, S2%solve(x=x, b=b, tol=tol, prob=prob), SUN_ERR_ARG_INCOMPATIBLE, 'SPTFQMR')
    end subroutine test_right_guess

    end module test_krylov
!*****************************************************************************************

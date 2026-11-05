!*****************************************************************************************
!>
!  Unit tests for optional features without a C example: a user linear
!  system function ([[CVodeSetLinSysFn]]), an alternative right-hand side
!  for the nonlinear system ([[CVodeSetNlsRhsFn]]), the linear solver
!  statistics and saved Jacobian getters, and the projection input checks.
!  Each feature is compared with an equivalent run through the internal
!  path, which must give identical results.

    module test_options

    use fvode
    use fvode_core
    use fvode_matrix, only: sun_matrix, dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

    private

    public :: collect_options

    integer,parameter :: neq = 3  !! number of equations (Robertson problem)

    type,extends(cvode_solver) :: robertson_solver
        !! the solver plus the data of the user linear system function and
        !! the counter of the alternative right-hand side
        real(wp),dimension(neq,neq) :: Jsaved = 0.0_wp  !! saved Jacobian (user linear system function)
        integer :: nls_f_calls = 0  !! calls of the alternative right-hand side
    end type robertson_solver

    contains
!*****************************************************************************************

    subroutine collect_options(testsuite)
    !! collect all tests in this suite
    type(unittest_type),allocatable,intent(out) :: testsuite(:)  !! the unit tests of this suite
    testsuite = [ new_unittest("lin_sys_fn",      test_lin_sys_fn), &
                  new_unittest("nls_rhs_fn",      test_nls_rhs_fn), &
                  new_unittest("lin_solve_stats", test_lin_solve_stats), &
                  new_unittest("proj_input",      test_proj_input) ]
    end subroutine collect_options

    function f(me, t, y, ydot) result(retval)
    !! Robertson chemical kinetics right-hand side
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: ydot  !! dy/dt
    integer :: retval  !! 0 on success
    ydot(1) = -0.04_wp * y(1) + 1.0e4_wp * y(2) * y(3)
    ydot(3) = 3.0e7_wp * y(2) * y(2)
    ydot(2) = -ydot(1) - ydot(3)
    retval = 0
    end function f

    function f_counted(me, t, y, ydot) result(retval)
    !! the same right-hand side, counting the calls
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: ydot  !! dy/dt
    integer :: retval  !! 0 on success
    select type (me)
    class is (robertson_solver)
        me%nls_f_calls = me%nls_f_calls + 1
    end select
    retval = f(me, t, y, ydot)
    end function f_counted

    pure subroutine robertson_jac(y, J)
    !! the analytic Jacobian
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(neq,neq),intent(out) :: J  !! df/dy
    J(1,1) = -0.04_wp
    J(1,2) = 1.0e4_wp * y(3)
    J(1,3) = 1.0e4_wp * y(2)
    J(3,1) = 0.0_wp
    J(3,2) = 6.0e7_wp * y(2)
    J(3,3) = 0.0_wp
    J(2,1) = -J(1,1) - J(3,1)
    J(2,2) = -J(1,2) - J(3,2)
    J(2,3) = -J(1,3) - J(3,3)
    end subroutine robertson_jac

    function jac(me, t, y, fy, Jm, tmp1, tmp2, tmp3) result(retval)
    !! user Jacobian function
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t,y)
    class(sun_matrix),intent(inout) :: Jm  !! the Jacobian
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors
    integer :: retval  !! 0 on success
    select type (Jm)
    type is (dense_matrix)
        call robertson_jac(y, Jm%data)
    end select
    retval = 0
    end function jac

    function linsys(me, t, y, fy, A, jok, jcur, gamma, tmp1, tmp2, tmp3) result(retval)
    !! user linear system function: A = I - gamma*J, reusing the saved J if
    !! jok (the same as the internal linear system function)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t,y)
    class(sun_matrix),intent(inout) :: A  !! the linear system matrix
    logical,intent(in) :: jok  !! may the saved J be reused?
    logical,intent(out) :: jcur  !! was J re-evaluated?
    real(wp),intent(in) :: gamma  !! the scalar gamma
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors
    integer :: retval  !! 0 on success
    integer :: i  !! diagonal index
    select type (me)
    class is (robertson_solver)
        if (jok) then
            jcur = .false.
        else
            call robertson_jac(y, me%Jsaved)
            jcur = .true.
        end if
        select type (A)
        type is (dense_matrix)
            A%data = -gamma * me%Jsaved
            do i = 1, neq
                A%data(i,i) = A%data(i,i) + 1.0_wp
            end do
        end select
    end select
    retval = 0
    end function linsys

    subroutine run(mode, y, nst, nje, ncalls)
    !! integrate the Robertson problem to t = 40 with BDF + dense Newton
    !! (mode 0: user Jacobian; 1: user linear system function; 2: user
    !! Jacobian and a counting nonlinear system right-hand side)
    integer,intent(in) :: mode  !! which variant
    real(wp),dimension(neq),intent(out) :: y  !! the solution at t = 40
    integer(li),intent(out) :: nst  !! number of steps
    integer(li),intent(out) :: nje  !! number of Jacobian evaluations
    integer,intent(out) :: ncalls  !! calls of the counting right-hand side
    type(robertson_solver) :: cv  !! the solver
    real(wp) :: t  !! return time
    integer :: flag  !! return flag
    y = [1.0_wp, 0.0_wp, 0.0_wp]
    call cv%create(CV_BDF, flag)
    call cv%init(f, 0.0_wp, y, flag)
    call cv%sv_tolerances(1.0e-4_wp, [1.0e-8_wp, 1.0e-14_wp, 1.0e-6_wp], flag)
    call cv%set_linear_solver(sun_linsol_dense(int(neq, ip)), dense_matrix(int(neq, ip), int(neq, ip)), flag)
    select case (mode)
    case (1)
        call cv%set_lin_sys_fn(linsys, flag)
    case (2)
        call cv%set_jac_fn(jac, flag)
        call cv%set_nls_rhs_fn(f_counted, flag)
    case default
        call cv%set_jac_fn(jac, flag)
    end select
    call cv%solve(40.0_wp, y, t, CV_NORMAL, flag)
    nst = cv%get_num_steps()
    nje = cv%get_num_jac_evals()
    ncalls = cv%nls_f_calls
    call cv%destroy()
    end subroutine run

!*****************************************************************************************
    subroutine test_lin_sys_fn(error)
    !! a user linear system function equivalent to the internal one gives
    !! identical results
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    real(wp),dimension(neq) :: y0, y1  !! solutions of the two runs
    integer(li) :: nst0, nst1, nje0, nje1  !! steps and Jacobian evaluations
    integer :: nc  !! (unused)
    call run(0, y0, nst0, nje0, nc)
    call run(1, y1, nst1, nje1, nc)
    call check(error, nst1 == nst0, 'same number of steps')
    if (allocated(error)) return
    call check(error, nje1 == nje0, 'same number of Jacobian evaluations')
    if (allocated(error)) return
    call check(error, all(y1 == y0), 'identical solution')
    end subroutine test_lin_sys_fn

    subroutine test_nls_rhs_fn(error)
    !! the alternative nonlinear system right-hand side is used, and an
    !! equivalent one gives identical results
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    real(wp),dimension(neq) :: y0, y2  !! solutions of the two runs
    integer(li) :: nst0, nst2, nje0, nje2  !! steps and Jacobian evaluations
    integer :: nc0, nc2  !! calls of the counting right-hand side
    call run(0, y0, nst0, nje0, nc0)
    call run(2, y2, nst2, nje2, nc2)
    call check(error, nc0 == 0 .and. nc2 > 0, 'alternative right-hand side called')
    if (allocated(error)) return
    call check(error, nst2 == nst0 .and. all(y2 == y0), 'identical results')
    end subroutine test_nls_rhs_fn

    subroutine test_lin_solve_stats(error)
    !! `get_lin_solve_stats` agrees with the individual getters, and
    !! `get_jac` returns the saved Jacobian
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(robertson_solver) :: cv  !! the solver
    real(wp),dimension(neq) :: y  !! solution
    real(wp) :: t  !! return time
    integer :: flag  !! return flag
    integer(li) :: nje, nfeLS, nli, ncfl, npe, nps, njts, njt  !! the statistics
    class(sun_matrix),allocatable :: J  !! the saved Jacobian
    y = [1.0_wp, 0.0_wp, 0.0_wp]
    call cv%create(CV_BDF, flag)
    call cv%init(f, 0.0_wp, y, flag)
    call cv%ss_tolerances(1.0e-4_wp, 1.0e-8_wp, flag)
    call cv%set_linear_solver(sun_linsol_dense(int(neq, ip)), dense_matrix(int(neq, ip), int(neq, ip)), flag)
    call cv%solve(1.0_wp, y, t, CV_NORMAL, flag)
    call cv%get_lin_solve_stats(nje, nfeLS, nli, ncfl, npe, nps, njts, njt)
    call check(error, nje == cv%get_num_jac_evals() .and. nje > 0, 'nje')
    if (allocated(error)) return
    call check(error, nfeLS == cv%get_num_lin_rhs_evals() .and. nfeLS == neq * nje, 'nfeLS')
    if (allocated(error)) return
    call check(error, nli == 0 .and. npe == 0 .and. nps == 0 .and. njts == 0 .and. njt == 0, &
               'no iterative statistics')
    if (allocated(error)) return
    call cv%get_jac(J, flag)
    call check(error, flag, CVLS_SUCCESS, 'get_jac flag')
    if (allocated(error)) return
    call check(error, allocated(J), 'saved Jacobian returned')
    if (allocated(error)) return
    select type (J)
    type is (dense_matrix)
        call check(error, abs(J%data(1,1) + 0.04_wp) < 1.0e-6_wp, 'J(1,1) ~ -0.04 (difference quotient)')
    class default
        call check(error, .false., 'dense Jacobian expected')
    end select
    call cv%destroy()
    end subroutine test_lin_solve_stats

    function proj_dummy(me, t, ycur, corr, epsProj, err) result(retval)
    !! a projection function that does nothing
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(inout) :: ycur  !! the solution
    real(wp),dimension(:),intent(out) :: corr  !! the correction (zero)
    real(wp),intent(in) :: epsProj  !! projection tolerance
    real(wp),dimension(:),intent(inout),optional :: err  !! the error estimate
    integer :: retval  !! 0 on success
    corr = 0.0_wp
    retval = 0
    end function proj_dummy

    subroutine test_proj_input(error)
    !! projection requires BDF, and its setters require the projection memory
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(neq) :: y  !! initial condition
    integer :: flag  !! return flag
    y = [1.0_wp, 0.0_wp, 0.0_wp]
    call cv%create(CV_ADAMS, flag)
    call cv%set_error_output(.false.)
    call cv%init(f, 0.0_wp, y, flag)
    call cv%set_proj_fn(proj_dummy, flag)
    call check(error, flag, CV_ILL_INPUT, 'projection with Adams')
    if (allocated(error)) return
    call cv%set_proj_frequency(1_li, flag)
    call check(error, flag, CV_PROJ_MEM_NULL, 'no projection memory')
    if (allocated(error)) return
    call cv%create(CV_BDF, flag)
    call cv%set_error_output(.false.)
    call cv%init(f, 0.0_wp, y, flag)
    call cv%set_proj_fn(proj_dummy, flag)
    call check(error, flag, CV_SUCCESS, 'projection with BDF')
    if (allocated(error)) return
    call cv%set_eps_proj(-1.0_wp, flag)
    call check(error, flag, CV_SUCCESS, 'eps_proj default')
    call cv%destroy()
    end subroutine test_proj_input

    end module test_options
!*****************************************************************************************

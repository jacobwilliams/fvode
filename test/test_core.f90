!*****************************************************************************************
!>
!  Unit tests for [[fvode_core]]: accuracy on a problem with a known
!  solution, interpolation with `get_dky`, and the input checks of the API.
!  (Step-by-step agreement with C CVODE is checked by the traces in
!  `reference/`.)

    module test_core

    use fvode
    use fvode_core
    use fvode_nonlinsol_fixedpoint
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

    private

    public :: collect_core

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Collect all tests in this suite.

    subroutine collect_core(testsuite)

    type(unittest_type),allocatable,intent(out) :: testsuite(:)  !! the unit tests of this suite

    testsuite = [ new_unittest("decay_adams_fp",  test_decay),       &
                  new_unittest("get_dky",         test_get_dky),     &
                  new_unittest("stop_time",       test_stop_time),   &
                  new_unittest("input_checks",    test_input_checks),&
                  new_unittest("flag_names",      test_flag_names)   ]

    end subroutine collect_core
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y' = -y` (two uncoupled components with different rates).

    function decay(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success

    ydot(1) = -y(1)
    ydot(2) = -0.5_wp * y(2)
    retval = 0

    end function decay
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set up a solver for [[decay]] with Adams + fixed point.

    subroutine setup(cv, y, error)

    type(cvode_solver),intent(inout) :: cv  !! the solver to set up
    real(wp),dimension(2),intent(inout) :: y  !! initial condition (set here)
    type(error_type),allocatable,intent(inout) :: error  !! error information (allocated on failure)

    integer :: flag  !! return flag of the last solver call

    y = [1.0_wp, 2.0_wp]
    call cv%create(CV_ADAMS, flag);                     call check(error, flag, CV_SUCCESS)
    if (allocated(error)) return
    call cv%init(decay, 0.0_wp, y, flag);               call check(error, flag, CV_SUCCESS)
    if (allocated(error)) return
    call cv%ss_tolerances(1.0e-8_wp, 1.0e-10_wp, flag); call check(error, flag, CV_SUCCESS)
    if (allocated(error)) return
    call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(2, 0), flag)
    call check(error, flag, CV_SUCCESS)

    end subroutine setup
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_decay(error)
    !! test the accuracy of the solution of [[decay]]

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y  !! solution vector
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag of the last solver call

    call setup(cv, y, error); if (allocated(error)) return

    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_SUCCESS);                 if (allocated(error)) return
    call check(error, t, 2.0_wp);                        if (allocated(error)) return
    call check(error, y(1), exp(-2.0_wp), thr=1.0e-6_wp, rel=.true.); if (allocated(error)) return
    call check(error, y(2), 2.0_wp*exp(-1.0_wp), thr=1.0e-6_wp, rel=.true.)

    end subroutine test_decay
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_get_dky(error)
    !! test interpolation with get_dky inside the last step

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y, dky  !! solution and interpolated derivative
    real(wp) :: t, tmid  !! time reached and a time inside the last step
    integer :: flag  !! return flag of the last solver call

    call setup(cv, y, error); if (allocated(error)) return

    ! take a few internal steps
    t = 0.0_wp
    do while (cv%get_num_steps() < 20)
        call cv%solve(10.0_wp, y, t, CV_ONE_STEP, flag)
        call check(error, flag, CV_SUCCESS); if (allocated(error)) return
    end do

    ! k = 0 at the current time is the returned solution
    call cv%get_dky(t, 0, dky, flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return
    call check(error, all(dky == y), 'dky(k=0) /= y at tn'); if (allocated(error)) return

    ! k = 0 and k = 1 inside the last step approximate y and y'
    tmid = t - 0.5_wp * cv%get_last_step()
    call cv%get_dky(tmid, 0, dky, flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return
    call check(error, dky(1), exp(-tmid), thr=1.0e-6_wp, rel=.true.); if (allocated(error)) return
    call cv%get_dky(tmid, 1, dky, flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return
    call check(error, dky(1), -exp(-tmid), thr=1.0e-5_wp, rel=.true.); if (allocated(error)) return

    ! outside the last step, and an illegal k
    call cv%set_error_output(.false.)
    call cv%get_dky(t + 1.0_wp, 0, dky, flag)
    call check(error, flag, CV_BAD_T);             if (allocated(error)) return
    call cv%get_dky(t, cv%get_current_order() + 1, dky, flag)
    call check(error, flag, CV_BAD_K)

    end subroutine test_get_dky
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_stop_time(error)
    !! test the stop time option

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y  !! solution vector
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag of the last solver call

    call setup(cv, y, error); if (allocated(error)) return

    call cv%set_stop_time(0.75_wp, flag)
    call check(error, flag, CV_SUCCESS);      if (allocated(error)) return
    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_TSTOP_RETURN); if (allocated(error)) return
    call check(error, t, 0.75_wp);            if (allocated(error)) return
    call check(error, y(1), exp(-0.75_wp), thr=1.0e-6_wp, rel=.true.); if (allocated(error)) return

    ! the stop time is cleared after it is reached: continue to tout
    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_SUCCESS);      if (allocated(error)) return
    call check(error, t, 2.0_wp)

    end subroutine test_stop_time
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_input_checks(error)
    !! test the input checks of the API (expected error flags)

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y  !! solution vector
    real(wp),dimension(3) :: y3  !! solution vector of the wrong size
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag of the last solver call

    ! bad lmm
    call cv%set_error_output(.false.)
    call cv%create(3, flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return

    ! calls before init
    call cv%create(CV_BDF, flag)
    call cv%set_error_output(.false.)
    call cv%solve(1.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_NO_MALLOC);    if (allocated(error)) return
    call cv%ss_tolerances(1.0e-6_wp, 1.0e-8_wp, flag)
    call check(error, flag, CV_NO_MALLOC);    if (allocated(error)) return

    ! no tolerances
    call cv%init(decay, 0.0_wp, [1.0_wp, 1.0_wp], flag)
    call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(2, 0), flag)
    call cv%solve(1.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return

    ! bad tolerances
    call cv%ss_tolerances(-1.0_wp, 1.0e-8_wp, flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return
    call cv%sv_tolerances(1.0e-6_wp, [1.0e-8_wp, -1.0_wp], flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return

    ! bad itask, wrong-size output, tout too close to t0
    call cv%ss_tolerances(1.0e-6_wp, 1.0e-8_wp, flag)
    call cv%solve(1.0_wp, y, t, 7, flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return
    call cv%solve(1.0_wp, y3, t, CV_NORMAL, flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return
    call cv%solve(0.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_TOO_CLOSE);    if (allocated(error)) return

    ! inconsistent step bounds
    call cv%set_max_step(1.0_wp, flag)
    call cv%set_min_step(2.0_wp, flag)
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return

    ! cannot increase the maximum order
    call cv%set_max_ord(6, flag)   ! BDF maximum is 5
    call check(error, flag, CV_ILL_INPUT)

    end subroutine test_input_checks
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_flag_names(error)
    !! test the return flag names

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    call check(error, CVodeGetReturnFlagName(CV_SUCCESS), 'CV_SUCCESS')
    if (allocated(error)) return
    call check(error, CVodeGetReturnFlagName(CV_TOO_MUCH_WORK), 'CV_TOO_MUCH_WORK')
    if (allocated(error)) return
    call check(error, CVodeGetReturnFlagName(12345), 'NONE')

    end subroutine test_flag_names
!*****************************************************************************************

    end module test_core
!*****************************************************************************************

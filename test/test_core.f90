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
                  new_unittest("flag_names",      test_flag_names),  &
                  new_unittest("root_location",   test_root_location), &
                  new_unittest("root_direction",  test_root_direction), &
                  new_unittest("root_at_t0",      test_root_at_t0),  &
                  new_unittest("constraints_input", test_constraints_input), &
                  new_unittest("stald_input",     test_stald_input), &
                  new_unittest("adams_max_order", test_adams_max_order), &
                  new_unittest("adams_high_order", test_adams_high_order)  ]

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

    ! cannot increase the BDF maximum order
    call cv%set_max_ord(6, flag)   ! BDF maximum is 5
    call check(error, flag, CV_ILL_INPUT);    if (allocated(error)) return
    call cv%set_max_ord(0, flag)
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

!*****************************************************************************************
!>
!  Root function for [[decay]]: g(1) = y(1) - 0.5 (a root at t = ln 2,
!  decreasing).

    function g_half(me, t, y, gout) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (not used)
    real(wp),intent(in) :: t  !! current time (not used)
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: gout  !! root function value
    integer :: retval  !! 0 (success)

    gout(1) = y(1) - 0.5_wp
    retval = 0

    end function g_half
!*****************************************************************************************

!*****************************************************************************************
!>
!  Root function for [[decay]]: g(1) = y(1) - 1, which is zero at t0.

    function g_at_t0(me, t, y, gout) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (not used)
    real(wp),intent(in) :: t  !! current time (not used)
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: gout  !! root function value
    integer :: retval  !! 0 (success)

    gout(1) = y(1) - 1.0_wp
    retval = 0

    end function g_at_t0
!*****************************************************************************************

!*****************************************************************************************
!>
!  A root of y1 - 0.5 is located at t = ln 2, decreasing, and the
!  integration then continues to tout.

    subroutine test_root_location(error)

    type(error_type),allocatable,intent(out) :: error  !! test error (allocated on failure)

    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(2) :: y  !! solution
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag
    integer,dimension(1) :: rootsfound  !! root information

    call setup(cv, y, error); if (allocated(error)) return
    call cv%root_init(1, g_half, flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return

    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_ROOT_RETURN);       if (allocated(error)) return
    call check(error, t, log(2.0_wp), thr=1.0e-6_wp, rel=.true.); if (allocated(error)) return
    call cv%get_root_info(rootsfound)
    call check(error, rootsfound(1), -1, 'root should be decreasing'); if (allocated(error)) return
    call check(error, cv%get_num_g_evals() > 0, 'no root function evaluations'); if (allocated(error)) return

    ! continue to tout: no more roots
    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return
    call check(error, t, 2.0_wp)

    end subroutine test_root_location
!*****************************************************************************************

!*****************************************************************************************
!>
!  With the direction restricted to increasing zero-crossings, the
!  (decreasing) root of y1 - 0.5 is not reported.

    subroutine test_root_direction(error)

    type(error_type),allocatable,intent(out) :: error  !! test error (allocated on failure)

    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(2) :: y  !! solution
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag

    call setup(cv, y, error); if (allocated(error)) return
    call cv%set_error_output(.false.)
    call cv%set_root_direction([1], flag)
    call check(error, flag, CV_ILL_INPUT, 'direction before root_init'); if (allocated(error)) return
    call cv%root_init(1, g_half, flag)
    call cv%set_root_direction([1], flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return

    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_SUCCESS, 'decreasing root should be ignored'); if (allocated(error)) return
    call check(error, t, 2.0_wp)

    end subroutine test_root_direction
!*****************************************************************************************

!*****************************************************************************************
!>
!  A root function that is zero at t0 is made inactive: no root is reported
!  at t0, and the integration reaches tout.

    subroutine test_root_at_t0(error)

    type(error_type),allocatable,intent(out) :: error  !! test error (allocated on failure)

    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(2) :: y  !! solution
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag

    call setup(cv, y, error); if (allocated(error)) return
    call cv%root_init(1, g_at_t0, flag)
    call cv%set_no_inactive_root_warn(flag)
    call cv%solve(1.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_SUCCESS, 'no root expected');  if (allocated(error)) return
    call check(error, t, 1.0_wp)

    end subroutine test_root_at_t0
!*****************************************************************************************

!*****************************************************************************************
!>
!  Input checks of `set_constraints`: illegal values, and y0 violating the
!  constraints.

    subroutine test_constraints_input(error)

    type(error_type),allocatable,intent(out) :: error  !! test error (allocated on failure)

    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(2) :: y  !! solution
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag

    call setup(cv, y, error); if (allocated(error)) return
    call cv%set_error_output(.false.)

    ! values outside {0, +-1, +-2}
    call cv%set_constraints([3.0_wp, 0.0_wp], flag)
    call check(error, flag, CV_ILL_INPUT, 'constraint 3 accepted'); if (allocated(error)) return
    call cv%set_constraints([0.0_wp, 0.0_wp], flag)
    call check(error, flag, CV_ILL_INPUT, 'all-zero constraints accepted'); if (allocated(error)) return

    ! y0 = (1, 2) violates y2 < 0
    call cv%set_constraints([0.0_wp, -2.0_wp], flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return
    call cv%solve(1.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_ILL_INPUT, 'y0 violating the constraints accepted'); if (allocated(error)) return

    ! satisfied constraints (y >= 0): the solution is unchanged, no failures
    call setup(cv, y, error); if (allocated(error)) return
    call cv%set_constraints([1.0_wp, 1.0_wp], flag)
    call cv%solve(2.0_wp, y, t, CV_NORMAL, flag)
    call check(error, flag, CV_SUCCESS);           if (allocated(error)) return
    call check(error, y(1), exp(-2.0_wp), thr=1.0e-6_wp, rel=.true.)

    end subroutine test_constraints_input
!*****************************************************************************************

!*****************************************************************************************
!>
!  Stability limit detection is only allowed with BDF.

    subroutine test_stald_input(error)

    type(error_type),allocatable,intent(out) :: error  !! test error (allocated on failure)

    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(2) :: y  !! solution
    integer :: flag  !! return flag

    call setup(cv, y, error); if (allocated(error)) return   ! (Adams)
    call cv%set_error_output(.false.)
    call cv%set_stab_lim_det(.true., flag)
    call check(error, flag, CV_ILL_INPUT, 'STALD accepted with Adams')

    end subroutine test_stald_input
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y1' = 20 t^19`, `y2' = 15 t^14` (solution `t^20`, `t^15`). The high
!  derivatives of the solution vanish, so high Adams orders pay off.

    function polynomial(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success

    ydot(1) = 20 * t**19
    ydot(2) = 15 * t**14
    retval = 0

    end function polynomial
!*****************************************************************************************

!*****************************************************************************************
!>
!  Harmonic oscillator `y1' = y2, y2' = -y1` (solution `cos t`, `-sin t`).

    function oscillator(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success

    ydot(1) =  y(2)
    ydot(2) = -y(1)
    retval = 0

    end function oscillator
!*****************************************************************************************

!*****************************************************************************************
!>
!  Integrate in one-step mode from `tstart` to `tout` (a
!  stop time), and update the largest order used.

    subroutine run_to(cv, tstart, tout, y, qused, error)

    type(cvode_solver),intent(inout) :: cv  !! the solver (already initialized)
    real(wp),intent(in) :: tstart  !! current time
    real(wp),intent(in) :: tout  !! final time
    real(wp),dimension(2),intent(inout) :: y  !! solution at tout
    integer,intent(inout) :: qused  !! largest order used so far
    type(error_type),allocatable,intent(inout) :: error  !! error information (allocated on failure)

    real(wp) :: t  !! time reached
    integer :: flag  !! return flag of the last solver call

    t = tstart
    call cv%set_stop_time(tout, flag); call check(error, flag, CV_SUCCESS)
    if (allocated(error)) return
    do while (t < tout)
        call cv%solve(tout, y, t, CV_ONE_STEP, flag)
        if (flag < 0) then
            call check(error, flag, CV_SUCCESS); return
        end if
        qused = max(qused, cv%get_last_order())
    end do

    end subroutine run_to
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_adams_max_order(error)
    !! the Adams maximum order can be decreased and increased, before the
    !! integration and during it (which enlarges the history array)

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    real(wp),parameter :: tf = 20.0_wp  !! final time
    real(wp),parameter :: tmid = 5.0_wp  !! time at which the maximum order is changed
    real(wp),parameter :: tol = 1.0e-9_wp  !! accuracy required of the solution

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y  !! solution vector
    real(wp),dimension(2) :: yref  !! solution with the default maximum order
    integer(li) :: nstref  !! number of steps with the default maximum order
    integer :: flag  !! return flag of the last solver call
    integer :: qused  !! largest order used
    integer :: icase  !! which test case

    do icase = 0, 3

        y = [1.0_wp, 0.0_wp]
        qused = 0
        call cv%create(CV_ADAMS, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return
        select case (icase)
        case (1)    ! decrease before init (zn is allocated for order 3)
            call cv%set_max_ord(3, flag)
        case (3)    ! increase before init
            call cv%set_max_ord(16, flag)
        end select
        call check(error, flag, CV_SUCCESS); if (allocated(error)) return
        call cv%init(oscillator, 0.0_wp, y, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return
        call cv%ss_tolerances(1.0e-13_wp, 1.0e-13_wp, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return
        call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(2, 0), flag)
        call check(error, flag, CV_SUCCESS); if (allocated(error)) return

        call run_to(cv, 0.0_wp, tmid, y, qused, error); if (allocated(error)) return

        select case (icase)
        case (1)    ! the order reached the cap; raise it (q == qmax when enlarged)
            call check(error, qused, 3); if (allocated(error)) return
            call cv%set_max_ord(8, flag); call check(error, flag, CV_SUCCESS)
            if (allocated(error)) return
        case (2)    ! raise the default maximum mid-run (q < qmax when enlarged)
            call cv%set_max_ord(16, flag); call check(error, flag, CV_SUCCESS)
            if (allocated(error)) return
        end select

        call run_to(cv, tmid, tf, y, qused, error); if (allocated(error)) return
        call check(error, y(1), cos(tf), thr=tol); if (allocated(error)) return
        call check(error, y(2), -sin(tf), thr=tol); if (allocated(error)) return

        select case (icase)
        case (0)
            yref = y
            nstref = cv%get_num_steps()
        case (1)
            call check(error, qused > 3 .and. qused <= 8, 'raised maximum order not used')
            if (allocated(error)) return
        case (2:3)
            ! orders above the default are not chosen for this problem in double
            ! precision, so the solution is the same as with the default maximum
            call check(error, cv%get_num_steps() == nstref, 'different number of steps')
            if (allocated(error)) return
            call check(error, y(1), yref(1), thr=10*epsilon(1.0_wp)); if (allocated(error)) return
            call check(error, y(2), yref(2), thr=10*epsilon(1.0_wp)); if (allocated(error)) return
        end select

    end do

    end subroutine test_adams_max_order
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_adams_high_order(error)
    !! with enough precision (a `-DREAL128` build), Adams orders above the
    !! default maximum of 12 are used once the maximum order is raised

    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)

    integer,parameter :: QDEF = 12  !! default Adams maximum order
    integer,parameter :: QBIG = 16  !! raised maximum order
    real(wp),parameter :: tf = 1.0_wp  !! final time
    real(wp),parameter :: tol = 1.0e-33_wp  !! relative and absolute tolerance

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y  !! solution vector
    integer :: flag  !! return flag of the last solver call
    integer :: qused  !! largest order used
    integer :: icase  !! 1: default maximum order, 2: raised

    if (precision(1.0_wp) < 30) return  ! needs quad precision

    do icase = 1, 2
        y = [0.0_wp, 0.0_wp]
        qused = 0
        call cv%create(CV_ADAMS, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return
        if (icase == 2) then
            call cv%set_max_ord(QBIG, flag); call check(error, flag, CV_SUCCESS)
            if (allocated(error)) return
        end if
        call cv%init(polynomial, 0.0_wp, y, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return
        call cv%ss_tolerances(tol, tol, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return
        call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(2, 0), flag)
        call check(error, flag, CV_SUCCESS); if (allocated(error)) return
        call cv%set_max_num_steps(100000_li, flag); call check(error, flag, CV_SUCCESS)
        if (allocated(error)) return

        call run_to(cv, 0.0_wp, tf, y, qused, error); if (allocated(error)) return
        call check(error, y(1), 1.0_wp, thr=1.0e-25_wp); if (allocated(error)) return
        call check(error, y(2), 1.0_wp, thr=1.0e-25_wp); if (allocated(error)) return

        if (icase == 1) then
            call check(error, qused, QDEF, 'the default maximum order was not reached')
        else
            call check(error, qused > QDEF .and. qused <= QBIG, &
                       'an order above the default maximum was not used')
        end if
        if (allocated(error)) return
    end do

    end subroutine test_adams_high_order
!*****************************************************************************************

    end module test_core
!*****************************************************************************************

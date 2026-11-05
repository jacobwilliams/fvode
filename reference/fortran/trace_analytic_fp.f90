!*****************************************************************************************
!>
!  Step trace: analytic problem (cv_analytic_fp), Adams + fixed point,
!  one row per CV_NORMAL output. Fortran counterpart of
!  `reference/trace_analytic_fp.c`.

    program trace_analytic_fp

    use fvode
    use fvode_core
    use fvode_nonlinsol_fixedpoint
    use trace_common

    implicit none

    real(wp),parameter :: lambda = -100.0_wp  !! stiffness parameter of the ODE
    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(1) :: y  !! solution vector
    real(wp) :: t, tout  !! time reached and next output time
    integer :: flag, iunit, i  !! return flag, output unit, and output counter

    y = 0.0_wp
    call cv%create(CV_ADAMS, flag)
    call cv%init(f, 0.0_wp, y, flag)
    call cv%ss_tolerances(1.0e-6_wp, 1.0e-10_wp, flag)
    call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(1, 0), flag)

    open(newunit=iunit, file='reference/out/fvode_analytic_fp.csv', status='replace')
    call trace_header(iunit, 1)
    t = 0.0_wp
    tout = 0.0_wp
    do i = 1, 10
        tout = tout + 1.0_wp
        call cv%solve(tout, y, t, CV_NORMAL, flag)
        call trace_row(iunit, cv, flag, t, y)
        if (flag < 0) exit
    end do
    close(iunit)

    contains

        function f(me, t, y, ydot) result(retval)
        !! right-hand side of the analytic problem
        class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
        real(wp),intent(in) :: t  !! current value of the independent variable
        real(wp),dimension(:),intent(in) :: y  !! current state vector
        real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
        integer :: retval  !! return value: 0 on success
        ydot(1) = lambda*y(1) + 1.0_wp/(1.0_wp + t*t) - lambda*atan(t)
        retval = 0
        end function f

    end program trace_analytic_fp
!*****************************************************************************************

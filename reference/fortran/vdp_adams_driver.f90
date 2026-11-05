!*****************************************************************************************
!>
!  Shared van der Pol driver (cvDirectDemo_ls problem 1): Adams + fixed-point
!  iteration, reltol = 0, abstol = 1e-6. Fortran counterpart of
!  `reference/vdp_adams_driver.h`.

    module vdp_adams_driver

    use fvode
    use fvode_core
    use fvode_nonlinsol_fixedpoint
    use trace_common

    implicit none

    private

    public :: run_vdp

    real(wp),parameter :: eta = 3.0_wp  !! van der Pol parameter
    real(wp),parameter :: t0 = 0.0_wp  !! initial time
    real(wp),parameter :: tfinal = 1.39283880203_wp + 3.0_wp * 2.214773875_wp  !! final time

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Run the problem and write the step trace to `outfile`.

    subroutine run_vdp(m, maxiters, beta, outfile)

    integer,intent(in) :: m           !! number of Anderson acceleration vectors (0: none)
    integer,intent(in) :: maxiters    !! max nonlinear iterations per step (0: default)
    real(wp),intent(in) :: beta       !! damping parameter (1: none)
    character(len=*),intent(in) :: outfile  !! file for the trace

    type(cvode_solver) :: cv  !! the CVODE solver
    real(wp),dimension(2) :: y  !! solution vector
    real(wp) :: t  !! time reached
    integer :: flag, iunit  !! return flag and output unit

    y = [2.0_wp, 0.0_wp]

    call cv%create(CV_ADAMS, flag);                    call check('create')
    call cv%init(f, t0, y, flag);                      call check('init')
    call cv%ss_tolerances(0.0_wp, 1.0e-6_wp, flag);    call check('ss_tolerances')
    call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(2, m, beta), flag)
    call check('set_nonlinear_solver')
    if (maxiters > 0) then
        call cv%set_max_nonlin_iters(maxiters, flag);  call check('set_max_nonlin_iters')
    end if
    call cv%set_max_num_steps(100000_li, flag);        call check('set_max_num_steps')

    open(newunit=iunit, file=outfile, status='replace')
    call trace_header(iunit, 2)
    t = t0
    do while (t < tfinal)
        call cv%solve(tfinal, y, t, CV_ONE_STEP, flag)
        call trace_row(iunit, cv, flag, t, y)
        if (flag < 0) exit
    end do
    close(iunit)

    contains

        subroutine check(name)
        !! stop if the last solver call returned an error
        character(len=*),intent(in) :: name  !! name of the call that was checked
        if (flag < 0) then
            write(*,'(a,i0)') name//' failed with flag ', flag
            error stop
        end if
        end subroutine check

    end subroutine run_vdp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Right-hand side of the van der Pol equation.

    function f(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success

    ydot(1) = y(2)
    ydot(2) = (1.0_wp - y(1)*y(1)) * eta * y(2) - y(1)
    retval = 0

    end function f
!*****************************************************************************************

    end module vdp_adams_driver
!*****************************************************************************************

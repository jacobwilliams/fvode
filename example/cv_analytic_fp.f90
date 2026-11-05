!*****************************************************************************************
!>
!  Analytic test problem solved with Adams + fixed-point iteration.
!
!  FVODE port of the SUNDIALS example `examples/cvode/F2003_serial/cv_analytic_fp_f2003.f90`
!  (same problem, settings and output format; compare with its `.out` file).
!
!  The ODE is
!
!      dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t),   y(0) = 0,
!
!  for t in [0, 10] with lambda = -100, whose solution is y(t) = atan(t).
!  The problem is mildly stiff, so the fixed-point iteration often fails to
!  converge and CVODE reduces the step size.
!
!  The solution printed at t = 1, ..., 10 matches the upstream output. The
!  solver statistics (number of steps etc.) depend on floating-point rounding,
!  because many convergence tests are decided near their threshold: they
!  differ between platforms and compilers (also for the C version; the
!  upstream `.out` file reports 2065 steps, the C library built here takes
!  2160, and this program typically takes about 2000-2100).

    program cv_analytic_fp

    use fvode                      ! kinds and constants
    use fvode_core                 ! the CVODE solver
    use fvode_nonlinsol_fixedpoint ! fixed-point nonlinear solver
    use iso_fortran_env, only: output_unit_ => output_unit

    implicit none

    real(wp),parameter :: lambda = -100.0_wp  !! stiffness parameter of the ODE
    real(wp),parameter :: tstart = 0.0_wp  !! initial time
    real(wp),parameter :: tend   = 10.0_wp  !! final time
    real(wp),parameter :: dtout  = 1.0_wp  !! output time interval
    real(wp),parameter :: rtol   = 1.0e-6_wp  !! relative tolerance
    real(wp),parameter :: atol   = 1.0e-10_wp  !! absolute tolerance

    type(cvode_solver) :: cvode  !! the CVODE solver
    real(wp),dimension(1) :: y  !! solution vector
    real(wp) :: tout, tcur  !! next output time and time reached by the solver
    integer :: ierr, nout, outstep  !! error flag, number of outputs, and output counter

    tcur = tstart
    tout = tstart
    nout = ceiling(tend/dtout)
    y = 0.0_wp

    ! create and initialize the solver
    call cvode%create(CV_ADAMS, ierr);                  call check('create')
    call cvode%init(rhs, tstart, y, ierr);              call check('init')
    call cvode%ss_tolerances(rtol, atol, ierr);         call check('ss_tolerances')
    call cvode%set_nonlinear_solver(sun_nonlinsol_fixedpoint(size(y), 0), ierr)
    call check('set_nonlinear_solver')

    print *, '   '
    print *, 'Finished initialization, starting time steps'
    print *, '   '
    print *, '       t           y        '
    print *, '----------------------------'
    print '(2x,2(es12.5,1x))', tcur, y(1)
    do outstep = 1, nout
        tout = min(tout + dtout, tend)
        call cvode%solve(tout, y, tcur, CV_NORMAL, ierr);  call check('solve')
        print '(2x,2(es12.5,1x))', tcur, y(1)
    end do

    call print_stats()

    print *, ' '
    print *, ' All Solver Stats:'
    call cvode%print_all_stats(output_unit_, SUN_OUTPUTFORMAT_TABLE, ierr)

    call cvode%destroy()

    contains

        function rhs(me, t, y, ydot) result(retval)
        !! right-hand side function
        class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
        real(wp),intent(in) :: t  !! current value of the independent variable
        real(wp),dimension(:),intent(in) :: y  !! current state vector
        real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
        integer :: retval  !! return value: 0 on success
        ydot(1) = lambda*y(1) + 1.0_wp/(1.0_wp + t*t) - lambda*atan(t)
        retval = 0
        end function rhs

        subroutine print_stats()
        !! same output as CVodeStats in the SUNDIALS example
        integer(li) :: nsteps, nfevals, nlinsetups, netfails, nniters, nncfails  !! number of steps, f evaluations, linear solver setups, error test failures, nonlinear iterations and nonlinear convergence failures
        integer :: qlast, qcur  !! method order of the last step and of the next step
        real(wp) :: hinused, hlast, hcur, tcur_  !! initial, last and next step size, and current internal time
        call cvode%get_integrator_stats(nsteps, nfevals, nlinsetups, netfails, &
                                        qlast, qcur, hinused, hlast, hcur, tcur_)
        call cvode%get_nonlin_solv_stats(nniters, nncfails)
        print *, ' '
        print *, ' General Solver Stats:'
        print '(4x,A,i9)',    'Total internal steps taken =', nsteps
        print '(4x,A,i9)',    'Total rhs function calls   =', nfevals
        print '(4x,A,i9)',    'Num lin solver setup calls =', nlinsetups
        print '(4x,A,i9)',    'Num error test failures    =', netfails
        print '(4x,A,i9)',    'Last method order          =', qlast
        print '(4x,A,i9)',    'Next method order          =', qcur
        print '(4x,A,es12.5)','First internal step size   =', hinused
        print '(4x,A,es12.5)','Last internal step size    =', hlast
        print '(4x,A,es12.5)','Next internal step size    =', hcur
        print '(4x,A,es12.5)','Current internal time      =', tcur_
        print '(4x,A,i9)',    'Num nonlinear solver iters =', nniters
        print '(4x,A,i9)',    'Num nonlinear solver fails =', nncfails
        end subroutine print_stats

        subroutine check(name)
        !! stop if the last solver call returned an error
        character(len=*),intent(in) :: name  !! name of the call that was checked
        if (ierr < 0) then
            print *, 'Error in '//name//', ierr = ', ierr, '; halting'
            error stop 1
        end if
        end subroutine check

    end program cv_analytic_fp
!*****************************************************************************************

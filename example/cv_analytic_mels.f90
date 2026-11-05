!*****************************************************************************************
!>
!  Scalar ODE with an analytical solution, BDF + Newton with a custom
!  "matrix-embedded" linear solver.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvAnalytic_mels.c`
!  (same problem, settings and output; compare with `cvAnalytic_mels.out`).
!
!  The problem is
!
!      dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
!
!  for t in [0, 10] with y(0) = 0 and lambda = -100; the solution is
!  `y = atan(t)`. The linear systems `(1 - gamma*lambda) x = b` of the
!  Newton iteration are solved by a user-defined linear solver of type
!  `SUNLINEARSOLVER_MATRIX_EMBEDDED`, which gets gamma from
!  [[CVodeGetNonlinearSystemData]].

    module cv_analytic_mels_module

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: sun_matrix
    use fvode_linsol, only: sun_linear_solver, ls_problem

    implicit none

    private

    type,extends(cvode_solver),public :: analytic_solver
        !! the CVODE solver plus the problem data (C: user data)
        real(wp) :: lambda = -100.0_wp  !! stiffness parameter
    end type analytic_solver

    type,extends(sun_linear_solver),public :: matrix_embedded_ls
        !! custom matrix-embedded linear solver (C: `MatrixEmbeddedLS`)
    contains
        procedure :: get_type => mels_type  !! `SUNLINEARSOLVER_MATRIX_EMBEDDED`
        procedure :: setup => mels_setup  !! nothing to set up (C: no setup operation)
        procedure :: solve => mels_solve  !! (1 - gamma*lambda) x = b
    end type matrix_embedded_ls

    public :: f

    contains

    function f(me, t, y, ydot) result(retval)
    !! right-hand side f(t,y)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: ydot  !! dy/dt
    integer :: retval  !! 0 on success
    real(wp) :: u  !! current solution value
    select type (me)
    class is (analytic_solver)
        u = y(1)
        ydot(1) = me%lambda * u + 1.0_wp / (1.0_wp + t * t) - me%lambda * atan(t)
    end select
    retval = 0
    end function f

    pure function mels_type(me) result(itype)
    !! the solver type
    class(matrix_embedded_ls),intent(in) :: me  !! the solver
    integer :: itype  !! `SUNLINEARSOLVER_MATRIX_EMBEDDED`
    itype = SUNLINEARSOLVER_MATRIX_EMBEDDED
    end function mels_type

    function mels_setup(me, A, prob) result(flag)
    !! nothing to set up
    class(matrix_embedded_ls),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used
    class(ls_problem),intent(inout),optional :: prob  !! not used
    integer :: flag  !! `SUN_SUCCESS`
    flag = SUN_SUCCESS
    end function mels_setup

    function mels_solve(me, A, x, b, tol, prob) result(flag)
    !! solve `(1 - gamma*lambda) x = b`, with gamma from the integrator
    class(matrix_embedded_ls),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used
    real(wp),dimension(:),intent(inout) :: x  !! the solution
    real(wp),dimension(:),intent(in) :: b  !! the right-hand side
    real(wp),intent(in) :: tol  !! tolerance (not used)
    class(ls_problem),intent(inout),optional :: prob  !! the integrator
    integer :: flag  !! `SUN_SUCCESS`, or -1 without the integrator
    real(wp) :: gamma  !! the current gamma
    flag = -1
    if (.not. present(prob)) return
    select type (prob)
    class is (analytic_solver)
        ! retrieve implicit system data from CVode
        call prob%get_nonlinear_system_data(gamma=gamma)
        ! perform linear solve: (1-gamma*lambda)*x = b
        x(1) = b(1) / (1 - gamma * prob%lambda)
        flag = SUN_SUCCESS
    end select
    end function mels_solve

    end module cv_analytic_mels_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_analytic_mels

    use fvode
    use cv_analytic_mels_module

    implicit none

    real(wp),parameter :: T0     = 0.0_wp    !! initial time
    real(wp),parameter :: Tf     = 10.0_wp   !! final time
    real(wp),parameter :: dTout  = 1.0_wp    !! time between outputs
    real(wp),parameter :: reltol = 1.0e-6_wp !! relative tolerance
    real(wp),parameter :: abstol = 1.0e-10_wp  !! absolute tolerance

    type(analytic_solver) :: cvode  !! the CVODE solver with the problem data
    type(matrix_embedded_ls) :: LS  !! the custom linear solver
    real(wp),dimension(1) :: y  !! solution vector
    real(wp) :: t, tout  !! current and output times
    integer :: retval  !! return flag

    ! initial diagnostics output
    write(*,'(/a)') 'Analytical ODE test problem:'
    write(*,'(a,i0)') '   lambda = ', nint(cvode%lambda)
    write(*,'(a,es7.1)') '   reltol = ', reltol
    write(*,'(a,es7.1/)') '   abstol = ', abstol

    ! specify initial condition
    y = 0.0_wp

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, T0, y, retval)
    call cvode%ss_tolerances(reltol, abstol, retval)

    ! attach the custom matrix-embedded linear solver
    call cvode%set_linear_solver(LS, flag=retval)
    if (retval < 0) then
        write(*,'(/a,i0/)') 'SUNDIALS_ERROR: CVodeSetLinearSolver() failed with flag = ', retval
        stop 1
    end if

    ! in loop, call CVode, print results, and test for error
    t    = T0
    tout = T0 + dTout
    write(*,'(a)') '        t           u'
    write(*,'(a)') '   ---------------------'
    do while (Tf - t > 1.0e-15_wp)
        call cvode%solve(tout, y, t, CV_NORMAL, retval)
        write(*,'(2x,f10.6,2x,f10.6)') t, y(1)
        if (retval >= 0) then
            ! successful solve: update time
            tout = tout + dTout
            tout = min(tout, Tf)
        else
            ! unsuccessful solve: break
            write(*,'(a)') 'Solver failure, stopping integration'
            exit
        end if
    end do
    write(*,'(a)') '   ---------------------'

    ! get/print some final statistics on how the solve progressed
    write(*,'(/a)') 'Final Solver Statistics:'
    write(*,'(a,i0)') '   Internal solver steps = ', cvode%get_num_steps()
    write(*,'(a,i0)') '   Total RHS evals = ', cvode%get_num_rhs_evals()
    write(*,'(a,i0)') '   Total linear solver setups = ', cvode%get_num_lin_solv_setups()
    write(*,'(a,i0)') '   Total RHS evals for setting up the linear system = ', &
        cvode%get_num_lin_rhs_evals()
    write(*,'(a,i0)') '   Total number of Jacobian evaluations = ', cvode%get_num_jac_evals()
    write(*,'(a,i0)') '   Total number of Newton iterations = ', cvode%get_num_nonlin_solv_iters()
    write(*,'(a,i0)') '   Total number of linear solver convergence failures = ', &
        cvode%get_num_nonlin_solv_conv_fails()
    write(*,'(a,i0/)') '   Total number of error test failures = ', cvode%get_num_err_test_fails()

    ! check the solution error
    retval = check_ans(y, t, reltol, abstol)

    call cvode%destroy()

    if (retval /= 0) stop 1

    contains

        function check_ans(y, t, rtol, atol) result(passfail)
        !! check the computed solution against atan(t)
        real(wp),dimension(1),intent(in) :: y  !! the computed solution
        real(wp),intent(in) :: t  !! the time
        real(wp),intent(in) :: rtol  !! relative tolerance
        real(wp),intent(in) :: atol  !! absolute tolerance
        integer :: passfail  !! 0: pass, 1: fail
        real(wp) :: ans, err, ewt  !! answer, error, and error weight
        ans = atan(t)
        ewt = 1.0_wp / (rtol * abs(ans) + atol)
        err = ewt * abs(y(1) - ans)
        passfail = merge(0, 1, err < 1.0_wp)
        if (passfail /= 0) write(*,'(/a,g0/)') 'SUNDIALS_WARNING: check_ans error=', err
        end function check_ans

    end program cv_analytic_mels
!*****************************************************************************************

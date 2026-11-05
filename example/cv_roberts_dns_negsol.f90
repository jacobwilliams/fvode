!*****************************************************************************************
!>
!  Robertson chemical kinetics problem with a recoverable right-hand-side
!  failure on negative solution components: BDF + Newton + dense linear
!  solver with the difference-quotient Jacobian.
!
!  FVODE port of the SUNDIALS example
!  `examples/cvode/serial/cvRoberts_dns_negsol.c` (same problem, settings and
!  output; compare with `cvRoberts_dns_negsol.out`).
!
!  The problem is integrated twice, out to t = 4e12. The first time negative
!  solution components are ignored (and the solution eventually blows up);
!  the second time the right-hand side function returns a recoverable error
!  (1) when a component is negative, so that CVODE retries the step with a
!  smaller step size, which keeps the solution nonnegative.

    module cv_roberts_dns_negsol_module

    use fvode
    use fvode_core, only: cvode_solver

    implicit none

    private

    type,extends(cvode_solver),public :: negsol_solver
        !! the CVODE solver plus the problem data (C: user_data)
        logical :: check_negative = .false.  !! return a recoverable error for negative y_i?
    end type negsol_solver

    public :: f

    contains

!*****************************************************************************************
!>
!  Right-hand side function f(t, y); returns 1 (recoverable error) for a
!  negative component if `check_negative` is set.

    function f(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (a negsol_solver)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: ydot  !! f(t, y)
    integer :: retval  !! 0: success, 1: negative component found

    retval = 0
    select type (me)
    class is (negsol_solver)
        if (me%check_negative .and. (y(1) < 0 .or. y(2) < 0 .or. y(3) < 0)) then
            retval = 1
            return
        end if
    end select

    ydot(1) = -0.04_wp * y(1) + 1.0e4_wp * y(2) * y(3)
    ydot(3) = 3.0e7_wp * y(2) * y(2)
    ydot(2) = -ydot(1) - ydot(3)

    end function f
!*****************************************************************************************

    end module cv_roberts_dns_negsol_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_roberts_dns_negsol

    use fvode
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use cv_roberts_dns_negsol_module

    implicit none

    integer,parameter  :: NEQ   = 3            !! number of equations
    real(wp),parameter :: Y1    = 1.0_wp       !! initial y1
    real(wp),parameter :: Y2    = 0.0_wp       !! initial y2
    real(wp),parameter :: Y3    = 0.0_wp       !! initial y3
    real(wp),parameter :: RTOL  = 1.0e-4_wp    !! scalar relative tolerance
    real(wp),parameter :: ATOL1 = 1.0e-7_wp    !! absolute tolerance of y1
    real(wp),parameter :: ATOL2 = 1.0e-13_wp   !! absolute tolerance of y2
    real(wp),parameter :: ATOL3 = 1.0e-5_wp    !! absolute tolerance of y3
    real(wp),parameter :: T0    = 0.0_wp       !! initial time
    real(wp),parameter :: T1    = 0.4_wp       !! first output time
    real(wp),parameter :: TMULT = 10.0_wp      !! output time factor
    integer,parameter  :: NOUT  = 14           !! number of output times

    type(negsol_solver) :: cvode  !! the solver
    real(wp),dimension(NEQ) :: y  !! solution vector
    real(wp) :: t, tout  !! time reached, next output time
    integer :: retval  !! return flag

    y = [Y1, Y2, Y3]

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, T0, y, retval)
    call cvode%sv_tolerances(RTOL, [ATOL1, ATOL2, ATOL3], retval)
    call cvode%set_linear_solver(sun_linsol_dense(int(NEQ,ip)), &
                                 dense_matrix(int(NEQ,ip), int(NEQ,ip)), retval)

    ! first run: ignore negative solution components
    write(*,'(a/)') 'Ignore negative solution components'
    cvode%check_negative = .false.
    call run()
    call print_final_stats()

    ! second run: intercept negative solution components
    write(*,'(a/)') 'Intercept negative solution components'
    cvode%check_negative = .true.
    y = [Y1, Y2, Y3]
    call cvode%reinit(T0, y, retval)
    call run()
    call print_final_stats()

    call cvode%destroy()

    contains

!*****************************************************************************************
!>
!  Integrate to the NOUT output times and print the solution.

        subroutine run()
        integer :: iout  !! output counter
        tout = T1
        do iout = 1, NOUT
            call cvode%solve(tout, y, t, CV_NORMAL, retval)
            write(*,'(a,es10.4,a,es14.6,2x,es14.6,2x,es14.6)') 'At t = ', t, '      y =', y(1), y(2), y(3)
            tout = tout * TMULT
        end do
        end subroutine run

!*****************************************************************************************
!>
!  Print the final statistics (as in the C example).

        subroutine print_final_stats()
        write(*,'(/a)') 'Final Statistics:'
        write(*,'(a)') 'nst = '//lj(cvode%get_num_steps())//' nfe = '//lj(cvode%get_num_rhs_evals())// &
            ' nsetups = '//lj(cvode%get_num_lin_solv_setups())//' nfeLS = '// &
            lj(cvode%get_num_lin_rhs_evals())//' nje = '//itoa(cvode%get_num_jac_evals())
        write(*,'(a/)') 'nni = '//lj(cvode%get_num_nonlin_solv_iters())//' nnf = '// &
            lj(cvode%get_num_nonlin_solv_conv_fails())//' netf = '// &
            lj(cvode%get_num_err_test_fails())//'    ncfn = '//lj(cvode%get_num_step_solve_fails())
        end subroutine print_final_stats

!*****************************************************************************************
!>
!  An integer as a string (C: %ld).

        function itoa(i) result(s)
        integer(li),intent(in) :: i  !! the integer
        character(len=:),allocatable :: s  !! its decimal representation
        character(len=24) :: buf  !! write buffer
        write(buf,'(i0)') i
        s = trim(buf)
        end function itoa

!*****************************************************************************************
!>
!  An integer left-justified in 6 characters (C: %-6ld).

        function lj(i) result(s)
        integer(li),intent(in) :: i  !! the integer
        character(len=6) :: s  !! the padded string
        s = itoa(i)
        end function lj

    end program cv_roberts_dns_negsol
!*****************************************************************************************

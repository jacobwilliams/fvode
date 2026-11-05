!*****************************************************************************************
!>
!  Robertson chemical kinetics problem with the inequality constraints
!  y_i >= 0: BDF + Newton + dense linear solver with a user-supplied
!  Jacobian, root finding, and constraints.
!
!  FVODE port of the SUNDIALS example
!  `examples/cvode/serial/cvRoberts_dns_constraints.c` (same problem, settings
!  and output; compare with `cvRoberts_dns_constraints.out`).
!
!  The problem is the stiff system
!
!      dy1/dt = -0.04*y1 + 1.0e4*y2*y3
!      dy2/dt = +0.04*y1 - 1.0e4*y2*y3 - 3.0e7*y2^2
!      dy3/dt = 3.0e7*y2^2
!
!  on 0 <= t <= 4.0e10, with y = (1, 0, 0) at t = 0. The roots of
!  `g1 = y1 - 1e-4` and `g2 = y3 - 0.01` are located during the integration.
!  The constraint y_i >= 0 is imposed on all components.

    program cv_roberts_dns_constraints

    use fvode
    use fvode_core
    use fvode_matrix, only: sun_matrix, dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use iso_fortran_env, only: error_unit

    implicit none

    integer,parameter  :: NEQ   = 3            !! number of equations
    real(wp),parameter :: Y1    = 1.0_wp       !! initial y components
    real(wp),parameter :: Y2    = 0.0_wp       !! initial y2
    real(wp),parameter :: Y3    = 0.0_wp       !! initial y3
    real(wp),parameter :: RTOL  = 1.0e-4_wp    !! scalar relative tolerance
    real(wp),parameter :: ATOL1 = 1.0e-6_wp    !! vector absolute tolerance components
    real(wp),parameter :: ATOL2 = 1.0e-11_wp   !! absolute tolerance of y2
    real(wp),parameter :: ATOL3 = 1.0e-5_wp    !! absolute tolerance of y3
    real(wp),parameter :: T0    = 0.0_wp       !! initial time
    real(wp),parameter :: T1    = 0.4_wp       !! first output time
    real(wp),parameter :: TMULT = 10.0_wp      !! output time factor
    integer,parameter  :: NOUT  = 12           !! number of output times

    type(cvode_solver) :: cvode  !! the solver
    real(wp),dimension(NEQ) :: y  !! solution vector
    real(wp),dimension(NEQ) :: abstol  !! absolute tolerances
    real(wp) :: t, tout  !! time reached, next output time
    integer :: retval, iout  !! return flag, output counter
    integer :: retvalr  !! return flag of check_ans
    integer,dimension(2) :: rootsfound  !! root information

    ! initial conditions and tolerances
    y = [Y1, Y2, Y3]
    abstol = [ATOL1, ATOL2, ATOL3]

    ! create the solver (BDF), initialize it, and set the tolerances
    call cvode%create(CV_BDF, retval);                 if (check(retval, 'create')) stop 1
    call cvode%init(f, T0, y, retval);                 if (check(retval, 'init')) stop 1
    call cvode%sv_tolerances(RTOL, abstol, retval);    if (check(retval, 'sv_tolerances')) stop 1

    ! root function g with 2 components
    call cvode%root_init(2, g, retval);                if (check(retval, 'root_init')) stop 1

    ! dense matrix and linear solver, user Jacobian
    call cvode%set_linear_solver(sun_linsol_dense(int(NEQ,ip)), dense_matrix(int(NEQ,ip), int(NEQ,ip)), retval)
    if (check(retval, 'set_linear_solver')) stop 1
    call cvode%set_jac_fn(jac, retval);                if (check(retval, 'set_jac_fn')) stop 1

    ! constraints: all components nonnegative
    call cvode%set_constraints([1.0_wp, 1.0_wp, 1.0_wp], retval)
    if (check(retval, 'set_constraints')) stop 1

    ! in loop, call solve, print results, and test for error.
    ! Break out of loop when NOUT preset output times have been reached.
    write(*,'(a)') ' '
    write(*,'(a/)') '3-species kinetics problem'

    iout = 0
    tout = T1
    do
        call cvode%solve(tout, y, t, CV_NORMAL, retval)
        write(*,'(a,es10.4,a,es14.6,2x,es14.6,2x,es14.6)') 'At t = ', t, '      y =', y(1), y(2), y(3)

        if (retval == CV_ROOT_RETURN) then
            call cvode%get_root_info(rootsfound)
            write(*,'(a,i3,1x,i3)') '    rootsfound[] = ', rootsfound(1), rootsfound(2)
        end if

        if (retval == CV_SUCCESS) then
            iout = iout + 1
            tout = tout * TMULT
        end if

        if (iout == NOUT) exit
    end do

    ! print some final statistics
    call print_final_stats()

    ! check the solution error
    retvalr = check_ans(y, RTOL, abstol)

    call cvode%destroy()

    contains

!*****************************************************************************************
!>
!  Right-hand side function f(t, y).

        function f(me, t, y, ydot) result(retval)
        class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        real(wp),intent(in) :: t  !! current time
        real(wp),dimension(:),intent(in) :: y  !! current solution
        real(wp),dimension(:),intent(out) :: ydot  !! f(t, y)
        integer :: retval  !! 0 (success)
        ydot(1) = -0.04_wp * y(1) + 1.0e4_wp * y(2) * y(3)
        ydot(3) = 3.0e7_wp * y(2) * y(2)
        ydot(2) = -ydot(1) - ydot(3)
        retval = 0
        end function f

!*****************************************************************************************
!>
!  Root function: g1 = y1 - 0.0001, g2 = y3 - 0.01.

        function g(me, t, y, gout) result(retval)
        class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        real(wp),intent(in) :: t  !! current time
        real(wp),dimension(:),intent(in) :: y  !! current solution
        real(wp),dimension(:),intent(out) :: gout  !! root function values
        integer :: retval  !! 0 (success)
        gout(1) = y(1) - 0.0001_wp
        gout(2) = y(3) - 0.01_wp
        retval = 0
        end function g

!*****************************************************************************************
!>
!  Jacobian J = df/dy.

        function jac(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
        class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        real(wp),intent(in) :: t  !! current time
        real(wp),dimension(:),intent(in) :: y  !! current solution
        real(wp),dimension(:),intent(in) :: fy  !! f(t, y) (not used)
        class(sun_matrix),intent(inout) :: J  !! the Jacobian (a dense matrix)
        real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors (not used)
        integer :: retval  !! 0 (success)
        select type (J)
        type is (dense_matrix)
            J%data(1,1) = -0.04_wp
            J%data(1,2) = 1.0e4_wp * y(3)
            J%data(1,3) = 1.0e4_wp * y(2)
            J%data(2,1) = 0.04_wp
            J%data(2,2) = -1.0e4_wp * y(3) - 6.0e7_wp * y(2)
            J%data(2,3) = -1.0e4_wp * y(2)
            J%data(3,1) = 0.0_wp
            J%data(3,2) = 6.0e7_wp * y(2)
            J%data(3,3) = 0.0_wp
        end select
        retval = 0
        end function jac

!*****************************************************************************************
!>
!  Report a failed call (negative flag) like `check_retval` in the C example.
!  Returns true on failure.

        logical function check(flag, name)
        integer,intent(in) :: flag  !! return flag
        character(len=*),intent(in) :: name  !! name of the called procedure
        check = flag < 0
        if (check) write(error_unit,'(/a,a,a,i0/)') &
            'SUNDIALS_ERROR: ', name, '() failed with retval = ', flag
        end function check

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
            lj(cvode%get_num_err_test_fails())//'    ncfn = '// &
            lj(cvode%get_num_step_solve_fails())//'  nge = '//itoa(cvode%get_num_g_evals())
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

!*****************************************************************************************
!>
!  Compare the solution at t = 4e10 with a reference solution (computed with
!  a tight tolerance). Returns 0 on pass, 1 on fail.

        function check_ans(y, rtol, atol) result(passfail)
        real(wp),dimension(:),intent(in) :: y  !! solution at t = 4e10
        real(wp),intent(in) :: rtol  !! relative tolerance
        real(wp),dimension(:),intent(in) :: atol  !! absolute tolerances
        integer :: passfail  !! 0: pass, 1: fail
        real(wp),dimension(NEQ) :: ref  !! reference solution
        real(wp),dimension(NEQ) :: ewt  !! error weight vector
        real(wp) :: err  !! WRMS error
        ref = [5.2083495894337328e-08_wp, 2.0833399429795671e-13_wp, 9.9999994791629776e-01_wp]
        ewt = rtol * abs(ref) + 10.0_wp * atol
        if (minval(ewt) <= 0.0_wp) then
            write(error_unit,'(/a/)') 'SUNDIALS_ERROR: check_ans failed - ewt <= 0'
            passfail = -1
            return
        end if
        ewt = 1.0_wp / ewt
        err = sqrt(sum(((y - ref) * ewt)**2) / NEQ)
        passfail = merge(0, 1, err < 1.0_wp)
        if (passfail /= 0) write(*,'(/a,g0/)') 'SUNDIALS_WARNING: check_ans error=', err
        end function check_ans

    end program cv_roberts_dns_constraints
!*****************************************************************************************

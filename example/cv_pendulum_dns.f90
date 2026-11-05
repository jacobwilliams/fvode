!*****************************************************************************************
!>
!  Simple pendulum in Cartesian coordinates, BDF + Newton + dense linear
!  solver, with and without projection onto the constraint manifold.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvPendulum_dns.c`
!  (same problem, settings and output, with the default command-line
!  options; compare with `cvPendulum_dns.out`).
!
!  The pendulum bob has mass 1 and hangs from the origin on a rod of
!  length 1:
!
!      x'  = vx
!      y'  = vy
!      vx' = -x * T
!      vy' = -y * T - g
!
!  with the constraints `x^2 + y^2 - 1 = 0` and `x*vx + y*vy = 0`, where T
!  is the tension in the rod and g is chosen so that the period is 2. The
!  initial condition is x = 1, y = 0, vx = 0, vy = 0. A reference solution
!  is computed from the angle formulation `theta'' = -g cos(theta)`. The
!  Cartesian system is solved to tf = 30 for five tolerances, with and
!  without projection, and the errors at tf, the constraint violation and
!  the integrator statistics are printed. The solutions are also written to
!  the files `cvPendulum_dns_*.txt`.

    module cv_pendulum_dns_module

    use fvode
    use fvode_core, only: cvode_solver

    implicit none

    private

    real(wp),parameter,public :: GRAV = 13.750371636040745654980191559621114395801712_wp  !! acceleration due to gravity (period 2)

    public :: f, fref, proj

    contains

    function fref(me, t, y, ydot) result(retval)
    !! right-hand side of the reference (angle) formulation
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! theta and theta'
    real(wp),dimension(:),intent(out) :: ydot  !! theta' and theta''
    integer :: retval  !! 0 on success
    ydot(1) = y(2)                  ! theta'
    ydot(2) = -GRAV * cos(y(1))     ! -g * cos(theta)
    retval = 0
    end function fref

    function f(me, t, y, ydot) result(retval)
    !! right-hand side of the Cartesian system
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! x, y, x', y'
    real(wp),dimension(:),intent(out) :: ydot  !! the time derivatives
    integer :: retval  !! 0 on success
    real(wp) :: x, yy, xd, yd, tmp  !! positions, velocities and the tension
    x  = y(1)
    yy = y(2)
    xd = y(3)
    yd = y(4)
    ! compute tension
    tmp = xd * xd + yd * yd - GRAV * yy
    ! compute RHS
    ydot(1) = xd
    ydot(2) = yd
    ydot(3) = -x * tmp
    ydot(4) = -yy * tmp - GRAV
    retval = 0
    end function f

    function proj(me, t, ycur, corr, epsProj, err) result(retval)
    !! project the positions to the unit circle and the velocities (and the
    !! error) to its tangent
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(inout) :: ycur  !! the solution (only the correction is returned)
    real(wp),dimension(:),intent(out) :: corr  !! the projection correction
    real(wp),intent(in) :: epsProj  !! projection tolerance (not used: exact projection)
    real(wp),dimension(:),intent(inout),optional :: err  !! the error estimate, projected
    integer :: retval  !! 0 on success
    real(wp) :: x, y, x_new, y_new  !! positions before and after projection
    real(wp) :: xd, yd, xd_new, yd_new  !! velocities before and after projection
    real(wp) :: e1, e2, e3, e4  !! error components
    real(wp) :: e1_new, e2_new, e3_new, e4_new  !! projected error components
    real(wp) :: R  !! distance from the origin

    ! extract current solution
    x  = ycur(1)
    y  = ycur(2)
    xd = ycur(3)
    yd = ycur(4)

    ! project positions
    R = sqrt(x * x + y * y)
    x_new = x / R
    y_new = y / R

    ! project velocities
    xd_new = xd * y_new * y_new - yd * x_new * y_new
    yd_new = -xd * x_new * y_new + yd * x_new * x_new

    ! return position and velocity corrections
    corr(1) = x_new - x
    corr(2) = y_new - y
    corr(3) = xd_new - xd
    corr(4) = yd_new - yd

    ! project error P * err
    if (present(err)) then
        e1 = err(1)
        e2 = err(2)
        e3 = err(3)
        e4 = err(4)
        e1_new = y_new * y_new * e1 - x_new * y_new * e2
        e2_new = -x_new * y_new * e1 + x_new * x_new * e2
        e3_new = y_new * y_new * e3 - x_new * y_new * e4
        e4_new = -x_new * y_new * e3 + x_new * x_new * e4
        err(1) = e1_new
        err(2) = e2_new
        err(3) = e3_new
        err(4) = e4_new
    end if

    retval = 0
    end function proj

    end module cv_pendulum_dns_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_pendulum_dns

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use cv_pendulum_dns_module

    implicit none

    integer,parameter :: nout = 1          !! number of outputs
    real(wp),parameter :: tf = 30.0_wp     !! final integration time
    logical,parameter :: projerr = .true.  !! enable error projection

    type(cvode_solver) :: cvode  !! the CVODE solver
    real(wp) :: rtol, atol  !! relative and absolute tolerances
    real(wp),dimension(4) :: yy0  !! initial condition
    real(wp),dimension(4) :: yref  !! reference solution at tf
    integer :: i, retval  !! run counter and return flag

    rtol = 1.0e-5_wp
    atol = 1.0e-5_wp

    ! compute reference solution
    retval = ref_sol(tf, yref, nout)
    if (retval /= 0) stop 1

    ! set the initial condition values
    yy0 = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]  ! x, y, xd, yd

    ! create and initialize CVODE
    call cvode%create(CV_BDF, retval)
    call cvode%init(f, 0.0_wp, yy0, retval)

    ! attach a dense matrix and linear solver
    call cvode%set_linear_solver(sun_linsol_dense(4_ip), dense_matrix(4_ip, 4_ip), retval)

    ! set a user-supplied projection function
    call cvode%set_proj_fn(proj, retval)
    if (check_retval(retval, 'CVodeSetProjFn')) stop 1

    ! set maximum number of steps between outputs
    call cvode%set_max_num_steps(50000_li, retval)

    ! compute the solution with various tolerances
    do i = 0, 4

        ! output tolerance and output header for this run
        write(*,'(//a,es8.2,a,es8.2)') 'rtol = ', rtol, ', atol = ', atol
        write(*,'(a)') 'Project    x         y         x''        y''     |     g      |    '// &
                       'nst     rhs eval    setups (J eval)  |   cf   ef'

        ! compute solution with projection
        retval = get_sol(rtol, atol, .true., projerr)
        if (retval /= 0) stop 1

        ! compute solution without projection
        retval = get_sol(rtol, atol, .false., .false.)
        if (retval /= 0) stop 1

        ! reduce tolerance for next run
        rtol = rtol / 10.0_wp
        atol = atol / 10.0_wp

    end do

    call cvode%destroy()

    contains

        function get_sol(rtol, atol, proj_on, projerr) result(retval)
        !! integrate the Cartesian system to tf and print the errors and stats
        real(wp),intent(in) :: rtol  !! relative tolerance
        real(wp),intent(in) :: atol  !! absolute tolerance
        logical,intent(in) :: proj_on  !! use projection?
        logical,intent(in) :: projerr  !! project the error estimate?
        integer :: retval  !! 0 on success
        real(wp),dimension(4) :: yy  !! solution vector
        real(wp) :: dtout, tout, t  !! output interval, output time and return time
        real(wp) :: x, y, xd, yd, g  !! errors and the constraint violation
        integer :: out, iunit  !! output counter and file unit
        character(len=100) :: outname  !! output file name

        ! enable or disable projection
        if (proj_on) then
            write(*,'(a)',advance='no') '  YES   '
            call cvode%set_proj_frequency(1_li, retval)
            if (check_retval(retval, 'CVodeSetProjFrequency')) return
            ! enable or disable error projection
            call cvode%set_proj_err_est(projerr, retval)
            if (check_retval(retval, 'CVodeSetProjErrEst')) return
        else
            call cvode%set_proj_frequency(0_li, retval)
            if (check_retval(retval, 'CVodeSetProjFrequency')) return
            write(*,'(a)',advance='no') '  NO    '
        end if

        ! copy initial condition into solution vector
        yy = yy0

        ! reinitialize CVODE for this run
        call cvode%reinit(0.0_wp, yy0, retval)
        if (check_retval(retval, 'CVodeReInit')) return

        ! set integration tolerances for this run
        call cvode%ss_tolerances(rtol, atol, retval)
        if (check_retval(retval, 'CVodeSStolerances')) return

        ! open output file
        if (proj_on) then
            outname = 'cvPendulum_dns_rtol_'//e2(rtol)//'_atol_'//e2(atol)//'_proj.txt'
        else
            outname = 'cvPendulum_dns_rtol_'//e2(rtol)//'_atol_'//e2(atol)//'.txt'
        end if
        open(newunit=iunit, file=trim(outname), status='replace', action='write')

        ! output initial condition
        call write_line(iunit, 0.0_wp, yy)

        ! integrate to tf and periodically output the solution
        dtout = tf / nout
        tout  = dtout
        do out = 0, nout - 1
            ! set stop time (do not interpolate output)
            call cvode%set_stop_time(tout, retval)
            if (check_retval(retval, 'CVodeSetStopTime')) then
                close(iunit)
                return
            end if
            ! integrate to tout
            call cvode%solve(tout, yy, t, CV_NORMAL, retval)
            if (check_retval(retval, 'CVode')) then
                close(iunit)
                return
            end if
            ! write output
            call write_line(iunit, t, yy)
            ! update output time
            if (out < nout - 1) then
                tout = tout + dtout
            else
                tout = tf
            end if
        end do

        close(iunit)

        ! compute the constraint violation
        g = abs(yy(1) * yy(1) + yy(2) * yy(2) - 1.0_wp)

        ! compute the absolute error compared to the reference solution
        yy = abs(yy - yref)
        x  = yy(1)
        y  = yy(2)
        xd = yy(3)
        yd = yy(4)

        ! output errors
        write(*,'(es8.2,2x,es8.2,2x,es8.2,2x,es8.2,a,es8.2,a)',advance='no') &
            x, y, xd, yd, '  |  ', g, '  |'

        ! output stats
        write(*,'(1x,i6,3x,i6,"+",a4,5x,i4," (",i3,")",5x,"|  ",i3,2x,i3)') &
            cvode%get_num_steps(), cvode%get_num_rhs_evals(), &
            lj4(cvode%get_num_lin_rhs_evals()), cvode%get_num_lin_solv_setups(), &
            cvode%get_num_jac_evals(), cvode%get_num_nonlin_solv_conv_fails(), &
            cvode%get_num_err_test_fails()

        retval = 0
        end function get_sol

        function ref_sol(tf, yref, nout) result(retval)
        !! integrate the reference (angle) system to tf to tolerance 1e-14 and
        !! return the Cartesian state
        real(wp),intent(in) :: tf  !! final time
        real(wp),dimension(4),intent(out) :: yref  !! reference solution at tf
        integer,intent(in) :: nout  !! number of outputs
        integer :: retval  !! 0 on success
        type(cvode_solver) :: cv  !! the CVODE solver for the reference system
        real(wp),dimension(2) :: yy  !! theta and theta'
        real(wp) :: dtout, tout, t, th, thd  !! output interval, output and return time, theta and theta'
        real(wp),parameter :: tol = 1.0e-14_wp  !! integration tolerance
        integer :: out, iunit  !! output counter and file unit

        yy = 0.0_wp  ! theta, theta'

        call cv%create(CV_BDF, retval)
        call cv%init(fref, 0.0_wp, yy, retval)
        call cv%ss_tolerances(tol, tol, retval)
        call cv%set_linear_solver(sun_linsol_dense(2_ip), dense_matrix(2_ip, 2_ip), retval)
        call cv%set_max_num_steps(100000_li, retval)
        call cv%set_stop_time(tf, retval)

        open(newunit=iunit, file='cvPendulum_dns_ref.txt', status='replace', action='write')

        ! output initial condition
        th  = yy(1)
        thd = yy(2)
        call write_line(iunit, 0.0_wp, [cos(th), sin(th), -thd * sin(th), thd * cos(th)])

        ! integrate to tf and periodically output the solution
        dtout = tf / nout
        tout  = dtout
        do out = 0, nout - 1
            ! set stop time (do not interpolate output)
            call cv%set_stop_time(tout, retval)
            if (check_retval(retval, 'CVodeSetStopTime')) exit
            ! integrate to tout
            call cv%solve(tf, yy, t, CV_NORMAL, retval)
            if (check_retval(retval, 'CVode')) exit
            ! write output
            th  = yy(1)
            thd = yy(2)
            call write_line(iunit, t, [cos(th), sin(th), -thd * sin(th), thd * cos(th)])
            ! update output time
            if (out < nout - 1) then
                tout = tout + dtout
            else
                tout = tf
            end if
        end do

        close(iunit)
        if (retval < 0) then
            call cv%destroy()
            return
        end if

        ! save reference solution
        yref = [cos(th), sin(th), -thd * sin(th), thd * cos(th)]

        call cv%destroy()
        retval = 0
        end function ref_sol

        subroutine write_line(iunit, t, y)
        !! write t and the four state components (C: `%24.16e` each)
        integer,intent(in) :: iunit  !! file unit
        real(wp),intent(in) :: t  !! time
        real(wp),dimension(4),intent(in) :: y  !! the state
        character(len=124) :: buf  !! the line
        write(buf,'(es24.16,4(1x,es24.16))') t, y
        write(iunit,'(a)') lower_e(trim(buf))
        end subroutine write_line

        function e2(x) result(s)
        !! C's `%03.2e`, e.g. `1.00e-05`
        real(wp),intent(in) :: x  !! the value
        character(len=:),allocatable :: s  !! the formatted value
        character(len=16) :: buf  !! write buffer
        write(buf,'(es8.2)') x
        s = lower_e(trim(adjustl(buf)))
        end function e2

        function lower_e(str) result(s)
        !! replace the exponent letter `E` by `e` (as printed by C)
        character(len=*),intent(in) :: str  !! the string
        character(len=len(str)) :: s  !! the string with lowercase exponents
        integer :: k  !! character index
        s = str
        do k = 1, len(s)
            if (s(k:k) == 'E') s(k:k) = 'e'
        end do
        end function lower_e

        function lj4(i) result(s)
        !! integer left-justified in 4 characters (C: `%-4ld`)
        integer(li),intent(in) :: i  !! the integer
        character(len=4) :: s  !! the formatted integer
        write(s,'(i0)') i
        s = adjustl(s)
        end function lj4

        logical function check_retval(retval, funcname)
        !! report a negative return flag (C: `check_retval(&retval, funcname, 1)`)
        integer,intent(in) :: retval  !! the return flag
        character(len=*),intent(in) :: funcname  !! name of the called function
        check_retval = retval < 0
        if (check_retval) write(*,'(/3a,i0/)') 'SUNDIALS_ERROR: ', funcname, '() failed with retval = ', retval
        end function check_retval

    end program cv_pendulum_dns
!*****************************************************************************************

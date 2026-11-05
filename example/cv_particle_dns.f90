!*****************************************************************************************
!>
!  Particle moving counterclockwise on the unit circle, BDF + Newton +
!  dense linear solver with a user Jacobian, with projection onto the
!  circle.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvParticle_dns.c`
!  (same problem, settings, command-line options and output; compare with
!  `cvParticle_dns.out`).
!
!  The system is
!
!      x' = -alpha * y
!      y' =  alpha * x
!
!  subject to the constraint `x^2 + y^2 - 1 = 0`, with x = 1, y = 0 at
!  t = 0 and the analytic solution `x = cos(alpha t)`, `y = sin(alpha t)`.
!
!  Command-line options: `--alpha <vel>`, `--orbits <n>`, `--rtol <rtol>`,
!  `--atol <atol>`, `--proj <1|0>`, `--projerr <1|0>`, `--nout <n>`
!  (outputs per orbit, written to `cvParticle_solution.txt` and
!  `cvParticle_error.txt`), `--tstop` (stop at the output times) and
!  `--help`.

    module cv_particle_dns_module

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: sun_matrix, dense_matrix

    implicit none

    private

    real(wp),parameter,public :: PI = 3.141592653589793238462643383279502884197169_wp  !! pi

    type,extends(cvode_solver),public :: particle_solver
        !! the CVODE solver plus the problem data (C: UserData)
        real(wp) :: alpha = 1.0_wp  !! particle velocity
        integer :: orbits = 100  !! number of orbits
        real(wp) :: torbit = 2.0_wp * PI  !! orbit time
        real(wp) :: rtol = 1.0e-4_wp  !! relative tolerance
        real(wp) :: atol = 1.0e-9_wp  !! absolute tolerance
        integer :: proj = 1  !! enable/disable solution projection
        integer :: projerr = 0  !! enable/disable error projection
        integer :: tstop_mode = 0  !! use tstop mode (stop at the output times)
        integer :: nout = 0  !! number of outputs per orbit
    end type particle_solver

    public :: f, jac, proj_fn

    contains

    function f(me, t, y, ydot) result(retval)
    !! right-hand side, y' = f(t,y)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! x and y
    real(wp),dimension(:),intent(out) :: ydot  !! x' and y'
    integer :: retval  !! 0 on success
    select type (me)
    class is (particle_solver)
        ydot(1) = -(me%alpha) * y(2)
        ydot(2) = (me%alpha) * y(1)
    end select
    retval = 0
    end function f

    function jac(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
    !! the Jacobian df/dy
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t,y)
    class(sun_matrix),intent(inout) :: J  !! the Jacobian (filled in here)
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors (not used)
    integer :: retval  !! 0 on success
    select type (me)
    class is (particle_solver)
        select type (J)
        type is (dense_matrix)
            J%data(1,1) = 0.0_wp
            J%data(2,1) = (me%alpha)
            J%data(1,2) = -(me%alpha)
            J%data(2,2) = 0.0_wp
        end select
    end select
    retval = 0
    end function jac

    function proj_fn(me, t, ycur, corr, epsProj, err) result(retval)
    !! project the solution (and the error) onto the unit circle
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(inout) :: ycur  !! the solution (only the correction is returned)
    real(wp),dimension(:),intent(out) :: corr  !! the projection correction
    real(wp),intent(in) :: epsProj  !! projection tolerance (not used: exact projection)
    real(wp),dimension(:),intent(inout),optional :: err  !! the error estimate, projected
    integer :: retval  !! 0 on success
    real(wp) :: x, y, xp, yp, r  !! solution, projected solution and radius
    real(wp) :: errxp, erryp  !! projected error
    x = ycur(1)
    y = ycur(2)
    ! project onto the unit circle
    r = sqrt(x * x + y * y)
    xp = x / r
    yp = y / r
    ! correction to the unprojected solution
    corr(1) = xp - x
    corr(2) = yp - y
    ! project the error
    if (present(err)) then
        errxp = err(1) * yp * yp - err(2) * xp * yp
        erryp = -err(1) * xp * yp + err(2) * xp * xp
        err(1) = errxp
        err(2) = erryp
    end if
    retval = 0
    end function proj_fn

    end module cv_particle_dns_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_particle_dns

    use fvode
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use cv_particle_dns_module

    implicit none

    type(particle_solver) :: cvode  !! the CVODE solver with the problem data
    real(wp),dimension(2) :: y  !! solution vector
    real(wp),dimension(2) :: e  !! error vector
    real(wp) :: t, dtout, tout, ec  !! current time, output spacing, next output time, constraint error
    integer :: retval, out, totalout  !! return flag, output counter, number of outputs
    integer :: yfid, efid  !! solution and error output files

    t = 0.0_wp
    ec = 0.0_wp

    ! initialize user data
    if (init_user_data() /= 0) stop 1

    ! set initial condition and error
    y = [1.0_wp, 0.0_wp]
    e = 0.0_wp

    ! create and initialize CVODE
    call cvode%create(CV_BDF, retval)
    call cvode%init(f, t, y, retval)
    call cvode%ss_tolerances(cvode%rtol, cvode%atol, retval)

    ! attach the matrix and linear solver, and a user Jacobian
    call cvode%set_linear_solver(sun_linsol_dense(2_ip), dense_matrix(2_ip, 2_ip), retval)
    call cvode%set_jac_fn(jac, retval)

    ! set a user-supplied projection function
    if (cvode%proj /= 0) then
        call cvode%set_proj_fn(proj_fn, retval)
        if (check_retval(retval, 'CVodeSetProjFn')) stop 1
        call cvode%set_proj_err_est(cvode%projerr /= 0, retval)
        if (check_retval(retval, 'CVodeSetProjErrEst')) stop 1
    end if

    ! set max steps between outputs
    call cvode%set_max_num_steps(100000_li, retval)

    ! output problem setup
    call print_user_data()

    ! output initial condition
    write(*,'(/a)') '     t            x              y             err x          err y       err constr'
    call write_screen(t, y, e, ec)

    if (cvode%nout > 0) then
        open(newunit=yfid, file='cvParticle_solution.txt', status='replace', action='write')
        open(newunit=efid, file='cvParticle_error.txt', status='replace', action='write')
        call write_files(t, y, e, ec)
    end if

    ! integrate in time and periodically output the solution and error
    if (cvode%nout > 0) then
        totalout = cvode%orbits * cvode%nout
        dtout    = cvode%torbit / cvode%nout
    else
        totalout = 1
        dtout    = cvode%torbit * cvode%orbits
    end if

    tout = dtout
    do out = 0, totalout - 1

        ! stop at output time (do not interpolate output)
        if (cvode%tstop_mode /= 0 .or. cvode%nout == 0) then
            call cvode%set_stop_time(tout, retval)
            if (check_retval(retval, 'CVodeSetStopTime')) stop 1
        end if

        ! advance in time
        call cvode%solve(tout, y, t, CV_NORMAL, retval)
        if (check_retval(retval, 'CVode')) exit

        ! output solution and error
        if (cvode%nout > 0) then
            call compute_error(t, y, e, ec)
            call write_files(t, y, e, ec)
        end if

        ! update output time
        if (out < totalout - 1) then
            tout = tout + dtout
        else
            tout = cvode%torbit * cvode%orbits
        end if

    end do

    ! close output files
    if (cvode%nout > 0) then
        close(yfid)
        close(efid)
    end if

    ! output final solution and error to screen
    call compute_error(t, y, e, ec)
    call write_screen(t, y, e, ec)

    ! print some final statistics
    write(*,'(/a)') 'Integration Statistics:'
    write(*,'(2a)') 'Number of steps taken = ', lj(cvode%get_num_steps())
    write(*,'(2a)') 'Number of function evaluations = ', lj(cvode%get_num_rhs_evals())
    write(*,'(2a)') 'Number of linear solver setups = ', lj(cvode%get_num_lin_solv_setups())
    write(*,'(2a)') 'Number of Jacobian evaluations = ', lj(cvode%get_num_jac_evals())
    write(*,'(2a)') 'Number of nonlinear solver iterations = ', lj(cvode%get_num_nonlin_solv_iters())
    write(*,'(2a)') 'Number of convergence failures = ', lj(cvode%get_num_nonlin_solv_conv_fails())
    write(*,'(2a)') 'Number of error test failures = ', lj(cvode%get_num_err_test_fails())

    call cvode%destroy()

    contains

        function init_user_data() result(retval)
        !! set the problem data from the defaults and the command line
        !! (C: `InitUserData`)
        integer :: retval  !! 0 on success, -1 for --help or an invalid input
        integer :: arg_idx  !! argument index
        character(len=256) :: arg  !! command-line argument
        retval = 0
        arg_idx = 1
        do while (arg_idx <= command_argument_count())
            call get_command_argument(arg_idx, arg)
            select case (trim(arg))
            case ('--alpha')
                cvode%alpha = real_arg(arg_idx + 1)
                cvode%torbit = (2.0_wp * PI) / cvode%alpha
                arg_idx = arg_idx + 2
            case ('--orbits')
                cvode%orbits = int_arg(arg_idx + 1)
                arg_idx = arg_idx + 2
            case ('--rtol')
                cvode%rtol = real_arg(arg_idx + 1)
                arg_idx = arg_idx + 2
            case ('--atol')
                cvode%atol = real_arg(arg_idx + 1)
                arg_idx = arg_idx + 2
            case ('--proj')
                cvode%proj = int_arg(arg_idx + 1)
                arg_idx = arg_idx + 2
            case ('--projerr')
                cvode%projerr = int_arg(arg_idx + 1)
                arg_idx = arg_idx + 2
            case ('--nout')
                cvode%nout = int_arg(arg_idx + 1)
                arg_idx = arg_idx + 2
            case ('--tstop')
                cvode%tstop_mode = 1
                arg_idx = arg_idx + 1
            case ('--help')
                call input_help()
                retval = -1
                return
            case default
                write(*,'(2a)',advance='no') 'ERROR: Invalid input ', trim(arg)
                call input_help()
                retval = -1
                return
            end select
        end do
        ! if projection is disabled then disable error projection
        if (cvode%proj == 0) cvode%projerr = 0
        end function init_user_data

        function real_arg(i) result(x)
        !! command-line argument i as a real
        integer,intent(in) :: i  !! argument index
        real(wp) :: x  !! the value
        character(len=256) :: arg  !! the argument
        call get_command_argument(i, arg)
        read(arg, *) x
        end function real_arg

        function int_arg(i) result(n)
        !! command-line argument i as an integer
        integer,intent(in) :: i  !! argument index
        integer :: n  !! the value
        character(len=256) :: arg  !! the argument
        call get_command_argument(i, arg)
        read(arg, *) n
        end function int_arg

        subroutine input_help()
        !! print the command-line options
        write(*,'(/a)') 'Command line options:'
        write(*,'(a)') '  --alpha <vel>      : particle velocity'
        write(*,'(a)') '  --orbits <orbits>  : number of orbits to perform'
        write(*,'(a)') '  --rtol <rtol>      : relative tolerance'
        write(*,'(a)') '  --atol <atol>      : absolute tolerance'
        write(*,'(a)') '  --proj <1 or 0>    : enable (1) / disable (0) projection'
        write(*,'(a)') '  --projerr <1 or 0> : enable (1) / disable (0) error projection'
        write(*,'(a)') '  --nout <nout>      : outputs per period'
        write(*,'(a)') '  --tstop            : stop at output time (do not interpolate)'
        end subroutine input_help

        subroutine print_user_data()
        !! print the problem setup
        write(*,'(/a)') 'Particle traveling on the unit circle example'
        write(*,'(a)') '---------------------------------------------'
        write(*,'(a,es10.4)') 'alpha      = ', cvode%alpha
        write(*,'(a,i0)') 'num orbits = ', cvode%orbits
        write(*,'(a)') '---------------------------------------------'
        write(*,'(2a)') 'rtol       = ', cg(cvode%rtol)
        write(*,'(2a)') 'atol       = ', cg(cvode%atol)
        write(*,'(a,i0)') 'proj sol   = ', cvode%proj
        write(*,'(a,i0)') 'proj err   = ', cvode%projerr
        write(*,'(a,i0)') 'nout       = ', cvode%nout
        write(*,'(a,i0)') 'tstop      = ', cvode%tstop_mode
        write(*,'(a)') '---------------------------------------------'
        end subroutine print_user_data

        subroutine compute_error(t, y, e, ec)
        !! the error in the solution and the constraint (C: `ComputeError`)
        real(wp),intent(in) :: t  !! current time
        real(wp),dimension(2),intent(in) :: y  !! the solution
        real(wp),dimension(2),intent(out) :: e  !! the solution error
        real(wp),intent(out) :: ec  !! the constraint error
        ! solution error
        e = [cos((cvode%alpha) * t), sin((cvode%alpha) * t)]
        e = y - e
        ! constraint error
        ec = y(1) * y(1) + y(2) * y(2) - 1.0_wp
        end subroutine compute_error

        subroutine write_screen(t, y, e, ec)
        !! output solution and error to screen
        real(wp),intent(in) :: t  !! current time
        real(wp),dimension(2),intent(in) :: y  !! the solution
        real(wp),dimension(2),intent(in) :: e  !! the solution error
        real(wp),intent(in) :: ec  !! the constraint error
        write(*,'(es10.4,5(1x,es14.6))') t, y(1), y(2), e(1), e(2), ec
        end subroutine write_screen

        subroutine write_files(t, y, e, ec)
        !! output solution and error to disk (C: `%24.16e`)
        real(wp),intent(in) :: t  !! current time
        real(wp),dimension(2),intent(in) :: y  !! the solution
        real(wp),dimension(2),intent(in) :: e  !! the solution error
        real(wp),intent(in) :: ec  !! the constraint error
        character(len=128) :: buf  !! a line
        write(buf,'(es24.16,2(1x,es24.16))') t, y(1), y(2)
        write(yfid,'(a)') lower_e(trim(buf))
        write(buf,'(es24.16,3(1x,es24.16))') t, e(1), e(2), ec
        write(efid,'(a)') lower_e(trim(buf))
        end subroutine write_files

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

        function cg(x) result(s)
        !! a real formatted like C's `%g`: 6 significant digits, fixed notation
        !! if the decimal exponent is in [-4, 6), otherwise exponential
        !! notation, trailing zeros removed
        real(wp),intent(in) :: x  !! value to format
        character(len=:),allocatable :: s  !! formatted value
        integer,parameter :: p = 6  !! number of significant digits
        character(len=64) :: buf  !! write buffer
        character(len=16) :: fmt  !! run-time format
        integer :: ie  !! position of the exponent letter in buf
        integer :: iexp  !! decimal exponent (after rounding to p digits)
        if (x == 0.0_wp) then
            s = '0'
            return
        end if
        write(buf, '(es16.5e3)') x
        buf = adjustl(buf)
        ie = index(buf, 'E')
        read(buf(ie+1:), *) iexp
        if (iexp < -4 .or. iexp >= p) then
            s = strip_zeros(buf(1:ie-1))
            write(buf, '(a,"e",sp,i0.2)') s, iexp
            s = trim(buf)
        else
            write(fmt, '(a,i0,a)') '(f0.', p - 1 - iexp, ')'
            write(buf, fmt) x
            buf = adjustl(buf)
            if (buf(1:1) == '.') buf = '0'//buf
            if (buf(1:2) == '-.') buf = '-0'//buf(2:)
            s = strip_zeros(trim(buf))
        end if
        end function cg

        function strip_zeros(str) result(s)
        !! remove trailing zeros after a decimal point (and the point itself)
        character(len=*),intent(in) :: str  !! the number
        character(len=:),allocatable :: s  !! the number without trailing zeros
        integer :: n  !! length of the result
        s = trim(str)
        if (index(s, '.') == 0) return
        n = len(s)
        do while (s(n:n) == '0')
            n = n - 1
        end do
        if (s(n:n) == '.') n = n - 1
        s = s(1:n)
        end function strip_zeros

        function lj(i) result(s)
        !! integer left-justified in 6 characters (C: `%-6ld`)
        integer(li),intent(in) :: i  !! the integer
        character(len=6) :: s  !! the formatted integer
        write(s,'(i0)') i
        s = adjustl(s)
        end function lj

        logical function check_retval(retval, funcname)
        !! report a negative return flag (C: `check_retval(&retval, funcname, 1)`)
        integer,intent(in) :: retval  !! the return flag
        character(len=*),intent(in) :: funcname  !! name of the called function
        check_retval = retval < 0
        if (check_retval) write(*,'(/3a,i0/)') 'SUNDIALS_ERROR: ', funcname, '() failed with retval = ', retval
        end function check_retval

    end program cv_particle_dns
!*****************************************************************************************

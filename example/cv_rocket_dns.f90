!*****************************************************************************************
!>
!  Accelerating rocket: root finding with a change of the root functions,
!  BDF + Newton + dense linear solver with a user-supplied Jacobian.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvRocket_dns.c`
!  (same problem, settings and output; compare with `cvRocket_dns.out`).
!
!  A rocket of mass `massr` with fuel of mass `massf0` is accelerated
!  upwards by a constant force `Force`, burning fuel at rate `brate`, with
!  drag and gravity:
!
!      dH/dt = v
!      dv/dt = a(t) - Drag*v - grav,   a = Force/(massr + massf0 - brate*t)
!
!  While the engine is on, the root functions are `massf0 - brate*t` (fuel
!  exhausted) and `H - Hcut` (cutoff height). At the first root the engine
!  is turned off, and the only root function becomes `v` (maximum height).
!  The integration stops when the rocket falls below H = 0.

    module cv_rocket_dns_module

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: sun_matrix, dense_matrix

    implicit none

    private

    real(wp),parameter,public :: Force  = 2200.0_wp  !! engine force
    real(wp),parameter,public :: massr  = 10.0_wp    !! mass of the rocket
    real(wp),parameter,public :: massf0 = 1.0_wp     !! initial mass of the fuel
    real(wp),parameter,public :: brate  = 0.1_wp     !! fuel burn rate
    real(wp),parameter,public :: Drag   = 0.3_wp     !! drag coefficient
    real(wp),parameter,public :: grav   = 32.0_wp    !! acceleration due to gravity
    real(wp),parameter,public :: Hcut   = 4000.0_wp  !! height of the engine cutoff

    type,extends(cvode_solver),public :: rocket_solver
        !! the CVODE solver plus the problem data (C: user_data)
        logical :: engine_on = .true.  !! is the engine on?
    end type rocket_solver

    public :: f, g, jac

    contains

!*****************************************************************************************
!>
!  Right-hand side: dH/dt = v, dv/dt = a - Drag*v - grav.

    function f(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (a rocket_solver)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! (H, v)
    real(wp),dimension(:),intent(out) :: ydot  !! (dH/dt, dv/dt)
    integer :: retval  !! 0 (success)

    real(wp) :: acc  !! engine acceleration

    select type (me)
    class is (rocket_solver)
        ydot(1) = y(2)
        if (me%engine_on) then
            acc = Force / (massr + massf0 - brate * t)
        else
            acc = 0.0_wp
        end if
        ydot(2) = acc - Drag * y(2) - grav
    end select
    retval = 0

    end function f
!*****************************************************************************************

!*****************************************************************************************
!>
!  Jacobian: J = [0 1; 0 -Drag].

    function jac(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution (not used)
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y) (not used)
    class(sun_matrix),intent(inout) :: J  !! the Jacobian (a dense matrix)
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors (not used)
    integer :: retval  !! 0 (success)

    select type (J)
    type is (dense_matrix)
        J%data(1,1) = 0.0_wp
        J%data(2,1) = 0.0_wp
        J%data(1,2) = 1.0_wp
        J%data(2,2) = -Drag
    end select
    retval = 0

    end function jac
!*****************************************************************************************

!*****************************************************************************************
!>
!  Root functions: engine on: (fuel left, H - Hcut); engine off: (v).

    function g(me, t, y, gout) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (a rocket_solver)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! (H, v)
    real(wp),dimension(:),intent(out) :: gout  !! root function values
    integer :: retval  !! 0 (success)

    select type (me)
    class is (rocket_solver)
        if (me%engine_on) then
            gout(1) = massf0 - brate * t
            gout(2) = y(1) - Hcut
        else
            gout(1) = y(2)
        end if
    end select
    retval = 0

    end function g
!*****************************************************************************************

    end module cv_rocket_dns_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_rocket_dns

    use fvode
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use cv_rocket_dns_module

    implicit none

    integer,parameter  :: NEQ   = 2          !! number of equations
    real(wp),parameter :: Y1    = 0.0_wp     !! initial height
    real(wp),parameter :: Y2    = 0.0_wp     !! initial velocity
    real(wp),parameter :: RTOL  = 1.0e-5_wp  !! relative tolerance
    real(wp),parameter :: ATOL1 = 1.0e-2_wp  !! absolute tolerance of H
    real(wp),parameter :: ATOL2 = 1.0e-1_wp  !! absolute tolerance of v
    real(wp),parameter :: T0    = 0.0_wp     !! initial time
    real(wp),parameter :: T1    = 1.0_wp     !! first output time
    real(wp),parameter :: TINC  = 1.0_wp     !! output time increment
    integer,parameter  :: NOUT  = 70         !! number of output times

    type(rocket_solver) :: cvode  !! the solver
    real(wp),dimension(NEQ) :: y  !! (H, v)
    real(wp) :: t, tout  !! time reached, next output time
    integer :: retval, iout  !! return flag, output counter
    integer :: numroot  !! number of root functions
    integer,dimension(2) :: rootsfound  !! root information

    y = [Y1, Y2]

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, T0, y, retval)
    call cvode%sv_tolerances(RTOL, [ATOL1, ATOL2], retval)
    call cvode%root_init(2, g, retval)
    call cvode%set_linear_solver(sun_linsol_dense(int(NEQ,ip)), &
                                 dense_matrix(int(NEQ,ip), int(NEQ,ip)), retval)
    call cvode%set_jac_fn(jac, retval)

    write(*,'(a)') ' '
    write(*,'(a/)') 'Accelerating rocket problem'

    iout = 0
    tout = T1
    cvode%engine_on = .true.
    numroot = 2
    do
        call cvode%solve(tout, y, t, CV_NORMAL, retval)
        write(*,'(a,es10.4,a,es14.6,2x,es14.6)') 'At t = ', t, '      y =', y(1), y(2)

        if (cvode%engine_on .and. (retval == CV_ROOT_RETURN)) then
            ! engine cutoff: one root function (v), restart the integration
            call cvode%get_root_info(rootsfound)
            call print_root_info()
            cvode%engine_on = .false.
            numroot = 1
            ! (as in the C example, these calls overwrite retval, so tout is
            ! advanced below even though this output was a root)
            call cvode%root_init(1, g, retval)
            call cvode%reinit(t, y, retval)
        else if ((.not. cvode%engine_on) .and. (retval == CV_ROOT_RETURN)) then
            call cvode%get_root_info(rootsfound)
            call print_root_info()
        end if

        if (retval == CV_SUCCESS) then
            iout = iout + 1
            tout = tout + TINC
        end if

        if (iout == NOUT) exit
        if (y(1) < 0.0_wp) exit
    end do

    write(*,'(/a)') 'Final Statistics:'
    write(*,'(a)') 'nst = '//lj(cvode%get_num_steps())//' nfe  = '//lj(cvode%get_num_rhs_evals())// &
        ' nsetups = '//lj(cvode%get_num_lin_solv_setups())//' nfeLS = '// &
        lj(cvode%get_num_lin_rhs_evals())//' nje =  '//itoa(cvode%get_num_jac_evals())
    write(*,'(a/a)') 'nni = '//lj(cvode%get_num_nonlin_solv_iters())//' ncfn = '// &
        lj(cvode%get_num_nonlin_solv_conv_fails())//' netf = '// &
        lj(cvode%get_num_err_test_fails())//' nge = '//itoa(cvode%get_num_g_evals()), ' '

    call cvode%destroy()

    contains

!*****************************************************************************************
!>
!  Print the root information (C: PrintRootInfo).

        subroutine print_root_info()
        if (numroot == 2) write(*,'(a,i3,1x,i3)') '    rootsfound[] = ', rootsfound(1), rootsfound(2)
        if (numroot == 1) write(*,'(a,i3)') '    rootsfound[] = ', rootsfound(1)
        end subroutine print_root_info

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

    end program cv_rocket_dns
!*****************************************************************************************

!*****************************************************************************************
!>
!  Integration across discontinuities with stop times and re-initialization:
!  BDF + Newton + dense linear solver with the difference-quotient Jacobian.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvDisc_dns.c`
!  (same problem, settings and output; compare with `cvDisc_dns.out`).
!
!  The problem `y' = -y` (RHS1) or `y' = -5y` (RHS2), `y(0) = 1`, is
!  integrated over [0, 2] in one-step mode, stopping exactly at t = 1 with
!  a stop time:
!
!  1. a discontinuity in the solution: y is reset to 1 at t = 1, and the
!     solver is re-initialized;
!  2. a discontinuity in the right-hand side (RHS1 -> RHS2 at t = 1),
!     treated explicitly by re-initializing at t = 1;
!  3. the same discontinuity in the right-hand side, without
!     re-initialization (CVODE deals with it by reducing the step size).

    module cv_disc_dns_module

    use fvode
    use fvode_core, only: cvode_solver

    implicit none

    private

    integer,parameter,public :: RHS1 = 1  !! right-hand side y' = -y
    integer,parameter,public :: RHS2 = 2  !! right-hand side y' = -5y

    type,extends(cvode_solver),public :: disc_solver
        !! the CVODE solver plus the problem data (C: user_data)
        integer :: flag = RHS1  !! which right-hand side to use: RHS1 or RHS2
    end type disc_solver

    public :: f

    contains

!*****************************************************************************************
!>
!  Right-hand side: `-y` (RHS1) or `-5y` (RHS2).

    function f(me, t, y, ydot) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (a disc_solver)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: ydot  !! f(t, y)
    integer :: retval  !! 0 (success)

    select type (me)
    class is (disc_solver)
        select case (me%flag)
        case (RHS1); ydot(1) = -y(1)
        case (RHS2); ydot(1) = -5.0_wp * y(1)
        end select
    end select
    retval = 0

    end function f
!*****************************************************************************************

    end module cv_disc_dns_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_disc_dns

    use fvode
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use cv_disc_dns_module

    implicit none

    integer,parameter  :: NEQ = 1             !! number of equations
    real(wp),parameter :: reltol = 1.0e-3_wp  !! relative tolerance
    real(wp),parameter :: abstol = 1.0e-4_wp  !! absolute tolerance
    real(wp),parameter :: t0 = 0.0_wp         !! initial time
    real(wp),parameter :: t1 = 1.0_wp         !! time of the discontinuity
    real(wp),parameter :: t2 = 2.0_wp         !! final time

    type(disc_solver) :: cvode  !! the solver
    real(wp),dimension(NEQ) :: y  !! solution
    real(wp) :: t  !! time reached
    integer :: retval  !! return flag
    integer(li) :: nst1, nst2, nst  !! steps before / after t1, total

    y(1) = 1.0_wp

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, t0, y, retval)
    call cvode%ss_tolerances(reltol, abstol, retval)
    call cvode%set_linear_solver(sun_linsol_dense(int(NEQ,ip)), &
                                 dense_matrix(int(NEQ,ip), int(NEQ,ip)), retval)

    !-------------------------------------------------------------------------
    ! discontinuity in the solution
    !-------------------------------------------------------------------------
    write(*,'(/a/)') 'Discontinuity in solution'

    ! first leg: integrate to t1 (stop exactly there)
    call cvode%set_stop_time(t1, retval)
    cvode%flag = RHS1
    t = t0
    call print_output()
    call integrate_to(t1)
    nst1 = cvode%get_num_steps()

    ! reset y, re-initialize, and integrate to t2
    y(1) = 1.0_wp
    call cvode%reinit(t1, y, retval)
    call cvode%set_stop_time(t2, retval)
    cvode%flag = RHS1
    t = t1
    call print_output()
    call integrate_to(t2)
    nst2 = cvode%get_num_steps()
    nst = nst1 + nst2
    call print_steps()

    !-------------------------------------------------------------------------
    ! discontinuity in the right-hand side: explicit treatment
    !-------------------------------------------------------------------------
    write(*,'(/a/)') 'Discontinuity in RHS: Case 1 - explicit treatment'

    y(1) = 1.0_wp
    call cvode%reinit(t0, y, retval)
    call cvode%set_stop_time(t1, retval)
    cvode%flag = RHS1
    t = t0
    call print_output()
    call integrate_to(t1)
    nst1 = cvode%get_num_steps()

    ! re-initialize at t1 and switch the right-hand side
    call cvode%reinit(t1, y, retval)
    call cvode%set_stop_time(t2, retval)
    cvode%flag = RHS2
    t = t1
    call print_output()
    call integrate_to(t2)
    nst2 = cvode%get_num_steps()
    nst = nst1 + nst2
    call print_steps()

    !-------------------------------------------------------------------------
    ! discontinuity in the right-hand side: let CVODE deal with it
    !-------------------------------------------------------------------------
    write(*,'(/a/)') 'Discontinuity in RHS: Case 2 - let CVODE deal with it'

    y(1) = 1.0_wp
    call cvode%reinit(t0, y, retval)
    call cvode%set_stop_time(t1, retval)
    cvode%flag = RHS1
    t = t0
    call print_output()
    call integrate_to(t1)
    nst1 = cvode%get_num_steps()

    ! switch the right-hand side without re-initializing
    call cvode%set_stop_time(t2, retval)
    cvode%flag = RHS2
    t = t1
    call print_output()
    call integrate_to(t2)
    nst = cvode%get_num_steps()
    nst2 = nst - nst1
    call print_steps()

    call cvode%destroy()

    contains

!*****************************************************************************************
!>
!  Integrate in one-step mode until tend, printing each step.

        subroutine integrate_to(tend)
        real(wp),intent(in) :: tend  !! end of the interval (also the stop time)
        do while (t < tend)
            call cvode%solve(tend, y, t, CV_ONE_STEP, retval)
            call print_output()
        end do
        end subroutine integrate_to

!*****************************************************************************************
!>
!  Print t and y (C: "%12.8e  %12.8e").

        subroutine print_output()
        write(*,'(es14.8,2x,es14.8)') t, y(1)
        end subroutine print_output

!*****************************************************************************************
!>
!  Print the number of steps on both sides of the discontinuity.

        subroutine print_steps()
        write(*,'(/a,i0,a,i0,a,i0)') 'Number of steps: ', nst1, ' + ', nst2, ' = ', nst
        end subroutine print_steps

    end program cv_disc_dns
!*****************************************************************************************

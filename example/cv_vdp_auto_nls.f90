!*****************************************************************************************
!>
!  Van der Pol oscillator (mu = 100), BDF with the nonlinear solver that
!  switches automatically between fixed-point and Newton iteration
!  ([[sun_nonlinsol_auto]]), and a dense linear solver with a user
!  Jacobian for the Newton phases.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvVdp_auto_nls.c`
!  (same problem, settings and output; compare with `cvVdp_auto_nls.out`).
!  The C example also passes its command-line arguments to
!  `CVodeSetOptions`; that is not ported.
!
!  The system is
!
!      y1' = y2
!      y2' = mu (1 - y1^2) y2 - y1
!
!  with y1(0) = 2, y2(0) = 0, for 0 <= t <= 250.

    module cv_vdp_auto_nls_module

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: sun_matrix, dense_matrix

    implicit none

    private

    type,extends(cvode_solver),public :: vdp_solver
        !! the CVODE solver plus the problem data (C: UserData)
        real(wp) :: mu = 100.0_wp  !! the Van der Pol parameter
    end type vdp_solver

    public :: f, jac

    contains

    function f(me, t, y, ydot) result(retval)
    !! right-hand side f(t,y)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: ydot  !! dy/dt
    integer :: retval  !! 0 on success
    select type (me)
    class is (vdp_solver)
        ydot(1) = y(2)
        ydot(2) = me%mu * (1.0_wp - y(1) * y(1)) * y(2) - y(1)
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
    class is (vdp_solver)
        select type (J)
        type is (dense_matrix)
            J%data(1,1) = 0.0_wp
            J%data(1,2) = 1.0_wp
            J%data(2,1) = -2.0_wp * me%mu * y(1) * y(2) - 1.0_wp
            J%data(2,2) = me%mu * (1.0_wp - y(1) * y(1))
        end select
    end select
    retval = 0
    end function jac

    end module cv_vdp_auto_nls_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_vdp_auto_nls

    use iso_fortran_env, only: output_unit
    use fvode
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use fvode_nonlinsol_auto, only: sun_nonlinsol_auto, SUNNONLINSOL_AUTO_NEWTON
    use cv_vdp_auto_nls_module

    implicit none

    integer,parameter  :: NEQ   = 2          !! number of equations
    real(wp),parameter :: T0    = 0.0_wp     !! initial time
    real(wp),parameter :: TF    = 250.0_wp   !! final time
    real(wp),parameter :: DTOUT = 10.0_wp    !! output interval
    real(wp),parameter :: y10   = 2.0_wp     !! initial y1
    real(wp),parameter :: y20   = 0.0_wp     !! initial y2
    real(wp),parameter :: reltol = 1.0e-4_wp !! relative tolerance
    real(wp),parameter :: abstol = 1.0e-4_wp !! absolute tolerance

    type(vdp_solver) :: cvode  !! the CVODE solver with the problem data
    real(wp),dimension(NEQ) :: y  !! solution vector
    real(wp) :: tout, tret  !! output and return times
    integer :: retval, iout, Nt  !! return flag, output counter, number of outputs
    integer(li) :: nfp, nnewt  !! total fixed-point and Newton iterations

    y = [y10, y20]

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, T0, y, retval)
    if (check_retval(retval, 'CVodeInit')) stop 1
    call cvode%ss_tolerances(reltol, abstol, retval)
    if (check_retval(retval, 'CVodeSStolerances')) stop 1
    call cvode%set_max_num_steps(10000_li, retval)
    if (check_retval(retval, 'CVodeSetMaxNumSteps')) stop 1

    ! create nonlinear solver (auto)
    call cvode%set_nonlinear_solver(sun_nonlinsol_auto(NEQ, 0, SUNNONLINSOL_AUTO_NEWTON), retval)
    if (check_retval(retval, 'CVodeSetNonlinearSolver')) stop 1

    ! provide dense linear solver and Jacobian for when Newton is active
    call cvode%set_linear_solver(sun_linsol_dense(int(NEQ, ip)), dense_matrix(int(NEQ, ip), int(NEQ, ip)), retval)
    if (check_retval(retval, 'CVodeSetLinearSolver')) stop 1
    call cvode%set_jac_fn(jac, retval)
    if (check_retval(retval, 'CVodeSetJacFn')) stop 1

    write(*,'(/a)') 'Van der Pol oscillator (CVODE):'
    write(*,'(a,f0.6,a,f0.6)') '    initial conditions: y1 = ', y10, ', y2 = ', y20
    write(*,'(a,f0.6)') '    mu = ', cvode%mu
    write(*,'(a,es8.2,a,es8.2/)') '    reltol = ', reltol, ', abstol = ', abstol
    write(*,'(a)') '        t           y1           y2'
    write(*,'(a)') '   -----------------------------------'
    write(*,'(2x,f10.6,2x,f10.6,2x,f10.6)') T0, y(1), y(2)

    Nt = ceiling(TF / DTOUT)
    tout = T0 + DTOUT
    do iout = 0, Nt - 1
        call cvode%solve(tout, y, tret, CV_NORMAL, retval)
        write(*,'(2x,f10.6,2x,f10.6,2x,f10.6)') tret, y(1), y(2)
        if (retval == CV_SUCCESS) then
            tout = tout + DTOUT
            tout = min(tout, TF)
        else
            write(*,'(a)') 'Solver failure, stopping integration'
            exit
        end if
    end do
    write(*,'(a)') '   -----------------------------------'

    call cvode%print_all_stats(output_unit, SUN_OUTPUTFORMAT_TABLE, retval)
    if (check_retval(retval, 'CVodePrintAllStats')) stop 1

    ! the totals of the attached (copied) auto solver
    select type (NLS => cvode%NLS)
    type is (sun_nonlinsol_auto)
        call NLS%get_total_num_iters_by_type(nfp, nnewt)
    end select
    write(*,'(a,i0,a,i0)') '   Auto nonlinear solver iteration totals: newton = ', nnewt, &
        ', fixed-point = ', nfp

    call cvode%destroy()

    contains

        logical function check_retval(retval, funcname)
        !! report a negative return flag
        integer,intent(in) :: retval  !! the return flag
        character(len=*),intent(in) :: funcname  !! name of the called function
        check_retval = retval < 0
        if (check_retval) write(*,'(3a,i0)') 'ERROR: ', funcname, '() returned ', retval
        end function check_retval

    end program cv_vdp_auto_nls
!*****************************************************************************************

!*****************************************************************************************
!>
!  Step trace: van der Pol (cvDirectDemo_ls problem 1), Adams + fixed-point
!  iteration, scalar tolerances (reltol = 0, abstol = 1e-6).
!  Fortran counterpart of `reference/trace_vdp_adams_fp.c`; writes
!  `reference/out/fvode_vdp_adams_fp.csv`.

    program trace_vdp_adams_fp

    use fvode
    use fvode_core
    use fvode_nonlinsol_fixedpoint
    use trace_common

    implicit none

    real(wp),parameter :: eta = 3.0_wp
    real(wp),parameter :: t0 = 0.0_wp
    real(wp),parameter :: tfinal = 1.39283880203_wp + 3.0_wp * 2.214773875_wp

    type(cvode_solver) :: cv
    real(wp),dimension(2) :: y
    real(wp) :: t
    integer :: flag, iunit

    y = [2.0_wp, 0.0_wp]

    call cv%create(CV_ADAMS, flag);                                   call check('create')
    call cv%init(f, t0, y, flag);                                     call check('init')
    call cv%ss_tolerances(0.0_wp, 1.0e-6_wp, flag);                   call check('ss_tolerances')
    call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(2, 0), flag); call check('set_nonlinear_solver')
    call cv%set_max_num_steps(100000_li, flag);                       call check('set_max_num_steps')

    open(newunit=iunit, file='reference/out/fvode_vdp_adams_fp.csv', status='replace')
    call trace_header(iunit, 2)
    t = t0
    do while (t < tfinal)
        call cv%solve(tfinal, y, t, CV_ONE_STEP, flag)
        call trace_row(iunit, cv, flag, t, y)
        if (flag < 0) exit
    end do
    close(iunit)

    contains

        function f(me, t, y, ydot) result(retval)
        class(cvode_solver),intent(inout) :: me
        real(wp),intent(in) :: t
        real(wp),dimension(:),intent(in) :: y
        real(wp),dimension(:),intent(out) :: ydot
        integer :: retval
        ydot(1) = y(2)
        ydot(2) = (1.0_wp - y(1)*y(1)) * eta * y(2) - y(1)
        retval = 0
        end function f

        subroutine check(name)
        character(len=*),intent(in) :: name
        if (flag < 0) then
            write(*,'(a,i0)') name//' failed with flag ', flag
            error stop
        end if
        end subroutine check

    end program trace_vdp_adams_fp
!*****************************************************************************************

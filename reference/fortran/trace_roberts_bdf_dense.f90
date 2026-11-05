!*****************************************************************************************
!>
!  Step trace: Robertson kinetics (cvRoberts_dns without root finding), BDF +
!  Newton + dense linear solver with the internal difference-quotient
!  Jacobian, vector absolute tolerances. Fortran counterpart of
!  `reference/trace_roberts_bdf_dense.c`.

    program trace_roberts_bdf_dense

    use fvode
    use fvode_core
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use trace_common

    implicit none

    real(wp),parameter :: t0 = 0.0_wp
    real(wp),parameter :: tfinal = 4.0e10_wp

    type(cvode_solver) :: cv
    real(wp),dimension(3) :: y
    real(wp) :: t
    integer :: flag, iunit

    y = [1.0_wp, 0.0_wp, 0.0_wp]

    call cv%create(CV_BDF, flag);                         call check('create')
    call cv%init(f, t0, y, flag);                         call check('init')
    call cv%sv_tolerances(1.0e-4_wp, [1.0e-8_wp, 1.0e-14_wp, 1.0e-6_wp], flag)
    call check('sv_tolerances')
    call cv%set_linear_solver(sun_linsol_dense(3_ip), dense_matrix(3_ip, 3_ip), flag)
    call check('set_linear_solver')
    call cv%set_max_num_steps(100000_li, flag);           call check('set_max_num_steps')

    open(newunit=iunit, file='reference/out/fvode_roberts_bdf_dense.csv', status='replace')
    call trace_header(iunit, 3)
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
        ydot(1) = -0.04_wp * y(1) + 1.0e4_wp * y(2) * y(3)
        ydot(3) = 3.0e7_wp * y(2) * y(2)
        ydot(2) = -ydot(1) - ydot(3)
        retval = 0
        end function f

        subroutine check(name)
        character(len=*),intent(in) :: name
        if (flag < 0) then
            write(*,'(a,i0)') name//' failed with flag ', flag
            error stop
        end if
        end subroutine check

    end program trace_roberts_bdf_dense
!*****************************************************************************************

!*****************************************************************************************
!>
!  Step trace: stability limit detection (STALD) on a damped fast oscillation
!  (eigenvalues -1 +- 100i) plus a slow decay (-0.01); BDF + Newton + dense
!  DQ Jacobian. Fortran counterpart of `reference/trace_osc_bdf_stald.c`.

    program trace_osc_bdf_stald

    use fvode
    use fvode_core
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use trace_common

    implicit none

    real(wp),parameter :: t0 = 0.0_wp      !! initial time
    real(wp),parameter :: tfinal = 500.0_wp  !! final time

    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(3) :: y  !! solution vector
    real(wp) :: t  !! time reached
    integer :: flag  !! return flag
    integer :: iunit  !! unit of the trace file

    y = [1.0_wp, 0.0_wp, 1.0_wp]

    call cv%create(CV_BDF, flag)
    call cv%init(f, t0, y, flag)
    call cv%ss_tolerances(1.0e-6_wp, 1.0e-10_wp, flag)
    call cv%set_linear_solver(sun_linsol_dense(3_ip), dense_matrix(3_ip, 3_ip), flag)
    call cv%set_stab_lim_det(.true., flag)
    call cv%set_max_num_steps(100000_li, flag)

    open(newunit=iunit, file='reference/out/fvode_osc_bdf_stald.csv', status='replace')
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
        !! right-hand side
        class(cvode_solver),intent(inout) :: me  !! the CVODE solver (not used)
        real(wp),intent(in) :: t  !! current time (not used)
        real(wp),dimension(:),intent(in) :: y  !! current solution
        real(wp),dimension(:),intent(out) :: ydot  !! f(t, y)
        integer :: retval  !! 0 (success)
        ydot(1) = -y(1) + 100.0_wp * y(2)
        ydot(2) = -100.0_wp * y(1) - y(2)
        ydot(3) = -0.01_wp * y(3)
        retval = 0
        end function f

    end program trace_osc_bdf_stald
!*****************************************************************************************

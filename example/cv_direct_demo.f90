!*****************************************************************************************
!>
!  Demonstration of the direct linear solvers: dense, band and diagonal,
!  with user-supplied and difference-quotient Jacobians.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvDirectDemo_ls.c`
!  (same problems, settings and output, except the workspace lengths, which
!  FVODE does not track; compare with `cvDirectDemo_ls.out`).
!
!  * Problem 1: van der Pol oscillator `x'' - 3(1 - x^2)x' + x = 0`,
!    `x(0) = 2`, `x'(0) = 0`, solved with fixed-point iteration, then Newton
!    with a dense matrix (user and difference-quotient Jacobian) and with the
!    diagonal approximation; each with Adams and BDF.
!  * Problem 2: `y' = A y`, A a banded lower triangular matrix from a 2-D
!    advection PDE on a 5 x 5 mesh (exact solution known); solved with
!    fixed-point iteration, then Newton with the diagonal approximation and a
!    band matrix (user and difference-quotient Jacobian); each with Adams
!    and BDF.

    program cv_direct_demo

    use fvode
    use fvode_core
    use fvode_matrix, only: sun_matrix, dense_matrix, band_matrix
    use fvode_linsol_direct, only: sun_linsol_dense, sun_linsol_band
    use fvode_nonlinsol_fixedpoint, only: sun_nonlinsol_fixedpoint
    use fvode_nonlinsol_newton, only: sun_nonlinsol_newton

    implicit none

    ! shared problem constants
    real(wp),parameter :: ATOL = 1.0e-6_wp  !! absolute tolerance
    real(wp),parameter :: RTOL = 0.0_wp  !! relative tolerance

    ! problem 1 constants
    integer,parameter  :: P1_NEQ = 2  !! problem 1: number of equations
    real(wp),parameter :: P1_ETA = 3.0_wp  !! problem 1: van der Pol parameter
    integer,parameter  :: P1_NOUT = 4  !! problem 1: number of output times
    real(wp),parameter :: P1_T0 = 0.0_wp  !! problem 1: initial time
    real(wp),parameter :: P1_T1 = 1.39283880203_wp  !! problem 1: first output time
    real(wp),parameter :: P1_DTOUT = 2.214773875_wp  !! problem 1: output time interval
    real(wp),parameter :: P1_TOL_FACTOR = 1.0e4_wp  !! problem 1: error check factor

    ! problem 2 constants
    integer,parameter  :: P2_MESHX = 5  !! problem 2: mesh points in x
    integer,parameter  :: P2_MESHY = 5  !! problem 2: mesh points in y
    integer,parameter  :: P2_NEQ = P2_MESHX * P2_MESHY  !! problem 2: number of equations
    real(wp),parameter :: P2_ALPH1 = 1.0_wp  !! problem 2: advection coefficient in x
    real(wp),parameter :: P2_ALPH2 = 1.0_wp  !! problem 2: advection coefficient in y
    integer,parameter  :: P2_NOUT = 5  !! problem 2: number of output times
    integer,parameter  :: P2_ML = 5  !! problem 2: lower bandwidth
    integer,parameter  :: P2_MU = 0  !! problem 2: upper bandwidth
    real(wp),parameter :: P2_T0 = 0.0_wp  !! problem 2: initial time
    real(wp),parameter :: P2_T1 = 0.01_wp  !! problem 2: first output time
    real(wp),parameter :: P2_TOUT_MULT = 10.0_wp  !! problem 2: output time multiplier
    real(wp),parameter :: P2_TOL_FACTOR = 1.0e3_wp  !! problem 2: error check factor

    ! linear solver options
    integer,parameter :: FUNC = 0, DENSE_USER = 1, DENSE_DQ = 2, DIAG = 3, &
                         BAND_USER = 4, BAND_DQ = 5  !! linear solver options: fixed point, dense (user/DQ Jacobian), diagonal, band (user/DQ Jacobian)

    integer :: nerr  !! total number of errors

    nerr = problem1()
    nerr = nerr + problem2()
    call print_err_info(nerr)

    contains

!*****************************************************************************************
    function problem1() result(nerr)
    !! run problem 1 with Adams and BDF and each solver option; returns the number of errors

    integer :: nerr  !! number of errors

    type(cvode_solver) :: cvode  !! the CVODE solver
    real(wp),dimension(P1_NEQ) :: y  !! solution vector
    real(wp) :: t, tout, ero, er  !! time reached, output time, maximum error ratio, and error ratio at an output
    integer :: lmm, miter, retval, iout  !! method, solver option, return flag, and output counter

    nerr = 0

    write(*,'(a)') 'Demonstration program for CVODE package - direct linear solvers'
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') 'Problem 1: Van der Pol oscillator'
    write(*,'(a)') ' xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0'
    write(*,'(a)') ' neq = 2,  reltol = 0,  abstol = 1e-06'

    do lmm = CV_ADAMS, CV_BDF
        call cvode%create(lmm, retval)
        do miter = FUNC, DIAG
            ero = 0.0_wp
            y = [2.0_wp, 0.0_wp]
            if (miter == FUNC) then
                call cvode%init(f1, P1_T0, y, retval)
                call cvode%ss_tolerances(RTOL, ATOL, retval)
            else
                call cvode%reinit(P1_T0, y, retval)
            end if
            call prepare_next_run(cvode, lmm, miter, P1_NEQ, 0, 0)
            write(*,'(/a)') '     t           x              xdot         qu     hu '
            tout = P1_T1
            do iout = 1, P1_NOUT
                call cvode%solve(tout, y, t, CV_NORMAL, retval)
                write(*,'(f10.5,4x,es12.5,3x,es12.5,3x,i2,4x,es10.4)') &
                    t, y(1), y(2), cvode%get_last_order(), cvode%get_last_step()
                if (retval /= CV_SUCCESS) then
                    nerr = nerr + 1
                    exit
                end if
                if (mod(iout, 2) == 0) then
                    er = abs(y(1)) / ATOL
                    if (er > ero) ero = er
                    if (er > P1_TOL_FACTOR) then
                        nerr = nerr + 1
                        call print_err_output(P1_TOL_FACTOR)
                    end if
                end if
                tout = tout + P1_DTOUT
            end do
            call print_final_stats(cvode, miter, ero)
        end do
        call cvode%destroy()
    end do

    end function problem1
!*****************************************************************************************

!*****************************************************************************************
    function f1(me, t, y, ydot) result(retval)
    !! right-hand side of problem 1 (van der Pol)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success
    ydot(1) = y(2)
    ydot(2) = (1.0_wp - y(1)**2) * P1_ETA * y(2) - y(1)
    retval = 0
    end function f1
!*****************************************************************************************

!*****************************************************************************************
    function jac1(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
    !! dense Jacobian of problem 1
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y, fy  !! current state vector and f(t,y)
    class(sun_matrix),intent(inout) :: J  !! Jacobian matrix df/dy (filled in here)
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors (not used)
    integer :: retval  !! return value: 0 on success
    select type (J)
    type is (dense_matrix)
        J%data(1,2) = 1.0_wp
        J%data(2,1) = -2.0_wp * P1_ETA * y(1) * y(2) - 1.0_wp
        J%data(2,2) = P1_ETA * (1.0_wp - y(1)**2)
    end select
    retval = 0
    end function jac1
!*****************************************************************************************

!*****************************************************************************************
    function problem2() result(nerr)
    !! run problem 2 with Adams and BDF and each solver option; returns the number of errors

    integer :: nerr  !! number of errors

    type(cvode_solver) :: cvode  !! the CVODE solver
    real(wp),dimension(P2_NEQ) :: y  !! solution vector
    real(wp) :: t, tout, er, erm, ero  !! time reached, output time, error ratio, maximum error at an output, and maximum error ratio
    integer :: lmm, miter, retval, iout  !! method, solver option, return flag, and output counter

    nerr = 0

    write(*,'(/a)') ''
    write(*,'(a)') '-------------------------------------------------------------'
    write(*,'(a)') '-------------------------------------------------------------'
    write(*,'(/a)') 'Problem 2: ydot = A * y, where A is a banded lower'
    write(*,'(a/)') 'triangular matrix derived from 2-D advection PDE'
    write(*,'(a)') ' neq = 25, ml = 5, mu = 0'
    write(*,'(a)') ' itol = CV_SS, reltol = 0, abstol = 1e-06'

    do lmm = CV_ADAMS, CV_BDF
        call cvode%create(lmm, retval)
        do miter = FUNC, BAND_DQ
            if (miter == DENSE_USER .or. miter == DENSE_DQ) cycle
            ero = 0.0_wp
            y = 0.0_wp
            y(1) = 1.0_wp
            if (miter == FUNC) then
                call cvode%init(f2, P2_T0, y, retval)
                call cvode%ss_tolerances(RTOL, ATOL, retval)
            else
                call cvode%reinit(P2_T0, y, retval)
            end if
            call prepare_next_run(cvode, lmm, miter, P2_NEQ, P2_MU, P2_ML)
            write(*,'(/a)') '      t        max.err      qu     hu '
            tout = P2_T1
            do iout = 1, P2_NOUT
                call cvode%solve(tout, y, t, CV_NORMAL, retval)
                erm = max_error(y, t)
                write(*,'(f10.3,2x,es12.4,3x,i2,3x,es12.4)') &
                    t, erm, cvode%get_last_order(), cvode%get_last_step()
                if (retval /= CV_SUCCESS) then
                    nerr = nerr + 1
                    exit
                end if
                er = erm / ATOL
                if (er > ero) ero = er
                if (er > P2_TOL_FACTOR) then
                    nerr = nerr + 1
                    call print_err_output(P2_TOL_FACTOR)
                end if
                tout = tout * P2_TOUT_MULT
            end do
            call print_final_stats(cvode, miter, ero)
        end do
        call cvode%destroy()
    end do

    end function problem2
!*****************************************************************************************

!*****************************************************************************************
    function f2(me, t, y, ydot) result(retval)
    !! excluding boundaries,
    !! ydot(i,j) = -2 y(i,j) + alpha1 * y(i-1,j) + alpha2 * y(i,j-1)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success
    integer :: i, j, k  !! mesh indices in x and y, and index in the solution vector
    real(wp) :: d  !! accumulated derivative of y(k)
    do j = 0, P2_MESHY - 1
        do i = 0, P2_MESHX - 1
            k = 1 + i + j * P2_MESHX
            d = -2.0_wp * y(k)
            if (i /= 0) d = d + P2_ALPH1 * y(k - 1)
            if (j /= 0) d = d + P2_ALPH2 * y(k - P2_MESHX)
            ydot(k) = d
        end do
    end do
    retval = 0
    end function f2
!*****************************************************************************************

!*****************************************************************************************
    function jac2(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
    !! band Jacobian of problem 2
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y, fy  !! current state vector and f(t,y)
    class(sun_matrix),intent(inout) :: J  !! Jacobian matrix df/dy (filled in here)
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors (not used)
    integer :: retval  !! return value: 0 on success
    integer(ip) :: i, jj, k  !! mesh indices in x and y, and index in the solution vector
    select type (J)
    type is (band_matrix)
        do jj = 0, P2_MESHY - 1
            do i = 0, P2_MESHX - 1
                k = 1 + i + jj * P2_MESHX
                call J%set(k, k, -2.0_wp)
                if (i /= P2_MESHX - 1) call J%set(k + 1, k, P2_ALPH1)
                if (jj /= P2_MESHY - 1) call J%set(k + P2_MESHX, k, P2_ALPH2)
            end do
        end do
    end select
    retval = 0
    end function jac2
!*****************************************************************************************

!*****************************************************************************************
    function max_error(y, t) result(maxerr)
    !! maximum error against the exact solution
    !! y(i,j) = t^(i+j) exp(-2t) / (i! j!)
    real(wp),dimension(:),intent(in) :: y  !! computed solution
    real(wp),intent(in) :: t  !! time of the solution
    real(wp) :: maxerr  !! maximum absolute error
    integer :: i, j, k  !! mesh indices in x and y, and index in the solution vector
    real(wp) :: er, ex, yt, ifact_inv, jfact_inv  !! error at a point, exp(-2t), exact solution at a point, and 1/i!, 1/j!
    maxerr = 0.0_wp
    if (t == 0.0_wp) return
    ex = 0.0_wp
    if (t <= 30.0_wp) ex = exp(-2.0_wp * t)
    jfact_inv = 1.0_wp
    do j = 0, P2_MESHY - 1
        ifact_inv = 1.0_wp
        do i = 0, P2_MESHX - 1
            k = 1 + i + j * P2_MESHX
            yt = t**(i + j) * ex * ifact_inv * jfact_inv
            er = abs(y(k) - yt)
            if (er > maxerr) maxerr = er
            ifact_inv = ifact_inv / (i + 1)
        end do
        jfact_inv = jfact_inv / (j + 1)
    end do
    end function max_error
!*****************************************************************************************

!*****************************************************************************************
    subroutine prepare_next_run(cvode, lmm, miter, neq, mu, ml)
    !! attach the nonlinear and linear solvers for the next run, and print its header
    type(cvode_solver),intent(inout) :: cvode  !! the CVODE solver
    integer,intent(in) :: lmm, miter, neq, mu, ml  !! method, solver option, number of equations, and upper/lower bandwidths
    integer :: retval  !! return flag
    integer(ip) :: n  !! number of equations (index kind)
    n = neq
    write(*,'(/a)') ''
    write(*,'(a)') '-------------------------------------------------------------'
    write(*,'(/a)',advance='no') 'Linear Multistep Method : '
    if (lmm == CV_ADAMS) then
        write(*,'(a)') 'ADAMS'
    else
        write(*,'(a)') 'BDF'
    end if
    write(*,'(a)',advance='no') 'Iteration               : '
    if (miter == FUNC) then
        write(*,'(a)') 'FIXEDPOINT'
        call cvode%set_nonlinear_solver(sun_nonlinsol_fixedpoint(neq, 0), retval)
    else
        write(*,'(a)') 'NEWTON'
        call cvode%set_nonlinear_solver(sun_nonlinsol_newton(neq), retval)
        write(*,'(a)',advance='no') 'Linear Solver           : '
        select case (miter)
        case (DENSE_USER)
            write(*,'(a)') 'Dense, User-Supplied Jacobian'
            call cvode%set_linear_solver(sun_linsol_dense(n), dense_matrix(n, n), retval)
            call cvode%set_jac_fn(jac1, retval)
        case (DENSE_DQ)
            write(*,'(a)') 'Dense, Difference Quotient Jacobian'
            call cvode%set_linear_solver(sun_linsol_dense(n), dense_matrix(n, n), retval)
            call cvode%set_jac_fn(flag=retval)
        case (DIAG)
            write(*,'(a)') 'Diagonal Jacobian'
            call cvode%diag(retval)
        case (BAND_USER)
            write(*,'(a)') 'Band, User-Supplied Jacobian'
            call cvode%set_linear_solver(sun_linsol_band(n), &
                                         band_matrix(n, int(mu,ip), int(ml,ip)), retval)
            call cvode%set_jac_fn(jac2, retval)
        case (BAND_DQ)
            write(*,'(a)') 'Band, Difference Quotient Jacobian'
            call cvode%set_linear_solver(sun_linsol_band(n), &
                                         band_matrix(n, int(mu,ip), int(ml,ip)), retval)
            call cvode%set_jac_fn(flag=retval)
        end select
    end if
    end subroutine prepare_next_run
!*****************************************************************************************

!*****************************************************************************************
    subroutine print_err_output(tol_factor)
    !! print the error-exceeds-tolerance message
    real(wp),intent(in) :: tol_factor  !! tolerance factor that was exceeded
    write(*,'(//a,g0,a//)') ' Error exceeds ', tol_factor, ' * tolerance '
    end subroutine print_err_output
!*****************************************************************************************

!*****************************************************************************************
    subroutine print_final_stats(cvode, miter, ero)
    !! print the statistics of one run
    type(cvode_solver),intent(in) :: cvode  !! the CVODE solver
    integer,intent(in) :: miter  !! solver option of the run
    real(wp),intent(in) :: ero  !! maximum error ratio of the run
    integer(li) :: nje, nfeLS  !! number of Jacobian evaluations and of f evaluations in the linear solver
    write(*,'(/a/)') ' Final statistics for this run:'
    write(*,'(a,i4,a)') ' Number of steps                          = ', cvode%get_num_steps(), ' '
    write(*,'(a,i4,a)') ' Number of f-s                            = ', cvode%get_num_rhs_evals(), ' '
    write(*,'(a,i4,a)') ' Number of setups                         = ', cvode%get_num_lin_solv_setups(), ' '
    write(*,'(a,i4,a)') ' Number of nonlinear iterations           = ', cvode%get_num_nonlin_solv_iters(), ' '
    write(*,'(a,i4,a)') ' Number of nonlinear convergence failures = ', &
                        cvode%get_num_nonlin_solv_conv_fails(), ' '
    write(*,'(a,i4,a/)') ' Number of error test failures            = ', cvode%get_num_err_test_fails(), ' '
    if (miter /= FUNC) then
        if (miter /= DIAG) then
            nje   = cvode%get_num_jac_evals()
            nfeLS = cvode%get_num_lin_rhs_evals()
        else
            nje   = cvode%get_num_lin_solv_setups()
            nfeLS = cvode%get_diag_num_rhs_evals()
        end if
        write(*,'(a,i4,a)') ' Number of Jacobian evaluations           = ', nje, ' '
        write(*,'(a,i4,a/)') ' Number of f evals. in linear solver      = ', nfeLS, ' '
    end if
    write(*,'(a,f0.3,a)') ' Error overrun = ', ero, ' '
    end subroutine print_final_stats
!*****************************************************************************************

!*****************************************************************************************
    subroutine print_err_info(nerr)
    !! print the total number of errors
    integer,intent(in) :: nerr  !! total number of errors
    write(*,'(/a)') ''
    write(*,'(a)') '-------------------------------------------------------------'
    write(*,'(a)') '-------------------------------------------------------------'
    write(*,'(/a,i0,a)') ' Number of errors encountered = ', nerr, ' '
    end subroutine print_err_info
!*****************************************************************************************

    end program cv_direct_demo
!*****************************************************************************************

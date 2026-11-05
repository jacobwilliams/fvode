!*****************************************************************************************
!>
!  Shared driver for the `CVodeResizeHistory` step traces: the Fortran
!  counterpart of `reference/resize_driver.h` (after the SUNDIALS unit test
!  `cv_test_resize_history.cpp`).

    module resize_driver

    use fvode
    use fvode_core
    use fvode_matrix, only: sun_matrix, band_matrix
    use fvode_linsol_direct, only: sun_linsol_band
    use fvode_nonlinsol_fixedpoint, only: sun_nonlinsol_fixedpoint
    use trace_common

    implicit none

    private

    public :: run_resize

    real(wp),parameter :: lambda = -1.0_wp  !! the decay rate

    contains
!*****************************************************************************************

    function ode_rhs(me, t, y, f) result(retval)
    !! `y' = lambda*(y - atan(t)) + 1/(1+t^2)` for every component
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: f  !! dy/dt
    integer :: retval  !! 0 on success
    retval = rhs(t, y, f)
    end function ode_rhs

    function rhs(t, y, f) result(retval)
    !! the right-hand side (without the solver)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(out) :: f  !! dy/dt
    integer :: retval  !! 0 on success
    real(wp) :: phi_t, phi_prime_t  !! atan(t) and its derivative
    integer :: i  !! component index
    phi_t = atan(t)
    phi_prime_t = 1.0_wp / (1.0_wp + t * t)
    do i = 1, size(y)
        f(i) = lambda * (y(i) - phi_t) + phi_prime_t
    end do
    retval = 0
    end function rhs

    function ode_rhs_jac(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
    !! the (diagonal) Jacobian
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t,y)
    class(sun_matrix),intent(inout) :: J  !! the Jacobian
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors
    integer :: retval  !! 0 on success
    select type (J)
    type is (band_matrix)
        J%data(0,:) = lambda
    end select
    retval = 0
    end function ode_rhs_jac

    subroutine run_resize(mode, method, newton, fname)
    !! integrate 60 steps, resizing the history after every step
    integer,intent(in) :: mode  !! 1: same size; 2: grow by one element each step
    integer,intent(in) :: method  !! `CV_ADAMS` or `CV_BDF`
    logical,intent(in) :: newton  !! Newton + band (true) or fixed point (false)
    character(len=*),intent(in) :: fname  !! output file
    type(cvode_solver) :: cv  !! the solver
    real(wp),dimension(:),allocatable :: y  !! solution vector
    real(wp),dimension(:,:),allocatable :: y_hist, f_hist  !! solution and right-hand side history
    real(wp),dimension(12) :: t_hist  !! time history
    real(wp) :: t_ret  !! return time
    integer :: hist_size, i, n_hist, flag, iunit, n  !! history size, step, used history, flag, unit, size

    allocate(y(1), source=0.0_wp)  ! atan(0)
    call cv%create(method, flag)
    call cv%init(ode_rhs, 0.0_wp, y, flag)
    call cv%ss_tolerances(1.0e-5_wp, 1.0e-8_wp, flag)
    call attach_solvers(1)

    hist_size = merge(5, 12, method == CV_BDF)
    t_hist = 0.0_wp
    allocate(y_hist(1, hist_size), f_hist(1, hist_size), source=0.0_wp)
    y_hist(:,1) = y
    flag = rhs(t_hist(1), y_hist(:,1), f_hist(:,1))

    open(newunit=iunit, file=fname, status='replace')
    call trace_header(iunit, 1)
    do i = 1, 60
        call cv%solve(10.0_wp, y, t_ret, CV_ONE_STEP, flag)
        call trace_row(iunit, cv, flag, t_ret, y(1:1))
        if (flag < 0) exit

        call save_history(t_ret, y, i)
        if (mode == 2) call resize_history()
        n_hist = merge(i + 1, hist_size, i < hist_size)
        call cv%resize_history(t_hist(1:n_hist), y_hist(:,1:n_hist), f_hist(:,1:n_hist), flag)
        if (flag < 0) then
            write(*,'(a,i0)') 'CVodeResizeHistory failed with ', flag
            exit
        end if
        if (mode == 2) then
            n = size(y_hist, 1)
            deallocate(y)
            allocate(y(n), source=0.0_wp)
            call attach_solvers(n)
        end if
    end do
    close(iunit)
    call cv%destroy()

    contains

        subroutine attach_solvers(n)
        !! attach the nonlinear (fixed point) or linear (band) solver
        integer,intent(in) :: n  !! problem size
        if (.not. newton) then
            call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(n, 2), flag)
            call cv%set_max_nonlin_iters(10, flag)
        else
            call cv%set_linear_solver(sun_linsol_band(int(n, ip)), &
                                      band_matrix(int(n, ip), 0_ip, 0_ip), flag)
            call cv%set_jac_fn(ode_rhs_jac, flag)
        end if
        end subroutine attach_solvers

        subroutine save_history(t_n, y_n, step)
        !! shift the history and insert the newest (t, y, f)
        real(wp),intent(in) :: t_n  !! time of the step
        real(wp),dimension(:),intent(in) :: y_n  !! solution of the step
        integer,intent(in) :: step  !! step number
        integer :: i_start, k  !! last shifted entry (0-based) and index
        i_start = hist_size - 2
        if (step < hist_size) i_start = step - 1
        do k = i_start, 0, -1
            t_hist(k + 2) = t_hist(k + 1)
        end do
        t_hist(1) = t_n
        do k = i_start, 0, -1
            y_hist(:,k + 2) = y_hist(:,k + 1)
            f_hist(:,k + 2) = f_hist(:,k + 1)
        end do
        y_hist(:,1) = y_n
        flag = rhs(t_n, y_n, f_hist(:,1))
        end subroutine save_history

        subroutine resize_history()
        !! grow the history vectors by one element (all set to the first)
        real(wp),dimension(:,:),allocatable :: y_new  !! the new solution history
        integer :: new_size, k  !! new size and history index
        new_size = size(y_hist, 1) + 1
        allocate(y_new(new_size, hist_size))
        do k = 1, hist_size
            y_new(:,k) = y_hist(1,k)
        end do
        call move_alloc(y_new, y_hist)
        deallocate(f_hist)
        allocate(f_hist(new_size, hist_size))
        do k = 1, hist_size
            flag = rhs(t_hist(k), y_hist(:,k), f_hist(:,k))
        end do
        end subroutine resize_history

    end subroutine run_resize

    end module resize_driver
!*****************************************************************************************

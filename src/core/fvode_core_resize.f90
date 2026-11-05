!*****************************************************************************************
!>
!  CVODE history resizing (SUNDIALS v7.9.0 `cvode_resize.c`).
!
!  [[CVodeResizeHistory]] changes the size of the ODE system between
!  steps (e.g. after a mesh change): the user gives the solution and
!  right-hand side history at the new size, and the Nordsieck array is
!  rebuilt from it by polynomial interpolation (Newton form of the f
!  history for Adams, Hermite form of the y history for BDF). The linear
!  solver (and other size-dependent objects such as vector tolerances
!  and constraints) must be set again by the user; the default Newton
!  solver is recreated at the new size.

    submodule (fvode_core) fvode_core_resize

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeResizeHistory`: resize the integrator to the size of `y_hist`,
!  and rebuild the history array from the given history:
!
!  * `t_hist(1:)`: the times, `t_hist(1)` is the current time `tn`;
!  * `y_hist(:,1:)`: the solution at these times (Adams: at least 2
!    columns; BDF: at least `min(q+1, qmax)`);
!  * `f_hist(:,1:)`: the right-hand side at these times (Adams: at least
!    `min(q+1, qmax)` columns; BDF: at least 2).

    module subroutine CVodeResizeHistory(me, t_hist, y_hist, f_hist, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: t_hist  !! times of the history
    real(wp),dimension(:,:),intent(in) :: y_hist  !! solution history (new size)
    real(wp),dimension(:,:),intent(in) :: f_hist  !! right-hand side history (new size)
    integer,intent(out) :: flag  !! `CV_SUCCESS`, `CV_ILL_INPUT` or `CV_MEM_FAIL`

    integer :: n_hist  !! number of history values needed for the current (next) order
    integer :: nnew  !! the new problem size
    integer :: wrk_space_size  !! number of work vectors
    integer :: i, j  !! indices
    real(wp),dimension(:,:),allocatable :: resize_wrk  !! work vectors (0-based columns)

    ! check that the input history is sufficient for the current (next) order
    n_hist = min(me%q + 1, me%qmax)
    if (size(t_hist) < n_hist) then
        ! (C has no size information and reads t_hist as needed)
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeResizeHistory', 'Insufficient time history')
        flag = CV_ILL_INPUT
        return
    end if
    if (me%lmm == CV_ADAMS) then
        if (size(y_hist, 2) < 2) then
            call cvProcessError(me, CV_ILL_INPUT, 'CVodeResizeHistory', 'Insufficient solution history')
            flag = CV_ILL_INPUT
            return
        end if
        if (size(f_hist, 2) < n_hist) then
            call cvProcessError(me, CV_ILL_INPUT, 'CVodeResizeHistory', &
                                'Insufficient right-hand side history')
            flag = CV_ILL_INPUT
            return
        end if
    else
        if (size(f_hist, 2) < 2) then
            call cvProcessError(me, CV_ILL_INPUT, 'CVodeResizeHistory', &
                                'Insufficient right-hand side history')
            flag = CV_ILL_INPUT
            return
        end if
        if (size(y_hist, 2) < n_hist) then
            call cvProcessError(me, CV_ILL_INPUT, 'CVodeResizeHistory', 'Insufficient solution history')
            flag = CV_ILL_INPUT
            return
        end if
    end if
    nnew = size(y_hist, 1)
    if (size(f_hist, 1) /= nnew) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeResizeHistory', &
                            'The solution and right-hand side histories differ in size')
        flag = CV_ILL_INPUT
        return
    end if

    !--------------- resize vectors

    call resize(me%ewt)
    call resize(me%acor)
    call resize(me%tempv)
    call resize(me%ftemp)
    call resize(me%vtemp1)
    call resize(me%vtemp2)
    call resize(me%vtemp3)
    call resize(me%y)

    ! user will need to set a new vector of absolute tolerances
    if (allocated(me%Vabstol)) call resize(me%Vabstol)

    ! user will need to set a new constraints vector
    if (allocated(me%constraints)) deallocate(me%constraints)

    deallocate(me%zn)
    allocate(me%zn(nnew, 0:me%qmax_alloc), source=ZERO)

    me%n = nnew

    !--------------- resize nonlinear solver

    if (allocated(me%NLS) .and. me%ownNLS) then
        call CVodeSetNonlinearSolverInternal(me, sun_nonlinsol_newton(nnew), flag)
        if (flag /= CV_SUCCESS) then
            call cvProcessError(me, CV_MEM_FAIL, 'CVodeResizeHistory', &
                                'Error attaching default Newton solver')
            flag = CV_MEM_FAIL
            return
        end if
        me%ownNLS = .true.
    end if

    !--------------- create workspace for resizing

    wrk_space_size = max(me%q, me%qprime)
    if (me%lmm == CV_BDF) wrk_space_size = wrk_space_size + 1
    allocate(resize_wrk(nnew, 0:wrk_space_size-1))

    ! construct Nordsieck array at the old time but with the new size to
    ! compute correction vector at the new state size
    if (me%q < me%qmax) then

        ! compute z_{n-1} with new history size
        if (me%lmm == CV_ADAMS) then
            call build_nordsieck_adams(t_hist(2:), y_hist(:,2), f_hist(:,2:), me%q)
        else
            call build_nordsieck_bdf(t_hist(2:), y_hist(:,2:), f_hist(:,2), me%q)
        end if

        ! get predicted value
        me%vtemp1 = me%zn(:,0)
        do j = 1, me%q
            me%vtemp1 = me%zn(:,j) + me%vtemp1
        end do

        ! resized correction
        me%zn(:,me%qmax) = y_hist(:,1) - me%vtemp1

    end if

    ! construct new Nordsieck array
    if (me%lmm == CV_ADAMS) then
        call build_nordsieck_adams(t_hist, y_hist(:,1), f_hist, me%qprime)
    else
        call build_nordsieck_bdf(t_hist, y_hist, f_hist(:,1), me%qprime)
    end if

    ! ensure internal time and step history match the input history
    me%tn = t_hist(1)
    do i = 1, n_hist - 1
        me%tau(i) = t_hist(i) - t_hist(i + 1)
    end do

    ! in the next step, perform initialization needed after a resize
    me%first_step_after_resize = .true.

    flag = CV_SUCCESS

    contains

        subroutine resize(v)
        !! reallocate a vector at the new size (C: `N_VDestroy` + `N_VClone`)
        real(wp),dimension(:),allocatable,intent(inout) :: v  !! the vector
        if (allocated(v)) deallocate(v)
        allocate(v(nnew), source=ZERO)
        end subroutine resize

        subroutine build_nordsieck_adams(t, y, f, order)
        !! `cvBuildNordsieckArrayAdams`: the Adams Nordsieck array from the
        !! f history (Newton polynomial) and the current y
        real(wp),dimension(0:),intent(in) :: t  !! times (0-based)
        real(wp),dimension(:),intent(in) :: y  !! the current solution
        real(wp),dimension(:,0:),intent(in) :: f  !! the f history (0-based columns)
        integer,intent(in) :: order  !! the order
        integer :: i, j  !! indices
        real(wp) :: delta_t  !! inverse of a time difference
        real(wp) :: scale  !! hscale**i / i!
        associate (zn => me%zn, wrk => resize_wrk)
        if (order > 1) then
            ! compute Newton polynomial coefficients interpolating f history
            do i = 0, order - 1
                wrk(:,i) = f(:,i)
            end do
            do i = 1, order - 1
                do j = order - 1, i, -1
                    ! divided difference
                    delta_t = ONE / (t(j - i) - t(j))
                    wrk(:,j) = delta_t * (wrk(:,j - 1) - wrk(:,j))
                end do
            end do
            ! compute derivatives of Newton polynomial of f history
            zn(:,1) = wrk(:,order - 1)
            do i = 2, order
                zn(:,i) = ZERO
            end do
            do i = order - 2, 0, -1
                do j = order - 1, 1, -1
                    zn(:,j + 1) = (t(0) - t(i)) * zn(:,j + 1) + real(j, wp) * zn(:,j)
                end do
                zn(:,1) = (t(0) - t(i)) * zn(:,1) + wrk(:,i)
            end do
        end if
        ! overwrite first two columns with input values
        zn(:,0) = y
        zn(:,1) = f(:,0)
        ! scale entries
        scale = ONE
        do i = 1, order
            scale = scale * (me%hscale / real(i, wp))
            zn(:,i) = scale * zn(:,i)
        end do
        end associate
        end subroutine build_nordsieck_adams

        subroutine build_nordsieck_bdf(t, y, f, order)
        !! `cvBuildNordsieckArrayBDF`: the BDF Nordsieck array from the y
        !! history and the current f (Hermite polynomial)
        real(wp),dimension(0:),intent(in) :: t  !! times (0-based)
        real(wp),dimension(:,0:),intent(in) :: y  !! the y history (0-based columns)
        real(wp),dimension(:),intent(in) :: f  !! the current right-hand side
        integer,intent(in) :: order  !! the order
        integer :: i, j  !! indices
        real(wp) :: delta_t  !! inverse of a time difference
        real(wp) :: scale  !! hscale**i / i!
        real(wp),dimension(0:BDF_Q_MAX) :: t_ext  !! times, with t(0) repeated for the derivative value
        associate (zn => me%zn, wrk => resize_wrk)
        if (order > 1) then
            ! setup extended array of times to incorporate derivative value
            t_ext(0) = t(0)
            do i = 1, order
                t_ext(i) = t(i - 1)
            end do
            ! compute Hermite polynomial coefficients interpolating y history and f
            wrk(:,0) = y(:,0)
            do i = 1, order
                wrk(:,i) = y(:,i - 1)
            end do
            do i = 1, order
                do j = order, i, -1
                    if (i == 1 .and. j == 1) then
                        ! replace with actual derivative value
                        wrk(:,j) = f
                    else
                        ! divided difference
                        delta_t = ONE / (t_ext(j - i) - t_ext(j))
                        wrk(:,j) = delta_t * (wrk(:,j - 1) - wrk(:,j))
                    end if
                end do
            end do
            ! compute derivatives of Hermite polynomial
            zn(:,0) = wrk(:,order)
            do i = 1, order
                zn(:,i) = ZERO
            end do
            do i = order - 1, 0, -1
                do j = order, 1, -1
                    zn(:,j) = (t_ext(0) - t_ext(i)) * zn(:,j) + real(j, wp) * zn(:,j - 1)
                end do
                zn(:,0) = (t_ext(0) - t_ext(i)) * zn(:,0) + wrk(:,i)
            end do
        end if
        ! overwrite first two columns with input values
        zn(:,0) = y(:,0)
        zn(:,1) = f
        ! scale entries
        scale = ONE
        do i = 1, order
            scale = scale * (me%hscale / real(i, wp))
            zn(:,i) = scale * zn(:,i)
        end do
        end associate
        end subroutine build_nordsieck_bdf

    end subroutine CVodeResizeHistory
!*****************************************************************************************

    end submodule fvode_core_resize
!*****************************************************************************************

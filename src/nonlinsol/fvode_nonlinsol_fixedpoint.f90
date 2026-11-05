!*****************************************************************************************
!>
!  Fixed-point nonlinear solver (SUNDIALS v7.9.0 `SUNNonlinSol_FixedPoint`).
!
!  Translated from `src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c`.
!
!  @note Anderson acceleration (`m > 0`) is not ported yet. Constructing the
!  solver with `m > 0` is allowed, but [[fixedpoint_initialize]] then returns
!  `SUN_ERR_NOT_IMPLEMENTED`.

    module fvode_nonlinsol_fixedpoint

    use fvode_kinds, only: wp, li
    use fvode_parameters, only: SUNNONLINEARSOLVER_FIXEDPOINT, SUN_SUCCESS, &
                                SUN_NLS_SWITCH, SUN_NLS_CONTINUE, SUN_NLS_CONV_RECVR, &
                                SUN_ERR_NOT_IMPLEMENTED, SUN_ERR_ARG_OUTOFRANGE
    use fvode_nonlinsol, only: sun_nonlinear_solver, nls_problem

    implicit none

    private

    real(wp),parameter :: ONE = 1.0_wp

    type,extends(sun_nonlinear_solver),public :: sun_nonlinsol_fixedpoint
        !! Fixed-point iteration (`SUNNonlinearSolverContent_FixedPoint`).
        private
        integer :: m = 0               !! number of acceleration vectors (0: none)
        logical :: damping = .false.   !! apply damping (Anderson only)
        real(wp) :: beta = ONE         !! damping parameter (Anderson only)
        integer :: curiter = 0         !! current iteration number in a solve
        integer :: maxiters = 3        !! maximum iterations per solve
        integer(li) :: niters = 0      !! total iterations in the last solve
        integer(li) :: nconvfails = 0  !! convergence failures in the last solve
        real(wp),dimension(:),allocatable :: yprev  !! previous iterate
        real(wp),dimension(:),allocatable :: gy     !! G(y)
        real(wp),dimension(:),allocatable :: delta  !! update (y - yprev)
    contains
        procedure :: get_type           => fixedpoint_get_type
        procedure :: initialize         => fixedpoint_initialize
        procedure :: solve              => fixedpoint_solve
        procedure :: set_max_iters      => fixedpoint_set_max_iters
        procedure :: get_num_iters      => fixedpoint_get_num_iters
        procedure :: get_cur_iter       => fixedpoint_get_cur_iter
        procedure :: get_num_conv_fails => fixedpoint_get_num_conv_fails
        procedure :: set_damping        => fixedpoint_set_damping
    end type sun_nonlinsol_fixedpoint

    interface sun_nonlinsol_fixedpoint
        module procedure fixedpoint_create
    end interface sun_nonlinsol_fixedpoint

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNNonlinSol_FixedPoint(y, m)`): `n` is the vector length
!  and `m` the number of Anderson acceleration vectors (0 for none).

    function fixedpoint_create(n, m) result(nls)

    integer,intent(in) :: n
    integer,intent(in) :: m
    type(sun_nonlinsol_fixedpoint) :: nls

    nls%m          = m
    nls%damping    = .false.
    nls%beta       = ONE
    nls%curiter    = 0
    nls%maxiters   = 3
    nls%niters     = 0
    nls%nconvfails = 0
    allocate(nls%yprev(n), nls%gy(n), nls%delta(n))

    end function fixedpoint_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetType_FixedPoint`

    function fixedpoint_get_type(me) result(itype)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me
    integer :: itype

    itype = SUNNONLINEARSOLVER_FIXEDPOINT

    end function fixedpoint_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolInitialize_FixedPoint`: reset the counters.

    function fixedpoint_initialize(me) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me
    integer :: flag

    if (me%m > 0) then
        flag = SUN_ERR_NOT_IMPLEMENTED   ! Anderson acceleration: not ported yet
        return
    end if

    me%niters     = 0
    me%nconvfails = 0
    flag = SUN_SUCCESS

    end function fixedpoint_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSolve_FixedPoint`: fixed-point iteration `ycor <- G(ycor)`
!  until the convergence test passes or `maxiters` is reached.

    function fixedpoint_solve(me, prob, y0, ycor, w, tol, call_lsetup) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me
    class(nls_problem),intent(inout) :: prob
    real(wp),dimension(:),intent(in) :: y0       !! unused (as in C)
    real(wp),dimension(:),intent(inout) :: ycor
    real(wp),dimension(:),intent(in) :: w
    real(wp),intent(in) :: tol
    logical,intent(in) :: call_lsetup                   !! unused (as in C)
    integer :: flag

    integer :: retval

    ! reset counter
    me%niters     = 0
    me%nconvfails = 0

    ! looping point for attempts at solution of the nonlinear system:
    ! evaluate fixed-point function (store in gy), perform the update,
    ! and test for convergence
    me%curiter = 0
    do while (me%curiter < me%maxiters)

        ! update previous solution guess
        me%yprev = ycor

        ! compute fixed-point iteration function, store in gy
        retval = prob%nls_sys(ycor, me%gy)
        if (retval /= 0) then
            flag = retval
            return
        end if

        ! perform fixed point update
        if (me%m == 0) then
            ycor = me%gy                          ! update solution with gy
        else
            flag = SUN_ERR_NOT_IMPLEMENTED        ! Anderson acceleration
            return
        end if

        ! increment nonlinear solver iteration counter
        me%niters = me%niters + 1

        ! compute change in solution
        me%delta = ycor - me%yprev

        ! test for convergence
        retval = prob%nls_convtest(me, ycor, me%delta, tol, w)

        ! return if successful
        if (retval == 0) then
            flag = SUN_SUCCESS
            return
        else if (retval == SUN_NLS_SWITCH) then
            flag = SUN_NLS_SWITCH
            return
        end if

        ! check if the iterations should continue; otherwise increment the
        ! convergence failure count and return error flag
        if (retval /= SUN_NLS_CONTINUE) then
            me%nconvfails = me%nconvfails + 1
            flag = retval
            return
        end if

        me%curiter = me%curiter + 1
    end do

    ! if we've reached this point, then we exhausted the iteration limit;
    ! increment the convergence failure count and return
    me%nconvfails = me%nconvfails + 1
    flag = SUN_NLS_CONV_RECVR

    end function fixedpoint_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSetMaxIters_FixedPoint` (`maxiters >= 1`)

    function fixedpoint_set_max_iters(me, maxiters) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me
    integer,intent(in) :: maxiters
    integer :: flag

    if (maxiters < 1) then
        flag = SUN_ERR_ARG_OUTOFRANGE
        return
    end if
    me%maxiters = maxiters
    flag = SUN_SUCCESS

    end function fixedpoint_set_max_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSetDamping_FixedPoint` (only used with Anderson acceleration)

    function fixedpoint_set_damping(me, beta) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me
    real(wp),intent(in) :: beta
    integer :: flag

    if (beta <= 0.0_wp) then
        flag = SUN_ERR_ARG_OUTOFRANGE
        return
    end if

    if (beta < ONE) then
        ! enable damping
        me%beta    = beta
        me%damping = .true.
    else
        ! disable damping
        me%beta    = ONE
        me%damping = .false.
    end if
    flag = SUN_SUCCESS

    end function fixedpoint_set_damping
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetNumIters_FixedPoint`

    function fixedpoint_get_num_iters(me) result(val)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me
    integer(li) :: val

    val = me%niters

    end function fixedpoint_get_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetCurIter_FixedPoint`

    function fixedpoint_get_cur_iter(me) result(val)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me
    integer :: val

    val = me%curiter

    end function fixedpoint_get_cur_iter
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetNumConvFails_FixedPoint`

    function fixedpoint_get_num_conv_fails(me) result(val)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me
    integer(li) :: val

    val = me%nconvfails

    end function fixedpoint_get_num_conv_fails
!*****************************************************************************************

    end module fvode_nonlinsol_fixedpoint
!*****************************************************************************************

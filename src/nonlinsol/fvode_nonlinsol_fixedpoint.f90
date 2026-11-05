!*****************************************************************************************
!>
!  Fixed-point nonlinear solver with optional Anderson acceleration
!  (SUNDIALS v7.9.0 `SUNNonlinSol_FixedPoint`).
!
!  Translated from `src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c`.
!
!  @note The solver object is copied when it is attached to CVODE, so set
!  options (`m`, damping) before attaching it.

    module fvode_nonlinsol_fixedpoint

    use fvode_kinds, only: wp, li
    use fvode_parameters, only: SUNNONLINEARSOLVER_FIXEDPOINT, SUN_SUCCESS, &
                                SUN_NLS_SWITCH, SUN_NLS_CONTINUE, SUN_NLS_CONV_RECVR, &
                                SUN_ERR_ARG_OUTOFRANGE
    use fvode_nonlinsol, only: sun_nonlinear_solver, nls_problem

    implicit none

    private

    real(wp),parameter :: ZERO = 0.0_wp  !! 0.0
    real(wp),parameter :: ONE  = 1.0_wp  !! 1.0

    type,extends(sun_nonlinear_solver),public :: sun_nonlinsol_fixedpoint
        !! Fixed-point iteration (`SUNNonlinearSolverContent_FixedPoint`).
        private
        integer :: m = 0               !! number of acceleration vectors (0: none)
        logical :: damping = .false.   !! apply damping in the acceleration
        real(wp) :: beta = ONE         !! damping parameter
        integer :: curiter = 0         !! current iteration number in a solve
        integer :: maxiters = 3        !! maximum iterations per solve
        integer(li) :: niters = 0      !! total iterations in the last solve
        integer(li) :: nconvfails = 0  !! convergence failures in the last solve
        real(wp),dimension(:),allocatable :: yprev  !! previous iterate
        real(wp),dimension(:),allocatable :: gy     !! G(y)
        real(wp),dimension(:),allocatable :: delta  !! update (y - yprev); also fv in AA
        ! Anderson acceleration (m > 0)
        integer,dimension(:),allocatable :: imap    !! index map, 0:m-1
        real(wp),dimension(:,:),allocatable :: R    !! QR factor R(0:m-1,0:m-1)
        real(wp),dimension(:),allocatable :: gamma  !! least-squares coefficients, 0:m-1
        real(wp),dimension(:),allocatable :: fold   !! previous residual
        real(wp),dimension(:),allocatable :: gold   !! previous G(y)
        real(wp),dimension(:,:),allocatable :: df   !! residual differences, (n,0:m-1)
        real(wp),dimension(:,:),allocatable :: dg   !! G differences, (n,0:m-1)
        real(wp),dimension(:,:),allocatable :: q    !! QR factor Q, (n,0:m-1)
    contains
        procedure :: get_type           => fixedpoint_get_type  !! `SUNNonlinSolGetType`
        procedure :: initialize         => fixedpoint_initialize  !! `SUNNonlinSolInitialize`
        procedure :: solve              => fixedpoint_solve  !! `SUNNonlinSolSolve`
        procedure :: set_max_iters      => fixedpoint_set_max_iters  !! `SUNNonlinSolSetMaxIters`
        procedure :: get_num_iters      => fixedpoint_get_num_iters  !! `SUNNonlinSolGetNumIters`
        procedure :: get_cur_iter       => fixedpoint_get_cur_iter  !! `SUNNonlinSolGetCurIter`
        procedure :: get_num_conv_fails => fixedpoint_get_num_conv_fails  !! `SUNNonlinSolGetNumConvFails`
        procedure :: set_damping        => fixedpoint_set_damping  !! `SUNNonlinSolSetDamping_FixedPoint`
        procedure,private :: anderson_accelerate  !! Anderson acceleration update
    end type sun_nonlinsol_fixedpoint

    interface sun_nonlinsol_fixedpoint
        module procedure fixedpoint_create
    end interface sun_nonlinsol_fixedpoint

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNNonlinSol_FixedPoint(y, m)` + `AllocateContent`):
!  `n` is the vector length, `m` the number of Anderson acceleration vectors
!  (0 for plain fixed-point iteration), and `beta` an optional damping
!  parameter (`SUNNonlinSolSetDamping_FixedPoint`).

    function fixedpoint_create(n, m, beta) result(nls)

    integer,intent(in) :: n  !! vector length
    integer,intent(in) :: m  !! number of Anderson acceleration vectors (0: none)
    real(wp),intent(in),optional :: beta  !! damping parameter, 0 < beta <= 1 (default 1: none)
    type(sun_nonlinsol_fixedpoint) :: nls  !! the new solver

    integer :: istat  !! status flag (unused)

    nls%m = max(m, 0)
    allocate(nls%yprev(n), nls%gy(n), nls%delta(n))

    ! allocate all m-dependent content
    if (nls%m > 0) then
        allocate(nls%fold(n), nls%gold(n))
        allocate(nls%imap(0:nls%m-1), nls%gamma(0:nls%m-1))
        allocate(nls%R(0:nls%m-1, 0:nls%m-1), source=ZERO)
        allocate(nls%df(n, 0:nls%m-1), nls%dg(n, 0:nls%m-1), nls%q(n, 0:nls%m-1))
    end if

    if (present(beta)) istat = nls%set_damping(beta)

    end function fixedpoint_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetType_FixedPoint`

    function fixedpoint_get_type(me) result(itype)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me  !! the fixed-point solver
    integer :: itype  !! solver type (SUNNONLINEARSOLVER_*)

    itype = SUNNONLINEARSOLVER_FIXEDPOINT

    end function fixedpoint_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolInitialize_FixedPoint`: reset the counters.

    function fixedpoint_initialize(me) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me  !! the fixed-point solver
    integer :: flag  !! return flag (0 on success)

    me%niters     = 0
    me%nconvfails = 0
    flag = SUN_SUCCESS

    end function fixedpoint_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSolve_FixedPoint`: fixed-point iteration `ycor <- G(ycor)`
!  (accelerated if `m > 0`) until the convergence test passes or `maxiters`
!  is reached.

    function fixedpoint_solve(me, prob, y0, ycor, w, tol, call_lsetup) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me  !! the fixed-point solver
    class(nls_problem),intent(inout),target :: prob  !! the problem (integrator) providing the system functions
    real(wp),dimension(:),intent(in) :: y0       !! unused (as in C)
    real(wp),dimension(:),intent(inout) :: ycor  !! correction: initial guess on input, solution on output
    real(wp),dimension(:),intent(in) :: w  !! error weights
    real(wp),intent(in) :: tol  !! convergence tolerance
    logical,intent(in) :: call_lsetup            !! unused (as in C)
    integer :: flag  !! return flag (0 on success)

    integer :: retval  !! return flag of the called function
    integer :: i  !! element index

    ! reset counter
    me%niters     = 0
    me%nconvfails = 0

    ! looping point for attempts at solution of the nonlinear system:
    ! evaluate fixed-point function (store in gy), perform the update,
    ! and test for convergence
    me%curiter = 0
    do while (me%curiter < me%maxiters)

        ! update previous solution guess
        call copy(ycor, me%yprev)

        ! compute fixed-point iteration function, store in gy
        retval = prob%nls_sys(ycor, me%gy)
        if (retval /= 0) then
            flag = retval
            return
        end if

        ! perform fixed point update, based on choice of acceleration or not
        if (me%m == 0) then
            call copy(me%gy, ycor)          ! basic fixed-point solver
        else
            call me%anderson_accelerate(me%gy, ycor, me%yprev, me%curiter)
        end if

        ! increment nonlinear solver iteration counter
        me%niters = me%niters + 1

        ! compute change in solution
        do i = 1, size(ycor)
            me%delta(i) = ycor(i) - me%yprev(i)
        end do

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

    contains

        pure subroutine copy(x, z)
        !! z = x with an explicit loop (C: N_VScale(ONE, x, z)); for the
        !! small vectors typical here this avoids a memmove call
        real(wp),dimension(:),intent(in) :: x  !! source
        real(wp),dimension(:),intent(out) :: z  !! destination
        integer :: j  !! element index
        do j = 1, size(x)
            z(j) = x(j)
        end do
        end subroutine copy

    end function fixedpoint_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `AndersonAccelerate`: compute the Anderson-accelerated iterate
!
!  `x = gval - sum_i gamma(i) dg(i) [- (1-beta) (fv - sum_i gamma(i) df(i))]`
!
!  where `fv = gval - xold` and `gamma` solves the least-squares problem
!  `min || fv - sum_i gamma(i) df(i) ||` using a QR factorization of the
!  residual differences `df`, updated by one column per iteration (modified
!  Gram-Schmidt) and, once `iter > m`, downdated with Givens rotations.
!
!  C stores `R` column-major as `R[j*m + i]`; here that is `R(i,j)`.

    subroutine anderson_accelerate(me, gval, x, xold, iter)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me  !! the fixed-point solver
    real(wp),dimension(:),intent(in) :: gval     !! G(xold)
    real(wp),dimension(:),intent(inout) :: x     !! new iterate (also used as work vector)
    real(wp),dimension(:),intent(in) :: xold     !! previous iterate
    integer,intent(in) :: iter                   !! iteration number (0-based)

    integer :: i_pt  !! index of the newest column in the circular df/dg storage
    integer :: i  !! loop index
    integer :: j  !! loop index
    integer :: k  !! position in the index map
    integer :: lAA  !! current dimension of the acceleration subspace, min(iter, m)
    integer :: maa  !! maximum dimension of the acceleration subspace (m)
    real(wp) :: a  !! Givens rotation input (first element)
    real(wp) :: b  !! Givens rotation input (second element)
    real(wp) :: rtemp  !! temporary (rotation result)
    real(wp) :: c  !! Givens rotation cosine
    real(wp) :: s  !! Givens rotation sine
    real(wp) :: onembeta  !! 1 - beta (damping)

    associate (fv => me%delta, R => me%R, Q => me%q, df => me%df, dg => me%dg, &
               ipt_map => me%imap, gamma => me%gamma, vtemp => x)

    maa = me%m

    ! reset ipt_map, the index map for the current subspace
    ipt_map = 0
    i_pt = iter - 1 - ((iter - 1) / maa) * maa

    ! update dg[i_pt], df[i_pt], fv, gold and fold
    fv = gval - xold
    if (iter > 0) then
        dg(:,i_pt) = gval - me%gold   ! dg_new = gval - gold
        df(:,i_pt) = fv - me%fold     ! df_new = fv - fold
    end if
    me%gold = gval
    me%fold = fv

    ! on the first iteration, just do basic fixed-point update
    if (iter == 0) then
        x = gval
        return
    end if

    ! update QR factorization
    if (iter == 1) then
        ! second iteration
        R(0,0) = sqrt(dot_product(df(:,i_pt), df(:,i_pt)))
        Q(:,i_pt) = (ONE / R(0,0)) * df(:,i_pt)
        ipt_map(0) = 0

    else if (iter <= maa) then
        ! another iteration before we've reached maa
        vtemp = df(:,i_pt)
        do j = 0, iter - 2
            ipt_map(j) = j
            R(j,iter-1) = dot_product(Q(:,j), vtemp)
            vtemp = vtemp - R(j,iter-1) * Q(:,j)
        end do
        R(iter-1,iter-1) = sqrt(dot_product(vtemp, vtemp))
        if (R(iter-1,iter-1) == ZERO) then
            Q(:,i_pt) = ZERO
        else
            Q(:,i_pt) = (ONE / R(iter-1,iter-1)) * vtemp
        end if
        ipt_map(iter-1) = iter - 1

    else
        ! we've filled the acceleration subspace, so start recycling:
        ! delete the left-most column vector from the QR factorization
        do i = 0, maa - 2
            a     = R(i,i+1)
            b     = R(i+1,i+1)
            rtemp = sqrt(a*a + b*b)
            c     = a / rtemp
            s     = b / rtemp
            R(i,i+1)   = rtemp
            R(i+1,i+1) = ZERO
            do j = i + 2, maa - 1
                a = R(i,j)
                b = R(i+1,j)
                rtemp    = c*a + s*b
                R(i+1,j) = -s*a + c*b
                R(i,j)   = rtemp
            end do
            vtemp    = c*Q(:,i) + s*Q(:,i+1)
            Q(:,i+1) = -s*Q(:,i) + c*Q(:,i+1)
            Q(:,i)   = vtemp
        end do

        ! shift R to the left by one column
        do i = 1, maa - 1
            R(0:maa-2,i-1) = R(0:maa-2,i)
        end do

        ! add the new df vector
        vtemp = df(:,i_pt)
        do j = 0, maa - 2
            R(j,maa-1) = dot_product(Q(:,j), vtemp)
            vtemp = vtemp - R(j,maa-1) * Q(:,j)
        end do
        R(maa-1,maa-1) = sqrt(dot_product(vtemp, vtemp))
        Q(:,maa-1) = (ONE / R(maa-1,maa-1)) * vtemp

        ! update the iteration map
        k = 0
        do i = i_pt + 1, maa - 1
            ipt_map(k) = i
            k = k + 1
        end do
        do i = 0, i_pt
            ipt_map(k) = i
            k = k + 1
        end do
    end if

    ! solve least squares problem and update solution
    lAA = min(iter, maa)
    do i = 0, lAA - 1
        gamma(i) = dot_product(fv, Q(:,i))
    end do

    ! back substitution with R, accumulating the new iterate
    x = gval
    do i = lAA - 1, 0, -1
        do j = i + 1, lAA - 1
            gamma(i) = gamma(i) - R(i,j) * gamma(j)
        end do
        if (gamma(i) /= ZERO) gamma(i) = gamma(i) / R(i,i)
        x = x - gamma(i) * dg(:,ipt_map(i))
    end do

    ! if enabled, apply damping
    if (me%damping) then
        onembeta = ONE - me%beta
        x = x - onembeta * fv
        do i = lAA - 1, 0, -1
            x = x + (onembeta * gamma(i)) * df(:,ipt_map(i))
        end do
    end if

    end associate

    end subroutine anderson_accelerate
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSetMaxIters_FixedPoint` (`maxiters >= 1`)

    function fixedpoint_set_max_iters(me, maxiters) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me  !! the fixed-point solver
    integer,intent(in) :: maxiters  !! maximum iterations per solve (>= 1)
    integer :: flag  !! return flag (0 on success)

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
!  `SUNNonlinSolSetDamping_FixedPoint`: damping parameter `0 < beta <= 1`
!  for Anderson acceleration (`beta = 1` disables damping).

    function fixedpoint_set_damping(me, beta) result(flag)

    class(sun_nonlinsol_fixedpoint),intent(inout) :: me  !! the fixed-point solver
    real(wp),intent(in) :: beta  !! damping parameter, 0 < beta <= 1 (1 disables damping)
    integer :: flag  !! return flag (0 on success)

    if (beta <= ZERO) then
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

    class(sun_nonlinsol_fixedpoint),intent(in) :: me  !! the fixed-point solver
    integer(li) :: val  !! iterations in the last solve

    val = me%niters

    end function fixedpoint_get_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetCurIter_FixedPoint`

    function fixedpoint_get_cur_iter(me) result(val)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me  !! the fixed-point solver
    integer :: val  !! current iteration number

    val = me%curiter

    end function fixedpoint_get_cur_iter
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetNumConvFails_FixedPoint`

    function fixedpoint_get_num_conv_fails(me) result(val)

    class(sun_nonlinsol_fixedpoint),intent(in) :: me  !! the fixed-point solver
    integer(li) :: val  !! convergence failures in the last solve

    val = me%nconvfails

    end function fixedpoint_get_num_conv_fails
!*****************************************************************************************

    end module fvode_nonlinsol_fixedpoint
!*****************************************************************************************

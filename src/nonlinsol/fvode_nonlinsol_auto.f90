!*****************************************************************************************
!>
!  Nonlinear solver that switches automatically between fixed-point and
!  Newton iteration (SUNDIALS v7.9.0 `SUNNonlinSol_Auto`), translated from
!  `src/sunnonlinsol/auto/sunnonlinsol_auto.c`.
!
!  The switching strategy is from Norsett and Thomsen, "Switching between
!  modified Newton and fix-point iteration for implicit ODE-solvers", BIT
!  26 (1986) 339-348: fixed point switches to Newton (after a failed
!  convergence test) when the convergence rate estimate of the integrator
!  reaches `fp_to_newt_threshold` (default 0.8), and Newton switches to
!  fixed point when its stiffness ratio `||F(y_{m+1})|| / ||delta_m||` is
!  below `newt_to_fp_threshold` (default 2.0). A switch is only allowed
!  after `fp_to_newt_delay` (1) or `newt_to_fp_delay` (10) solves with
!  the current solver.
!
!  In C the auto solver installs its own convergence test in the two
!  sub-solvers. Here the sub-solvers are given a proxy problem
!  ([[auto_proxy]]), which forwards everything to the integrator, selects
!  the residual or fixed-point form of the system function, and wraps the
!  integrator's convergence test with the switching logic.

    module fvode_nonlinsol_auto

    use fvode_kinds, only: wp, li
    use fvode_parameters, only: SUNNONLINEARSOLVER_HYBRID, SUN_SUCCESS, SUN_NLS_SWITCH, &
                                SUN_ERR_ARG_CORRUPT
    use fvode_nonlinsol, only: sun_nonlinear_solver, nls_problem
    use fvode_nonlinsol_newton, only: sun_nonlinsol_newton
    use fvode_nonlinsol_fixedpoint, only: sun_nonlinsol_fixedpoint

    implicit none

    private

    integer,parameter,public :: SUNNONLINSOL_AUTO_FIXEDPOINT = 0  !! `SUNNonlinSolAutoType`: fixed-point iteration
    integer,parameter,public :: SUNNONLINSOL_AUTO_NEWTON     = 1  !! `SUNNonlinSolAutoType`: Newton iteration

    integer,parameter :: SUNNLS_AUTO_MAX_SOLVE_ATTEMPTS = 3  !! maximum number of switches within one solve (+1)
    real(wp),parameter :: SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_THRESHOLD = 2.0_wp  !! default Newton -> fixed point threshold on the stiffness ratio
    real(wp),parameter :: SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_THRESHOLD = 0.8_wp  !! default fixed point -> Newton threshold on the convergence rate
    integer(li),parameter :: SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_DELAY = 10  !! default Newton -> fixed point delay (solves)
    integer(li),parameter :: SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_DELAY = 1   !! default fixed point -> Newton delay (solves)

    type,extends(nls_problem) :: auto_proxy
        !! The problem seen by the two sub-solvers: forwards to the
        !! integrator, and holds the switching state (C: the content of the
        !! auto solver and `SUNNonlinSolAutoConvTestData`).
        class(nls_problem),pointer :: prob => null()  !! the integrator (set for the duration of a solve)
        integer :: active_solver_type = SUNNONLINSOL_AUTO_NEWTON  !! the active sub-solver
        real(wp) :: fp_to_newt_threshold = SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_THRESHOLD  !! switch to Newton when crate >= this
        real(wp) :: newt_to_fp_threshold = SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_THRESHOLD  !! switch to fixed point when stiffr < this
        integer(li) :: fp_to_newt_delay = SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_DELAY  !! solves before a switch to Newton is allowed
        integer(li) :: newt_to_fp_delay = SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_DELAY  !! solves before a switch to fixed point is allowed
        integer(li) :: num_solves_since_switch = 0  !! number of solves since the last switch
        integer(li) :: switch_count = 0  !! total number of switches
    contains
        procedure :: nls_sys => proxy_sys  !! residual or fixed-point function, for the active sub-solver
        procedure :: nls_lsetup => proxy_lsetup  !! forwarded to the integrator
        procedure :: nls_lsolve => proxy_lsolve  !! forwarded to the integrator
        procedure :: nls_convtest => proxy_convtest  !! the integrator's test plus the switching logic
        procedure :: nls_get_update_norm => proxy_get_update_norm  !! forwarded to the integrator
        procedure :: nls_get_conv_rate => proxy_get_conv_rate  !! forwarded to the integrator
        procedure :: ls_atimes => proxy_atimes  !! forwarded to the integrator
        procedure :: ls_psetup => proxy_psetup  !! forwarded to the integrator
        procedure :: ls_psolve => proxy_psolve  !! forwarded to the integrator
    end type auto_proxy

    type,extends(sun_nonlinear_solver),public :: sun_nonlinsol_auto
        !! Automatic fixed-point/Newton switching solver (`SUNNonlinearSolverContent_Auto`).
        type(sun_nonlinsol_fixedpoint) :: fp_solver  !! the fixed-point sub-solver
        type(sun_nonlinsol_newton) :: newton_solver  !! the Newton sub-solver (computes the stiffness ratio)
        type(auto_proxy) :: proxy  !! the problem given to the sub-solvers, with the switching state
        integer(li) :: num_iters = 0  !! iterations in the last solve
        integer(li) :: num_conv_fails = 0  !! convergence failures in the last solve
        integer(li) :: fp_num_iters_total = 0  !! total fixed-point iterations
        integer(li) :: newton_num_iters_total = 0  !! total Newton iterations
        integer(li) :: fp_num_conv_fails_total = 0  !! total fixed-point convergence failures
        integer(li) :: newton_num_conv_fails_total = 0  !! total Newton convergence failures
    contains
        procedure :: get_type           => auto_get_type  !! `SUNNonlinSolGetType_Auto`: hybrid
        procedure :: initialize         => auto_initialize  !! `SUNNonlinSolInitialize_Auto`
        procedure :: solve              => auto_solve  !! `SUNNonlinSolSolve_Auto`
        procedure :: set_max_iters      => auto_set_max_iters  !! `SUNNonlinSolSetMaxIters_Auto`
        procedure :: get_num_iters      => auto_get_num_iters  !! `SUNNonlinSolGetNumIters_Auto`
        procedure :: get_cur_iter       => auto_get_cur_iter  !! `SUNNonlinSolGetCurIter_Auto`
        procedure :: get_num_conv_fails => auto_get_num_conv_fails  !! `SUNNonlinSolGetNumConvFails_Auto`
        procedure :: set_switching_parameters => auto_set_switching_parameters  !! `SUNNonlinSolSetSwitchingParameters_Auto`
        procedure :: get_active_solver_type => auto_get_active_solver_type  !! `SUNNonlinSolGetActiveSolverType_Auto`
        procedure :: get_switch_count => auto_get_switch_count  !! `SUNNonlinSolGetSwitchCount_Auto`
        procedure :: get_total_num_iters_by_type => auto_get_total_num_iters_by_type  !! `SUNNonlinSolGetTotalNumItersByType_Auto`
        procedure :: get_total_num_conv_fails_by_type => auto_get_total_num_conv_fails_by_type  !! `SUNNonlinSolGetTotalNumConvFailsByType_Auto`
    end type sun_nonlinsol_auto

    interface sun_nonlinsol_auto
        module procedure auto_create
    end interface sun_nonlinsol_auto

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNNonlinSol_Auto(y, m, initial_solver_type)`): `n` is
!  the vector length and `m` the Anderson acceleration depth of the
!  fixed-point sub-solver.

    function auto_create(n, m, initial_solver_type) result(nls)

    integer,intent(in) :: n  !! vector length
    integer,intent(in) :: m  !! Anderson acceleration depth of the fixed-point solver (0: none)
    integer,intent(in) :: initial_solver_type  !! `SUNNONLINSOL_AUTO_FIXEDPOINT` or `SUNNONLINSOL_AUTO_NEWTON`
    type(sun_nonlinsol_auto) :: nls  !! the solver

    integer :: flag  !! return flag

    nls%proxy%active_solver_type = initial_solver_type
    nls%fp_solver = sun_nonlinsol_fixedpoint(n, m)
    nls%newton_solver = sun_nonlinsol_newton(n)
    flag = nls%newton_solver%set_compute_stiffness_ratio(.true.)

    end function auto_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetType_Auto`

    function auto_get_type(me) result(itype)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer :: itype  !! `SUNNONLINEARSOLVER_HYBRID`

    itype = SUNNONLINEARSOLVER_HYBRID

    end function auto_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolInitialize_Auto`: initialize both sub-solvers (the linear
!  solver hooks go to the Newton solver, as in C `SetLSetupFn_Auto` and
!  `SetLSolveFn_Auto`).

    function auto_initialize(me) result(flag)

    class(sun_nonlinsol_auto),intent(inout) :: me  !! the solver
    integer :: flag  !! `SUN_SUCCESS` or the flag of a sub-solver

    me%newton_solver%has_lsetup = me%has_lsetup
    me%newton_solver%has_lsolve = me%has_lsolve

    flag = me%fp_solver%initialize()
    if (flag /= SUN_SUCCESS) return
    flag = me%newton_solver%initialize()

    end function auto_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSolve_Auto`: solve with the active sub-solver; if it asks
!  to switch (`SUN_NLS_SWITCH`), continue with the other one from the
!  current correction.

    function auto_solve(me, prob, y0, ycor, w, tol, call_lsetup) result(flag)

    class(sun_nonlinsol_auto),intent(inout) :: me  !! the solver
    class(nls_problem),intent(inout),target :: prob  !! the problem (integrator)
    real(wp),dimension(:),intent(in) :: y0  !! predicted solution
    real(wp),dimension(:),intent(inout) :: ycor  !! correction: initial guess on input, solution on output
    real(wp),dimension(:),intent(in) :: w  !! error weights
    real(wp),intent(in) :: tol  !! convergence tolerance
    logical,intent(in) :: call_lsetup  !! set up the linear solver before iterating
    integer :: flag  !! return flag of the sub-solver

    integer :: attempts  !! solve attempts (switches + 1)
    integer :: solve_solver_type  !! the sub-solver of this attempt
    integer(li) :: iters  !! iterations of the sub-solver
    integer(li) :: nconvfails  !! convergence failures of the sub-solver

    me%num_iters      = 0
    me%num_conv_fails = 0

    me%newton_solver%has_lsetup = me%has_lsetup
    me%newton_solver%has_lsolve = me%has_lsolve
    me%proxy%prob => prob

    do attempts = 0, SUNNLS_AUTO_MAX_SOLVE_ATTEMPTS - 1

        solve_solver_type = me%proxy%active_solver_type

        if (solve_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT) then
            flag = me%fp_solver%solve(me%proxy, y0, ycor, w, tol, call_lsetup)
            iters = me%fp_solver%get_num_iters()
            nconvfails = me%fp_solver%get_num_conv_fails()
            me%fp_num_iters_total = me%fp_num_iters_total + iters
            me%fp_num_conv_fails_total = me%fp_num_conv_fails_total + nconvfails
        else
            flag = me%newton_solver%solve(me%proxy, y0, ycor, w, tol, call_lsetup)
            iters = me%newton_solver%get_num_iters()
            nconvfails = me%newton_solver%get_num_conv_fails()
            me%newton_num_iters_total = me%newton_num_iters_total + iters
            me%newton_num_conv_fails_total = me%newton_num_conv_fails_total + nconvfails
        end if
        me%num_iters = me%num_iters + iters
        me%num_conv_fails = me%num_conv_fails + nconvfails

        if (flag == SUN_NLS_SWITCH) cycle

        me%proxy%num_solves_since_switch = me%proxy%num_solves_since_switch + 1
        me%proxy%prob => null()
        return

    end do

    me%proxy%prob => null()
    flag = SUN_SUCCESS

    end function auto_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSetMaxIters_Auto`: set the maximum iterations of both
!  sub-solvers.

    function auto_set_max_iters(me, maxiters) result(flag)

    class(sun_nonlinsol_auto),intent(inout) :: me  !! the solver
    integer,intent(in) :: maxiters  !! maximum iterations per solve
    integer :: flag  !! `SUN_SUCCESS` or the flag of a sub-solver

    flag = me%fp_solver%set_max_iters(maxiters)
    if (flag /= SUN_SUCCESS) return
    flag = me%newton_solver%set_max_iters(maxiters)

    end function auto_set_max_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetNumIters_Auto`: iterations in the last solve (both
!  sub-solvers).

    function auto_get_num_iters(me) result(val)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer(li) :: val  !! iterations in the last solve

    val = me%num_iters

    end function auto_get_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetCurIter_Auto`: the current iteration of the active
!  sub-solver.

    function auto_get_cur_iter(me) result(val)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer :: val  !! current iteration number

    if (me%proxy%active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT) then
        val = me%fp_solver%get_cur_iter()
    else
        val = me%newton_solver%get_cur_iter()
    end if

    end function auto_get_cur_iter
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetNumConvFails_Auto`: convergence failures in the last
!  solve (both sub-solvers).

    function auto_get_num_conv_fails(me) result(val)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer(li) :: val  !! convergence failures in the last solve

    val = me%num_conv_fails

    end function auto_get_num_conv_fails
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSetSwitchingParameters_Auto` (a negative value restores
!  the default of that parameter).

    function auto_set_switching_parameters(me, newt_to_fp_threshold, newt_to_fp_delay, &
                                           fp_to_newt_threshold, fp_to_newt_delay) result(flag)

    class(sun_nonlinsol_auto),intent(inout) :: me  !! the solver
    real(wp),intent(in) :: newt_to_fp_threshold  !! switch to fixed point when the stiffness ratio is below this
    integer(li),intent(in) :: newt_to_fp_delay  !! solves before a switch to fixed point is allowed
    real(wp),intent(in) :: fp_to_newt_threshold  !! switch to Newton when the convergence rate is at least this
    integer(li),intent(in) :: fp_to_newt_delay  !! solves before a switch to Newton is allowed
    integer :: flag  !! `SUN_SUCCESS`

    me%proxy%newt_to_fp_threshold = newt_to_fp_threshold
    if (newt_to_fp_threshold < 0.0_wp) me%proxy%newt_to_fp_threshold = SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_THRESHOLD
    me%proxy%newt_to_fp_delay = newt_to_fp_delay
    if (newt_to_fp_delay < 0) me%proxy%newt_to_fp_delay = SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_DELAY
    me%proxy%fp_to_newt_threshold = fp_to_newt_threshold
    if (fp_to_newt_threshold < 0.0_wp) me%proxy%fp_to_newt_threshold = SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_THRESHOLD
    me%proxy%fp_to_newt_delay = fp_to_newt_delay
    if (fp_to_newt_delay < 0) me%proxy%fp_to_newt_delay = SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_DELAY
    flag = SUN_SUCCESS

    end function auto_set_switching_parameters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetActiveSolverType_Auto`

    pure function auto_get_active_solver_type(me) result(val)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer :: val  !! `SUNNONLINSOL_AUTO_FIXEDPOINT` or `SUNNONLINSOL_AUTO_NEWTON`

    val = me%proxy%active_solver_type

    end function auto_get_active_solver_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetSwitchCount_Auto`

    pure function auto_get_switch_count(me) result(val)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer(li) :: val  !! total number of switches

    val = me%proxy%switch_count

    end function auto_get_switch_count
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetTotalNumItersByType_Auto`

    subroutine auto_get_total_num_iters_by_type(me, fp_iters, newt_iters)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer(li),intent(out) :: fp_iters  !! total fixed-point iterations
    integer(li),intent(out) :: newt_iters  !! total Newton iterations

    fp_iters   = me%fp_num_iters_total
    newt_iters = me%newton_num_iters_total

    end subroutine auto_get_total_num_iters_by_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolGetTotalNumConvFailsByType_Auto`

    subroutine auto_get_total_num_conv_fails_by_type(me, fp_nconvfails, newt_nconvfails)

    class(sun_nonlinsol_auto),intent(in) :: me  !! the solver
    integer(li),intent(out) :: fp_nconvfails  !! total fixed-point convergence failures
    integer(li),intent(out) :: newt_nconvfails  !! total Newton convergence failures

    fp_nconvfails   = me%fp_num_conv_fails_total
    newt_nconvfails = me%newton_num_conv_fails_total

    end subroutine auto_get_total_num_conv_fails_by_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  System function for the active sub-solver: the fixed-point form for
!  the fixed-point solver, the residual form for Newton (C:
!  `SUNNonlinSolSetSysFns_Auto`).

    function proxy_sys(me, ycor, res) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    real(wp),dimension(:),intent(in) :: ycor  !! current correction
    real(wp),dimension(:),intent(inout) :: res  !! system function value
    integer :: flag  !! return flag of the integrator

    if (me%active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT) then
        flag = me%prob%nls_fp_function(ycor, res)
    else
        flag = me%prob%nls_residual(ycor, res)
    end if

    end function proxy_sys
!*****************************************************************************************

!*****************************************************************************************
!>
!  Linear solver setup: forwarded to the integrator.

    function proxy_lsetup(me, jbad, jcur) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    logical,intent(in) :: jbad  !! the Jacobian is suspected to be out of date
    logical,intent(out) :: jcur  !! the Jacobian was re-evaluated
    integer :: flag  !! return flag of the integrator

    flag = me%prob%nls_lsetup(jbad, jcur)

    end function proxy_lsetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  Linear solve: forwarded to the integrator.

    function proxy_lsolve(me, b) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    real(wp),dimension(:),intent(inout) :: b  !! right-hand side on input, solution on output
    integer :: flag  !! return flag of the integrator

    flag = me%prob%nls_lsolve(b)

    end function proxy_lsolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolConvTest_Auto`: the integrator's convergence test, then
!  the switching decision (Norsett and Thomsen). Returns `SUN_NLS_SWITCH`
!  to continue the solve with the other sub-solver.

    function proxy_convtest(me, nls, ycor, del, tol, ewt) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    class(sun_nonlinear_solver),intent(in) :: nls  !! the active sub-solver
    real(wp),dimension(:),intent(in) :: ycor  !! current correction
    real(wp),dimension(:),intent(in) :: del  !! last update
    real(wp),intent(in) :: tol  !! convergence tolerance
    real(wp),dimension(:),intent(in) :: ewt  !! error weights
    integer :: flag  !! return flag

    integer :: retval  !! return flag of the integrator's test
    integer :: curiter  !! current iteration of the sub-solver
    integer :: ierr  !! return flag of the convergence rate
    real(wp) :: crate  !! convergence rate estimate of the integrator
    real(wp) :: stiffr  !! stiffness ratio of the Newton solver
    logical :: diverging  !! is fixed point diverging?
    logical :: contraction  !! is the fixed-point map a contraction?
    logical :: dont_delay  !! have enough solves been done since the last switch?

    retval = me%prob%nls_convtest(nls, ycor, del, tol, ewt)

    ! return early if convergence test error is unrecoverable
    if (retval < 0) then
        flag = retval
        return
    end if

    if (me%active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT) then

        ! if the integrator-provided convergence test passed, exit with
        ! success and don't consider switching
        if (retval == SUN_SUCCESS) then
            flag = SUN_SUCCESS
            return
        end if

        ! the convergence rate estimate is not meaningful on the first
        ! nonlinear iteration
        curiter = nls%get_cur_iter()
        if (curiter == 0) then
            flag = retval
            return
        end if

        ! get the convergence rate from the integrator
        ierr = me%prob%nls_get_conv_rate(crate)
        if (ierr /= SUN_SUCCESS) then
            flag = ierr
            return
        end if

        diverging  = crate >= me%fp_to_newt_threshold
        dont_delay = me%num_solves_since_switch >= me%fp_to_newt_delay
        if (diverging .and. dont_delay) then
            me%num_solves_since_switch = 0
            me%active_solver_type      = SUNNONLINSOL_AUTO_NEWTON
            me%switch_count            = me%switch_count + 1
            ! continue the solve with Newton
            flag = SUN_NLS_SWITCH
            return
        end if

    else

        ! Newton is active: check if we should switch to fixed point
        ! regardless of whether the convergence test passed
        select type (nls)
        type is (sun_nonlinsol_newton)
            stiffr = nls%get_stiffness_ratio()
        class default
            flag = SUN_ERR_ARG_CORRUPT
            return
        end select

        contraction = stiffr < me%newt_to_fp_threshold
        dont_delay  = me%num_solves_since_switch >= me%newt_to_fp_delay
        if (contraction .and. dont_delay) then
            me%num_solves_since_switch = 0
            me%active_solver_type      = SUNNONLINSOL_AUTO_FIXEDPOINT
            me%switch_count            = me%switch_count + 1
            ! if the test passed, stop now (fixed point is used from the
            ! next solve); otherwise continue this solve with fixed point
            if (retval == SUN_SUCCESS) then
                flag = SUN_SUCCESS
            else
                flag = SUN_NLS_SWITCH
            end if
            return
        end if

    end if

    flag = retval

    end function proxy_convtest
!*****************************************************************************************

!*****************************************************************************************
!>
!  Update norm: forwarded to the integrator.

    function proxy_get_update_norm(me, delnrm) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    real(wp),intent(out) :: delnrm  !! norm of the last update
    integer :: flag  !! return flag of the integrator

    flag = me%prob%nls_get_update_norm(delnrm)

    end function proxy_get_update_norm
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convergence rate: forwarded to the integrator.

    function proxy_get_conv_rate(me, crate) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    real(wp),intent(out) :: crate  !! estimated convergence rate
    integer :: flag  !! return flag of the integrator

    flag = me%prob%nls_get_conv_rate(crate)

    end function proxy_get_conv_rate
!*****************************************************************************************

!*****************************************************************************************
!>
!  A*v for iterative linear solvers: forwarded to the integrator.

    function proxy_atimes(me, v, z) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
    real(wp),dimension(:),intent(out) :: z  !! the product
    integer :: flag  !! return flag of the integrator

    flag = me%prob%ls_atimes(v, z)

    end function proxy_atimes
!*****************************************************************************************

!*****************************************************************************************
!>
!  Preconditioner setup: forwarded to the integrator.

    function proxy_psetup(me) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    integer :: flag  !! return flag of the integrator

    flag = me%prob%ls_psetup()

    end function proxy_psetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  Preconditioner solve: forwarded to the integrator.

    function proxy_psolve(me, r, z, tol, lr) result(flag)

    class(auto_proxy),intent(inout) :: me  !! the proxy
    real(wp),dimension(:),intent(in) :: r  !! the right-hand side
    real(wp),dimension(:),intent(out) :: z  !! the solution
    real(wp),intent(in) :: tol  !! tolerance
    integer,intent(in) :: lr  !! left or right preconditioner
    integer :: flag  !! return flag of the integrator

    flag = me%prob%ls_psolve(r, z, tol, lr)

    end function proxy_psolve
!*****************************************************************************************

    end module fvode_nonlinsol_auto
!*****************************************************************************************

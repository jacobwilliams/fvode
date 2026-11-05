!*****************************************************************************************
!>
!  Newton nonlinear solver (SUNDIALS v7.9.0 `SUNNonlinSol_Newton`).
!
!  Translated from `src/sunnonlinsol/newton/sunnonlinsol_newton.c`,
!  including the optional stiffness-ratio computation used by the "auto"
!  switching solver ([[sun_nonlinsol_auto]]).

    module fvode_nonlinsol_newton

    use fvode_kinds, only: wp, li
    use fvode_parameters, only: SUNNONLINEARSOLVER_ROOTFIND, SUN_SUCCESS, &
                                SUN_NLS_SWITCH, SUN_NLS_CONTINUE, SUN_NLS_CONV_RECVR, &
                                SUN_ERR_ARG_OUTOFRANGE, SUN_ERR_ARG_CORRUPT, &
                                SUN_ERR_NOT_IMPLEMENTED
    use fvode_nonlinsol, only: sun_nonlinear_solver, nls_problem
    use fvode_nvector, only: n_vwrmsnorm

    implicit none

    private

    type,extends(sun_nonlinear_solver),public :: sun_nonlinsol_newton
        !! Newton iteration (`SUNNonlinearSolverContent_Newton`).
        private
        logical :: jcur = .false.      !! Jacobian status, current = true / stale = false
        integer :: curiter = 0         !! current iteration number in a solve
        integer :: maxiters = 3        !! maximum iterations per solve
        integer(li) :: niters = 0      !! total iterations in the last solve
        integer(li) :: nconvfails = 0  !! convergence failures in the last solve
        real(wp),dimension(:),allocatable :: delta  !! Newton update vector
        logical :: compute_stiffr = .false.  !! compute the stiffness ratio?
        real(wp) :: stiffr = 0.0_wp  !! the last stiffness ratio ||F(y_{m+1})|| / ||delta_m||
        real(wp) :: delnrm = 0.0_wp  !! norm of the last Newton update
    contains
        procedure :: get_type           => newton_get_type  !! `SUNNonlinSolGetType`
        procedure :: initialize         => newton_initialize  !! `SUNNonlinSolInitialize`
        procedure :: solve              => newton_solve  !! `SUNNonlinSolSolve`
        procedure :: set_max_iters      => newton_set_max_iters  !! `SUNNonlinSolSetMaxIters`
        procedure :: get_num_iters      => newton_get_num_iters  !! `SUNNonlinSolGetNumIters`
        procedure :: get_cur_iter       => newton_get_cur_iter  !! `SUNNonlinSolGetCurIter`
        procedure :: get_num_conv_fails => newton_get_num_conv_fails  !! `SUNNonlinSolGetNumConvFails`
        procedure :: set_compute_stiffness_ratio => newton_set_compute_stiffr  !! `SUNNonlinSolSetComputeStiffnessRatio_Newton`
        procedure :: get_stiffness_ratio => newton_get_stiffr  !! `SUNNonlinSolGetStiffnessRatio_Newton`
    end type sun_nonlinsol_newton

    interface sun_nonlinsol_newton
        module procedure newton_create
    end interface sun_nonlinsol_newton

    contains
!*****************************************************************************************

    function newton_create(n) result(nls)
    !! `SUNNonlinSol_Newton(y)` for vectors of length n.
    integer,intent(in) :: n  !! vector length
    type(sun_nonlinsol_newton) :: nls  !! the new solver
    allocate(nls%delta(n))
    end function newton_create

    function newton_get_type(me) result(itype)
    !! `SUNNonlinSolGetType_Newton`: a root-finding solver.
    class(sun_nonlinsol_newton),intent(in) :: me  !! the Newton solver
    integer :: itype  !! solver type (SUNNONLINEARSOLVER_*)
    itype = SUNNONLINEARSOLVER_ROOTFIND
    end function newton_get_type

!*****************************************************************************************
!>
!  `SUNNonlinSolInitialize_Newton`: reset counters. Fails if no linear solve
!  function is attached (C asserts this in debug builds).

    function newton_initialize(me) result(flag)

    class(sun_nonlinsol_newton),intent(inout) :: me  !! the Newton solver
    integer :: flag  !! return flag (0 on success)

    if (.not. me%has_lsolve) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if

    me%niters     = 0
    me%nconvfails = 0
    me%stiffr     = 0.0_wp
    me%jcur       = .false.
    flag = SUN_SUCCESS

    end function newton_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNNonlinSolSolve_Newton`: Newton iteration for the correction `ycor`.
!
!  Each iteration solves `J delta = -F(ycor)` and updates `ycor += delta`.
!  On a recoverable failure with a stale Jacobian, the linear solver is set
!  up again (with `jbad = .true.`) and the iteration restarts from zero.

    function newton_solve(me, prob, y0, ycor, w, tol, call_lsetup) result(flag)

    class(sun_nonlinsol_newton),intent(inout) :: me  !! the Newton solver
    class(nls_problem),intent(inout),target :: prob  !! the problem (integrator) providing the system functions
    real(wp),dimension(:),intent(in) :: y0      !! unused (as in C)
    real(wp),dimension(:),intent(inout) :: ycor  !! correction: initial guess on input, solution on output
    real(wp),dimension(:),intent(in) :: w  !! error weights
    real(wp),intent(in) :: tol  !! convergence tolerance
    logical,intent(in) :: call_lsetup  !! set up the linear solver before iterating
    integer :: flag  !! return flag (0 on success)

    integer :: retval  !! return flag of the called function
    integer :: ierr  !! return flag of the update norm
    logical :: jbad  !! the Jacobian is suspected to be out of date
    logical :: callLSetup  !! set up the linear solver in this pass
    real(wp) :: delnrm  !! norm of the last update (stiffness ratio)
    real(wp) :: resnrm  !! norm of the new residual (stiffness ratio)

    if (call_lsetup .and. .not. me%has_lsetup) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if

    callLSetup = call_lsetup

    ! assume the Jacobian is good
    jbad = .false.

    ! initialize counter
    me%niters     = 0
    me%nconvfails = 0

    ! looping point for attempts at solution of the nonlinear system:
    !   evaluate the nonlinear residual function (store in delta),
    !   setup the linear solver if necessary,
    !   perform Newton iteration
    setup_loop: do

        ! compute the nonlinear residual, store in delta
        me%curiter = 0
        retval = prob%nls_sys(ycor, me%delta)
        if (retval /= SUN_SUCCESS) exit setup_loop

        ! if indicated, setup the linear system
        if (callLSetup) then
            retval = prob%nls_lsetup(jbad, me%jcur)
            if (retval /= SUN_SUCCESS) exit setup_loop
        end if

        ! looping point for Newton iteration. Break out on any error.
        newton_loop: do

            ! increment nonlinear solver iteration counter
            me%niters = me%niters + 1

            ! compute the negative of the residual for the linear system rhs
            me%delta = -me%delta

            ! solve the linear system to get Newton update delta
            retval = prob%nls_lsolve(me%delta)
            if (retval /= SUN_SUCCESS) exit newton_loop

            ! update the Newton iterate
            ycor = ycor + me%delta

            ! test for convergence
            retval = prob%nls_convtest(me, ycor, me%delta, tol, w)

            ! update here to keep curiter correct for the convergence test
            me%curiter = me%curiter + 1

            ! retrieve the update norm if it is available, otherwise compute it
            ierr = prob%nls_get_update_norm(me%delnrm)
            if (ierr == SUN_ERR_NOT_IMPLEMENTED) then
                me%delnrm = n_vwrmsnorm(me%delta, w)
            else if (ierr /= SUN_SUCCESS) then
                flag = ierr
                return
            end if

            ! if successful update Jacobian status and return
            if (retval == SUN_SUCCESS) then
                me%jcur = .false.
                flag = SUN_SUCCESS
                return
            else if (retval == SUN_NLS_SWITCH) then
                flag = SUN_NLS_SWITCH
                return
            end if

            ! check if the iteration should continue; otherwise exit Newton loop
            if (retval /= SUN_NLS_CONTINUE) exit newton_loop

            ! not yet converged. Increment curiter and test for max allowed.
            if (me%curiter >= me%maxiters) then
                retval = SUN_NLS_CONV_RECVR
                exit newton_loop
            end if

            ! compute the nonlinear residual, store in delta
            retval = prob%nls_sys(ycor, me%delta)
            if (retval /= SUN_SUCCESS) exit newton_loop

            if (me%compute_stiffr) then
                delnrm = me%delnrm
                resnrm = n_vwrmsnorm(me%delta, w)   ! C: norm_fn (cvNlsNorm) or N_VWrmsNorm
                ! Norsett's switching metric compares the next residual to
                ! the previous Newton update norm
                if (delnrm > 0.0_wp) then
                    me%stiffr = resnrm / delnrm
                else
                    me%stiffr = 0.0_wp
                end if
            end if

        end do newton_loop

        ! If there is a recoverable convergence failure and the Jacobian-related
        ! data appears not to be current, increment the convergence failure count,
        ! reset the initial correction to zero, and loop again with a call to
        ! lsetup in which jbad is TRUE. Otherwise break out and return.
        if ((retval > 0) .and. (.not. me%jcur) .and. me%has_lsetup) then
            me%nconvfails = me%nconvfails + 1
            callLSetup = .true.
            jbad = .true.
            ycor = 0.0_wp
            cycle setup_loop
        else
            exit setup_loop
        end if

    end do setup_loop

    ! increment number of convergence failures
    me%nconvfails = me%nconvfails + 1

    ! all errors go here
    flag = retval

    end function newton_solve
!*****************************************************************************************

    function newton_set_max_iters(me, maxiters) result(flag)
    !! `SUNNonlinSolSetMaxIters_Newton` (`maxiters >= 1`)
    class(sun_nonlinsol_newton),intent(inout) :: me  !! the Newton solver
    integer,intent(in) :: maxiters  !! maximum iterations per solve (>= 1)
    integer :: flag  !! return flag (0 on success)
    if (maxiters < 1) then
        flag = SUN_ERR_ARG_OUTOFRANGE
        return
    end if
    me%maxiters = maxiters
    flag = SUN_SUCCESS
    end function newton_set_max_iters

    function newton_get_num_iters(me) result(val)
    !! `SUNNonlinSolGetNumIters_Newton`
    class(sun_nonlinsol_newton),intent(in) :: me  !! the Newton solver
    integer(li) :: val  !! iterations in the last solve
    val = me%niters
    end function newton_get_num_iters

    function newton_get_cur_iter(me) result(val)
    !! `SUNNonlinSolGetCurIter_Newton`
    class(sun_nonlinsol_newton),intent(in) :: me  !! the Newton solver
    integer :: val  !! current iteration number
    val = me%curiter
    end function newton_get_cur_iter

    function newton_get_num_conv_fails(me) result(val)
    !! `SUNNonlinSolGetNumConvFails_Newton`
    class(sun_nonlinsol_newton),intent(in) :: me  !! the Newton solver
    integer(li) :: val  !! convergence failures in the last solve
    val = me%nconvfails
    end function newton_get_num_conv_fails

    function newton_set_compute_stiffr(me, onoff) result(flag)
    !! `SUNNonlinSolSetComputeStiffnessRatio_Newton`: compute the stiffness
    !! ratio `||F(y_{m+1})|| / ||delta_m||` in each iteration?
    class(sun_nonlinsol_newton),intent(inout) :: me  !! the Newton solver
    logical,intent(in) :: onoff  !! compute the ratio?
    integer :: flag  !! `SUN_SUCCESS`
    me%compute_stiffr = onoff
    if (.not. onoff) me%stiffr = 0.0_wp
    flag = SUN_SUCCESS
    end function newton_set_compute_stiffr

    pure function newton_get_stiffr(me) result(val)
    !! `SUNNonlinSolGetStiffnessRatio_Newton`: the last stiffness ratio
    class(sun_nonlinsol_newton),intent(in) :: me  !! the Newton solver
    real(wp) :: val  !! the stiffness ratio
    val = me%stiffr
    end function newton_get_stiffr

    end module fvode_nonlinsol_newton
!*****************************************************************************************

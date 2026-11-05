!*****************************************************************************************
!>
!  Scaled, preconditioned Bi-CGStab linear solver (SUNDIALS v7.9.0
!  `SUNLinSol_SPBCGS`), translated from
!  `src/sunlinsol/spbcgs/sunlinsol_spbcgs.c`.
!
!  Solves `A x = b` with at most `maxl` Bi-CGStab iterations (no restarts),
!  with optional left and/or right preconditioning and scaling: the method
!  is applied to `(Sb P1^{-1} A P2^{-1} Sx^{-1}) (Sx P2 x) = Sb P1^{-1} b`,
!  where `Sb = s1` and `Sx = s2`. Right preconditioning requires a zero
!  initial guess.

    module fvode_linsol_spbcgs

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUNLINEARSOLVER_ITERATIVE, SUN_SUCCESS, SUN_ERR_ARG_CORRUPT, &
                                SUN_ERR_ARG_INCOMPATIBLE, &
                                SUN_PREC_NONE, SUN_PREC_LEFT, SUN_PREC_RIGHT, SUN_PREC_BOTH, &
                                SUNSPBCGS_MAXL_DEFAULT, &
                                SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, &
                                SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, &
                                SUNLS_PSET_FAIL_UNREC, SUNLS_PSET_FAIL_REC, &
                                SUNLS_RES_REDUCED, SUNLS_CONV_FAIL
    use fvode_matrix, only: sun_matrix
    use fvode_linsol, only: sun_linear_solver, ls_problem

    implicit none

    private

    real(wp),parameter :: ZERO = 0.0_wp  !! 0
    real(wp),parameter :: ONE  = 1.0_wp  !! 1

    type,extends(sun_linear_solver),public :: sun_linsol_spbcgs
        !! Bi-CGStab solver (`SUNLinearSolverContent_SPBCGS`).
        integer :: maxl = SUNSPBCGS_MAXL_DEFAULT  !! maximum number of iterations
        integer :: pretype = SUN_PREC_NONE  !! preconditioning type (`SUN_PREC_*`)
        integer :: numiters = 0  !! number of iterations in the last solve
        real(wp) :: resnorm = ZERO  !! final (scaled) residual norm of the last solve
        real(wp),dimension(:),allocatable :: r_star  !! the initial (scaled, preconditioned) residual
        real(wp),dimension(:),allocatable :: r  !! the current residual
        real(wp),dimension(:),allocatable :: p  !! the search direction
        real(wp),dimension(:),allocatable :: q  !! r - alpha*Ap
        real(wp),dimension(:),allocatable :: u  !! A-tilde q
        real(wp),dimension(:),allocatable :: Ap  !! A-tilde p
        real(wp),dimension(:),allocatable :: vtemp  !! temporary vector
    contains
        procedure :: get_type   => spbcgs_get_type  !! `SUNLinSolGetType_SPBCGS`: iterative solver
        procedure :: initialize => spbcgs_initialize  !! `SUNLinSolInitialize_SPBCGS`: check options
        procedure :: setup      => spbcgs_setup  !! `SUNLinSolSetup_SPBCGS`: preconditioner setup
        procedure :: solve      => spbcgs_solve  !! `SUNLinSolSolve_SPBCGS`: Bi-CGStab iteration
        procedure :: num_iters  => spbcgs_num_iters  !! `SUNLinSolNumIters_SPBCGS`
        procedure :: res_norm   => spbcgs_res_norm  !! `SUNLinSolResNorm_SPBCGS`
        procedure :: supports_atimes => spbcgs_supports  !! SPBCGS uses an A*v function
        procedure :: supports_preconditioner => spbcgs_supports  !! SPBCGS supports preconditioning
        procedure :: supports_scaling_vectors => spbcgs_supports  !! SPBCGS supports scaling vectors
        procedure :: set_prec_type => spbcgs_set_prec_type  !! `SUNLinSol_SPBCGSSetPrecType`
        procedure :: set_maxl => spbcgs_set_maxl  !! `SUNLinSol_SPBCGSSetMaxl`
    end type sun_linsol_spbcgs

    interface sun_linsol_spbcgs
        module procedure spbcgs_create
    end interface sun_linsol_spbcgs

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNLinSol_SPBCGS(y, pretype, maxl)`): `n` is the vector
!  length. An illegal `pretype` means no preconditioning; `maxl <= 0`
!  means the default (5).

    function spbcgs_create(n, pretype, maxl) result(S)

    integer,intent(in) :: n  !! vector length
    integer,intent(in) :: pretype  !! preconditioning type
    integer,intent(in) :: maxl  !! maximum number of iterations
    type(sun_linsol_spbcgs) :: S  !! the solver

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        S%pretype = SUN_PREC_NONE
    else
        S%pretype = pretype
    end if
    S%maxl = maxl
    if (maxl <= 0) S%maxl = SUNSPBCGS_MAXL_DEFAULT

    allocate(S%r_star(n), S%r(n), S%p(n), S%q(n), S%u(n), S%Ap(n), S%vtemp(n))

    end function spbcgs_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver type: `SUNLINEARSOLVER_ITERATIVE`.

    pure function spbcgs_get_type(me) result(itype)

    class(sun_linsol_spbcgs),intent(in) :: me  !! the solver
    integer :: itype  !! `SUNLINEARSOLVER_ITERATIVE`

    itype = SUNLINEARSOLVER_ITERATIVE

    end function spbcgs_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  SPBCGS provides the A*v, preconditioner and scaling-vector operations.

    pure function spbcgs_supports(me) result(supported)

    class(sun_linsol_spbcgs),intent(in) :: me  !! the solver
    logical :: supported  !! always true

    supported = .true.

    end function spbcgs_supports
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPBCGSSetPrecType`

    function spbcgs_set_prec_type(me, pretype) result(flag)

    class(sun_linsol_spbcgs),intent(inout) :: me  !! the solver
    integer,intent(in) :: pretype  !! `SUN_PREC_NONE`, `_LEFT`, `_RIGHT` or `_BOTH`
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_CORRUPT`

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if
    me%pretype = pretype
    flag = SUN_SUCCESS

    end function spbcgs_set_prec_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPBCGSSetMaxl` (<= 0 restores the default, 5).

    function spbcgs_set_maxl(me, maxl) result(flag)

    class(sun_linsol_spbcgs),intent(inout) :: me  !! the solver
    integer,intent(in) :: maxl  !! maximum number of iterations
    integer :: flag  !! `SUN_SUCCESS`

    me%maxl = maxl
    if (maxl <= 0) me%maxl = SUNSPBCGS_MAXL_DEFAULT
    flag = SUN_SUCCESS

    end function spbcgs_set_maxl
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolInitialize_SPBCGS`: check the options (no memory to allocate).

    function spbcgs_initialize(me) result(flag)

    class(sun_linsol_spbcgs),intent(inout) :: me  !! the solver
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_CORRUPT`

    if (me%maxl <= 0) me%maxl = SUNSPBCGS_MAXL_DEFAULT
    if ((me%pretype /= SUN_PREC_LEFT) .and. (me%pretype /= SUN_PREC_RIGHT) .and. &
        (me%pretype /= SUN_PREC_BOTH)) me%pretype = SUN_PREC_NONE

    if (.not. me%has_atimes .or. &
        (me%pretype /= SUN_PREC_NONE .and. .not. me%has_psolve)) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if

    flag = SUN_SUCCESS

    end function spbcgs_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetup_SPBCGS`: call the preconditioner setup function (if any).

    function spbcgs_setup(me, A, prob) result(flag)

    class(sun_linsol_spbcgs),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used (matrix-free)
    class(ls_problem),intent(inout),optional :: prob  !! the problem (preconditioner setup)
    integer :: flag  !! `SUN_SUCCESS` or a `SUNLS_PSET_FAIL_*` flag

    integer :: status  !! return value of the preconditioner setup

    if (me%has_psetup .and. present(prob)) then
        status = prob%ls_psetup()
        if (status /= 0) then
            me%last_flag = merge(SUNLS_PSET_FAIL_UNREC, SUNLS_PSET_FAIL_REC, status < 0)
            flag = me%last_flag
            return
        end if
    end if

    me%last_flag = SUN_SUCCESS
    flag = SUN_SUCCESS

    end function spbcgs_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSolve_SPBCGS`: solve `A x = b` to the scaled residual tolerance
!  `delta`. Returns `SUN_SUCCESS` (converged), `SUNLS_RES_REDUCED` (not
!  converged, but the residual was reduced), `SUNLS_CONV_FAIL`, or a
!  failure flag of the A*v / preconditioner operations.

    function spbcgs_solve(me, A, x, b, tol, prob) result(flag)

    class(sun_linsol_spbcgs),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used (matrix-free)
    real(wp),dimension(:),intent(inout) :: x  !! initial guess (unless zeroguess) / solution
    real(wp),dimension(:),intent(in) :: b  !! right-hand side
    real(wp),intent(in) :: tol  !! tolerance delta on the scaled residual 2-norm
    class(ls_problem),intent(inout),optional :: prob  !! the problem (A*v, preconditioner)
    integer :: flag  !! return flag (see above)

    real(wp) :: alpha  !! step length along p
    real(wp) :: beta  !! update factor of the search direction
    real(wp) :: omega  !! step length along q
    real(wp) :: omega_denom  !! <u,u>
    real(wp) :: beta_num  !! <r_new, r_star>
    real(wp) :: beta_denom  !! <r_old, r_star>
    real(wp) :: r_norm  !! initial residual norm
    real(wp) :: rho  !! current residual norm
    real(wp) :: delta  !! tolerance
    logical :: preOnLeft, preOnRight  !! left / right preconditioning?
    logical :: scale_x, scale_b  !! scaling vectors in use?
    logical :: converged  !! has the iteration converged?
    integer :: l  !! iteration index
    integer :: l_max  !! maximum number of iterations
    integer :: status  !! return value of A*v / psolve

    delta = tol

    associate (r_star => me%r_star, r => me%r, p => me%p, q => me%q, u => me%u, &
               Ap => me%Ap, vtemp => me%vtemp)

    ! initialize counters and convergence flag
    l_max = me%maxl
    me%numiters = 0
    converged = .false.

    ! set booleans for internal solver options
    preOnLeft  = (me%pretype == SUN_PREC_LEFT) .or. (me%pretype == SUN_PREC_BOTH)
    preOnRight = (me%pretype == SUN_PREC_RIGHT) .or. (me%pretype == SUN_PREC_BOTH)
    scale_x = allocated(me%s2)
    scale_b = allocated(me%s1)

    ! check for unsupported use case
    if (preOnRight .and. .not. me%zeroguess) then
        flag = fail(SUN_ERR_ARG_INCOMPATIBLE)
        return
    end if

    ! check if Atimes function has been set, and psolve if preconditioning
    if (.not. present(prob) .or. .not. me%has_atimes .or. &
        ((preOnLeft .or. preOnRight) .and. .not. me%has_psolve)) then
        flag = fail(SUN_ERR_ARG_CORRUPT)
        return
    end if

    ! set r_star to initial (unscaled) residual r_0 = b - A*x_0
    if (me%zeroguess) then
        r_star = b
    else
        status = prob%ls_atimes(x, r_star)
        if (status /= 0) then
            flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
            return
        end if
        r_star = b - r_star
    end if

    ! apply left preconditioner and b-scaling to r_star = r_0
    if (preOnLeft) then
        status = prob%ls_psolve(r_star, r, delta, SUN_PREC_LEFT)
        if (status /= 0) then
            flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
            return
        end if
    else
        r = r_star
    end if

    if (scale_b) then
        r_star = me%s1 * r
    else
        r_star = r
    end if

    ! initialize beta_denom to the dot product of r0 with r0
    beta_denom = dot_product(r_star, r_star)

    ! set r_norm to L2 norm of r_star = sb P1_inv r_0, and return if small
    r_norm = sqrt(beta_denom)
    rho = r_norm
    me%resnorm = r_norm
    if (r_norm <= delta) then
        flag = fail(SUN_SUCCESS)
        return
    end if

    ! copy r_star to r and p
    r = r_star
    p = r_star

    ! set x = sx x if non-zero guess
    if (scale_x .and. .not. me%zeroguess) x = me%s2 * x

    ! begin main iteration loop
    do l = 0, l_max - 1

        me%numiters = me%numiters + 1

        ! generate Ap = A-tilde p, where A-tilde = sb P1_inv A P2_inv sx_inv

        ! apply x-scaling: vtemp = sx_inv p
        if (scale_x) then
            vtemp = p / me%s2
        else
            vtemp = p
        end if

        ! apply right preconditioner: vtemp = P2_inv sx_inv p
        if (preOnRight) then
            Ap = vtemp
            status = prob%ls_psolve(Ap, vtemp, delta, SUN_PREC_RIGHT)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
        end if

        ! apply A: Ap = A P2_inv sx_inv p
        status = prob%ls_atimes(vtemp, Ap)
        if (status /= 0) then
            flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
            return
        end if

        ! apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p
        if (preOnLeft) then
            status = prob%ls_psolve(Ap, vtemp, delta, SUN_PREC_LEFT)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
        else
            vtemp = Ap
        end if

        ! apply b-scaling: Ap = sb P1_inv A P2_inv sx_inv p
        if (scale_b) then
            Ap = me%s1 * vtemp
        else
            Ap = vtemp
        end if

        ! calculate alpha = <r,r_star>/<Ap,r_star>
        alpha = dot_product(Ap, r_star)
        alpha = beta_denom / alpha

        ! update q = r - alpha*Ap = r - alpha*(sb P1_inv A P2_inv sx_inv p)
        q = r - alpha * Ap

        ! generate u = A-tilde q

        ! apply x-scaling: vtemp = sx_inv q
        if (scale_x) then
            vtemp = q / me%s2
        else
            vtemp = q
        end if

        ! apply right preconditioner: vtemp = P2_inv sx_inv q
        if (preOnRight) then
            u = vtemp
            status = prob%ls_psolve(u, vtemp, delta, SUN_PREC_RIGHT)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
        end if

        ! apply A: u = A P2_inv sx_inv u
        status = prob%ls_atimes(vtemp, u)
        if (status /= 0) then
            flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
            return
        end if

        ! apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p
        if (preOnLeft) then
            status = prob%ls_psolve(u, vtemp, delta, SUN_PREC_LEFT)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
        else
            vtemp = u
        end if

        ! apply b-scaling: u = sb P1_inv A P2_inv sx_inv u
        if (scale_b) then
            u = me%s1 * vtemp
        else
            u = vtemp
        end if

        ! calculate omega = <u,q>/<u,u>
        omega_denom = dot_product(u, u)
        if (omega_denom == ZERO) omega_denom = ONE
        omega = dot_product(u, q)
        omega = omega / omega_denom

        ! update x = x + alpha*p + omega*q
        if (l == 0 .and. me%zeroguess) then
            x = alpha * p + omega * q
        else
            x = x + alpha * p + omega * q
        end if

        ! update the residual r = q - omega*u
        r = q - omega * u

        ! set rho = norm(r) and check convergence
        rho = sqrt(dot_product(r, r))
        me%resnorm = rho
        if (rho <= delta) then
            converged = .true.
            exit
        end if

        ! not yet converged, continue iteration
        ! update beta = <rnew,r_star> / <rold,r_start> * alpha / omega
        beta_num = dot_product(r, r_star)
        beta = ((beta_num / beta_denom) * (alpha / omega))

        ! update p = r + beta*(p - omega*Ap) = beta*p - beta*omega*Ap + r
        ! (C: the coefficient of Ap is -alpha*(beta_num/beta_denom))
        p = beta * p + (-alpha * (beta_num / beta_denom)) * Ap + r

        ! update beta_denom for next iteration
        beta_denom = beta_num

    end do

    ! main loop finished
    if (converged .or. (rho < r_norm)) then

        ! apply the x-scaling and right preconditioner: x = P2_inv sx_inv x
        if (scale_x) x = x / me%s2
        if (preOnRight) then
            status = prob%ls_psolve(x, vtemp, delta, SUN_PREC_RIGHT)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
            x = vtemp
        end if

        if (converged) then
            flag = fail(SUN_SUCCESS)
        else
            flag = fail(SUNLS_RES_REDUCED)
        end if

    else
        flag = fail(SUNLS_CONV_FAIL)
    end if

    end associate

    contains

        function fail(f) result(ff)
        !! return: reset zeroguess and record the flag
        integer,intent(in) :: f  !! the return flag
        integer :: ff  !! the same flag
        me%zeroguess = .false.
        me%last_flag = f
        ff = f
        end function fail

    end function spbcgs_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolNumIters_SPBCGS`

    pure function spbcgs_num_iters(me) result(n)

    class(sun_linsol_spbcgs),intent(in) :: me  !! the solver
    integer :: n  !! number of iterations in the last solve

    n = me%numiters

    end function spbcgs_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolResNorm_SPBCGS`

    pure function spbcgs_res_norm(me) result(r)

    class(sun_linsol_spbcgs),intent(in) :: me  !! the solver
    real(wp) :: r  !! final scaled residual norm of the last solve

    r = me%resnorm

    end function spbcgs_res_norm
!*****************************************************************************************

    end module fvode_linsol_spbcgs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Scaled, preconditioned transpose-free QMR linear solver (SUNDIALS
!  v7.9.0 `SUNLinSol_SPTFQMR`), translated from
!  `src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c`.
!
!  Solves `A x = b` with at most `maxl` TFQMR iterations (no restarts),
!  with optional left and/or right preconditioning and scaling: the method
!  is applied to `(Sb P1^{-1} A P2^{-1} Sx^{-1}) (Sx P2 x) = Sb P1^{-1} b`,
!  where `Sb = s1` and `Sx = s2`. Right preconditioning requires a zero
!  initial guess.

    module fvode_linsol_sptfqmr

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUNLINEARSOLVER_ITERATIVE, SUN_SUCCESS, SUN_ERR_ARG_CORRUPT, &
                                SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_ARG_OUTOFRANGE, &
                                SUN_PREC_NONE, SUN_PREC_LEFT, SUN_PREC_RIGHT, SUN_PREC_BOTH, &
                                SUNSPTFQMR_MAXL_DEFAULT, &
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

    type,extends(sun_linear_solver),public :: sun_linsol_sptfqmr
        !! TFQMR solver (`SUNLinearSolverContent_SPTFQMR`).
        integer :: maxl = SUNSPTFQMR_MAXL_DEFAULT  !! maximum number of iterations
        integer :: pretype = SUN_PREC_NONE  !! preconditioning type (`SUN_PREC_*`)
        integer :: numiters = 0  !! number of iterations in the last solve
        real(wp) :: resnorm = ZERO  !! final (scaled) residual norm of the last solve
        real(wp),dimension(:),allocatable :: r_star  !! the initial (scaled, preconditioned) residual
        real(wp),dimension(:),allocatable :: q  !! u - alpha*v
        real(wp),dimension(:),allocatable :: d  !! the update direction of x
        real(wp),dimension(:),allocatable :: v  !! A-tilde p
        real(wp),dimension(:),allocatable :: p  !! the search direction
        real(wp),dimension(:,:),allocatable :: r  !! the residuals r(:,0:1) of the previous and current iteration
        real(wp),dimension(:),allocatable :: u  !! r(:,1) + beta*q
        real(wp),dimension(:),allocatable :: vtemp1  !! temporary vector
        real(wp),dimension(:),allocatable :: vtemp2  !! temporary vector
        real(wp),dimension(:),allocatable :: vtemp3  !! the scaled, preconditioned right-hand side
    contains
        procedure :: get_type   => sptfqmr_get_type  !! `SUNLinSolGetType_SPTFQMR`: iterative solver
        procedure :: initialize => sptfqmr_initialize  !! `SUNLinSolInitialize_SPTFQMR`: check options
        procedure :: setup      => sptfqmr_setup  !! `SUNLinSolSetup_SPTFQMR`: preconditioner setup
        procedure :: solve      => sptfqmr_solve  !! `SUNLinSolSolve_SPTFQMR`: TFQMR iteration
        procedure :: num_iters  => sptfqmr_num_iters  !! `SUNLinSolNumIters_SPTFQMR`
        procedure :: res_norm   => sptfqmr_res_norm  !! `SUNLinSolResNorm_SPTFQMR`
        procedure :: supports_atimes => sptfqmr_supports  !! SPTFQMR uses an A*v function
        procedure :: supports_preconditioner => sptfqmr_supports  !! SPTFQMR supports preconditioning
        procedure :: supports_scaling_vectors => sptfqmr_supports  !! SPTFQMR supports scaling vectors
        procedure :: set_prec_type => sptfqmr_set_prec_type  !! `SUNLinSol_SPTFQMRSetPrecType`
        procedure :: set_maxl => sptfqmr_set_maxl  !! `SUNLinSol_SPTFQMRSetMaxl`
    end type sun_linsol_sptfqmr

    interface sun_linsol_sptfqmr
        module procedure sptfqmr_create
    end interface sun_linsol_sptfqmr

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNLinSol_SPTFQMR(y, pretype, maxl)`): `n` is the vector
!  length. An illegal `pretype` means no preconditioning; `maxl <= 0`
!  means the default (5).

    function sptfqmr_create(n, pretype, maxl) result(S)

    integer,intent(in) :: n  !! vector length
    integer,intent(in) :: pretype  !! preconditioning type
    integer,intent(in) :: maxl  !! maximum number of iterations
    type(sun_linsol_sptfqmr) :: S  !! the solver

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        S%pretype = SUN_PREC_NONE
    else
        S%pretype = pretype
    end if
    S%maxl = maxl
    if (maxl <= 0) S%maxl = SUNSPTFQMR_MAXL_DEFAULT

    allocate(S%r_star(n), S%q(n), S%d(n), S%v(n), S%p(n), S%r(n,0:1), S%u(n), &
             S%vtemp1(n), S%vtemp2(n), S%vtemp3(n))

    end function sptfqmr_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver type: `SUNLINEARSOLVER_ITERATIVE`.

    pure function sptfqmr_get_type(me) result(itype)

    class(sun_linsol_sptfqmr),intent(in) :: me  !! the solver
    integer :: itype  !! `SUNLINEARSOLVER_ITERATIVE`

    itype = SUNLINEARSOLVER_ITERATIVE

    end function sptfqmr_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  SPTFQMR provides the A*v, preconditioner and scaling-vector operations.

    pure function sptfqmr_supports(me) result(supported)

    class(sun_linsol_sptfqmr),intent(in) :: me  !! the solver
    logical :: supported  !! always true

    supported = .true.

    end function sptfqmr_supports
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPTFQMRSetPrecType`

    function sptfqmr_set_prec_type(me, pretype) result(flag)

    class(sun_linsol_sptfqmr),intent(inout) :: me  !! the solver
    integer,intent(in) :: pretype  !! `SUN_PREC_NONE`, `_LEFT`, `_RIGHT` or `_BOTH`
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_OUTOFRANGE`

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        flag = SUN_ERR_ARG_OUTOFRANGE
        return
    end if
    me%pretype = pretype
    flag = SUN_SUCCESS

    end function sptfqmr_set_prec_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPTFQMRSetMaxl` (<= 0 restores the default, 5).

    function sptfqmr_set_maxl(me, maxl) result(flag)

    class(sun_linsol_sptfqmr),intent(inout) :: me  !! the solver
    integer,intent(in) :: maxl  !! maximum number of iterations
    integer :: flag  !! `SUN_SUCCESS`

    me%maxl = maxl
    if (maxl <= 0) me%maxl = SUNSPTFQMR_MAXL_DEFAULT
    flag = SUN_SUCCESS

    end function sptfqmr_set_maxl
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolInitialize_SPTFQMR`: check the options (no memory to allocate).

    function sptfqmr_initialize(me) result(flag)

    class(sun_linsol_sptfqmr),intent(inout) :: me  !! the solver
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_CORRUPT`

    if (me%maxl <= 0) me%maxl = SUNSPTFQMR_MAXL_DEFAULT
    if ((me%pretype /= SUN_PREC_LEFT) .and. (me%pretype /= SUN_PREC_RIGHT) .and. &
        (me%pretype /= SUN_PREC_BOTH)) me%pretype = SUN_PREC_NONE

    if (.not. me%has_atimes .or. &
        (me%pretype /= SUN_PREC_NONE .and. .not. me%has_psolve)) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if

    flag = SUN_SUCCESS

    end function sptfqmr_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetup_SPTFQMR`: call the preconditioner setup function (if any).

    function sptfqmr_setup(me, A, prob) result(flag)

    class(sun_linsol_sptfqmr),intent(inout) :: me  !! the solver
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

    end function sptfqmr_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSolve_SPTFQMR`: solve `A x = b` to the scaled residual
!  tolerance `delta`. Returns `SUN_SUCCESS` (converged), `SUNLS_RES_REDUCED`
!  (not converged, but the residual was reduced), `SUNLS_CONV_FAIL`, or a
!  failure flag of the A*v / preconditioner operations.

    function sptfqmr_solve(me, A, x, b, tol, prob) result(flag)

    class(sun_linsol_sptfqmr),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used (matrix-free)
    real(wp),dimension(:),intent(inout) :: x  !! initial guess (unless zeroguess) / solution
    real(wp),dimension(:),intent(in) :: b  !! right-hand side
    real(wp),intent(in) :: tol  !! tolerance delta on the scaled residual 2-norm
    class(ls_problem),intent(inout),optional :: prob  !! the problem (A*v, preconditioner)
    integer :: flag  !! return flag (see above)

    real(wp) :: alpha  !! rho(0)/sigma
    real(wp) :: tau  !! QMR quasi-residual norm
    real(wp) :: eta  !! step length along d
    real(wp) :: beta  !! rho(1)/rho(0)
    real(wp) :: c  !! cosine of the QMR rotation
    real(wp) :: sigma  !! <r_star, v>
    real(wp) :: v_bar  !! omega/tau
    real(wp) :: omega  !! the QMR weight
    real(wp),dimension(0:1) :: rho  !! <r_star, r> for the previous and current iteration
    real(wp) :: r_init_norm  !! initial residual norm
    real(wp) :: r_curr_norm  !! current residual norm (or its estimate)
    real(wp) :: temp_val  !! ||r(:,1)||, saved for the second inner iteration
    real(wp) :: delta  !! tolerance
    logical :: preOnLeft, preOnRight  !! left / right preconditioning?
    logical :: scale_x, scale_b  !! scaling vectors in use?
    logical :: converged  !! has the iteration converged?
    logical :: b_ok  !! has the scaled, preconditioned b been computed (in vtemp3)?
    integer :: n, m  !! outer and inner iteration indices
    integer :: l_max  !! maximum number of iterations
    integer :: status  !! return value of A*v / psolve

    delta = tol

    associate (r_star => me%r_star, q => me%q, d => me%d, v => me%v, p => me%p, &
               r => me%r, u => me%u, vtemp1 => me%vtemp1, vtemp2 => me%vtemp2, &
               vtemp3 => me%vtemp3)

    ! initialize counters and convergence flag
    l_max = me%maxl
    me%numiters = 0
    temp_val = -ONE
    r_curr_norm = -ONE
    converged = .false.
    b_ok = .false.

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

    ! set r_star to initial (unscaled) residual r_star = r_0 = b - A*x_0
    ! (if x == 0 then just set residual to b and continue)
    if (me%zeroguess) then
        r_star = b
    else
        status = prob%ls_atimes(x, r_star)
        if (status /= 0) then
            flag = fail(atimes_fail(status))
            return
        end if
        r_star = b - r_star
    end if

    ! apply left preconditioner and b-scaling to r_star (or really just r_0)
    if (preOnLeft) then
        status = prob%ls_psolve(r_star, vtemp1, delta, SUN_PREC_LEFT)
        if (status /= 0) then
            flag = fail(psolve_fail(status))
            return
        end if
    else
        vtemp1 = r_star
    end if

    if (scale_b) then
        r_star = me%s1 * vtemp1
    else
        r_star = vtemp1
    end if

    ! initialize rho(0)
    rho(0) = dot_product(r_star, r_star)

    ! compute norm of initial residual (r_0) to see if we really need
    ! to do anything
    r_init_norm = sqrt(rho(0))
    me%resnorm = r_init_norm
    if (r_init_norm <= delta) then
        flag = fail(SUN_SUCCESS)
        return
    end if

    ! set v = A*r_0 (preconditioned and scaled)
    if (scale_x) then
        vtemp1 = r_star / me%s2
    else
        vtemp1 = r_star
    end if

    if (preOnRight) then
        v = vtemp1
        status = prob%ls_psolve(v, vtemp1, delta, SUN_PREC_RIGHT)
        if (status /= 0) then
            flag = fail(psolve_fail(status))
            return
        end if
    end if

    status = prob%ls_atimes(vtemp1, v)
    if (status /= 0) then
        flag = fail(atimes_fail(status))
        return
    end if

    if (preOnLeft) then
        status = prob%ls_psolve(v, vtemp1, delta, SUN_PREC_LEFT)
        if (status /= 0) then
            flag = fail(psolve_fail(status))
            return
        end if
    else
        vtemp1 = v
    end if

    if (scale_b) then
        v = me%s1 * vtemp1
    else
        v = vtemp1
    end if

    ! initialize remaining variables
    r(:,0) = r_star
    u = r_star
    p = r_star
    d = ZERO

    ! set x = sx x if non-zero guess
    if (scale_x .and. .not. me%zeroguess) x = me%s2 * x

    tau = r_init_norm
    v_bar = ZERO
    eta = ZERO

    ! start outer loop
    do n = 0, l_max - 1

        ! increment linear iteration counter
        me%numiters = me%numiters + 1

        ! sigma = r_star^T*v
        sigma = dot_product(r_star, v)

        ! alpha = rho(0)/sigma
        alpha = rho(0) / sigma

        ! q = u-alpha*v
        q = u - alpha * v

        ! r(1) = r(0)-alpha*A*(u+q)
        r(:,1) = u + q
        if (scale_x) r(:,1) = r(:,1) / me%s2

        if (preOnRight) then
            vtemp1 = r(:,1)
            status = prob%ls_psolve(vtemp1, r(:,1), delta, SUN_PREC_RIGHT)
            if (status /= 0) then
                flag = fail(psolve_fail(status))
                return
            end if
        end if

        status = prob%ls_atimes(r(:,1), vtemp1)
        if (status /= 0) then
            flag = fail(atimes_fail(status))
            return
        end if

        if (preOnLeft) then
            status = prob%ls_psolve(vtemp1, r(:,1), delta, SUN_PREC_LEFT)
            if (status /= 0) then
                flag = fail(psolve_fail(status))
                return
            end if
        else
            r(:,1) = vtemp1
        end if

        if (scale_b) then
            vtemp1 = me%s1 * r(:,1)
        else
            vtemp1 = r(:,1)
        end if

        r(:,1) = r(:,0) - alpha * vtemp1

        ! start inner loop
        do m = 0, 1

            ! d = [*]+(v_bar^2*eta/alpha)*d, where [*] = u if m == 0, and q
            ! if m == 1 (temp_val saves computations if the inner loop is
            ! executed twice)
            if (m == 0) then
                temp_val = sqrt(dot_product(r(:,1), r(:,1)))
                omega = dot_product(r(:,0), r(:,0))
                omega = sqrt(sqrt(omega) * temp_val)
                d = u + (v_bar**2 * eta / alpha) * d
            else
                omega = temp_val
                d = q + (v_bar**2 * eta / alpha) * d
            end if

            ! v_bar = omega/tau
            v_bar = omega / tau

            ! c = (1+v_bar^2)^(-1/2)
            c = ONE / sqrt(ONE + v_bar**2)

            ! tau = tau*v_bar*c
            tau = tau * v_bar * c

            ! eta = c^2*alpha
            eta = c**2 * alpha

            ! x = x+eta*d
            if (n == 0 .and. m == 0 .and. me%zeroguess) then
                x = eta * d
            else
                x = x + eta * d
            end if

            ! check for convergence, using the approximation to the norm of
            ! the residual if possible
            r_curr_norm = tau * sqrt(real(m + 1, wp))
            me%resnorm = r_curr_norm

            ! exit inner loop if iteration has converged based upon
            ! approximation to norm of current residual
            if (r_curr_norm <= delta) then
                converged = .true.
                exit
            end if

            ! decide if actual norm of residual vector should be computed:
            ! (1) if r_curr_norm > delta, check if the actual residual norm
            !     is OK; (2) if r_curr_norm >= r_init_norm and m == 1 and
            !     n == l_max, compute the actual residual norm to see if the
            !     iteration can be saved; (3) the scaled and preconditioned
            !     b is computed only once and saved in vtemp3
            if ((r_curr_norm > delta) .or. &
                (r_curr_norm >= r_init_norm .and. m == 1 .and. n == l_max)) then

                ! compute norm of residual ||b-A*x||_2 (preconditioned and scaled)
                if (scale_x) then
                    vtemp1 = x / me%s2
                else
                    vtemp1 = x
                end if

                if (preOnRight) then
                    status = prob%ls_psolve(vtemp1, vtemp2, delta, SUN_PREC_RIGHT)
                    if (status /= 0) then
                        flag = fail(psolve_fail(status))
                        return
                    end if
                    vtemp1 = vtemp2
                end if

                status = prob%ls_atimes(vtemp1, vtemp2)
                if (status /= 0) then
                    flag = fail(atimes_fail(status))
                    return
                end if

                if (preOnLeft) then
                    status = prob%ls_psolve(vtemp2, vtemp1, delta, SUN_PREC_LEFT)
                    if (status /= 0) then
                        flag = fail(psolve_fail(status))
                        return
                    end if
                else
                    vtemp1 = vtemp2
                end if

                if (scale_b) then
                    vtemp2 = me%s1 * vtemp1
                else
                    vtemp2 = vtemp1
                end if

                ! only precondition and scale b once (result saved for reuse)
                if (.not. b_ok) then
                    b_ok = .true.
                    if (preOnLeft) then
                        status = prob%ls_psolve(b, vtemp3, delta, SUN_PREC_LEFT)
                        if (status /= 0) then
                            flag = fail(psolve_fail(status))
                            return
                        end if
                    else
                        vtemp3 = b
                    end if
                    if (scale_b) vtemp3 = me%s1 * vtemp3
                end if

                vtemp1 = vtemp3 - vtemp2
                r_curr_norm = sqrt(dot_product(vtemp1, vtemp1))
                me%resnorm = r_curr_norm

                ! exit inner loop if we have converged
                if (r_curr_norm <= delta) then
                    converged = .true.
                    exit
                end if

            end if

        end do  ! end inner loop

        ! if converged, then exit outer loop as well
        if (converged) exit

        ! rho(1) = r_star^T*r(1)
        rho(1) = dot_product(r_star, r(:,1))

        ! beta = rho(1)/rho(0)
        beta = rho(1) / rho(0)

        ! u = r(1)+beta*q
        u = r(:,1) + beta * q

        ! p = u+beta*(q+beta*p) = beta*beta*p + beta*q + u
        p = beta**2 * p + beta * q + u

        ! v = A*p
        if (scale_x) then
            vtemp1 = p / me%s2
        else
            vtemp1 = p
        end if

        if (preOnRight) then
            v = vtemp1
            status = prob%ls_psolve(v, vtemp1, delta, SUN_PREC_RIGHT)
            if (status /= 0) then
                flag = fail(psolve_fail(status))
                return
            end if
        end if

        status = prob%ls_atimes(vtemp1, v)
        if (status /= 0) then
            flag = fail(atimes_fail(status))
            return
        end if

        if (preOnLeft) then
            status = prob%ls_psolve(v, vtemp1, delta, SUN_PREC_LEFT)
            if (status /= 0) then
                flag = fail(psolve_fail(status))
                return
            end if
        else
            vtemp1 = v
        end if

        if (scale_b) then
            v = me%s1 * vtemp1
        else
            v = vtemp1
        end if

        ! shift variable values
        r(:,0) = r(:,1)
        rho(0) = rho(1)

    end do  ! end outer loop

    ! determine return value: if iteration converged or residual was
    ! reduced, then return current iterate (x)
    if (converged .or. (r_curr_norm < r_init_norm)) then

        if (scale_x) x = x / me%s2

        if (preOnRight) then
            status = prob%ls_psolve(x, vtemp1, delta, SUN_PREC_RIGHT)
            if (status /= 0) then
                flag = fail(psolve_fail(status))
                return
            end if
            x = vtemp1
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

        pure function atimes_fail(status) result(f)
        !! the flag of a failed A*v product
        integer,intent(in) :: status  !! the A*v return value
        integer :: f  !! `SUNLS_ATIMES_FAIL_UNREC` or `SUNLS_ATIMES_FAIL_REC`
        f = merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0)
        end function atimes_fail

        pure function psolve_fail(status) result(f)
        !! the flag of a failed preconditioner solve
        integer,intent(in) :: status  !! the psolve return value
        integer :: f  !! `SUNLS_PSOLVE_FAIL_UNREC` or `SUNLS_PSOLVE_FAIL_REC`
        f = merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0)
        end function psolve_fail

        function fail(f) result(ff)
        !! return: reset zeroguess and record the flag
        integer,intent(in) :: f  !! the return flag
        integer :: ff  !! the same flag
        me%zeroguess = .false.
        me%last_flag = f
        ff = f
        end function fail

    end function sptfqmr_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolNumIters_SPTFQMR`

    pure function sptfqmr_num_iters(me) result(n)

    class(sun_linsol_sptfqmr),intent(in) :: me  !! the solver
    integer :: n  !! number of iterations in the last solve

    n = me%numiters

    end function sptfqmr_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolResNorm_SPTFQMR`

    pure function sptfqmr_res_norm(me) result(r)

    class(sun_linsol_sptfqmr),intent(in) :: me  !! the solver
    real(wp) :: r  !! final scaled residual norm of the last solve

    r = me%resnorm

    end function sptfqmr_res_norm
!*****************************************************************************************

    end module fvode_linsol_sptfqmr
!*****************************************************************************************

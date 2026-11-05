!*****************************************************************************************
!>
!  Preconditioned conjugate gradient linear solver (SUNDIALS v7.9.0
!  `SUNLinSol_PCG`), translated from `src/sunlinsol/pcg/sunlinsol_pcg.c`.
!
!  For symmetric positive definite `A` (and preconditioner P): at most
!  `maxl` CG iterations, with the residual measured in the 2-norm scaled by
!  the first scaling vector `s1` (the second one is not used). Any
!  preconditioning type other than none means the preconditioner is
!  applied (as `SUN_PREC_LEFT` in the psolve call).

    module fvode_linsol_pcg

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUNLINEARSOLVER_ITERATIVE, SUN_SUCCESS, SUN_ERR_ARG_CORRUPT, &
                                SUN_ERR_ARG_OUTOFRANGE, &
                                SUN_PREC_NONE, SUN_PREC_LEFT, SUN_PREC_RIGHT, SUN_PREC_BOTH, &
                                SUNPCG_MAXL_DEFAULT, &
                                SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, &
                                SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, &
                                SUNLS_PSET_FAIL_UNREC, SUNLS_PSET_FAIL_REC, &
                                SUNLS_RES_REDUCED, SUNLS_CONV_FAIL
    use fvode_matrix, only: sun_matrix
    use fvode_linsol, only: sun_linear_solver, ls_problem

    implicit none

    private

    real(wp),parameter :: ZERO = 0.0_wp  !! 0

    type,extends(sun_linear_solver),public :: sun_linsol_pcg
        !! PCG solver (`SUNLinearSolverContent_PCG`).
        integer :: maxl = SUNPCG_MAXL_DEFAULT  !! maximum number of iterations
        integer :: pretype = SUN_PREC_NONE  !! preconditioning type (`SUN_PREC_*`)
        integer :: numiters = 0  !! number of iterations in the last solve
        real(wp) :: resnorm = ZERO  !! final (scaled) residual norm of the last solve
        real(wp),dimension(:),allocatable :: r  !! the residual
        real(wp),dimension(:),allocatable :: p  !! the search direction
        real(wp),dimension(:),allocatable :: z  !! the preconditioned residual
        real(wp),dimension(:),allocatable :: Ap  !! A*p (also the scaled residual)
    contains
        procedure :: get_type   => pcg_get_type  !! `SUNLinSolGetType_PCG`: iterative solver
        procedure :: initialize => pcg_initialize  !! `SUNLinSolInitialize_PCG`: check options
        procedure :: setup      => pcg_setup  !! `SUNLinSolSetup_PCG`: preconditioner setup
        procedure :: solve      => pcg_solve  !! `SUNLinSolSolve_PCG`: CG iteration
        procedure :: num_iters  => pcg_num_iters  !! `SUNLinSolNumIters_PCG`
        procedure :: res_norm   => pcg_res_norm  !! `SUNLinSolResNorm_PCG`
        procedure :: supports_atimes => pcg_supports  !! PCG uses an A*v function
        procedure :: supports_preconditioner => pcg_supports  !! PCG supports preconditioning
        procedure :: supports_scaling_vectors => pcg_supports  !! PCG supports (the first) scaling vector
        procedure :: set_prec_type => pcg_set_prec_type  !! `SUNLinSol_PCGSetPrecType`
        procedure :: set_maxl => pcg_set_maxl  !! `SUNLinSol_PCGSetMaxl`
    end type sun_linsol_pcg

    interface sun_linsol_pcg
        module procedure pcg_create
    end interface sun_linsol_pcg

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNLinSol_PCG(y, pretype, maxl)`): `n` is the vector
!  length. An illegal `pretype` means no preconditioning; `maxl <= 0`
!  means the default (5).

    function pcg_create(n, pretype, maxl) result(S)

    integer,intent(in) :: n  !! vector length
    integer,intent(in) :: pretype  !! preconditioning type
    integer,intent(in) :: maxl  !! maximum number of iterations
    type(sun_linsol_pcg) :: S  !! the solver

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        S%pretype = SUN_PREC_NONE
    else
        S%pretype = pretype
    end if
    S%maxl = maxl
    if (maxl <= 0) S%maxl = SUNPCG_MAXL_DEFAULT

    allocate(S%r(n), S%p(n), S%z(n), S%Ap(n))

    end function pcg_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver type: `SUNLINEARSOLVER_ITERATIVE`.

    pure function pcg_get_type(me) result(itype)

    class(sun_linsol_pcg),intent(in) :: me  !! the solver
    integer :: itype  !! `SUNLINEARSOLVER_ITERATIVE`

    itype = SUNLINEARSOLVER_ITERATIVE

    end function pcg_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  PCG provides the A*v, preconditioner and scaling-vector operations.

    pure function pcg_supports(me) result(supported)

    class(sun_linsol_pcg),intent(in) :: me  !! the solver
    logical :: supported  !! always true

    supported = .true.

    end function pcg_supports
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_PCGSetPrecType`

    function pcg_set_prec_type(me, pretype) result(flag)

    class(sun_linsol_pcg),intent(inout) :: me  !! the solver
    integer,intent(in) :: pretype  !! `SUN_PREC_NONE`, `_LEFT`, `_RIGHT` or `_BOTH`
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_OUTOFRANGE`

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        flag = SUN_ERR_ARG_OUTOFRANGE
        return
    end if
    me%pretype = pretype
    flag = SUN_SUCCESS

    end function pcg_set_prec_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_PCGSetMaxl` (<= 0 restores the default, 5).

    function pcg_set_maxl(me, maxl) result(flag)

    class(sun_linsol_pcg),intent(inout) :: me  !! the solver
    integer,intent(in) :: maxl  !! maximum number of iterations
    integer :: flag  !! `SUN_SUCCESS`

    me%maxl = maxl
    if (maxl <= 0) me%maxl = SUNPCG_MAXL_DEFAULT
    flag = SUN_SUCCESS

    end function pcg_set_maxl
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolInitialize_PCG`: check the options (no memory to allocate).

    function pcg_initialize(me) result(flag)

    class(sun_linsol_pcg),intent(inout) :: me  !! the solver
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_CORRUPT`

    if (me%maxl <= 0) me%maxl = SUNPCG_MAXL_DEFAULT
    if ((me%pretype /= SUN_PREC_LEFT) .and. (me%pretype /= SUN_PREC_RIGHT) .and. &
        (me%pretype /= SUN_PREC_BOTH)) me%pretype = SUN_PREC_NONE

    if (.not. me%has_atimes .or. &
        (me%pretype /= SUN_PREC_NONE .and. .not. me%has_psolve)) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if

    flag = SUN_SUCCESS

    end function pcg_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetup_PCG`: call the preconditioner setup function (if any).

    function pcg_setup(me, A, prob) result(flag)

    class(sun_linsol_pcg),intent(inout) :: me  !! the solver
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

    end function pcg_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSolve_PCG`: solve `A x = b` to the scaled residual tolerance
!  `delta`. Returns `SUN_SUCCESS` (converged), `SUNLS_RES_REDUCED` (not
!  converged, but the residual was reduced), `SUNLS_CONV_FAIL`, or a
!  failure flag of the A*v / preconditioner operations.

    function pcg_solve(me, A, x, b, tol, prob) result(flag)

    class(sun_linsol_pcg),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used (matrix-free)
    real(wp),dimension(:),intent(inout) :: x  !! initial guess (unless zeroguess) / solution
    real(wp),dimension(:),intent(in) :: b  !! right-hand side
    real(wp),intent(in) :: tol  !! tolerance delta on the scaled residual 2-norm
    class(ls_problem),intent(inout),optional :: prob  !! the problem (A*v, preconditioner)
    integer :: flag  !! return flag (see above)

    real(wp) :: alpha  !! step length along p
    real(wp) :: beta  !! update factor of the search direction
    real(wp) :: r0_norm  !! initial residual norm
    real(wp) :: rho  !! current residual norm
    real(wp) :: rz  !! <r,z>
    real(wp) :: rz_old  !! <r,z> of the previous iteration
    real(wp) :: delta  !! tolerance
    logical :: UsePrec  !! is a preconditioner applied?
    logical :: UseScaling  !! is the scaling vector in use?
    logical :: converged  !! has the iteration converged?
    integer :: l  !! iteration index
    integer :: l_max  !! maximum number of iterations
    integer :: status  !! return value of A*v / psolve

    delta = tol

    associate (r => me%r, p => me%p, z => me%z, Ap => me%Ap)

    ! initialize counters and convergence flag
    l_max = me%maxl
    me%numiters = 0
    converged = .false.

    ! set booleans for internal solver options
    UsePrec = (me%pretype == SUN_PREC_BOTH) .or. (me%pretype == SUN_PREC_LEFT) .or. &
              (me%pretype == SUN_PREC_RIGHT)
    UseScaling = allocated(me%s1)

    ! check if Atimes function has been set, and psolve if preconditioning
    if (.not. present(prob) .or. .not. me%has_atimes .or. &
        (UsePrec .and. .not. me%has_psolve)) then
        flag = fail(SUN_ERR_ARG_CORRUPT)
        return
    end if

    ! set r to initial residual r_0 = b - A*x_0
    if (me%zeroguess) then
        r = b
    else
        status = prob%ls_atimes(x, r)
        if (status /= 0) then
            flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
            return
        end if
        r = b - r
    end if

    ! set rho to scaled L2 norm of r, and return if small
    if (UseScaling) then
        Ap = r * me%s1
    else
        Ap = r
    end if
    rho = sqrt(dot_product(Ap, Ap))
    r0_norm = rho
    me%resnorm = rho
    if (rho <= delta) then
        flag = fail(SUN_SUCCESS)
        return
    end if

    ! apply preconditioner and b-scaling to r = r_0
    if (UsePrec) then
        status = prob%ls_psolve(r, z, delta, SUN_PREC_LEFT)  ! z = P^{-1}r
        if (status /= 0) then
            flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
            return
        end if
    else
        z = r
    end if

    ! initialize rz to <r,z>
    rz = dot_product(r, z)

    ! copy z to p
    p = z

    ! begin main iteration loop
    do l = 0, l_max - 1

        ! increment counter
        me%numiters = me%numiters + 1

        ! generate Ap = A*p
        status = prob%ls_atimes(p, Ap)
        if (status /= 0) then
            flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
            return
        end if

        ! calculate alpha = <r,z> / <Ap,p>
        alpha = dot_product(Ap, p)
        alpha = rz / alpha

        ! update x = x + alpha*p
        if (l == 0 .and. me%zeroguess) then
            x = alpha * p
        else
            x = x + alpha * p
        end if

        ! update r = r - alpha*Ap
        r = r - alpha * Ap

        ! set rho and check convergence
        if (UseScaling) then
            Ap = r * me%s1
        else
            Ap = r
        end if
        rho = sqrt(dot_product(Ap, Ap))
        me%resnorm = rho
        if (rho <= delta) then
            converged = .true.
            exit
        end if

        ! exit early on last iteration
        if (l == l_max - 1) exit

        ! apply preconditioner: z = P^{-1}*r
        if (UsePrec) then
            status = prob%ls_psolve(r, z, delta, SUN_PREC_LEFT)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
        else
            z = r
        end if

        ! update rz
        rz_old = rz
        rz = dot_product(r, z)

        ! calculate beta = <r,z> / <r_old,z_old>
        beta = rz / rz_old

        ! update p = z + beta*p
        p = z + beta * p

    end do

    ! main loop finished, return with result
    if (converged) then
        flag = fail(SUN_SUCCESS)
    else if (rho < r0_norm) then
        flag = fail(SUNLS_RES_REDUCED)
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

    end function pcg_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolNumIters_PCG`

    pure function pcg_num_iters(me) result(n)

    class(sun_linsol_pcg),intent(in) :: me  !! the solver
    integer :: n  !! number of iterations in the last solve

    n = me%numiters

    end function pcg_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolResNorm_PCG`

    pure function pcg_res_norm(me) result(r)

    class(sun_linsol_pcg),intent(in) :: me  !! the solver
    real(wp) :: r  !! final scaled residual norm of the last solve

    r = me%resnorm

    end function pcg_res_norm
!*****************************************************************************************

    end module fvode_linsol_pcg
!*****************************************************************************************

!*****************************************************************************************
!>
!  Scaled, preconditioned GMRES linear solver (SUNDIALS v7.9.0
!  `SUNLinSol_SPGMR`), translated from
!  `src/sunlinsol/spgmr/sunlinsol_spgmr.c`.
!
!  Solves `A x = b` with GMRES(maxl), restarted up to `max_restarts` times,
!  with optional left and/or right preconditioning and scaling: the
!  method is applied to `(S1 P1^{-1} A P2^{-1} S2^{-1}) (S2 P2 x) = S1 P1^{-1} b`
!  and the residual is measured in the scaled 2-norm.

    module fvode_linsol_spgmr

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUNLINEARSOLVER_ITERATIVE, SUN_SUCCESS, SUN_ERR_ARG_CORRUPT, &
                                SUN_PREC_NONE, SUN_PREC_LEFT, SUN_PREC_RIGHT, SUN_PREC_BOTH, &
                                SUN_MODIFIED_GS, SUN_CLASSICAL_GS, SUNSPGMR_MAXL_DEFAULT, &
                                SUNSPGMR_MAXRS_DEFAULT, SUNSPGMR_GSTYPE_DEFAULT, &
                                SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, &
                                SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, &
                                SUNLS_PSET_FAIL_UNREC, SUNLS_PSET_FAIL_REC, &
                                SUNLS_QRFACT_FAIL, SUNLS_QRSOL_FAIL, SUNLS_RES_REDUCED, &
                                SUNLS_CONV_FAIL
    use fvode_matrix, only: sun_matrix
    use fvode_linsol, only: sun_linear_solver, ls_problem
    use fvode_iterative, only: modified_gs, classical_gs, qr_fact, qr_sol

    implicit none

    private

    real(wp),parameter :: ZERO = 0.0_wp  !! 0
    real(wp),parameter :: ONE  = 1.0_wp  !! 1

    type,extends(sun_linear_solver),public :: sun_linsol_spgmr
        !! GMRES solver (`SUNLinearSolverContent_SPGMR`).
        integer :: maxl = SUNSPGMR_MAXL_DEFAULT  !! maximum Krylov subspace dimension
        integer :: pretype = SUN_PREC_NONE  !! preconditioning type (`SUN_PREC_*`)
        integer :: gstype = SUNSPGMR_GSTYPE_DEFAULT  !! Gram-Schmidt type (`SUN_MODIFIED_GS` / `SUN_CLASSICAL_GS`)
        integer :: max_restarts = SUNSPGMR_MAXRS_DEFAULT  !! maximum number of restarts
        integer :: numiters = 0  !! number of iterations in the last solve
        real(wp) :: resnorm = ZERO  !! final (scaled) residual norm of the last solve
        real(wp),dimension(:,:),allocatable :: V  !! Krylov basis vectors V(:,0:maxl)
        real(wp),dimension(:,:),allocatable :: Hes  !! Hessenberg matrix Hes(0:maxl,0:maxl-1)
        real(wp),dimension(:),allocatable :: givens  !! Givens rotations, 0:2*maxl-1
        real(wp),dimension(:),allocatable :: xcor  !! correction to the initial guess
        real(wp),dimension(:),allocatable :: yg  !! least-squares solution, 0:maxl
        real(wp),dimension(:),allocatable :: vtemp  !! temporary vector
        real(wp),dimension(:),allocatable :: cv  !! scratch for classical Gram-Schmidt, 0:maxl
    contains
        procedure :: get_type   => spgmr_get_type  !! `SUNLinSolGetType_SPGMR`: iterative solver
        procedure :: initialize => spgmr_initialize  !! `SUNLinSolInitialize_SPGMR`: check options, allocate workspace
        procedure :: setup      => spgmr_setup  !! `SUNLinSolSetup_SPGMR`: preconditioner setup
        procedure :: solve      => spgmr_solve  !! `SUNLinSolSolve_SPGMR`: GMRES iteration
        procedure :: num_iters  => spgmr_num_iters  !! `SUNLinSolNumIters_SPGMR`
        procedure :: res_norm   => spgmr_res_norm  !! `SUNLinSolResNorm_SPGMR`
        procedure :: supports_atimes => spgmr_supports  !! SPGMR uses an A*v function
        procedure :: supports_preconditioner => spgmr_supports  !! SPGMR supports preconditioning
        procedure :: supports_scaling_vectors => spgmr_supports  !! SPGMR supports scaling vectors
        procedure :: set_prec_type => spgmr_set_prec_type  !! `SUNLinSol_SPGMRSetPrecType`
        procedure :: set_gs_type => spgmr_set_gs_type  !! `SUNLinSol_SPGMRSetGSType`
        procedure :: set_max_restarts => spgmr_set_max_restarts  !! `SUNLinSol_SPGMRSetMaxRestarts`
    end type sun_linsol_spgmr

    interface sun_linsol_spgmr
        module procedure spgmr_create
    end interface sun_linsol_spgmr

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor (`SUNLinSol_SPGMR(y, pretype, maxl)`): `n` is the vector
!  length. An illegal `pretype` means no preconditioning; `maxl <= 0`
!  means the default (5).

    function spgmr_create(n, pretype, maxl) result(S)

    integer,intent(in) :: n  !! vector length
    integer,intent(in) :: pretype  !! preconditioning type
    integer,intent(in) :: maxl  !! maximum Krylov subspace dimension
    type(sun_linsol_spgmr) :: S  !! the solver

    if ((pretype /= SUN_PREC_NONE) .and. (pretype /= SUN_PREC_LEFT) .and. &
        (pretype /= SUN_PREC_RIGHT) .and. (pretype /= SUN_PREC_BOTH)) then
        S%pretype = SUN_PREC_NONE
    else
        S%pretype = pretype
    end if
    S%maxl = maxl
    if (maxl <= 0) S%maxl = SUNSPGMR_MAXL_DEFAULT

    allocate(S%xcor(n), S%vtemp(n))

    end function spgmr_create
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver type: `SUNLINEARSOLVER_ITERATIVE`.

    pure function spgmr_get_type(me) result(itype)

    class(sun_linsol_spgmr),intent(in) :: me  !! the solver
    integer :: itype  !! `SUNLINEARSOLVER_ITERATIVE`

    itype = SUNLINEARSOLVER_ITERATIVE

    end function spgmr_get_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  SPGMR provides the A*v, preconditioner and scaling-vector operations.

    pure function spgmr_supports(me) result(supported)

    class(sun_linsol_spgmr),intent(in) :: me  !! the solver
    logical :: supported  !! always true

    supported = .true.

    end function spgmr_supports
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPGMRSetPrecType`

    function spgmr_set_prec_type(me, pretype) result(flag)

    class(sun_linsol_spgmr),intent(inout) :: me  !! the solver
    integer,intent(in) :: pretype  !! `SUN_PREC_NONE`, `_LEFT`, `_RIGHT` or `_BOTH`
    integer :: flag  !! `SUN_SUCCESS`

    me%pretype = pretype
    flag = SUN_SUCCESS

    end function spgmr_set_prec_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPGMRSetGSType`

    function spgmr_set_gs_type(me, gstype) result(flag)

    class(sun_linsol_spgmr),intent(inout) :: me  !! the solver
    integer,intent(in) :: gstype  !! `SUN_MODIFIED_GS` or `SUN_CLASSICAL_GS`
    integer :: flag  !! `SUN_SUCCESS`

    me%gstype = gstype
    flag = SUN_SUCCESS

    end function spgmr_set_gs_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSol_SPGMRSetMaxRestarts` (< 0 restores the default, 0).

    function spgmr_set_max_restarts(me, maxrs) result(flag)

    class(sun_linsol_spgmr),intent(inout) :: me  !! the solver
    integer,intent(in) :: maxrs  !! maximum number of restarts
    integer :: flag  !! `SUN_SUCCESS`

    me%max_restarts = maxrs
    if (maxrs < 0) me%max_restarts = SUNSPGMR_MAXRS_DEFAULT
    flag = SUN_SUCCESS

    end function spgmr_set_max_restarts
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolInitialize_SPGMR`: check the options and allocate the
!  Krylov workspace.

    function spgmr_initialize(me) result(flag)

    class(sun_linsol_spgmr),intent(inout) :: me  !! the solver
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_CORRUPT`

    integer :: n  !! vector length

    if (me%max_restarts < 0) me%max_restarts = SUNSPGMR_MAXRS_DEFAULT
    if ((me%pretype /= SUN_PREC_LEFT) .and. (me%pretype /= SUN_PREC_RIGHT) .and. &
        (me%pretype /= SUN_PREC_BOTH)) me%pretype = SUN_PREC_NONE

    if (.not. me%has_atimes .or. &
        (me%pretype /= SUN_PREC_NONE .and. .not. me%has_psolve)) then
        flag = SUN_ERR_ARG_CORRUPT
        return
    end if

    n = size(me%vtemp)
    if (.not. allocated(me%V))      allocate(me%V(n, 0:me%maxl))
    if (.not. allocated(me%Hes))    allocate(me%Hes(0:me%maxl, 0:me%maxl-1))
    if (.not. allocated(me%givens)) allocate(me%givens(0:2*me%maxl-1))
    if (.not. allocated(me%yg))     allocate(me%yg(0:me%maxl))
    if (.not. allocated(me%cv))     allocate(me%cv(0:me%maxl))

    me%last_flag = SUN_SUCCESS
    flag = SUN_SUCCESS

    end function spgmr_initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetup_SPGMR`: call the preconditioner setup function (if any).

    function spgmr_setup(me, A, prob) result(flag)

    class(sun_linsol_spgmr),intent(inout) :: me  !! the solver
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

    end function spgmr_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSolve_SPGMR`: solve `A x = b` to the scaled residual tolerance
!  `delta`. Returns `SUN_SUCCESS` (converged), `SUNLS_RES_REDUCED` (not
!  converged, but the residual was reduced), `SUNLS_CONV_FAIL`, or a
!  failure flag of the A*v / preconditioner / QR operations.

    function spgmr_solve(me, A, x, b, tol, prob) result(flag)

    class(sun_linsol_spgmr),intent(inout) :: me  !! the solver
    class(sun_matrix),intent(inout),optional :: A  !! not used (matrix-free)
    real(wp),dimension(:),intent(inout) :: x  !! initial guess (unless zeroguess) / solution
    real(wp),dimension(:),intent(in) :: b  !! right-hand side
    real(wp),intent(in) :: tol  !! tolerance delta on the scaled residual 2-norm
    class(ls_problem),intent(inout),optional :: prob  !! the problem (A*v, preconditioner)
    integer :: flag  !! return flag (see above)

    real(wp) :: beta  !! initial residual norm
    real(wp) :: rotation_product  !! product of the Givens sines
    real(wp) :: r_norm  !! residual norm at the start of a cycle
    real(wp) :: s_product  !! product of Givens sines (restart)
    real(wp) :: rho  !! current residual norm estimate
    real(wp) :: delta  !! tolerance
    logical :: preOnLeft, preOnRight  !! left / right preconditioning?
    logical :: scale1, scale2  !! scaling vectors in use?
    logical :: converged  !! has the iteration converged?
    integer :: i, k, l  !! indices
    integer :: l_plus_1  !! l + 1
    integer :: l_max  !! maximum Krylov dimension
    integer :: krydim  !! current Krylov dimension
    integer :: ntries  !! restart counter
    integer :: status  !! return value of A*v / psolve
    real(wp) :: hnorm  !! norm of the new Krylov vector (Hes(l+1,l))

    delta = tol

    associate (V => me%V, Hes => me%Hes, givens => me%givens, xcor => me%xcor, &
               yg => me%yg, vtemp => me%vtemp)

    if (.not. present(prob)) then
        flag = fail(SUN_ERR_ARG_CORRUPT)
        return
    end if

    ! initialize some variables
    l_plus_1 = 0
    krydim = 0
    l_max = me%maxl
    me%numiters = 0
    converged = .false.

    ! set booleans for preconditioning and scaling
    preOnLeft  = (me%pretype == SUN_PREC_LEFT) .or. (me%pretype == SUN_PREC_BOTH)
    preOnRight = (me%pretype == SUN_PREC_RIGHT) .or. (me%pretype == SUN_PREC_BOTH)
    scale1 = allocated(me%s1)
    scale2 = allocated(me%s2)

    ! set vtemp and V[0] to initial (unscaled) residual r_0 = b - A*x_0
    if (me%zeroguess) then
        vtemp = b
    else
        status = prob%ls_atimes(x, vtemp)
        if (status /= 0) then
            flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
            return
        end if
        vtemp = b - vtemp
    end if
    V(:,0) = vtemp

    ! apply left preconditioner and left scaling to V[0] = r_0
    if (preOnLeft) then
        status = prob%ls_psolve(V(:,0), vtemp, delta, SUN_PREC_LEFT)
        if (status /= 0) then
            flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
            return
        end if
    else
        vtemp = V(:,0)
    end if

    if (scale1) then
        V(:,0) = me%s1 * vtemp
    else
        V(:,0) = vtemp
    end if

    ! set r_norm = beta to L2 norm of V[0] = s1 P1_inv r_0, and
    ! return if small
    r_norm = sqrt(dot_product(V(:,0), V(:,0)))
    beta = r_norm
    me%resnorm = r_norm
    if (r_norm <= delta) then
        flag = done(SUN_SUCCESS)
        return
    end if

    ! initialize rho to avoid compiler warning message
    rho = beta

    ! set xcor = 0
    xcor = ZERO

    ! begin outer iterations: up to (max_restarts + 1) attempts
    do ntries = 0, me%max_restarts

        ! initialize the Hessenberg matrix Hes and Givens rotation
        ! product. Normalize the initial vector V[0]
        Hes = ZERO
        rotation_product = ONE
        V(:,0) = (ONE / r_norm) * V(:,0)

        ! inner loop: generate Krylov sequence and Arnoldi basis
        do l = 0, l_max - 1

            me%numiters = me%numiters + 1

            krydim = l + 1
            l_plus_1 = l + 1

            ! generate A-tilde V[l], where A-tilde = s1 P1_inv A P2_inv s2_inv

            ! apply right scaling: vtemp = s2_inv V[l]
            if (scale2) then
                vtemp = V(:,l) / me%s2
            else
                vtemp = V(:,l)
            end if

            ! apply right preconditioner: vtemp = P2_inv s2_inv V[l]
            if (preOnRight) then
                V(:,l_plus_1) = vtemp
                status = prob%ls_psolve(V(:,l_plus_1), vtemp, delta, SUN_PREC_RIGHT)
                if (status /= 0) then
                    flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                    return
                end if
            end if

            ! apply A: V[l+1] = A P2_inv s2_inv V[l]
            status = prob%ls_atimes(vtemp, V(:,l_plus_1))
            if (status /= 0) then
                flag = fail(merge(SUNLS_ATIMES_FAIL_UNREC, SUNLS_ATIMES_FAIL_REC, status < 0))
                return
            end if

            ! apply left preconditioning: vtemp = P1_inv A P2_inv s2_inv V[l]
            if (preOnLeft) then
                status = prob%ls_psolve(V(:,l_plus_1), vtemp, delta, SUN_PREC_LEFT)
                if (status /= 0) then
                    flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                    return
                end if
            else
                vtemp = V(:,l_plus_1)
            end if

            ! apply left scaling: V[l+1] = s1 P1_inv A P2_inv s2_inv V[l]
            if (scale1) then
                V(:,l_plus_1) = me%s1 * vtemp
            else
                V(:,l_plus_1) = vtemp
            end if

            ! orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde
            if (me%gstype == SUN_CLASSICAL_GS) then
                call classical_gs(V, Hes, l_plus_1, l_max, hnorm, me%cv)
            else
                call modified_gs(V, Hes, l_plus_1, l_max, hnorm)
            end if
            Hes(l_plus_1,l) = hnorm

            ! update the QR factorization of Hes
            if (qr_fact(krydim, Hes, givens, l) /= 0) then
                flag = fail(SUNLS_QRFACT_FAIL)
                return
            end if

            ! update residual norm estimate; break if convergence test passes
            rotation_product = rotation_product * givens(2*l+1)
            rho = abs(rotation_product * r_norm)
            me%resnorm = rho
            if (rho <= delta) then
                converged = .true.
                exit
            end if

            ! normalize V[l+1] with norm value from the Gram-Schmidt routine
            V(:,l_plus_1) = (ONE / Hes(l_plus_1,l)) * V(:,l_plus_1)

        end do

        ! inner loop is done. Compute the new correction vector xcor

        ! construct g, then solve for y
        yg(0) = r_norm
        yg(1:krydim) = ZERO
        if (qr_sol(krydim, Hes, givens, yg) /= 0) then
            flag = fail(SUNLS_QRSOL_FAIL)
            return
        end if

        ! add correction vector V_l y to xcor
        do k = 0, krydim - 1
            xcor = xcor + yg(k) * V(:,k)
        end do

        ! if converged, construct the final solution vector x and return
        if (converged) then
            call finish_solution(status)
            if (status /= 0) then
                flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
                return
            end if
            flag = done(SUN_SUCCESS)
            return
        end if

        ! not yet converged; if allowed, prepare for restart
        if (ntries == me%max_restarts) exit

        ! construct last column of Q in yg
        s_product = ONE
        do i = krydim, 1, -1
            yg(i) = s_product * givens(2*i-2)
            s_product = s_product * givens(2*i-1)
        end do
        yg(0) = s_product

        ! scale r_norm and yg
        r_norm = r_norm * s_product
        yg(0:krydim) = r_norm * yg(0:krydim)
        r_norm = abs(r_norm)

        ! multiply yg by V_(krydim+1) to get last residual vector; restart
        V(:,0) = yg(0) * V(:,0)
        do k = 1, krydim
            V(:,0) = V(:,0) + yg(k) * V(:,k)
        end do

    end do

    ! failed to converge, even after allowed restarts.
    ! If the residual norm was reduced below its initial value, compute
    ! and return x anyway. Otherwise return failure flag.
    if (rho < beta) then
        call finish_solution(status)
        if (status /= 0) then
            flag = fail(merge(SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_FAIL_REC, status < 0))
            return
        end if
        flag = done(SUNLS_RES_REDUCED)
        return
    end if

    flag = done(SUNLS_CONV_FAIL)

    end associate

    contains

        subroutine finish_solution(status)
        !! x = x + P2_inv s2_inv xcor (or x = P2_inv s2_inv xcor for a zero guess)
        integer,intent(out) :: status  !! return value of psolve (0 if not called)
        status = 0
        ! apply right scaling and right preconditioning: vtemp = P2_inv s2_inv xcor
        if (scale2) me%xcor = me%xcor / me%s2
        if (preOnRight) then
            status = prob%ls_psolve(me%xcor, me%vtemp, delta, SUN_PREC_RIGHT)
            if (status /= 0) return
        else
            me%vtemp = me%xcor
        end if
        ! add vtemp to initial x to get final solution x, and return
        if (me%zeroguess) then
            x = me%vtemp
        else
            x = x + me%vtemp
        end if
        end subroutine finish_solution

        function fail(f) result(ff)
        !! a failure: reset zeroguess and record the flag
        integer,intent(in) :: f  !! the failure flag
        integer :: ff  !! the same flag
        me%zeroguess = .false.
        me%last_flag = f
        ff = f
        end function fail

        function done(f) result(ff)
        !! normal return: reset zeroguess and record the flag
        integer,intent(in) :: f  !! the return flag
        integer :: ff  !! the same flag
        me%zeroguess = .false.
        me%last_flag = f
        ff = f
        end function done

    end function spgmr_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolNumIters_SPGMR`

    pure function spgmr_num_iters(me) result(n)

    class(sun_linsol_spgmr),intent(in) :: me  !! the solver
    integer :: n  !! number of iterations in the last solve

    n = me%numiters

    end function spgmr_num_iters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolResNorm_SPGMR`

    pure function spgmr_res_norm(me) result(r)

    class(sun_linsol_spgmr),intent(in) :: me  !! the solver
    real(wp) :: r  !! final scaled residual norm of the last solve

    r = me%resnorm

    end function spgmr_res_norm
!*****************************************************************************************

    end module fvode_linsol_spgmr
!*****************************************************************************************

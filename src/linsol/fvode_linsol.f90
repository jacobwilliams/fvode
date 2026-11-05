!*****************************************************************************************
!>
!  Generic linear solver interface (SUNDIALS v7.9.0 `SUNLinearSolver`).
!
!  Translated from `include/sundials/sundials_linearsolver.h` and
!  `src/sundials/sundials_linearsolver.c`.
!
!  Iterative (Krylov) solvers need the matrix-vector product `A*v` and,
!  optionally, a preconditioner. In C these are function pointers plus
!  `void*` data given to `SUNLinSolSetATimes` / `SUNLinSolSetPreconditioner`.
!  Here the integrator implements the abstract type [[ls_problem]], which is
!  passed to `setup` and `solve`; the flags `has_atimes`, `has_psetup` and
!  `has_psolve` record which of the functions are available (C: the
!  pointers are not NULL).

    module fvode_linsol

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUN_SUCCESS, SUN_ERR_NOT_IMPLEMENTED
    use fvode_matrix, only: sun_matrix

    implicit none

    private

    type,abstract,public :: ls_problem
        !! The problem side of an iterative linear solve: the integrator
        !! implements this (C: `SUNATimesFn`, `SUNPSetupFn`, `SUNPSolveFn`).
    contains
        procedure(ls_atimes_fn),deferred :: ls_atimes  !! `SUNATimesFn`: z = A*v
        procedure(ls_psetup_fn),deferred :: ls_psetup  !! `SUNPSetupFn`: set up the preconditioner
        procedure(ls_psolve_fn),deferred :: ls_psolve  !! `SUNPSolveFn`: solve P z = r
    end type ls_problem

    type,abstract,public :: sun_linear_solver
        !! Base type for linear solvers (`SUNLinearSolver`).
        integer :: last_flag = 0          !! last error flag (`SUNLinSolLastFlag`)
        logical :: has_atimes = .false.   !! an A*v function is attached (iterative solvers)
        logical :: has_psetup = .false.   !! a preconditioner setup function is attached
        logical :: has_psolve = .false.   !! a preconditioner solve function is attached
        logical :: zeroguess = .false.    !! the initial guess is zero (`SUNLinSolSetZeroGuess`)
        real(wp),dimension(:),allocatable :: s1  !! first scaling vector (allocated: in use)
        real(wp),dimension(:),allocatable :: s2  !! second scaling vector (allocated: in use)
    contains
        procedure(ls_gettype_fn),deferred :: get_type   !! `SUNLinSolGetType`
        procedure(ls_setup_fn),deferred   :: setup      !! `SUNLinSolSetup`
        procedure(ls_solve_fn),deferred   :: solve      !! `SUNLinSolSolve`
        procedure :: initialize => ls_initialize_default !! `SUNLinSolInitialize`
        procedure :: supports_atimes => ls_supports_none  !! C: `ops->setatimes != NULL`
        procedure :: supports_preconditioner => ls_supports_none  !! C: `ops->setpreconditioner != NULL`
        procedure :: supports_scaling_vectors => ls_supports_none  !! C: `ops->setscalingvectors != NULL`
        procedure :: set_scaling_vectors => ls_set_scaling_vectors  !! `SUNLinSolSetScalingVectors`
        procedure :: set_zero_guess => ls_set_zero_guess  !! `SUNLinSolSetZeroGuess`
        procedure :: num_iters => ls_num_iters_default  !! `SUNLinSolNumIters` (0 for direct solvers)
        procedure :: res_norm => ls_res_norm_default    !! `SUNLinSolResNorm` (0 for direct solvers)
    end type sun_linear_solver

    abstract interface

        function ls_atimes_fn(me, v, z) result(flag)
            !! Compute the product z = A*v.
            import :: ls_problem, wp
            implicit none
            class(ls_problem),intent(inout) :: me  !! the problem (integrator)
            real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
            real(wp),dimension(:),intent(out) :: z  !! the product A*v
            integer :: flag  !! 0 on success, > 0 recoverable failure, < 0 unrecoverable
        end function ls_atimes_fn

        function ls_psetup_fn(me) result(flag)
            !! Set up the preconditioner.
            import :: ls_problem
            implicit none
            class(ls_problem),intent(inout) :: me  !! the problem (integrator)
            integer :: flag  !! 0 on success, > 0 recoverable failure, < 0 unrecoverable
        end function ls_psetup_fn

        function ls_psolve_fn(me, r, z, tol, lr) result(flag)
            !! Solve the preconditioner system P z = r.
            import :: ls_problem, wp
            implicit none
            class(ls_problem),intent(inout) :: me  !! the problem (integrator)
            real(wp),dimension(:),intent(in) :: r  !! the right-hand side
            real(wp),dimension(:),intent(out) :: z  !! the solution
            real(wp),intent(in) :: tol  !! tolerance for an iterative preconditioner solve
            integer,intent(in) :: lr  !! `SUN_PREC_LEFT` or `SUN_PREC_RIGHT`
            integer :: flag  !! 0 on success, > 0 recoverable failure, < 0 unrecoverable
        end function ls_psolve_fn

        pure function ls_gettype_fn(me) result(itype)
            !! Linear solver type (`SUNLinSolGetType`).
            import :: sun_linear_solver
            implicit none
            class(sun_linear_solver),intent(in) :: me  !! the linear solver
            integer :: itype   !! `SUNLINEARSOLVER_DIRECT`, `_ITERATIVE`, ...
        end function ls_gettype_fn

        function ls_setup_fn(me, A, prob) result(flag)
            !! Prepare to solve with A (e.g. LU factorization, or preconditioner setup).
            import :: sun_linear_solver, sun_matrix, ls_problem
            implicit none
            class(sun_linear_solver),intent(inout) :: me  !! the linear solver
            class(sun_matrix),intent(inout),optional :: A  !! the matrix A (direct solvers; may be overwritten by its factors)
            class(ls_problem),intent(inout),optional :: prob  !! the problem (iterative solvers)
            integer :: flag  !! `SUN_SUCCESS`, a recoverable failure (> 0) or an unrecoverable one (< 0)
        end function ls_setup_fn

        function ls_solve_fn(me, A, x, b, tol, prob) result(flag)
            !! Solve A x = b (b is not modified).
            import :: sun_linear_solver, sun_matrix, ls_problem, wp
            implicit none
            class(sun_linear_solver),intent(inout) :: me  !! the linear solver
            class(sun_matrix),intent(inout),optional :: A  !! the matrix A, as prepared by `setup` (direct solvers)
            real(wp),dimension(:),intent(inout) :: x  !! on input an initial guess (ignored if `zeroguess`), on output the solution
            real(wp),dimension(:),intent(in) :: b  !! the right-hand side b
            real(wp),intent(in) :: tol  !! tolerance for iterative solvers (ignored by direct solvers)
            class(ls_problem),intent(inout),optional :: prob  !! the problem (iterative solvers)
            integer :: flag  !! `SUN_SUCCESS`, a recoverable failure (> 0) or an unrecoverable one (< 0)
        end function ls_solve_fn

    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Default `initialize` (C: optional op; success if absent).

    function ls_initialize_default(me) result(flag)

    class(sun_linear_solver),intent(inout) :: me  !! the linear solver
    integer :: flag  !! always `SUN_SUCCESS`

    me%last_flag = SUN_SUCCESS
    flag = SUN_SUCCESS

    end function ls_initialize_default
!*****************************************************************************************

!*****************************************************************************************
!>
!  Default for the `supports_*` queries: the op is not provided.

    pure function ls_supports_none(me) result(supported)

    class(sun_linear_solver),intent(in) :: me  !! the linear solver
    logical :: supported  !! always false

    supported = .false.

    end function ls_supports_none
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetScalingVectors`: set the scaling vectors (copies are kept;
!  an absent vector means no scaling). Fails for solvers without scaling.

    function ls_set_scaling_vectors(me, s1, s2) result(flag)

    class(sun_linear_solver),intent(inout) :: me  !! the linear solver
    real(wp),dimension(:),intent(in),optional :: s1  !! first scaling vector
    real(wp),dimension(:),intent(in),optional :: s2  !! second scaling vector
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_NOT_IMPLEMENTED`

    if (.not. me%supports_scaling_vectors()) then
        flag = SUN_ERR_NOT_IMPLEMENTED
        return
    end if

    if (present(s1)) then
        me%s1 = s1
    else if (allocated(me%s1)) then
        deallocate(me%s1)
    end if
    if (present(s2)) then
        me%s2 = s2
    else if (allocated(me%s2)) then
        deallocate(me%s2)
    end if
    flag = SUN_SUCCESS

    end function ls_set_scaling_vectors
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetZeroGuess`: the next solve starts from x = 0.

    function ls_set_zero_guess(me, onoff) result(flag)

    class(sun_linear_solver),intent(inout) :: me  !! the linear solver
    logical,intent(in) :: onoff  !! is the initial guess zero?
    integer :: flag  !! `SUN_SUCCESS`

    me%zeroguess = onoff
    flag = SUN_SUCCESS

    end function ls_set_zero_guess
!*****************************************************************************************

!*****************************************************************************************
!>
!  Default `SUNLinSolNumIters`: 0 (direct solvers).

    pure function ls_num_iters_default(me) result(n)

    class(sun_linear_solver),intent(in) :: me  !! the linear solver
    integer :: n  !! number of iterations in the last solve

    n = 0

    end function ls_num_iters_default
!*****************************************************************************************

!*****************************************************************************************
!>
!  Default `SUNLinSolResNorm`: 0 (direct solvers).

    pure function ls_res_norm_default(me) result(r)

    class(sun_linear_solver),intent(in) :: me  !! the linear solver
    real(wp) :: r  !! final residual norm of the last solve

    r = 0.0_wp

    end function ls_res_norm_default
!*****************************************************************************************

    end module fvode_linsol
!*****************************************************************************************

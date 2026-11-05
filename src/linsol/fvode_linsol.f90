!*****************************************************************************************
!>
!  Generic linear solver interface (SUNDIALS v7.9.0 `SUNLinearSolver`).
!
!  Translated from `include/sundials/sundials_linearsolver.h` and
!  `src/sundials/sundials_linearsolver.c`. Only the operations used by the
!  direct solvers are here so far; the iterative ones (`setatimes`,
!  `setpreconditioner`, `setscalingvectors`, `numiters`, `resnorm`) come with
!  the Krylov solvers.

    module fvode_linsol

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUN_SUCCESS
    use fvode_matrix, only: sun_matrix

    implicit none

    private

    type,abstract,public :: sun_linear_solver
        !! Base type for linear solvers (`SUNLinearSolver`).
        integer :: last_flag = 0   !! last error flag (`SUNLinSolLastFlag`)
    contains
        procedure(ls_gettype_fn),deferred :: get_type   !! `SUNLinSolGetType`
        procedure(ls_setup_fn),deferred   :: setup      !! `SUNLinSolSetup`
        procedure(ls_solve_fn),deferred   :: solve      !! `SUNLinSolSolve`
        procedure :: initialize => ls_initialize_default !! `SUNLinSolInitialize`
    end type sun_linear_solver

    abstract interface

        pure function ls_gettype_fn(me) result(itype)
            !! Linear solver type (`SUNLinSolGetType`).
            import :: sun_linear_solver
            implicit none
            class(sun_linear_solver),intent(in) :: me  !! the linear solver
            integer :: itype   !! `SUNLINEARSOLVER_DIRECT`, `_ITERATIVE`, ...
        end function ls_gettype_fn

        function ls_setup_fn(me, A) result(flag)
            !! Prepare to solve with A (e.g. LU factorization).
            import :: sun_linear_solver, sun_matrix
            implicit none
            class(sun_linear_solver),intent(inout) :: me  !! the linear solver
            class(sun_matrix),intent(inout) :: A  !! the matrix A (overwritten by its factors)
            integer :: flag  !! `SUN_SUCCESS`, > 0: recoverable failure, < 0: unrecoverable
        end function ls_setup_fn

        function ls_solve_fn(me, A, x, b, tol) result(flag)
            !! Solve A x = b (b is not modified).
            import :: sun_linear_solver, sun_matrix, wp
            implicit none
            class(sun_linear_solver),intent(inout) :: me  !! the linear solver
            class(sun_matrix),intent(inout) :: A  !! the matrix A (as prepared by `setup`)
            real(wp),dimension(:),intent(inout) :: x  !! in: initial guess (iterative only); out: x
            real(wp),dimension(:),intent(in) :: b  !! the right-hand side b
            real(wp),intent(in) :: tol  !! tolerance (iterative solvers only)
            integer :: flag  !! `SUN_SUCCESS`, > 0: recoverable failure, < 0: unrecoverable
        end function ls_solve_fn

    end interface

    contains
!*****************************************************************************************

    function ls_initialize_default(me) result(flag)
    !! Default `initialize` (C: optional op; success if absent).
    class(sun_linear_solver),intent(inout) :: me  !! the linear solver
    integer :: flag  !! always `SUN_SUCCESS`
    me%last_flag = SUN_SUCCESS
    flag = SUN_SUCCESS
    end function ls_initialize_default

    end module fvode_linsol
!*****************************************************************************************

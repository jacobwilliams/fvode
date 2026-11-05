!*****************************************************************************************
!>
!  Generic nonlinear solver interface (SUNDIALS v7.9.0 `SUNNonlinearSolver`).
!
!  In C, an integrator attaches function pointers plus a `void*` memory block to
!  the nonlinear solver (`SUNNonlinSolSetSysFn`, `SetLSetupFn`, `SetLSolveFn`,
!  `SetConvTestFn`). Here the integrator instead implements the abstract type
!  [[nls_problem]], and the solver calls its deferred bindings. Whether the
!  linear-solver hooks exist (C: function pointer is NULL or not) is recorded in
!  the solver's `has_lsetup` / `has_lsolve` flags, set by the integrator
!  (`cvNlsInit`).
!
!  Translated from `include/sundials/sundials_nonlinearsolver.h` and
!  `src/sundials/sundials_nonlinearsolver.c`.

    module fvode_nonlinsol

    use fvode_kinds, only: wp, li

    implicit none

    private

    type,abstract,public :: nls_problem
        !! The problem side of a nonlinear solve: the integrator implements this.
    contains
        procedure(nls_sys_fn),deferred      :: nls_sys       !! `SUNNonlinSolSysFn`
        procedure(nls_lsetup_fn),deferred   :: nls_lsetup    !! `SUNNonlinSolLSetupFn`
        procedure(nls_lsolve_fn),deferred   :: nls_lsolve    !! `SUNNonlinSolLSolveFn`
        procedure(nls_convtest_fn),deferred :: nls_convtest  !! `SUNNonlinSolConvTestFn`
    end type nls_problem

    type,abstract,public :: sun_nonlinear_solver
        !! Base type for nonlinear solvers (`SUNNonlinearSolver`).
        logical :: has_lsetup = .false.  !! an lsetup function is attached
        logical :: has_lsolve = .false.  !! an lsolve function is attached
    contains
        procedure(nls_gettype_fn),deferred    :: get_type        !! `SUNNonlinSolGetType`
        procedure(nls_initialize_fn),deferred :: initialize      !! `SUNNonlinSolInitialize`
        procedure(nls_solve_fn),deferred      :: solve           !! `SUNNonlinSolSolve`
        procedure(nls_setmaxiters_fn),deferred :: set_max_iters  !! `SUNNonlinSolSetMaxIters`
        procedure(nls_getlong_fn),deferred    :: get_num_iters   !! `SUNNonlinSolGetNumIters`
        procedure(nls_getint_fn),deferred     :: get_cur_iter    !! `SUNNonlinSolGetCurIter`
        procedure(nls_getlong_fn),deferred    :: get_num_conv_fails !! `SUNNonlinSolGetNumConvFails`
        procedure :: setup => nls_setup_default                  !! `SUNNonlinSolSetup` (optional op)
    end type sun_nonlinear_solver

    abstract interface

        function nls_sys_fn(me, ycor, res) result(flag)
            !! Evaluate the nonlinear system at the correction `ycor`.
            import :: nls_problem, wp
            implicit none
            class(nls_problem),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: ycor
            real(wp),dimension(:),intent(inout) :: res
            integer :: flag
        end function nls_sys_fn

        function nls_lsetup_fn(me, jbad, jcur) result(flag)
            !! Set up the linear solver (Newton only).
            import :: nls_problem
            implicit none
            class(nls_problem),intent(inout) :: me
            logical,intent(in) :: jbad
            logical,intent(out) :: jcur
            integer :: flag
        end function nls_lsetup_fn

        function nls_lsolve_fn(me, b) result(flag)
            !! Solve the linear system, overwriting `b` (Newton only).
            import :: nls_problem, wp
            implicit none
            class(nls_problem),intent(inout) :: me
            real(wp),dimension(:),intent(inout) :: b
            integer :: flag
        end function nls_lsolve_fn

        function nls_convtest_fn(me, nls, ycor, del, tol, ewt) result(flag)
            !! Convergence test: returns 0 (converged), `SUN_NLS_CONTINUE`,
            !! `SUN_NLS_CONV_RECVR`, or a negative failure code.
            import :: nls_problem, sun_nonlinear_solver, wp
            implicit none
            class(nls_problem),intent(inout) :: me
            class(sun_nonlinear_solver),intent(in) :: nls
            real(wp),dimension(:),intent(in) :: ycor
            real(wp),dimension(:),intent(in) :: del
            real(wp),intent(in) :: tol
            real(wp),dimension(:),intent(in) :: ewt
            integer :: flag
        end function nls_convtest_fn

        function nls_gettype_fn(me) result(itype)
            import :: sun_nonlinear_solver
            implicit none
            class(sun_nonlinear_solver),intent(in) :: me
            integer :: itype   !! `SUNNONLINEARSOLVER_ROOTFIND` / `_FIXEDPOINT` / `_HYBRID`
        end function nls_gettype_fn

        function nls_initialize_fn(me) result(flag)
            !! Initialize the solver (reset counters). Work space is allocated
            !! by the constructor, as in C.
            import :: sun_nonlinear_solver
            implicit none
            class(sun_nonlinear_solver),intent(inout) :: me
            integer :: flag
        end function nls_initialize_fn

        function nls_solve_fn(me, prob, y0, ycor, w, tol, call_lsetup) result(flag)
            !! Solve the nonlinear system, starting from (and returning) the
            !! correction `ycor` relative to the predictor `y0`.
            import :: sun_nonlinear_solver, nls_problem, wp
            implicit none
            class(sun_nonlinear_solver),intent(inout) :: me
            class(nls_problem),intent(inout) :: prob
            real(wp),dimension(:),intent(in) :: y0
            real(wp),dimension(:),intent(inout) :: ycor
            real(wp),dimension(:),intent(in) :: w
            real(wp),intent(in) :: tol
            logical,intent(in) :: call_lsetup
            integer :: flag
        end function nls_solve_fn

        function nls_setmaxiters_fn(me, maxiters) result(flag)
            import :: sun_nonlinear_solver
            implicit none
            class(sun_nonlinear_solver),intent(inout) :: me
            integer,intent(in) :: maxiters
            integer :: flag
        end function nls_setmaxiters_fn

        function nls_getlong_fn(me) result(val)
            import :: sun_nonlinear_solver, li
            implicit none
            class(sun_nonlinear_solver),intent(in) :: me
            integer(li) :: val
        end function nls_getlong_fn

        function nls_getint_fn(me) result(val)
            import :: sun_nonlinear_solver
            implicit none
            class(sun_nonlinear_solver),intent(in) :: me
            integer :: val
        end function nls_getint_fn

    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Default `setup` (C: `ops->setup == NULL`, so `SUNNonlinSolSetup` returns
!  success without doing anything).

    function nls_setup_default(me, prob, y) result(flag)

    class(sun_nonlinear_solver),intent(inout) :: me
    class(nls_problem),intent(inout) :: prob
    real(wp),dimension(:),intent(in) :: y
    integer :: flag

    flag = 0

    end function nls_setup_default
!*****************************************************************************************

    end module fvode_nonlinsol
!*****************************************************************************************

!*****************************************************************************************
!>
!  The CVODE integrator (SUNDIALS v7.9.0): main solver type and time stepping.
!
!  Translated from `src/cvode/cvode.c`, `cvode_impl.h` and `cvode_nls.c`. Each
!  procedure names the C function it was translated from. Statement order and
!  floating-point expressions follow the C source, so that the step-by-step
!  behavior (step sizes, orders, counters) is identical.
!
!### Naming
!
!  `cv_mem->cv_X` becomes `me%X`. Exceptions, because Fortran is not
!  case-sensitive:
!
!  * `cv_L` (= q+1) is `me%LL` (it would clash with the coefficient array `cv_l`)
!  * the local `M[3]` in `cvSetAdams` is `mm` (clashes with `m`)
!
!  Arrays keep the C index ranges: `zn(:,0:qmax)`, `l(0:L_MAX-1)`,
!  `tau(0:L_MAX)` (`tau(1..q+1)` used), `tq(1:NUM_TESTS)`.
!
!### Not yet ported
!
!  Linear solvers (and so the default Newton iteration), root finding,
!  inequality constraints, projection, stability limit detection, history
!  resizing and Anderson acceleration. The places where the C code handles
!  them are marked `! [not yet ported]`.

    module fvode_core

    use iso_fortran_env, only: error_unit
    use fvode_kinds, only: wp, ip, li
    use fvode_parameters
    use fvode_internal_parameters
    use fvode_nvector, only: n_vwrmsnorm
    use fvode_nonlinsol, only: sun_nonlinear_solver, nls_problem

    implicit none

    private

    type,extends(nls_problem),public :: cvode_solver
        !! CVODE solver (`CVodeMemRec`).
        !!
        !! Pass problem data to the callbacks by extending this type.
        private

        real(wp) :: uround = SUN_UNIT_ROUNDOFF  !! machine unit roundoff

        !---- problem specification
        procedure(cv_rhs_fn),pointer,nopass :: f => null()  !! y' = f(t,y)
        integer :: lmm = CV_BDF        !! CV_ADAMS or CV_BDF
        integer :: itol = CV_NN        !! CV_SS, CV_SV, CV_WF, CV_NN
        real(wp) :: reltol = ZERO      !! relative tolerance
        real(wp) :: Sabstol = ZERO     !! scalar absolute tolerance
        real(wp),dimension(:),allocatable :: Vabstol  !! vector absolute tolerance
        logical :: atolmin0 = .true.   !! min(abstol) == 0
        logical :: user_efun = .false. !! user-supplied ewt function
        procedure(cv_ewt_fn),pointer,nopass :: efun_user => null()  !! user ewt function

        !---- Nordsieck history array
        real(wp),dimension(:,:),allocatable :: zn  !! zn(:,j) = h^j/j! * y^(j), j=0..qmax

        !---- vectors of length N
        integer(ip) :: n = 0                          !! problem size
        real(wp),dimension(:),allocatable :: ewt      !! error weight vector
        real(wp),dimension(:),allocatable :: y        !! work vector (C: user's yout)
        real(wp),dimension(:),allocatable :: acor     !! accumulated correction / local error
        real(wp),dimension(:),allocatable :: tempv    !! temporary storage
        real(wp),dimension(:),allocatable :: ftemp    !! temporary storage
        real(wp),dimension(:),allocatable :: vtemp1   !! temporary storage
        real(wp),dimension(:),allocatable :: vtemp2   !! temporary storage
        real(wp),dimension(:),allocatable :: vtemp3   !! temporary storage

        !---- tstop information
        logical :: tstopset = .false.
        logical :: tstopinterp = .false.
        real(wp) :: tstop = ZERO

        !---- step data
        integer :: q = 0          !! current order
        integer :: qprime = 0     !! order to be used on the next step
        integer :: next_q = 0     !! order to be used on the next step
        integer :: qwait = 0      !! steps to wait before considering a change in q
        integer :: LL = 0         !! L = q + 1  (C: cv_L)
        real(wp) :: hin = ZERO    !! initial step size
        real(wp) :: h = ZERO      !! current step size
        real(wp) :: hprime = ZERO !! step size to be used on the next step
        real(wp) :: next_h = ZERO !! step size to be used on the next step
        real(wp) :: eta = ZERO    !! eta = hprime / h
        real(wp) :: hscale = ZERO !! value of h used in zn
        real(wp) :: tn = ZERO     !! current internal value of t
        real(wp) :: tretlast = ZERO !! last value of t returned by solve
        real(wp),dimension(0:L_MAX) :: tau = ZERO     !! previous step sizes, tau(1..q+1)
        real(wp),dimension(NUM_TESTS) :: tq = ZERO    !! test quantities, tq(1..5)
        real(wp),dimension(0:L_MAX-1) :: l = ZERO     !! coefficients of l(x)
        real(wp) :: rl1 = ZERO     !! 1/l(1)
        real(wp) :: gamma = ZERO   !! gamma = h * rl1
        real(wp) :: gammap = ZERO  !! gamma at the last setup call
        real(wp) :: gamrat = ZERO  !! gamma / gammap
        real(wp) :: crate = ZERO   !! estimated corrector convergence rate
        real(wp) :: delp = ZERO    !! norm of previous nonlinear solver update
        real(wp) :: delnrm = ZERO  !! norm of current nonlinear solver update
        real(wp) :: acnrm = ZERO   !! |acor|
        logical :: acnrmcur = .false. !! is |acor| current?
        real(wp) :: nlscoef = CORTES  !! coefficient in nonlinear convergence test

        !---- limits
        integer :: qmax = 0                    !! q <= qmax
        integer(li) :: mxstep = MXSTEP_DEFAULT !! max internal steps per call
        integer :: mxhnil = MXHNIL_DEFAULT     !! max "t + h == t" warnings
        integer :: maxnef = MXNEF              !! max error test failures
        integer :: maxncf = MXNCF              !! max nonlinear convergence failures
        real(wp) :: hmin = HMIN_DEFAULT        !! |h| >= hmin
        real(wp) :: hmax_inv = HMAX_INV_DEFAULT !! |h| <= 1/hmax_inv
        real(wp) :: etamax = ZERO              !! eta <= etamax
        real(wp) :: eta_min_fx = ETA_MIN_FX_DEFAULT !! eta_min_fx < eta < eta_max_fx: keep h
        real(wp) :: eta_max_fx = ETA_MAX_FX_DEFAULT
        real(wp) :: eta_max_fs = ETA_MAX_FS_DEFAULT !! eta <= eta_max_fs on the first step
        real(wp) :: eta_max_es = ETA_MAX_ES_DEFAULT !! eta <= eta_max_es on early steps
        real(wp) :: eta_max_gs = ETA_MAX_GS_DEFAULT !! eta <= eta_max_gs on a general step
        real(wp) :: eta_min = ETA_MIN_DEFAULT       !! eta >= eta_min on a general step
        real(wp) :: eta_min_ef = ETA_MIN_EF_DEFAULT !! eta >= eta_min_ef after an error test failure
        real(wp) :: eta_max_ef = ETA_MAX_EF_DEFAULT !! eta on multiple error test failures
        real(wp) :: eta_cf = ETA_CF_DEFAULT         !! eta on a nonlinear convergence failure
        integer(li) :: small_nst = SMALL_NST_DEFAULT !! nst <= small_nst: use eta_max_es
        integer :: small_nef = SMALL_NEF_DEFAULT     !! nef >= small_nef: use eta_max_ef

        !---- counters
        integer(li) :: nst = 0     !! internal steps taken
        integer(li) :: nfe = 0     !! f evaluations
        integer(li) :: ncfn = 0    !! corrector convergence failures
        integer(li) :: nni = 0     !! nonlinear iterations
        integer(li) :: nnf = 0     !! nonlinear convergence failures
        integer(li) :: netf = 0    !! error test failures
        integer(li) :: nsetups = 0 !! linear solver setup calls
        integer :: nhnil = 0       !! "t + h == t" warnings issued

        !---- step size ratios
        real(wp) :: etaqm1 = ZERO  !! ratio of new to old h for order q-1
        real(wp) :: etaq = ZERO    !! ratio of new to old h for order q
        real(wp) :: etaqp1 = ZERO  !! ratio of new to old h for order q+1

        !---- nonlinear solver
        class(sun_nonlinear_solver),allocatable :: NLS  !! nonlinear solver
        integer :: convfail = CV_NO_FAILURES            !! for the linear solver setup

        !---- linear solver  [not yet ported]
        logical :: lsetup_exists = .false.  !! C: cv_lsetup != NULL
        integer(li) :: msbp = MSBP_DEFAULT  !! max steps between lsetup calls
        real(wp) :: dgmax_lsetup = DGMAX_LSETUP_DEFAULT !! gamma ratio threshold for lsetup

        !---- saved values
        integer :: qu = 0            !! last successful q
        integer(li) :: nstlp = 0     !! step number of last setup call
        real(wp) :: h0u = ZERO       !! actual initial step size
        real(wp) :: hu = ZERO        !! last successful h
        real(wp) :: saved_tq5 = ZERO !! saved value of tq(5)
        logical :: jcur = .false.    !! is Jacobian info current?
        real(wp) :: tolsf = ONE      !! tolerance scale factor
        integer :: qmax_alloc = 0    !! qmax used when allocating zn
        integer :: indx_acor = 0     !! index of the zn vector with saved acor

        !---- flags
        logical :: MallocDone = .false.  !! init has been called

        !---- stability limit detection  [not yet ported]
        logical :: sldeton = .false.
        integer :: nscon = 0         !! counter for STALD method
        integer(li) :: nor = 0       !! number of order reductions
        real(wp),dimension(6,4) :: ssdat = ZERO

        !---- root finding  [not yet ported]
        integer :: nrtfn = 0
        integer :: irfnd = 0
        integer(li) :: nge = 0

        !---- error reporting (replaces the SUNContext error handler)
        logical :: print_errors = .true.     !! print error and warning messages
        integer :: errunit = error_unit      !! unit for messages

    contains

        ! user-callable
        procedure,public :: create              => CVodeCreate
        procedure,public :: init                => CVodeInit
        procedure,public :: reinit              => CVodeReInit
        procedure,public :: ss_tolerances       => CVodeSStolerances
        procedure,public :: sv_tolerances       => CVodeSVtolerances
        procedure,public :: wf_tolerances       => CVodeWFtolerances
        procedure,public :: set_nonlinear_solver => CVodeSetNonlinearSolver
        procedure,public :: solve               => CVode
        procedure,public :: get_dky             => CVodeGetDky
        procedure,public :: destroy             => CVodeFree

        ! optional inputs (cvode_io.c)
        procedure,public :: set_error_output    => CVodeSetErrorOutput
        procedure,public :: set_max_ord         => CVodeSetMaxOrd
        procedure,public :: set_max_num_steps   => CVodeSetMaxNumSteps
        procedure,public :: set_max_hnil_warns  => CVodeSetMaxHnilWarns
        procedure,public :: set_init_step       => CVodeSetInitStep
        procedure,public :: set_min_step        => CVodeSetMinStep
        procedure,public :: set_max_step        => CVodeSetMaxStep
        procedure,public :: set_stop_time       => CVodeSetStopTime
        procedure,public :: set_interpolate_stop_time => CVodeSetInterpolateStopTime
        procedure,public :: clear_stop_time     => CVodeClearStopTime
        procedure,public :: set_max_err_test_fails => CVodeSetMaxErrTestFails
        procedure,public :: set_max_conv_fails  => CVodeSetMaxConvFails
        procedure,public :: set_nonlin_conv_coef => CVodeSetNonlinConvCoef

        ! optional outputs (cvode_io.c)
        procedure,public :: get_num_steps       => CVodeGetNumSteps
        procedure,public :: get_num_rhs_evals   => CVodeGetNumRhsEvals
        procedure,public :: get_num_lin_solv_setups => CVodeGetNumLinSolvSetups
        procedure,public :: get_num_err_test_fails => CVodeGetNumErrTestFails
        procedure,public :: get_last_order      => CVodeGetLastOrder
        procedure,public :: get_current_order   => CVodeGetCurrentOrder
        procedure,public :: get_actual_init_step => CVodeGetActualInitStep
        procedure,public :: get_last_step       => CVodeGetLastStep
        procedure,public :: get_current_step    => CVodeGetCurrentStep
        procedure,public :: get_current_time    => CVodeGetCurrentTime
        procedure,public :: get_tol_scale_factor => CVodeGetTolScaleFactor
        procedure,public :: get_num_nonlin_solv_iters => CVodeGetNumNonlinSolvIters
        procedure,public :: get_num_nonlin_solv_conv_fails => CVodeGetNumNonlinSolvConvFails
        procedure,public :: get_num_step_solve_fails => CVodeGetNumStepSolveFails

        ! nonlinear solver interface (cvode_nls.c)
        procedure :: nls_sys      => cvNlsSys
        procedure :: nls_lsetup   => cvNlsLSetup
        procedure :: nls_lsolve   => cvNlsLSolve
        procedure :: nls_convtest => cvNlsConvTest

    end type cvode_solver

    abstract interface

        function cv_rhs_fn(me, t, y, ydot) result(retval)
            !! Right-hand side function `ydot = f(t, y)` (`CVRhsFn`).
            !! Return 0 on success, > 0 for a recoverable error, < 0 for an
            !! unrecoverable error.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the solver (extend it to carry data)
            real(wp),intent(in) :: t
            real(wp),dimension(:),intent(in) :: y
            real(wp),dimension(:),intent(out) :: ydot
            integer :: retval
        end function cv_rhs_fn

        function cv_ewt_fn(me, y, ewt) result(retval)
            !! User error-weight function (`CVEwtFn`). Return 0 on success,
            !! nonzero on failure.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: y
            real(wp),dimension(:),intent(out) :: ewt
            integer :: retval
        end function cv_ewt_fn

    end interface

    public :: cv_rhs_fn, cv_ewt_fn

    contains
!*****************************************************************************************

!=========================================================================================
!  User-callable functions (cvode.c)
!=========================================================================================

!*****************************************************************************************
!>
!  `CVodeCreate`: set the linear multistep method and the default values of
!  all optional inputs.

    subroutine CVodeCreate(me, lmm, flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: lmm   !! `CV_ADAMS` or `CV_BDF`
    integer,intent(out) :: flag

    integer :: maxord

    if ((lmm /= CV_ADAMS) .and. (lmm /= CV_BDF)) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeCreate', &
            'Illegal value for lmm. The legal values are CV_ADAMS and CV_BDF.')
        flag = CV_ILL_INPUT
        return
    end if

    call CVodeFree(me)     ! start from a clean object (C: memset to 0)

    maxord = merge(ADAMS_Q_MAX, BDF_Q_MAX, lmm == CV_ADAMS)

    ! copy input parameters into cv_mem
    me%lmm = lmm

    ! set uround
    me%uround = SUN_UNIT_ROUNDOFF

    ! set default values for integrator optional inputs
    me%f            => null()
    me%itol         = CV_NN
    me%atolmin0     = .true.
    me%user_efun    = .false.
    me%efun_user    => null()
    me%qmax         = maxord
    me%mxstep       = MXSTEP_DEFAULT
    me%mxhnil       = MXHNIL_DEFAULT
    me%sldeton      = .false.
    me%hin          = ZERO
    me%hmin         = HMIN_DEFAULT
    me%hmax_inv     = HMAX_INV_DEFAULT
    me%eta_min_fx   = ETA_MIN_FX_DEFAULT
    me%eta_max_fx   = ETA_MAX_FX_DEFAULT
    me%eta_max_fs   = ETA_MAX_FS_DEFAULT
    me%eta_max_es   = ETA_MAX_ES_DEFAULT
    me%eta_max_gs   = ETA_MAX_GS_DEFAULT
    me%eta_min      = ETA_MIN_DEFAULT
    me%eta_min_ef   = ETA_MIN_EF_DEFAULT
    me%eta_max_ef   = ETA_MAX_EF_DEFAULT
    me%eta_cf       = ETA_CF_DEFAULT
    me%small_nst    = SMALL_NST_DEFAULT
    me%small_nef    = SMALL_NEF_DEFAULT
    me%tstopset     = .false.
    me%tstopinterp  = .false.
    me%maxnef       = MXNEF
    me%maxncf       = MXNCF
    me%nlscoef      = CORTES
    me%msbp         = MSBP_DEFAULT
    me%dgmax_lsetup = DGMAX_LSETUP_DEFAULT
    me%convfail     = CV_NO_FAILURES

    ! [not yet ported] constraints, root finding, projection, resize

    ! set the saved value qmax_alloc
    me%qmax_alloc = maxord

    ! no mallocs have been done yet
    me%MallocDone = .false.

    ! the nonlinear solver is set by init (C) / set_nonlinear_solver
    me%lsetup_exists = .false.

    flag = CV_SUCCESS

    end subroutine CVodeCreate
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeInit`: allocate and initialize memory for the problem
!  `y' = f(t, y)`, `y(t0) = y0`.
!
!  @note In C, `CVodeInit` also attaches a default Newton solver. Newton is
!  not ported yet: attach a solver with [[CVodeSetNonlinearSolver]] after
!  calling this.

    subroutine CVodeInit(me, f, t0, y0, flag)

    class(cvode_solver),intent(inout) :: me
    procedure(cv_rhs_fn) :: f
    real(wp),intent(in) :: t0
    real(wp),dimension(:),intent(in) :: y0
    integer,intent(out) :: flag

    integer :: i, k

    if (me%qmax == 0) then
        call cvProcessError(me, CV_MEM_NULL, 'CVodeInit', 'cvode_mem = NULL illegal.')
        flag = CV_MEM_NULL   ! create has not been called
        return
    end if

    ! allocate the vectors (using y0 as a template)
    call cvAllocVectors(me, size(y0, kind=ip))

    ! copy the input parameters into CVODE state
    me%f  => f
    me%tn = t0

    ! initialize zn[0] in the history array
    me%zn(:,0) = y0

    ! [not yet ported] create and attach the default Newton solver

    ! initialize all the counters and other optional output values
    me%q      = 1
    me%LL     = 2
    me%qwait  = me%LL
    me%etamax = me%eta_max_fs

    me%qu    = 0
    me%hu    = ZERO
    me%tolsf = ONE

    ! set the linear solver addresses to NULL
    me%lsetup_exists = .false.

    ! initialize all the counters
    me%nst     = 0
    me%nfe     = 0
    me%ncfn    = 0
    me%netf    = 0
    me%nni     = 0
    me%nnf     = 0
    me%nsetups = 0
    me%nhnil   = 0
    me%nstlp   = 0
    me%nscon   = 0
    me%nge     = 0

    me%irfnd = 0

    ! initialize other integrator optional outputs
    me%h0u    = ZERO
    me%next_h = ZERO
    me%next_q = 0

    ! initialize stability limit detection data
    me%nor = 0
    do i = 1, 5
        do k = 1, 3
            me%ssdat(i,k) = ZERO
        end do
    end do

    ! problem has been successfully initialized
    me%MallocDone = .true.

    flag = CV_SUCCESS

    end subroutine CVodeInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeReInit`: re-initialize for a new problem of the same size.

    subroutine CVodeReInit(me, t0, y0, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t0
    real(wp),dimension(:),intent(in) :: y0
    integer,intent(out) :: flag

    integer :: i, k

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeReInit', 'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    if (size(y0, kind=ip) /= me%n) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeReInit', 'y0 has the wrong size.')
        flag = CV_ILL_INPUT
        return
    end if

    ! copy the input parameters into CVODE state
    me%tn = t0

    ! set step parameters
    me%q      = 1
    me%LL     = 2
    me%qwait  = me%LL
    me%etamax = me%eta_max_fs

    me%qu    = 0
    me%hu    = ZERO
    me%tolsf = ONE

    ! initialize zn[0] in the history array
    me%zn(:,0) = y0

    ! initialize all the counters
    me%nst     = 0
    me%nfe     = 0
    me%ncfn    = 0
    me%netf    = 0
    me%nni     = 0
    me%nnf     = 0
    me%nsetups = 0
    me%nhnil   = 0
    me%nstlp   = 0
    me%nscon   = 0
    me%nge     = 0

    me%irfnd = 0

    ! [not yet ported] constraint counters, linear solver reinit

    ! initialize other integrator optional outputs
    me%h0u    = ZERO
    me%next_h = ZERO
    me%next_q = 0

    ! initialize stability limit detection data
    me%nor = 0
    do i = 1, 5
        do k = 1, 3
            me%ssdat(i,k) = ZERO
        end do
    end do

    flag = CV_SUCCESS

    end subroutine CVodeReInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSStolerances`: scalar relative and absolute tolerances.

    subroutine CVodeSStolerances(me, reltol, abstol, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: reltol
    real(wp),intent(in) :: abstol
    integer,intent(out) :: flag

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeSStolerances', 'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    ! check inputs
    if (reltol < ZERO) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSStolerances', 'reltol < 0 illegal.')
        flag = CV_ILL_INPUT
        return
    end if
    if (abstol < ZERO) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSStolerances', &
                            'abstol has negative component(s) (illegal).')
        flag = CV_ILL_INPUT
        return
    end if

    ! copy tolerances into memory
    me%reltol   = reltol
    me%Sabstol  = abstol
    me%atolmin0 = (abstol == ZERO)

    me%itol = CV_SS

    me%user_efun = .false.
    me%efun_user => null()

    flag = CV_SUCCESS

    end subroutine CVodeSStolerances
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSVtolerances`: scalar relative and vector absolute tolerances.

    subroutine CVodeSVtolerances(me, reltol, abstol, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: reltol
    real(wp),dimension(:),intent(in) :: abstol
    integer,intent(out) :: flag

    real(wp) :: atolmin

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeSVtolerances', 'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    ! check inputs
    if (reltol < ZERO) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSVtolerances', 'reltol < 0 illegal.')
        flag = CV_ILL_INPUT
        return
    end if
    if (size(abstol, kind=ip) /= me%n) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSVtolerances', 'abstol has the wrong size.')
        flag = CV_ILL_INPUT
        return
    end if

    atolmin = minval(abstol)
    if (atolmin < ZERO) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSVtolerances', &
                            'abstol has negative component(s) (illegal).')
        flag = CV_ILL_INPUT
        return
    end if

    ! copy tolerances into memory
    if (.not. allocated(me%Vabstol)) allocate(me%Vabstol(me%n))

    me%reltol = reltol
    me%Vabstol = abstol
    me%atolmin0 = (atolmin == ZERO)

    me%itol = CV_SV

    me%user_efun = .false.
    me%efun_user => null()

    flag = CV_SUCCESS

    end subroutine CVodeSVtolerances
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeWFtolerances`: user-supplied error weight function.

    subroutine CVodeWFtolerances(me, efun, flag)

    class(cvode_solver),intent(inout) :: me
    procedure(cv_ewt_fn) :: efun
    integer,intent(out) :: flag

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeWFtolerances', 'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    me%itol = CV_WF

    me%user_efun = .true.
    me%efun_user => efun

    flag = CV_SUCCESS

    end subroutine CVodeWFtolerances
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVode`: integrate towards `tout` (`itask = CV_NORMAL`) or take one
!  internal step (`itask = CV_ONE_STEP`). On return `yout` holds the solution
!  at `tret`.
!
!  The return flag is `CV_SUCCESS`, `CV_TSTOP_RETURN`, `CV_ROOT_RETURN`, or
!  a negative error code.

    subroutine CVode(me, tout, yout, tret, itask, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: tout
    real(wp),dimension(:),intent(inout) :: yout
    real(wp),intent(out) :: tret
    integer,intent(in) :: itask
    integer,intent(out) :: flag

    integer(li) :: nstloc
    integer :: retval, hflag, kflag, istate, ier
    integer :: ewtsetOK
    real(wp) :: troundoff, tout_hin, rh, nrm

    tret = me%tretlast

    !------------------------------------------------------------------
    ! check for legal inputs in all cases
    !------------------------------------------------------------------

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVode', 'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    if (size(yout, kind=ip) /= me%n) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVode', 'yout has the wrong size.')
        flag = CV_ILL_INPUT
        return
    end if

    if ((itask /= CV_NORMAL) .and. (itask /= CV_ONE_STEP)) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVode', 'Illegal value for itask.')
        flag = CV_ILL_INPUT
        return
    end if

    !------------------------------------------------------------------
    ! initializations and checks to be done on the first call only
    !------------------------------------------------------------------

    if (me%nst == 0) then

        me%tretlast = me%tn
        tret = me%tretlast

        ! check inputs for corectness
        ier = cvInitialSetup(me, tout)
        if (ier /= CV_SUCCESS) then
            flag = ier
            return
        end if

        ! call f at (t0,y0), set zn[1] = y'(t0)
        retval = me%f(me, me%tn, me%zn(:,0), me%zn(:,1))
        me%nfe = me%nfe + 1
        if (retval < 0) then
            call cvProcessError(me, CV_RHSFUNC_FAIL, 'CVode', &
                'At '//t_str(me%tn)//', the right-hand side routine failed in an unrecoverable manner.')
            flag = CV_RHSFUNC_FAIL
            return
        end if
        if (retval > 0) then
            call cvProcessError(me, CV_FIRST_RHSFUNC_ERR, 'CVode', &
                'The right-hand side routine failed at the first call.')
            flag = CV_FIRST_RHSFUNC_ERR
            return
        end if

        ! test input tstop for legality
        if (me%tstopset) then
            if ((me%tstop - me%tn) * (tout - me%tn) <= ZERO) then
                call cvProcessError(me, CV_ILL_INPUT, 'CVode', &
                    'The value tstop = '//g_str(me%tstop)//' is behind current '// &
                    t_str(me%tn)//' in the direction of integration.')
                flag = CV_ILL_INPUT
                return
            end if
        end if

        ! set initial h (from H0 or cvHin)
        me%h = me%hin
        if ((me%h /= ZERO) .and. ((tout - me%tn) * me%h < ZERO)) then
            call cvProcessError(me, CV_ILL_INPUT, 'CVode', 'h0 and tout - t0 inconsistent.')
            flag = CV_ILL_INPUT
            return
        end if
        if (me%h == ZERO) then
            tout_hin = tout
            if (me%tstopset .and. (tout - me%tn) * (tout - me%tstop) > ZERO) then
                tout_hin = me%tstop
            end if
            hflag = cvHin(me, tout_hin)
            if (hflag /= CV_SUCCESS) then
                istate = cvHandleFailure(me, hflag)
                flag = istate
                return
            end if
        end if

        ! enforce hmax and hmin
        rh = abs(me%h) * me%hmax_inv
        if (rh > ONE) me%h = me%h / rh
        if (abs(me%h) < me%hmin) me%h = me%h * (me%hmin / abs(me%h))

        ! check for approach to tstop
        if (me%tstopset) then
            if ((me%tn + me%h - me%tstop) * me%h > ZERO) then
                me%h = (me%tstop - me%tn) * (ONE - FOUR * me%uround)
            end if
        end if

        ! scale zn[1] by h
        me%hscale = me%h
        me%h0u    = me%h
        me%hprime = me%h

        me%zn(:,1) = me%h*me%zn(:,1)

        ! [not yet ported] root finding: cvRcheck1

    end if ! end of first call block

    !------------------------------------------------------------------
    ! at following steps, perform stop tests
    !------------------------------------------------------------------

    if (me%nst > 0) then

        ! estimate an infinitesimal time interval to be used as
        ! a roundoff for time quantities (based on current time
        ! and step size)
        troundoff = FUZZ_FACTOR * me%uround * (abs(me%tn) + abs(me%h))

        ! [not yet ported] root finding: cvRcheck2 / cvRcheck3

        ! in CV_NORMAL mode, test if tout was reached
        ! (tstop checks come first)
        if (me%tstopset) then

            ! test for tn at tstop
            if (abs(me%tn - me%tstop) <= troundoff) then
                ! ensure tout >= tstop, otherwise check for tout return below
                if ((tout - me%tstop) * me%h >= ZERO .or. &
                    abs(tout - me%tstop) <= troundoff) then
                    if (me%tstopinterp) then
                        call CVodeGetDky(me, me%tstop, 0, yout, ier)
                        if (ier /= CV_SUCCESS) then
                            call cvProcessError(me, CV_ILL_INPUT, 'CVode', &
                                'Trouble interpolating at tstop = '//g_str(me%tstop)// &
                                '. tstop too far back in direction of integration')
                            flag = CV_ILL_INPUT
                            return
                        end if
                    else
                        yout = me%zn(:,0)
                    end if
                    me%tretlast = me%tstop
                    tret = me%tretlast
                    me%tstopset = .false.
                    flag = CV_TSTOP_RETURN
                    return
                end if

            ! if next step would overtake tstop, adjust stepsize
            else if ((me%tn + me%hprime - me%tstop) * me%h > ZERO) then
                me%hprime = (me%tstop - me%tn) * (ONE - FOUR * me%uround)
                me%eta = me%hprime / me%h
            end if

        end if

        ! in CV_NORMAL mode, test if tout was reached
        if ((itask == CV_NORMAL) .and. ((me%tn - tout) * me%h >= ZERO)) then
            me%tretlast = tout
            tret = me%tretlast
            call CVodeGetDky(me, tout, 0, yout, ier)
            if (ier /= CV_SUCCESS) then
                call cvProcessError(me, CV_ILL_INPUT, 'CVode', &
                    'Trouble interpolating at tout = '//g_str(tout)// &
                    '. tout too far back in direction of integration')
                flag = CV_ILL_INPUT
                return
            end if
            flag = CV_SUCCESS
            return
        end if

        ! in CV_ONE_STEP mode, test if tn was returned
        if (itask == CV_ONE_STEP .and. abs(me%tn - me%tretlast) > troundoff) then
            me%tretlast = me%tn
            tret = me%tretlast
            yout = me%zn(:,0)
            flag = CV_SUCCESS
            return
        end if

    end if ! end stopping tests block

    !------------------------------------------------------------------
    ! looping point for internal steps
    !------------------------------------------------------------------

    nstloc = 0
    do

        me%next_h = me%h
        me%next_q = me%q

        ! reset and check ewt
        if (me%nst > 0) then
            ewtsetOK = cvEwtSet(me, me%zn(:,0), me%ewt)
            if (ewtsetOK /= 0) then
                if (me%itol == CV_WF) then
                    call cvProcessError(me, CV_ILL_INPUT, 'CVode', &
                        'At '//t_str(me%tn)//', the user-provide EwtSet function failed.')
                else
                    call cvProcessError(me, CV_ILL_INPUT, 'CVode', &
                        'At '//t_str(me%tn)//', a component of ewt has become <= 0.')
                end if
                istate = CV_ILL_INPUT
                me%tretlast = me%tn
                tret = me%tretlast
                yout = me%zn(:,0)
                exit
            end if
        end if

        ! check for too many steps
        if ((me%mxstep > 0) .and. (nstloc >= me%mxstep)) then
            call cvProcessError(me, CV_TOO_MUCH_WORK, 'CVode', &
                'At '//t_str(me%tn)//', mxstep steps taken before reaching tout.')
            istate = CV_TOO_MUCH_WORK
            me%tretlast = me%tn
            tret = me%tretlast
            yout = me%zn(:,0)
            exit
        end if

        ! check for too much accuracy requested
        nrm = n_vwrmsnorm(me%zn(:,0), me%ewt)
        me%tolsf = me%uround * nrm
        if (me%tolsf > ONE) then
            call cvProcessError(me, CV_TOO_MUCH_ACC, 'CVode', &
                'At '//t_str(me%tn)//', too much accuracy requested.')
            istate = CV_TOO_MUCH_ACC
            me%tretlast = me%tn
            tret = me%tretlast
            yout = me%zn(:,0)
            me%tolsf = me%tolsf * TWO
            exit
        else
            me%tolsf = ONE
        end if

        ! check for h below roundoff level in tn
        if (me%tn + me%h == me%tn) then
            me%nhnil = me%nhnil + 1
            if (me%nhnil <= me%mxhnil) then
                call cvProcessError(me, CV_WARNING, 'CVode', &
                    'Internal '//th_str(me%tn, me%h)// &
                    ' are such that t + h = t on the next step. The solver will continue anyway.')
            end if
            if (me%nhnil == me%mxhnil) then
                call cvProcessError(me, CV_WARNING, 'CVode', &
                    'The above warning has been issued mxhnil times and will not be issued '// &
                    'again for this problem.')
            end if
        end if

        ! call cvStep to take a step
        kflag = cvStep(me)

        ! process failed step cases, and exit loop
        if (kflag /= CV_SUCCESS) then
            istate = cvHandleFailure(me, kflag)
            me%tretlast = me%tn
            tret = me%tretlast
            yout = me%zn(:,0)
            exit
        end if

        nstloc = nstloc + 1

        ! if tstop is set and was reached, reset tn = tstop
        if (me%tstopset) then
            troundoff = FUZZ_FACTOR * me%uround * (abs(me%tn) + abs(me%h))
            if (abs(me%tn - me%tstop) <= troundoff) me%tn = me%tstop
        end if

        ! [not yet ported] root finding: cvRcheck3

        ! in CV_NORMAL mode, test if tout was reached
        ! (tstop checks come first)
        if (me%tstopset) then

            ! test for tn at tstop
            troundoff = FUZZ_FACTOR * me%uround * (abs(me%tn) + abs(me%h))
            if (abs(me%tn - me%tstop) <= troundoff) then
                ! ensure tout >= tstop, otherwise check for tout return below
                if ((tout - me%tstop) * me%h >= ZERO .or. &
                    abs(tout - me%tstop) <= troundoff) then
                    if (me%tstopinterp) then
                        call CVodeGetDky(me, me%tstop, 0, yout, ier)
                    else
                        yout = me%zn(:,0)
                    end if
                    me%tretlast = me%tstop
                    tret = me%tretlast
                    me%tstopset = .false.
                    istate = CV_TSTOP_RETURN
                    exit
                end if

            ! if next step would overtake tstop, adjust stepsize
            else if ((me%tn + me%hprime - me%tstop) * me%h > ZERO) then
                me%hprime = (me%tstop - me%tn) * (ONE - FOUR * me%uround)
                me%eta = me%hprime / me%h
            end if

        end if

        ! in CV_NORMAL mode, test if tout was reached
        if ((itask == CV_NORMAL) .and. (me%tn - tout) * me%h >= ZERO) then
            istate = CV_SUCCESS
            me%tretlast = tout
            tret = me%tretlast
            call CVodeGetDky(me, tout, 0, yout, ier)
            me%next_q = me%qprime
            me%next_h = me%hprime
            exit
        end if

        ! in CV_ONE_STEP mode, return now
        if (itask == CV_ONE_STEP) then
            istate = CV_SUCCESS
            me%tretlast = me%tn
            tret = me%tretlast
            yout = me%zn(:,0)
            me%next_q = me%qprime
            me%next_h = me%hprime
            exit
        end if

    end do ! end looping for internal steps

    flag = istate

    end subroutine CVode
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeGetDky`: the k-th derivative of the interpolating polynomial at `t`
!  (`tn - hu <= t <= tn`, `0 <= k <= q`).

    subroutine CVodeGetDky(me, t, k, dky, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    integer,intent(in) :: k
    real(wp),dimension(:),intent(inout) :: dky
    integer,intent(out) :: flag

    real(wp) :: s, r, c
    real(wp) :: tfuzz, tp, tn1
    integer :: i, j

    if (size(dky, kind=ip) /= me%n) then
        call cvProcessError(me, CV_BAD_DKY, 'CVodeGetDky', 'dky = NULL illegal.')
        flag = CV_BAD_DKY
        return
    end if

    if ((k < 0) .or. (k > me%q)) then
        call cvProcessError(me, CV_BAD_K, 'CVodeGetDky', 'Illegal value for k.')
        flag = CV_BAD_K
        return
    end if

    ! allow for some slack
    tfuzz = FUZZ_FACTOR * me%uround * (abs(me%tn) + abs(me%hu))
    if (me%hu < ZERO) tfuzz = -tfuzz
    tp  = me%tn - me%hu - tfuzz
    tn1 = me%tn + tfuzz
    if ((t - tp) * (t - tn1) > ZERO) then
        call cvProcessError(me, CV_BAD_T, 'CVodeGetDky', &
            'Illegal value for t. t = '//g_str(t)//' is not between tcur - hu = '// &
            g_str(me%tn - me%hu)//' and tcur = '//g_str(me%tn)//'.')
        flag = CV_BAD_T
        return
    end if

    ! sum the differentiated interpolating polynomial:
    ! dky = sum_{j=k}^{q} [j!/(j-k)!] s^(j-k) zn(:,j)
    s = (t - me%tn) / me%h
    dky = ZERO
    do j = me%q, k, -1
        c = ONE
        do i = j, j - k + 1, -1
            c = c * i
        end do
        c = c * s**(j - k)
        dky = dky + c * me%zn(:,j)
    end do

    if (k == 0) then
        flag = CV_SUCCESS
        return
    end if
    r = (me%h)**(-k)
    dky = r*dky

    flag = CV_SUCCESS

    end subroutine CVodeGetDky
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeFree`: free all memory (also called by `create`).

    subroutine CVodeFree(me)

    class(cvode_solver),intent(inout) :: me

    call cvFreeVectors(me)
    if (allocated(me%NLS)) deallocate(me%NLS)
    me%MallocDone = .false.
    me%f => null()
    me%efun_user => null()

    end subroutine CVodeFree
!*****************************************************************************************

!=========================================================================================
!  Private helper functions (cvode.c)
!=========================================================================================

!*****************************************************************************************
!>
!  `cvAllocVectors`: allocate the CVODE vectors ewt, acor, tempv, ftemp,
!  vtemp1-3 and zn[0..qmax], all of length `n`.

    subroutine cvAllocVectors(me, n)

    class(cvode_solver),intent(inout) :: me
    integer(ip),intent(in) :: n

    call cvFreeVectors(me)

    me%n = n
    allocate(me%ewt(n), me%acor(n), me%tempv(n), me%ftemp(n), &
             me%vtemp1(n), me%vtemp2(n), me%vtemp3(n), me%y(n))
    allocate(me%zn(n, 0:me%qmax))

    ! update solver workspace lengths (not tracked) and
    ! store the value of qmax used here
    me%qmax_alloc = me%qmax

    end subroutine cvAllocVectors
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvFreeVectors`

    subroutine cvFreeVectors(me)

    class(cvode_solver),intent(inout) :: me

    if (allocated(me%ewt))     deallocate(me%ewt)
    if (allocated(me%acor))    deallocate(me%acor)
    if (allocated(me%tempv))   deallocate(me%tempv)
    if (allocated(me%ftemp))   deallocate(me%ftemp)
    if (allocated(me%vtemp1))  deallocate(me%vtemp1)
    if (allocated(me%vtemp2))  deallocate(me%vtemp2)
    if (allocated(me%vtemp3))  deallocate(me%vtemp3)
    if (allocated(me%y))       deallocate(me%y)
    if (allocated(me%zn))      deallocate(me%zn)
    if (allocated(me%Vabstol)) deallocate(me%Vabstol)
    me%n = 0

    end subroutine cvFreeVectors
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvInitialSetup`: checks done at the first call to `CVode`; computes the
!  initial error weights and initializes the nonlinear solver.

    function cvInitialSetup(me, tout) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: tout
    integer :: flag

    integer :: ier
    real(wp) :: tdist, tround

    ! check that tout is not too close to the initial time
    tdist  = abs(tout - me%tn)
    tround = me%uround * max(abs(me%tn), abs(tout))

    if (tdist == ZERO .or. tdist < TWO * tround) then
        call cvProcessError(me, CV_TOO_CLOSE, 'cvInitialSetup', &
                            'tout too close to t0 to start integration.')
        flag = CV_TOO_CLOSE
        return
    end if

    ! did the user specify tolerances?
    if (me%itol == CV_NN) then
        call cvProcessError(me, CV_ILL_INPUT, 'cvInitialSetup', &
                            'No integration tolerances have been specified.')
        flag = CV_ILL_INPUT
        return
    end if

    ! [not yet ported] check constraints against y0

    ! load initial error weights
    ier = cvEwtSet(me, me%zn(:,0), me%ewt)
    if (ier /= 0) then
        if (me%itol == CV_WF) then
            call cvProcessError(me, CV_ILL_INPUT, 'cvInitialSetup', &
                                'The user-provide EwtSet function failed.')
        else
            call cvProcessError(me, CV_ILL_INPUT, 'cvInitialSetup', &
                                'Initial ewt has component(s) equal to zero (illegal).')
        end if
        flag = CV_ILL_INPUT
        return
    end if

    ! [not yet ported] call linit function (if it exists)

    ! initialize the nonlinear solver (must occur after linear solver is
    ! initialized) so that lsetup and lsolve pointer have been set
    ier = cvNlsInit(me)
    if (ier /= 0) then
        call cvProcessError(me, CV_NLS_INIT_FAIL, 'cvInitialSetup', &
                            "The nonlinear solver's init routine failed.")
        flag = CV_NLS_INIT_FAIL
        return
    end if

    ! [not yet ported] projection initialization

    flag = CV_SUCCESS

    end function cvInitialSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvHin`: compute a tentative initial step size h0.
!
!  The step size is chosen so that the local error of a first-order step is
!  about the requested tolerance, based on an estimate of the second
!  derivative `y''` from difference quotients. It is bounded below by
!  `HLB_FACTOR * tround` and above by [[cvUpperBoundH0]].

    function cvHin(me, tout) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: tout
    integer :: flag

    integer :: retval, sgn, count1, count2
    real(wp) :: tdiff, tdist, tround, hlb, hub
    real(wp) :: hg, hgs, hs, hnew, hrat, h0, yddnrm
    logical :: hgOK

    ! if tout is too close to tn, give up
    tdiff  = tout - me%tn
    sgn    = merge(1, -1, tdiff > ZERO)
    tdist  = abs(tdiff)
    tround = me%uround * max(abs(me%tn), abs(tout))

    ! set lower and upper bounds on h0, and take geometric mean as first
    ! trial value. Exit with this value if the bounds cross each other.
    hlb = HLB_FACTOR * tround
    hub = cvUpperBoundH0(me, tdist)

    hg = sqrt(hlb * hub)

    if (hub < hlb) then
        if (sgn == -1) then
            me%h = -hg
        else
            me%h = hg
        end if
        flag = CV_SUCCESS
        return
    end if

    ! outer loop
    hs = hg  ! safeguard against 'uninitialized variable' warning

    hnew = hs
    do count1 = 1, MAX_ITERS

        ! attempts to estimate ydd
        hgOK = .false.

        do count2 = 1, MAX_ITERS
            hgs    = hg * sgn
            retval = cvYddNorm(me, hgs, yddnrm)
            ! if a RHS function failed unrecoverably, give up
            if (retval < 0) then
                flag = CV_RHSFUNC_FAIL
                return
            end if
            ! if successful, we can use ydd
            if (retval == CV_SUCCESS) then
                hgOK = .true.
                exit
            end if
            ! a RHS function failed recoverably; cut step size and test it again
            hg = hg * POINT2
        end do

        ! if a RHS function failed recoverably MAX_ITERS times
        if (.not. hgOK) then
            ! exit if this is the first or second pass. No recovery possible
            if (count1 <= 2) then
                flag = CV_REPTD_RHSFUNC_ERR
                return
            end if
            ! we have a fall-back option. The value hs is a previous hnew which
            ! passed through f(). Use it and break
            hnew = hs
            exit
        end if

        ! the proposed step size is feasible. Save it.
        hs = hg

        ! propose new step size
        if (yddnrm * hub * hub > TWO) then
            hnew = sqrt(TWO / yddnrm)
        else
            hnew = sqrt(hg * hub)
        end if

        ! if last pass, stop now with hnew
        if (count1 == MAX_ITERS) exit

        hrat = hnew / hg

        ! accept hnew if it does not differ from hg by more than a factor of 2
        if ((hrat > HALF) .and. (hrat < TWO)) exit

        ! after one pass, if ydd seems to be bad, use fall-back value
        if ((count1 > 1) .and. (hrat > TWO)) then
            hnew = hg
            exit
        end if

        ! send this value back through f()
        hg = hnew

    end do

    ! apply bounds, bias factor, and attach sign
    h0 = H_BIAS * hnew
    if (h0 < hlb) h0 = hlb
    if (h0 > hub) h0 = hub
    if (sgn == -1) h0 = -h0
    me%h = h0

    flag = CV_SUCCESS

    end function cvHin
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvUpperBoundH0`: upper bound on the initial step size, based on
!  `tdist = |tout - t0|` and the initial values of y and y'.

    function cvUpperBoundH0(me, tdist) result(hub)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: tdist
    real(wp) :: hub

    real(wp) :: hub_inv
    integer :: ier

    ! bound based on |y0|/|y0'| -- allow at most an increase of
    ! HUB_FACTOR in y0 (based on a forward Euler step). The weight
    ! factor is used as a safeguard against zero components in y0.
    ! temp1 = tempv, temp2 = acor
    me%acor = abs(me%zn(:,0))
    ier = cvEwtSet(me, me%zn(:,0), me%tempv)
    me%tempv = ONE / me%tempv
    me%tempv = HUB_FACTOR*me%acor + me%tempv

    me%acor = abs(me%zn(:,1))

    me%tempv = me%acor / me%tempv
    hub_inv = maxval(abs(me%tempv))

    ! bound based on tdist -- with upper bound
    hub = HUB_FACTOR * tdist

    ! use the smaller of the two
    if (hub * hub_inv > ONE) hub = ONE / hub_inv

    end function cvUpperBoundH0
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvYddNorm`: estimate the WRMS norm of y'' from
!  `[f(t0+hg, y0+hg*y0') - y0'] / hg`.

    function cvYddNorm(me, hg, yddnrm) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: hg
    real(wp),intent(out) :: yddnrm
    integer :: flag

    integer :: retval

    me%y = hg*me%zn(:,1) + me%zn(:,0)
    retval = me%f(me, me%tn + hg, me%y, me%tempv)
    me%nfe = me%nfe + 1
    if (retval < 0) then
        flag = CV_RHSFUNC_FAIL
        return
    end if
    if (retval > 0) then
        flag = RHSFUNC_RECVR
        return
    end if

    me%tempv = (ONE / hg)*me%tempv + (-ONE / hg)*me%zn(:,1)

    yddnrm = n_vwrmsnorm(me%tempv, me%ewt)

    flag = CV_SUCCESS

    end function cvYddNorm
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvStep`: take one internal step.
!
!  Predict, solve the nonlinear system, run the error test, and either
!  complete the step or retry with a smaller step / lower order. Returns
!  `CV_SUCCESS` or a failure flag for [[cvHandleFailure]].

    function cvStep(me) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer :: flag

    real(wp) :: saved_t  ! time to restore to if a failure occurs
    real(wp) :: dsm      ! local truncation error estimate
    integer :: nflag, kflag ! nonlinear solver flags
    integer :: eflag     ! error test return flag

    ! are we computing sensitivities with the simultaneous approach?
    integer :: ncf  ! corrector failures in this step attempt
    integer :: nef  ! error test failures in this step attempt

    ncf = 0
    nef = 0

    ! if the step size has changed, update the history array
    if ((me%nst > 0) .and. (me%hprime /= me%h)) call cvAdjustParams(me)

    ! [not yet ported] check if projection is needed

    ! looping point for attempts to take a step
    saved_t = me%tn
    nflag   = FIRST_CALL

    do

        call cvPredict(me)
        call cvSet(me)

        nflag = cvNls(me, nflag)
        kflag = cvHandleNFlag(me, nflag, saved_t, ncf)

        ! go back in loop if we need to predict again (nflag=PREV_CONV_FAIL)
        if (kflag == PREDICT_AGAIN) cycle

        ! return if nonlinear solve failed and recovery not possible
        if (kflag /= DO_ERROR_TEST) then
            flag = kflag
            return
        end if

        ! [not yet ported] check if constraints are satisfied
        ! [not yet ported] perform projection

        ! perform error test (nflag=CV_SUCCESS)
        eflag = cvDoErrorTest(me, nflag, saved_t, nef, dsm)

        ! go back in loop if we need to predict again (nflag=PREV_ERR_FAIL)
        if (eflag == TRY_AGAIN) cycle

        ! return if error test failed and recovery not possible
        if (eflag /= CV_SUCCESS) then
            flag = eflag
            return
        end if

        ! error test passed (eflag=CV_SUCCESS), break from loop
        exit

    end do

    ! nonlinear system solve and error test were both successful;
    ! update data, and consider change of step and/or order
    call cvCompleteStep(me)

    call cvPrepareNextStep(me, dsm)

    ! if Stablilty Limit Detection is turned on, call stability limit
    ! detection routine for possible order reduction
    ! [not yet ported] if (me%sldeton) call cvBDFStab(me)

    if (me%nst <= me%small_nst) then
        me%etamax = me%eta_max_es
    else
        me%etamax = me%eta_max_gs
    end if

    ! finally, we rescale the acor array to be the estimated local error vector
    me%acor = me%tq(2)*me%acor

    flag = CV_SUCCESS

    end function cvStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAdjustParams`: called when a change in step size was decided upon;
!  adjusts the history array (order change, then rescaling).

    subroutine cvAdjustParams(me)

    class(cvode_solver),intent(inout) :: me

    if (me%qprime /= me%q) then
        ! [not yet ported] history resize: skip order adjustment after a resize
        call cvAdjustOrder(me, me%qprime - me%q)
        me%q     = me%qprime
        me%LL    = me%q + 1
        me%qwait = me%LL
    end if
    call cvRescale(me)

    end subroutine cvAdjustParams
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAdjustOrder`: adjust the history array on a change of order q by
!  `deltaq` (+1 or -1).

    subroutine cvAdjustOrder(me, deltaq)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: deltaq

    if ((me%q == 2) .and. (deltaq /= 1)) return

    select case (me%lmm)
    case (CV_ADAMS); call cvAdjustAdams(me, deltaq)
    case (CV_BDF);   call cvAdjustBDF(me, deltaq)
    end select

    end subroutine cvAdjustOrder
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAdjustAdams`: adjust the history array on a change of order for the
!  Adams method.

    subroutine cvAdjustAdams(me, deltaq)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: deltaq

    integer :: i, j
    real(wp) :: xi, hsum

    ! on an order increase, set new column of zn to zero and return
    if (deltaq == 1) then
        me%zn(:,me%LL) = ZERO
        return
    end if

    ! on an order decrease, each zn[j] is adjusted by a multiple of zn[q].
    ! The coeffs. in the adjustment are the coeffs. of the polynomial:
    !        x
    ! q * INT { u * ( u + xi_1 ) * ... * ( u + xi_{q-2} ) } du
    !        0
    ! where xi_j = [t_n - t_(n-j)]/h => xi_0 = 0

    do i = 0, me%qmax
        me%l(i) = ZERO
    end do
    me%l(1) = ONE
    hsum = ZERO
    do j = 1, me%q - 2
        hsum = hsum + me%tau(j)
        xi = hsum / me%hscale
        do i = j + 1, 1, -1
            me%l(i) = me%l(i) * xi + me%l(i - 1)
        end do
    end do

    do j = 1, me%q - 2
        me%l(j + 1) = me%q * (me%l(j) / (j + 1))
    end do

    ! zn(:,j) = zn(:,j) - l(j) * zn(:,q), j = 2..q-1
    do j = 2, me%q - 1
        me%zn(:,j) = me%zn(:,j) - me%l(j) * me%zn(:,me%q)
    end do

    end subroutine cvAdjustAdams
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAdjustBDF`: adjust the history array on a change of order for BDF.

    subroutine cvAdjustBDF(me, deltaq)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: deltaq

    select case (deltaq)
    case (1);  call cvIncreaseBDF(me)
    case (-1); call cvDecreaseBDF(me)
    end select

    end subroutine cvAdjustBDF
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvIncreaseBDF`: adjust the history array on an increase in the order q
!  for BDF. A new column zn[q+1] is set equal to a multiple of the saved
!  vector (= acor) in zn[indx_acor]. Then each zn[j] is adjusted by a
!  multiple of zn[q+1].

    subroutine cvIncreaseBDF(me)

    class(cvode_solver),intent(inout) :: me

    real(wp) :: alpha0, alpha1, prod, xi, xiold, hsum, A1
    integer :: i, j

    do i = 0, me%qmax
        me%l(i) = ZERO
    end do
    me%l(2) = ONE
    alpha1  = ONE
    prod    = ONE
    xiold   = ONE
    alpha0  = -ONE
    hsum    = me%hscale
    if (me%q > 1) then
        do j = 1, me%q - 1
            hsum = hsum + me%tau(j + 1)
            xi = hsum / me%hscale
            prod = prod * xi
            alpha0 = alpha0 - ONE / (j + 1)
            alpha1 = alpha1 + ONE / xi
            do i = j + 2, 2, -1
                me%l(i) = me%l(i) * xiold + me%l(i - 1)
            end do
            xiold = xi
        end do
    end if
    A1 = (-alpha0 - alpha1) / prod
    me%zn(:,me%LL) = A1*me%zn(:,me%indx_acor)

    ! zn(:,j) = zn(:,j) + l(j) * zn(:,L), j = 2..q
    do j = 2, me%q
        me%zn(:,j) = me%zn(:,j) + me%l(j) * me%zn(:,me%LL)
    end do

    end subroutine cvIncreaseBDF
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvDecreaseBDF`: adjust the history array on a decrease in the order q
!  for BDF. Each zn[j] is adjusted by a multiple of zn[q].

    subroutine cvDecreaseBDF(me)

    class(cvode_solver),intent(inout) :: me

    real(wp) :: hsum, xi
    integer :: i, j

    do i = 0, me%qmax
        me%l(i) = ZERO
    end do
    me%l(2) = ONE
    hsum = ZERO
    do j = 1, me%q - 2
        hsum = hsum + me%tau(j)
        xi = hsum / me%hscale
        do i = j + 2, 2, -1
            me%l(i) = me%l(i) * xi + me%l(i - 1)
        end do
    end do

    ! zn(:,j) = zn(:,j) - l(j) * zn(:,q), j = 2..q-1
    do j = 2, me%q - 1
        me%zn(:,j) = me%zn(:,j) - me%l(j) * me%zn(:,me%q)
    end do

    end subroutine cvDecreaseBDF
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvRescale`: rescale the Nordsieck array by multiplying the j-th column
!  zn[j] by eta^j, j = 1, ..., q. Then the value of h is rescaled by eta,
!  and hscale is reset to h.

    subroutine cvRescale(me)

    class(cvode_solver),intent(inout) :: me

    integer :: j
    real(wp) :: factor

    ! zn(:,j) = eta**j * zn(:,j), j = 1..q
    factor = ONE
    do j = 1, me%q
        factor = factor * me%eta
        me%zn(:,j) = factor * me%zn(:,j)
    end do

    me%h      = me%hscale * me%eta
    me%next_h = me%h
    me%hscale = me%h
    me%nscon  = 0

    end subroutine cvRescale
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvPredict`: advance tn by the tentative step size h, and compute the
!  predicted array z_n(0), which is overwritten on zn. The prediction of zn
!  is done by repeated additions. If tstop is enabled, it is possible for
!  tn + h to be past tstop by roundoff, and in that case, we reset tn
!  (after incrementing by h) to tstop.

    subroutine cvPredict(me)

    class(cvode_solver),intent(inout) :: me

    integer :: j, k

    me%tn = me%tn + me%h
    if (me%tstopset) then
        if ((me%tn - me%tstop) * me%h > ZERO) me%tn = me%tstop
    end if

    do k = 1, me%q
        do j = me%q, k, -1
            me%zn(:,j-1) = me%zn(:,j-1) + me%zn(:,j)
        end do
    end do

    end subroutine cvPredict
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvSet`: compute the coefficients l and tq used in the current step,
!  then rl1 = 1/l(1), gamma = h*rl1 and gamrat.

    subroutine cvSet(me)

    class(cvode_solver),intent(inout) :: me

    select case (me%lmm)
    case (CV_ADAMS); call cvSetAdams(me)
    case (CV_BDF);   call cvSetBDF(me)
    end select
    me%rl1   = ONE / me%l(1)
    me%gamma = me%h * me%rl1
    if (me%nst == 0) me%gammap = me%gamma
    if (me%nst > 0) then
        me%gamrat = me%gamma / me%gammap
    else
        me%gamrat = ONE  ! protect x / x != 1.0
    end if

    end subroutine cvSet
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvSetAdams`: compute the coefficients l and tq for the Adams method.
!
!  The tq values are the error test quantities `tq(2)`, the quantities used
!  to decide on an order change `tq(1)` (q-1) and `tq(3)` (q+1), and
!  `tq(4)` for the nonlinear convergence test; `tq(5)` is used when
!  increasing the order.

    subroutine cvSetAdams(me)

    class(cvode_solver),intent(inout) :: me

    real(wp),dimension(0:L_MAX-1) :: m
    real(wp),dimension(0:2) :: mm   ! C: M[3]
    real(wp) :: hsum

    if (me%q == 1) then
        me%l(0)  = ONE
        me%l(1)  = ONE
        me%tq(1) = ONE
        me%tq(5) = ONE
        me%tq(2) = HALF
        me%tq(3) = ONE / TWELVE
        me%tq(4) = me%nlscoef / me%tq(2)  ! = 0.1 / tq[2]
        return
    end if

    hsum = cvAdamsStart(me, m)

    mm(0) = cvAltSum(me%q - 1, m, 1)
    mm(1) = cvAltSum(me%q - 1, m, 2)

    call cvAdamsFinish(me, m, mm, hsum)

    end subroutine cvSetAdams
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAdamsStart`: generate in m[] the coefficients of the product polynomial
!  needed for the Adams l and tq coefficients for q > 1.

    function cvAdamsStart(me, m) result(hsum)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(0:L_MAX-1),intent(inout) :: m
    real(wp) :: hsum

    real(wp) :: xi_inv, sum
    integer :: i, j

    hsum = me%h
    m(0) = ONE
    do i = 1, me%q
        m(i) = ZERO
    end do
    do j = 1, me%q - 1
        if ((j == me%q - 1) .and. (me%qwait == 1)) then
            sum = cvAltSum(me%q - 2, m, 2)
            me%tq(1) = me%q * sum / m(me%q - 2)
        end if
        xi_inv = me%h / hsum
        do i = j, 1, -1
            m(i) = m(i) + m(i - 1) * xi_inv
        end do
        hsum = hsum + me%tau(j)
        ! the m[i] are coefficients of product(1 to j) (1 + x/xi_i)
    end do

    end function cvAdamsStart
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAdamsFinish`: complete the calculation of the Adams l and tq.

    subroutine cvAdamsFinish(me, m, mm, hsum)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(0:L_MAX-1),intent(inout) :: m
    real(wp),dimension(0:2),intent(inout) :: mm
    real(wp),intent(in) :: hsum

    integer :: i
    real(wp) :: M0_inv, xi, xi_inv

    M0_inv = ONE / mm(0)

    me%l(0) = ONE
    do i = 1, me%q
        me%l(i) = M0_inv * (m(i - 1) / i)
    end do
    xi     = hsum / me%h
    xi_inv = ONE / xi

    me%tq(2) = mm(1) * M0_inv / xi
    me%tq(5) = xi / me%l(me%q)

    if (me%qwait == 1) then
        do i = me%q, 1, -1
            m(i) = m(i) + m(i - 1) * xi_inv
        end do
        mm(2) = cvAltSum(me%q, m, 2)
        me%tq(3) = mm(2) * M0_inv / me%LL
    end if

    me%tq(4) = me%nlscoef / me%tq(2)

    end subroutine cvAdamsFinish
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAltSum`: returns the value of the alternating sum
!  `sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ]`.
!  If iend < 0 then cvAltSum returns 0.

    pure function cvAltSum(iend, a, k) result(sum)

    integer,intent(in) :: iend
    real(wp),dimension(0:),intent(in) :: a
    integer,intent(in) :: k
    real(wp) :: sum

    integer :: i, sgn

    if (iend < 0) then
        sum = ZERO
        return
    end if

    sum = ZERO
    sgn = 1
    do i = 0, iend
        sum = sum + sgn * (a(i) / (i + k))
        sgn = -sgn
    end do

    end function cvAltSum
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvSetBDF`: compute the coefficients l and tq for the case lmm == CV_BDF.
!  Calls cvSetTqBDF to set the test quantity array tq.
!
!  The components of the array l are the coefficients of a polynomial
!  Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
!                                  q-1
!  (d/dx) Lambda(x) = c * PRODUCT (1 + x / xi_i) , where
!                                  i=1
!  Lambda(-1) = 0, Lambda(0) = 1, and c is a normalization factor.
!  Here xi_i = [t_n - t_(n-i)] / h.

    subroutine cvSetBDF(me)

    class(cvode_solver),intent(inout) :: me

    real(wp) :: alpha0, alpha0_hat, xi_inv, xistar_inv, hsum
    integer :: i, j

    me%l(0)    = ONE
    me%l(1)    = ONE
    xi_inv     = ONE
    xistar_inv = ONE
    do i = 2, me%q
        me%l(i) = ZERO
    end do
    alpha0     = -ONE
    alpha0_hat = -ONE
    hsum       = me%h

    ! [not yet ported] projection coefficients proj_p

    if (me%q > 1) then
        do j = 2, me%q - 1
            hsum = hsum + me%tau(j - 1)
            xi_inv = me%h / hsum
            alpha0 = alpha0 - ONE / j
            do i = j, 1, -1
                me%l(i) = me%l(i) + me%l(i - 1) * xi_inv
            end do
            ! the l[i] are coefficients of product(1 to j) (1 + x/xi_i)
        end do

        ! j = q
        alpha0     = alpha0 - ONE / me%q
        xistar_inv = -me%l(1) - alpha0
        hsum       = hsum + me%tau(me%q - 1)
        xi_inv     = me%h / hsum
        alpha0_hat = -me%l(1) - xi_inv

        do i = me%q, 1, -1
            me%l(i) = me%l(i) + me%l(i - 1) * xistar_inv
        end do
    end if

    call cvSetTqBDF(me, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv)

    end subroutine cvSetBDF
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvSetTqBDF`: set the test quantity array tq for BDF.

    subroutine cvSetTqBDF(me, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: hsum
    real(wp),intent(in) :: alpha0
    real(wp),intent(in) :: alpha0_hat
    real(wp),intent(in) :: xi_inv
    real(wp),intent(in) :: xistar_inv

    real(wp) :: A1, A2, A3, A4, A5, A6
    real(wp) :: C, Cpinv, Cppinv
    real(wp) :: hsum_, xi_inv_   ! C modifies its by-value copies

    hsum_   = hsum
    xi_inv_ = xi_inv

    A1 = ONE - alpha0_hat + alpha0
    A2 = ONE + me%q * A1
    me%tq(2) = abs(A1 / (alpha0 * A2))
    me%tq(5) = abs(A2 * xistar_inv / (me%l(me%q) * xi_inv_))
    if (me%qwait == 1) then
        if (me%q > 1) then
            C        = xistar_inv / me%l(me%q)
            A3       = alpha0 + ONE / me%q
            A4       = alpha0_hat + xi_inv_
            Cpinv    = (ONE - A4 + A3) / A3
            me%tq(1) = abs(C * Cpinv)
        else
            me%tq(1) = ONE
        end if
        hsum_    = hsum_ + me%tau(me%q)
        xi_inv_  = me%h / hsum_
        A5       = alpha0 - (ONE / (me%q + 1))
        A6       = alpha0_hat - xi_inv_
        Cppinv   = (ONE - A6 + A5) / A2
        me%tq(3) = abs(Cppinv / (xi_inv_ * (me%q + 2) * A5))
    end if
    me%tq(4) = me%nlscoef / me%tq(2)

    end subroutine cvSetTqBDF
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNls`: attempt to solve the nonlinear system using the attached
!  nonlinear solver.

    function cvNls(me, nflag) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: nflag
    integer :: flag

    logical :: callSetup
    integer(li) :: nni_inc, nnf_inc

    flag    = CV_SUCCESS
    nni_inc = 0
    nnf_inc = 0

    ! decide whether or not to call setup routine (if one exists) and
    ! set flag convfail (input to lsetup for its evaluation decision)
    if (me%lsetup_exists) then
        if ((nflag == FIRST_CALL) .or. (nflag == PREV_ERR_FAIL)) then
            me%convfail = CV_NO_FAILURES
        else
            me%convfail = CV_FAIL_OTHER
        end if

        callSetup = (nflag == PREV_CONV_FAIL) .or. (nflag == PREV_ERR_FAIL) .or. &
                    (me%nst == 0) .or. &
                    (me%nst >= me%nstlp + me%msbp) .or. &
                    (abs(me%gamrat - ONE) > me%dgmax_lsetup)
    else
        me%crate  = ONE
        callSetup = .false.
    end if

    ! initial guess for the correction to the predictor
    me%acor = ZERO

    ! call nonlinear solver setup if it exists
    flag = me%NLS%setup(me, me%acor)
    if (flag < 0) then
        flag = CV_NLS_SETUP_FAIL
        return
    end if
    if (flag > 0) then
        flag = SUN_NLS_CONV_RECVR
        return
    end if

    ! solve the nonlinear system
    flag = me%NLS%solve(me, me%zn(:,0), me%acor, me%ewt, me%tq(4), callSetup)

    ! increment counters
    nni_inc = me%NLS%get_num_iters()
    me%nni = me%nni + nni_inc

    nnf_inc = me%NLS%get_num_conv_fails()
    me%nnf = me%nnf + nnf_inc

    ! if the solve failed return
    if (flag /= SUN_SUCCESS) return

    ! solve successful

    ! update the state based on the final correction from the nonlinear solver
    me%y = me%zn(:,0) + me%acor

    ! compute acnrm if is was not already done by the nonlinear solver
    if (.not. me%acnrmcur) me%acnrm = n_vwrmsnorm(me%acor, me%ewt)

    ! update Jacobian status
    me%jcur = .false.

    end function cvNls
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvHandleNFlag`: take action on the return value nflag = cvNls, and
!  return `DO_ERROR_TEST`, `PREDICT_AGAIN` or a failure flag.

    function cvHandleNFlag(me, nflag, saved_t, ncf) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(inout) :: nflag
    real(wp),intent(in) :: saved_t
    integer,intent(inout) :: ncf
    integer :: flag

    if (nflag == CV_SUCCESS) then
        flag = DO_ERROR_TEST
        return
    end if

    ! the nonlinear solve failed
    me%ncfn = me%ncfn + 1
    call cvRestore(me, saved_t)

    ! return if failed unrecoverably
    if (nflag < 0) then
        if (nflag == CV_LSETUP_FAIL) then
            flag = CV_LSETUP_FAIL
        else if (nflag == CV_LSOLVE_FAIL) then
            flag = CV_LSOLVE_FAIL
        else if (nflag == CV_RHSFUNC_FAIL) then
            flag = CV_RHSFUNC_FAIL
        else
            flag = CV_NLS_FAIL
        end if
        return
    end if

    ! at this point, a recoverable error occurred
    ncf = ncf + 1
    me%etamax = ONE

    ! if we had maxncf failures or |h| = hmin,
    ! return CV_CONV_FAILURE, CV_REPTD_RHSFUNC_ERR
    if ((abs(me%h) <= me%hmin * ONEPSM) .or. (ncf == me%maxncf)) then
        if (nflag == SUN_NLS_CONV_RECVR) then
            flag = CV_CONV_FAILURE
            return
        end if
        if (nflag == RHSFUNC_RECVR) then
            flag = CV_REPTD_RHSFUNC_ERR
            return
        end if
    end if

    ! reduce step size; return to reattempt the step
    me%eta = max(me%eta_cf, me%hmin / abs(me%h))
    nflag  = PREV_CONV_FAIL
    call cvRescale(me)

    flag = PREDICT_AGAIN

    end function cvHandleNFlag
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvRestore`: restore the value of tn to saved_t and undo the
!  prediction. After execution of cvRestore, the Nordsieck array zn has the
!  same values as before the call to cvPredict.

    subroutine cvRestore(me, saved_t)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: saved_t

    integer :: j, k

    me%tn = saved_t
    do k = 1, me%q
        do j = me%q, k, -1
            me%zn(:,j-1) = me%zn(:,j-1) - me%zn(:,j)
        end do
    end do

    end subroutine cvRestore
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvDoErrorTest`: perform the local error test.
!
!  The weighted local error norm dsm is loaded into dsm, and the test
!  dsm ?<= 1 is made. If the test passes, `CV_SUCCESS` is returned.
!  Otherwise the step is restored and `TRY_AGAIN` is returned after reducing
!  h (and possibly q), or `CV_ERR_FAILURE` if recovery is not possible.

    function cvDoErrorTest(me, nflag, saved_t, nef, dsm) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(inout) :: nflag
    real(wp),intent(in) :: saved_t
    integer,intent(inout) :: nef
    real(wp),intent(out) :: dsm
    integer :: flag

    integer :: retval

    dsm = me%acnrm * me%tq(2)

    ! if est. local error norm dsm passes test, return CV_SUCCESS
    if (dsm <= ONE) then
        flag = CV_SUCCESS
        return
    end if

    ! test failed; increment counters, set nflag, and restore zn array
    nef     = nef + 1
    me%netf = me%netf + 1
    nflag   = PREV_ERR_FAIL
    call cvRestore(me, saved_t)

    ! at maxnef failures or |h| = hmin, return CV_ERR_FAILURE
    if ((abs(me%h) <= me%hmin * ONEPSM) .or. (nef == me%maxnef)) then
        flag = CV_ERR_FAILURE
        return
    end if

    ! set etamax = 1 to prevent step size increase at end of this step
    me%etamax = ONE

    ! set h ratio eta from dsm, rescale, and return for retry of step
    if (nef <= MXNEF1) then
        me%eta = ONE / ((BIAS2 * dsm)**(ONE / me%LL) + ADDON)
        me%eta = max(me%eta_min_ef, max(me%eta, me%hmin / abs(me%h)))
        if (nef >= me%small_nef) me%eta = min(me%eta, me%eta_max_ef)

        call cvRescale(me)

        flag = TRY_AGAIN
        return
    end if

    ! after MXNEF1 failures, force an order reduction and retry step
    if (me%q > 1) then
        me%eta = max(me%eta_min_ef, me%hmin / abs(me%h))
        call cvAdjustOrder(me, -1)
        me%LL    = me%q
        me%q     = me%q - 1
        me%qwait = me%LL
        call cvRescale(me)

        flag = TRY_AGAIN
        return
    end if

    ! if already at order 1, restart: reload zn, zn[1] from scratch
    me%eta = max(me%eta_min_ef, me%hmin / abs(me%h))
    me%h      = me%h * me%eta
    me%next_h = me%h
    me%hscale = me%h
    me%qwait  = LONG_WAIT
    me%nscon  = 0

    retval = me%f(me, me%tn, me%zn(:,0), me%tempv)
    me%nfe = me%nfe + 1
    if (retval < 0) then
        flag = CV_RHSFUNC_FAIL
        return
    end if
    if (retval > 0) then
        flag = CV_UNREC_RHSFUNC_ERR
        return
    end if

    me%zn(:,1) = me%h*me%tempv

    flag = TRY_AGAIN

    end function cvDoErrorTest
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvCompleteStep`: perform various update operations when the solution to
!  the nonlinear system has passed the local error test: increment the step
!  counter nst, record the values hu and qu, update the tau array, and apply
!  the corrections to the zn array. The tau[i] are the last q values of h,
!  with tau[1] the most recent. The counter qwait is decremented, and if
!  qwait == 1 (and q < qmax) we save acor and tq[5] for a possible order
!  increase.

    subroutine cvCompleteStep(me)

    class(cvode_solver),intent(inout) :: me

    integer :: i, j

    me%nst   = me%nst + 1
    me%nscon = me%nscon + 1
    me%hu    = me%h
    me%qu    = me%q

    do i = me%q, 2, -1
        me%tau(i) = me%tau(i - 1)
    end do
    if ((me%q == 1) .and. (me%nst > 1)) me%tau(2) = me%tau(1)
    me%tau(1) = me%h

    ! apply correction to column j of zn: l_j * Delta_n
    do j = 0, me%q
        me%zn(:,j) = me%zn(:,j) + me%l(j) * me%acor
    end do

    ! [not yet ported] apply the projection correction

    me%qwait = me%qwait - 1
    if ((me%qwait == 1) .and. (me%q /= me%qmax)) then
        me%zn(:,me%qmax) = me%acor
        me%saved_tq5 = me%tq(5)
        me%indx_acor = me%qmax
    end if

    end subroutine cvCompleteStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvPrepareNextStep`: handle the setting of stepsize and order for the
!  next step -- hprime and qprime. Along with hprime, it sets the ratio
!  eta = hprime/h. It also updates other state variables related to a
!  change of step size or order.

    subroutine cvPrepareNextStep(me, dsm)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: dsm

    ! if etamax = 1, defer step size or order changes
    if (me%etamax == ONE) then
        me%qwait  = max(me%qwait, 2)
        me%qprime = me%q
        me%hprime = me%h
        me%eta    = ONE
    else

        ! etaq is the ratio of new to old h at the current order
        me%etaq = ONE / ((BIAS2 * dsm)**(ONE / me%LL) + ADDON)

        ! if no order change, adjust eta and acor in cvSetEta and return
        if (me%qwait /= 0) then
            me%eta    = me%etaq
            me%qprime = me%q
            call cvSetEta(me)
        else
            ! if qwait = 0, consider an order change. etaqm1 and etaqp1 are
            ! the ratios of new to old h at orders q-1 and q+1, respectively.
            ! cvChooseEta selects the largest; cvSetEta adjusts eta and acor
            me%qwait  = 2
            me%etaqm1 = cvComputeEtaqm1(me)
            me%etaqp1 = cvComputeEtaqp1(me)
            call cvChooseEta(me)
            call cvSetEta(me)
        end if

    end if

    end subroutine cvPrepareNextStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvSetEta`: adjust the value of eta according to the various heuristic
!  limits and the optional input hmax.

    subroutine cvSetEta(me)

    class(cvode_solver),intent(inout) :: me

    ! if eta below the threshold THRESH, reject a change of step size
    if ((me%eta > me%eta_min_fx) .and. (me%eta < me%eta_max_fx)) then
        ! keep the current step size
        me%eta    = ONE
        me%hprime = me%h
    else
        if (me%eta >= me%eta_max_fx) then
            ! increase the step size, limit by etamax and hmax
            me%eta = min(me%eta, me%etamax)
            me%eta = me%eta / max(ONE, abs(me%h) * me%hmax_inv * me%eta)
        else
            ! reduce the step size, limit by eta_min and hmin
            me%eta = max(me%eta, me%eta_min)
            me%eta = max(me%eta, me%hmin / abs(me%h))
        end if

        ! set hprime
        me%hprime = me%h * me%eta
        if (me%qprime < me%q) me%nscon = 0
    end if

    end subroutine cvSetEta
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvComputeEtaqm1`: compute and return the value of etaqm1 for a possible
!  decrease in order by 1.

    function cvComputeEtaqm1(me) result(etaqm1)

    class(cvode_solver),intent(inout) :: me
    real(wp) :: etaqm1

    real(wp) :: ddn

    me%etaqm1 = ZERO
    if (me%q > 1) then
        ddn = n_vwrmsnorm(me%zn(:,me%q), me%ewt) * me%tq(1)
        me%etaqm1 = ONE / ((BIAS1 * ddn)**(ONE / me%q) + ADDON)
    end if
    etaqm1 = me%etaqm1

    end function cvComputeEtaqm1
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvComputeEtaqp1`: compute and return the value of etaqp1 for a possible
!  increase in order by 1.

    function cvComputeEtaqp1(me) result(etaqp1)

    class(cvode_solver),intent(inout) :: me
    real(wp) :: etaqp1

    real(wp) :: dup, cquot

    me%etaqp1 = ZERO
    if (me%q /= me%qmax) then
        if (me%saved_tq5 == ZERO) then
            etaqp1 = me%etaqp1
            return
        end if
        cquot = (me%tq(5) / me%saved_tq5) * (me%h / me%tau(2))**(me%LL)
        me%tempv = (-cquot)*me%zn(:,me%qmax) + me%acor
        dup = n_vwrmsnorm(me%tempv, me%ewt) * me%tq(3)
        me%etaqp1 = ONE / ((BIAS3 * dup)**(ONE / (me%LL + 1)) + ADDON)
    end if
    etaqp1 = me%etaqp1

    end function cvComputeEtaqp1
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvChooseEta`: given etaqm1, etaq, etaqp1 (the values of eta for qprime =
!  q - 1, q, or q + 1, respectively), choose the maximum eta value, set eta
!  to that value, and set qprime to the corresponding value of q. If there
!  is a tie, the preference order is to (1) keep the same order, then (2)
!  decrease the order, and finally (3) increase the order. If the maximum
!  eta value is within the threshold, keep the current step size and order.
!  For BDF, when an order increase is chosen, acor is saved in zn[qmax] for
!  the order increase in cvIncreaseBDF.

    subroutine cvChooseEta(me)

    class(cvode_solver),intent(inout) :: me

    real(wp) :: etam

    etam = max(me%etaqm1, max(me%etaq, me%etaqp1))

    if ((etam > me%eta_min_fx) .and. (etam < me%eta_max_fx)) then
        me%eta    = ONE
        me%qprime = me%q
    else
        if (etam == me%etaq) then
            me%eta    = me%etaq
            me%qprime = me%q
        else if (etam == me%etaqm1) then
            me%eta    = me%etaqm1
            me%qprime = me%q - 1
        else
            me%eta    = me%etaqp1
            me%qprime = me%q + 1

            if (me%lmm == CV_BDF) then
                ! store Delta_n in zn[qmax] to be used in order increase
                !
                ! This happens at the last step of order q before an increase
                ! to order q+1, so it represents Delta_n in the ELTE at q+1
                me%zn(:,me%qmax) = me%acor
            end if
        end if
    end if

    end subroutine cvChooseEta
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvHandleFailure`: print an error message for a failed step and return
!  the flag.

    function cvHandleFailure(me, flag) result(istate)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: flag
    integer :: istate

    character(len=*),parameter :: fn = 'CVode'

    ! set vector of absolute weighted local errors
    ! (not done here: CVODE computes it in CVodeGetErrWeights / GetEstLocalErrors)

    select case (flag)
    case (CV_ERR_FAILURE)
        call cvProcessError(me, flag, fn, 'At '//th_str(me%tn, me%h)// &
            ', the error test failed repeatedly or with |h| = hmin.')
    case (CV_CONV_FAILURE)
        call cvProcessError(me, flag, fn, 'At '//th_str(me%tn, me%h)// &
            ', the corrector convergence test failed repeatedly or with |h| = hmin.')
    case (CV_LSETUP_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the setup routine failed in an unrecoverable manner.')
    case (CV_LSOLVE_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the solve routine failed in an unrecoverable manner.')
    case (CV_RHSFUNC_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the right-hand side routine failed in an unrecoverable manner.')
    case (CV_UNREC_RHSFUNC_ERR)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the right-hand side failed in a recoverable manner, but no recovery is possible.')
    case (CV_REPTD_RHSFUNC_ERR)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ' repeated recoverable right-hand side function errors.')
    case (CV_RTFUNC_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the rootfinding routine failed in an unrecoverable manner.')
    case (CV_TOO_CLOSE)
        call cvProcessError(me, flag, fn, 'tout too close to t0 to start integration.')
    case (CV_MEM_NULL)
        call cvProcessError(me, flag, fn, 'cvode_mem = NULL illegal.')
    case (SUN_ERR_ARG_CORRUPT)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the nonlinear solver was passed a NULL input.')
    case (CV_NLS_SETUP_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the nonlinear solver setup failed unrecoverably.')
    case (CV_CONSTR_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', unable to satisfy inequality constraints.')
    case (CV_NLS_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ', the nonlinear solver failed in an unrecoverable manner.')
    case (CV_PROJ_MEM_NULL)
        call cvProcessError(me, flag, fn, 'proj_mem = NULL illegal.')
    case (CV_PROJFUNC_FAIL)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ' the projection function failed with an unrecoverable error.')
    case (CV_REPTD_PROJFUNC_ERR)
        call cvProcessError(me, flag, fn, 'At '//t_str(me%tn)// &
            ' the projection function had repeated recoverable errors.')
    case default
        ! this return should never happen
        call cvProcessError(me, CV_UNRECOGNIZED_ERR, fn, &
            'CVODE encountered an unrecognized error. Please report this to the '// &
            'SUNDIALS developers at sundials-users@llnl.gov')
        istate = CV_UNRECOGNIZED_ERR
        return
    end select

    istate = flag

    end function cvHandleFailure
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvEwtSet`: set the error weight vector ewt according to the tolerance
!  type (scalar/vector), or by calling the user function. Returns 0 if
!  successful, nonzero if a component of ewt would be <= 0.
!
!  (In C this is the default `cv_efun`, called with `e_data = cv_mem`.)

    function cvEwtSet(me, ycur, weight) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: ycur
    real(wp),dimension(:),intent(inout) :: weight
    integer :: flag

    flag = 0

    if (me%user_efun) then
        flag = me%efun_user(me, ycur, weight)
        return
    end if

    select case (me%itol)
    case (CV_SS); flag = cvEwtSetSS(me, ycur, weight)
    case (CV_SV); flag = cvEwtSetSV(me, ycur, weight)
    end select

    end function cvEwtSet
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvEwtSetSS`: ewt(i) = 1 / (reltol * |ycur(i)| + abstol)

    function cvEwtSetSS(me, ycur, weight) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: ycur
    real(wp),dimension(:),intent(inout) :: weight
    integer :: flag

    me%tempv = abs(ycur)
    me%tempv = me%reltol*me%tempv
    me%tempv = me%tempv + me%Sabstol
    if (me%atolmin0) then
        if (minval(me%tempv) <= ZERO) then
            flag = -1
            return
        end if
    end if
    weight = ONE / me%tempv

    flag = 0

    end function cvEwtSetSS
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvEwtSetSV`: ewt(i) = 1 / (reltol * |ycur(i)| + abstol(i))

    function cvEwtSetSV(me, ycur, weight) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: ycur
    real(wp),dimension(:),intent(inout) :: weight
    integer :: flag

    me%tempv = abs(ycur)
    me%tempv = me%reltol*me%tempv + me%Vabstol
    if (me%atolmin0) then
        if (minval(me%tempv) <= ZERO) then
            flag = -1
            return
        end if
    end if
    weight = ONE / me%tempv

    flag = 0

    end function cvEwtSetSV
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvProcessError`: report an error or warning (replaces the SUNContext
!  error handler and logger). Messages go to `errunit` unless printing is
!  turned off with [[CVodeSetErrorOutput]].

    subroutine cvProcessError(me, error_code, func, msg)

    class(cvode_solver),intent(in) :: me
    integer,intent(in) :: error_code
    character(len=*),intent(in) :: func
    character(len=*),intent(in) :: msg

    if (.not. me%print_errors) return

    if (error_code == CV_WARNING) then
        write(me%errunit, '(a)') '[WARNING]['//func//'] '//msg
    else
        write(me%errunit, '(a,i0,a)') '[ERROR]['//func//'] ('//'flag ', error_code, ') '//msg
    end if

    end subroutine cvProcessError
!*****************************************************************************************


!=========================================================================================
!  Nonlinear solver interface (cvode_nls.c)
!=========================================================================================

!*****************************************************************************************
!>
!  `CVodeSetNonlinearSolver`: attach a nonlinear solver. The solver object is
!  copied into the CVODE object.

    subroutine CVodeSetNonlinearSolver(me, NLS, flag)

    class(cvode_solver),intent(inout) :: me
    class(sun_nonlinear_solver),intent(in) :: NLS
    integer,intent(out) :: flag

    integer :: retval, itype

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeSetNonlinearSolver', &
                            'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    itype = NLS%get_type()
    if (itype /= SUNNONLINEARSOLVER_ROOTFIND .and. &
        itype /= SUNNONLINEARSOLVER_FIXEDPOINT) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetNonlinearSolver', &
                            'Invalid nonlinear solver type')
        flag = CV_ILL_INPUT
        return
    end if

    ! free any existing nonlinear solver and attach the new one
    if (allocated(me%NLS)) deallocate(me%NLS)
    allocate(me%NLS, source=NLS)

    ! (the system function is selected in cvNlsSys from the solver type)

    ! set the nonlinear solver's maximum number of iterations
    retval = me%NLS%set_max_iters(NLS_MAXCOR)
    if (retval /= CV_SUCCESS) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetNonlinearSolver', &
                            'Setting maximum number of nonlinear iterations failed')
        flag = CV_ILL_INPUT
        return
    end if

    ! reset the acnrmcur flag
    me%acnrmcur = .false.

    flag = CV_SUCCESS

    end subroutine CVodeSetNonlinearSolver
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsInit`: record whether the linear solver functions exist, then
!  initialize the nonlinear solver.

    function cvNlsInit(me) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer :: flag

    integer :: retval

    if (.not. allocated(me%NLS)) then
        ! [not yet ported] C attaches a default Newton solver in CVodeInit
        call cvProcessError(me, CV_NLS_INIT_FAIL, 'cvNlsInit', &
            'No nonlinear solver attached (the default Newton solver is not ported yet).')
        flag = CV_NLS_INIT_FAIL
        return
    end if

    ! set the linear solver setup / solve wrappers (if they exist)
    me%NLS%has_lsetup = me%lsetup_exists
    me%NLS%has_lsolve = me%lsetup_exists

    ! initialize nonlinear solver
    retval = me%NLS%initialize()
    if (retval /= CV_SUCCESS) then
        call cvProcessError(me, CV_NLS_INIT_FAIL, 'cvNlsInit', &
                            "The nonlinear solver's init routine failed.")
        flag = CV_NLS_INIT_FAIL
        return
    end if

    flag = CV_SUCCESS

    end function cvNlsInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  The nonlinear system function: `cvNlsResidual` for root-finding solvers
!  (Newton), `cvNlsFPFunction` for fixed-point solvers. (C attaches one of
!  them in `CVodeSetNonlinearSolver` according to the solver type.)

    function cvNlsSys(me, ycor, res) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: ycor
    real(wp),dimension(:),intent(inout) :: res
    integer :: flag

    if (me%NLS%get_type() == SUNNONLINEARSOLVER_FIXEDPOINT) then
        flag = cvNlsFPFunction(me, ycor, res)
    else
        flag = cvNlsResidual(me, ycor, res)
    end if

    end function cvNlsSys
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsResidual`: `res = rl1*zn[1] + ycor - gamma*f(t, zn[0] + ycor)`

    function cvNlsResidual(me, ycor, res) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: ycor
    real(wp),dimension(:),intent(inout) :: res
    integer :: flag

    integer :: retval

    ! update the state based on the current correction
    me%y = me%zn(:,0) + ycor

    ! evaluate the rhs function
    retval = me%f(me, me%tn, me%y, me%ftemp)
    me%nfe = me%nfe + 1
    if (retval < 0) then
        flag = CV_RHSFUNC_FAIL
        return
    end if
    if (retval > 0) then
        flag = RHSFUNC_RECVR
        return
    end if

    ! compute the resiudal
    res = me%rl1*me%zn(:,1) + ycor
    res = (-me%gamma)*me%ftemp + res

    flag = CV_SUCCESS

    end function cvNlsResidual
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsFPFunction`: `res = rl1*(h*f(t, zn[0] + ycor) - zn[1])`

    function cvNlsFPFunction(me, ycor, res) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: ycor
    real(wp),dimension(:),intent(inout) :: res
    integer :: flag

    integer :: retval

    ! update the state based on the current correction
    me%y = me%zn(:,0) + ycor

    ! evaluate the rhs function
    retval = me%f(me, me%tn, me%y, res)
    me%nfe = me%nfe + 1
    if (retval < 0) then
        flag = CV_RHSFUNC_FAIL
        return
    end if
    if (retval > 0) then
        flag = RHSFUNC_RECVR
        return
    end if

    res = me%h*res - me%zn(:,1)
    res = me%rl1*res

    flag = CV_SUCCESS

    end function cvNlsFPFunction
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsLSetup`  [not yet ported: linear solvers]

    function cvNlsLSetup(me, jbad, jcur) result(flag)

    class(cvode_solver),intent(inout) :: me
    logical,intent(in) :: jbad
    logical,intent(out) :: jcur
    integer :: flag

    ! if the nonlinear solver marked the Jacobian as bad update convfail
    if (jbad) me%convfail = CV_FAIL_BAD_J

    jcur = me%jcur
    flag = CV_LSETUP_FAIL   ! no linear solver yet

    end function cvNlsLSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsLSolve`  [not yet ported: linear solvers]

    function cvNlsLSolve(me, b) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(inout) :: b
    integer :: flag

    if (size(b) < 0) flag = 0   ! (unused)
    flag = CV_LSOLVE_FAIL   ! no linear solver yet

    end function cvNlsLSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsConvTest`: the CVODE nonlinear convergence test, based on the
!  estimated convergence rate `crate` and the norm of the update `del`.

    function cvNlsConvTest(me, nls, ycor, del, tol, ewt) result(flag)

    class(cvode_solver),intent(inout) :: me
    class(sun_nonlinear_solver),intent(in) :: nls
    real(wp),dimension(:),intent(in) :: ycor
    real(wp),dimension(:),intent(in) :: del
    real(wp),intent(in) :: tol
    real(wp),dimension(:),intent(in) :: ewt
    integer :: flag

    integer :: m
    real(wp) :: dcon

    ! compute the norm of the correction
    me%delnrm = n_vwrmsnorm(del, ewt)   ! cvNlsNorm

    ! get the current nonlinear solver iteration count
    m = nls%get_cur_iter()

    ! test for convergence, first calculating the convergence rate constant,
    ! then the norm of the correction
    if (m > 0) me%crate = max(CRDOWN * me%crate, me%delnrm / me%delp)
    dcon = me%delnrm * min(ONE, me%crate) / tol

    if (dcon <= ONE) then
        if (m == 0) then
            me%acnrm = me%delnrm
        else
            me%acnrm = n_vwrmsnorm(ycor, ewt)
        end if
        me%acnrmcur = .true.
        flag = CV_SUCCESS  ! nonlinear system was solved successfully
        return
    end if

    ! check if the iteration seems to be diverging
    if ((m >= 1) .and. (me%delnrm > RDIV * me%delp)) then
        flag = SUN_NLS_CONV_RECVR
        return
    end if

    ! save norm of correction and loop again
    me%delp = me%delnrm

    ! not yet converged
    flag = SUN_NLS_CONTINUE

    end function cvNlsConvTest
!*****************************************************************************************

!=========================================================================================
!  Optional inputs (cvode_io.c)
!=========================================================================================

!*****************************************************************************************
!>
!  Turn error and warning messages on/off, and select the output unit
!  (replaces the SUNContext error handler / logger configuration).

    subroutine CVodeSetErrorOutput(me, print_errors, unit)

    class(cvode_solver),intent(inout) :: me
    logical,intent(in) :: print_errors
    integer,intent(in),optional :: unit

    me%print_errors = print_errors
    if (present(unit)) me%errunit = unit

    end subroutine CVodeSetErrorOutput
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxOrd`: maximum method order (may only be decreased).

    subroutine CVodeSetMaxOrd(me, maxord, flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: maxord
    integer,intent(out) :: flag

    integer :: qmax_alloc

    if (maxord <= 0) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetMaxOrd', 'maxord <= 0 illegal.')
        flag = CV_ILL_INPUT
        return
    end if

    ! cannot increase maximum order beyond the value that
    ! was used when allocating memory
    qmax_alloc = me%qmax_alloc

    if (maxord > qmax_alloc) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetMaxOrd', &
                            'Illegal attempt to increase maximum method order.')
        flag = CV_ILL_INPUT
        return
    end if

    me%qmax = maxord

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxOrd
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxNumSteps`: maximum number of internal steps per call
!  (0 means the default; a negative value disables the test).

    subroutine CVodeSetMaxNumSteps(me, mxsteps, flag)

    class(cvode_solver),intent(inout) :: me
    integer(li),intent(in) :: mxsteps
    integer,intent(out) :: flag

    ! passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test.
    if (mxsteps == 0) then
        me%mxstep = MXSTEP_DEFAULT
    else
        me%mxstep = mxsteps
    end if

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxNumSteps
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxHnilWarns`: maximum number of "t + h = t" warnings.

    subroutine CVodeSetMaxHnilWarns(me, mxhnil, flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: mxhnil
    integer,intent(out) :: flag

    me%mxhnil = mxhnil

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxHnilWarns
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetInitStep`: initial step size (0 means estimate it).

    subroutine CVodeSetInitStep(me, hin, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: hin
    integer,intent(out) :: flag

    me%hin = hin

    ! if an initial step size was set, take it into account at the first step
    flag = CV_SUCCESS

    end subroutine CVodeSetInitStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMinStep`: minimum absolute step size (0 means no limit).

    subroutine CVodeSetMinStep(me, hmin, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: hmin
    integer,intent(out) :: flag

    if (hmin < ZERO) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetMinStep', 'hmin < 0 illegal.')
        flag = CV_ILL_INPUT
        return
    end if

    ! passing 0 sets hmin = zero
    if (hmin == ZERO) then
        me%hmin = HMIN_DEFAULT
        flag = CV_SUCCESS
        return
    end if

    if (hmin * me%hmax_inv > ONE) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetMinStep', &
                            'Inconsistent step size limits: hmin > hmax.')
        flag = CV_ILL_INPUT
        return
    end if

    me%hmin = hmin

    flag = CV_SUCCESS

    end subroutine CVodeSetMinStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxStep`: maximum absolute step size (0 means no limit).

    subroutine CVodeSetMaxStep(me, hmax, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: hmax
    integer,intent(out) :: flag

    real(wp) :: hmax_inv

    if (hmax < ZERO) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetMaxStep', 'hmax < 0 illegal.')
        flag = CV_ILL_INPUT
        return
    end if

    ! passing 0 sets hmax = infinity
    if (hmax == ZERO) then
        me%hmax_inv = HMAX_INV_DEFAULT
        flag = CV_SUCCESS
        return
    end if

    hmax_inv = ONE / hmax
    if (hmax_inv * me%hmin > ONE) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetMaxStep', &
                            'Inconsistent step size limits: hmin > hmax.')
        flag = CV_ILL_INPUT
        return
    end if

    me%hmax_inv = hmax_inv

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetStopTime`: value of the independent variable past which the
!  solution is not to proceed.

    subroutine CVodeSetStopTime(me, tstop, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: tstop
    integer,intent(out) :: flag

    ! If CVode was called at least once, test if tstop is legal
    ! (i.e. if it was not already passed).
    ! If CVodeSetStopTime is called before the first call to CVode,
    ! tstop will be checked in CVode.
    if (me%nst > 0) then
        if ((tstop - me%tn) * me%h < ZERO) then
            call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetStopTime', &
                'The value tstop = '//g_str(tstop)//' is behind current '// &
                t_str(me%tn)//' in the direction of integration.')
            flag = CV_ILL_INPUT
            return
        end if
    end if

    me%tstop    = tstop
    me%tstopset = .true.

    flag = CV_SUCCESS

    end subroutine CVodeSetStopTime
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetInterpolateStopTime`: return the interpolated solution at tstop
!  (instead of the internal solution near tstop).

    subroutine CVodeSetInterpolateStopTime(me, interp, flag)

    class(cvode_solver),intent(inout) :: me
    logical,intent(in) :: interp
    integer,intent(out) :: flag

    me%tstopinterp = interp

    flag = CV_SUCCESS

    end subroutine CVodeSetInterpolateStopTime
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeClearStopTime`: disable the stop time.

    subroutine CVodeClearStopTime(me, flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(out) :: flag

    me%tstopset = .false.

    flag = CV_SUCCESS

    end subroutine CVodeClearStopTime
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxErrTestFails`: maximum number of error test failures per step
!  (<= 0 restores the default).

    subroutine CVodeSetMaxErrTestFails(me, maxnef, flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: maxnef
    integer,intent(out) :: flag

    if (maxnef <= 0) then
        me%maxnef = MXNEF
    else
        me%maxnef = maxnef
    end if

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxErrTestFails
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxConvFails`: maximum number of nonlinear convergence failures
!  per step (<= 0 restores the default).

    subroutine CVodeSetMaxConvFails(me, maxncf, flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: maxncf
    integer,intent(out) :: flag

    if (maxncf <= 0) then
        me%maxncf = MXNCF
    else
        me%maxncf = maxncf
    end if

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxConvFails
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetNonlinConvCoef`: coefficient in the nonlinear convergence test
!  (<= 0 restores the default).

    subroutine CVodeSetNonlinConvCoef(me, nlscoef, flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: nlscoef
    integer,intent(out) :: flag

    if (nlscoef <= ZERO) then
        me%nlscoef = CORTES
    else
        me%nlscoef = nlscoef
    end if

    flag = CV_SUCCESS

    end subroutine CVodeSetNonlinConvCoef
!*****************************************************************************************

!=========================================================================================
!  Optional outputs (cvode_io.c)
!=========================================================================================

    pure integer(li) function CVodeGetNumSteps(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%nst
    end function CVodeGetNumSteps

    pure integer(li) function CVodeGetNumRhsEvals(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%nfe
    end function CVodeGetNumRhsEvals

    pure integer(li) function CVodeGetNumLinSolvSetups(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%nsetups
    end function CVodeGetNumLinSolvSetups

    pure integer(li) function CVodeGetNumErrTestFails(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%netf
    end function CVodeGetNumErrTestFails

    pure integer function CVodeGetLastOrder(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%qu
    end function CVodeGetLastOrder

    pure integer function CVodeGetCurrentOrder(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%next_q
    end function CVodeGetCurrentOrder

    pure real(wp) function CVodeGetActualInitStep(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%h0u
    end function CVodeGetActualInitStep

    pure real(wp) function CVodeGetLastStep(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%hu
    end function CVodeGetLastStep

    pure real(wp) function CVodeGetCurrentStep(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%next_h
    end function CVodeGetCurrentStep

    pure real(wp) function CVodeGetCurrentTime(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%tn
    end function CVodeGetCurrentTime

    pure real(wp) function CVodeGetTolScaleFactor(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%tolsf
    end function CVodeGetTolScaleFactor

    pure integer(li) function CVodeGetNumNonlinSolvIters(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%nni
    end function CVodeGetNumNonlinSolvIters

    pure integer(li) function CVodeGetNumNonlinSolvConvFails(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%nnf
    end function CVodeGetNumNonlinSolvConvFails

    pure integer(li) function CVodeGetNumStepSolveFails(me) result(v)
    class(cvode_solver),intent(in) :: me
    v = me%ncfn
    end function CVodeGetNumStepSolveFails

!=========================================================================================
!  Message formatting helpers (SUN_FORMAT_G is "%.15g")
!=========================================================================================

    function g_str(x) result(s)
    !! A real formatted with 15 significant digits.
    real(wp),intent(in) :: x
    character(len=:),allocatable :: s
    character(len=32) :: buf
    write(buf, '(es22.15)') x
    s = trim(adjustl(buf))
    end function g_str

    function t_str(t) result(s)
    !! `MSG_TIME`: "t = <t>"
    real(wp),intent(in) :: t
    character(len=:),allocatable :: s
    s = 't = '//g_str(t)
    end function t_str

    function th_str(t, h) result(s)
    !! `MSG_TIME_H`: "t = <t> and h = <h>"
    real(wp),intent(in) :: t
    real(wp),intent(in) :: h
    character(len=:),allocatable :: s
    s = 't = '//g_str(t)//' and h = '//g_str(h)
    end function th_str

!*****************************************************************************************
    end module fvode_core
!*****************************************************************************************

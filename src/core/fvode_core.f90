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
!  The linear solver interface (`cvode_ls.c`) is in the submodule
!  `fvode_core_ls`.
!
!### Not yet ported
!
!  Projection and history resizing. Root finding is in the submodule `fvode_core_root`, stability
!  limit detection in `fvode_core_stald`, the band preconditioner in
!  `fvode_core_bandpre`.
!  The places where the C code handles them are marked `! [not yet ported]`.

    module fvode_core

    use iso_fortran_env, only: error_unit
    use fvode_kinds, only: wp, ip, li
    use fvode_parameters
    use fvode_internal_parameters
    use fvode_nvector, only: n_vwrmsnorm, n_vinvtest, n_vconstrmask, n_vminquotient
    use fvode_nonlinsol, only: sun_nonlinear_solver, nls_problem
    use fvode_nonlinsol_newton, only: sun_nonlinsol_newton
    use fvode_matrix, only: sun_matrix, dense_matrix, band_matrix, SUNMATRIX_DENSE, SUNMATRIX_BAND
    use fvode_linsol, only: sun_linear_solver
    use fvode_linsol_direct, only: sun_linsol_band

    implicit none

    private

    ! which linear solver interface is attached (C: cv_linit, cv_lsetup, ... pointers)
    integer,parameter :: LS_NONE    = 0  !! none
    integer,parameter :: LS_GENERIC = 1  !! SUNLinearSolver interface (cvode_ls.c)
    integer,parameter :: LS_DIAG    = 2  !! diagonal approximate Jacobian (cvode_diag.c)

    ! which preconditioner is attached (C: cvls_mem->pset, psolve, P_data)
    integer,parameter :: PRE_NONE    = 0  !! none
    integer,parameter :: PRE_USER    = 1  !! user psetup/psolve functions
    integer,parameter :: PRE_BANDPRE = 2  !! band preconditioner (cvode_bandpre.c)

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
        logical :: tstopset = .false.  !! is a stop time set?
        logical :: tstopinterp = .false.  !! interpolate the solution at tstop?
        real(wp) :: tstop = ZERO  !! stop time

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
        real(wp) :: eta_max_fx = ETA_MAX_FX_DEFAULT  !! upper bound of the 'keep h' range
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

        !---- linear solver (C: cv_linit, cv_lsetup, ... and cv_lmem)
        integer :: ls_kind = LS_NONE        !! which interface: LS_NONE, LS_GENERIC, LS_DIAG
        logical :: lsetup_exists = .false.  !! C: cv_lsetup != NULL
        logical :: lsolve_exists = .false.  !! C: cv_lsolve != NULL
        integer(li) :: msbp = MSBP_DEFAULT  !! max steps between lsetup calls
        real(wp) :: dgmax_lsetup = DGMAX_LSETUP_DEFAULT !! gamma ratio threshold for lsetup

        !---- linear solver interface memory (CVLsMemRec, cvode_ls_impl.h)
        logical :: ls_iterative = .false.    !! is the solver iterative?
        logical :: ls_matrixbased = .false.  !! is a matrix structure used?
        logical :: jacDQ = .false.           !! use the internal DQ Jacobian
        procedure(cv_jac_fn),pointer,nopass :: jac_user => null()  !! user Jacobian
        logical :: jbad = .true.             !! heuristic suggestion for pset
        real(wp) :: dgmax_jbad = CVLS_DGMAX  !! |gamma/gammap-1| < dgmax_jbad: J is bad
        logical :: scalesol = .false.        !! scale the solution for a change in gamma
        real(wp) :: eplifac = CVLS_EPLIN     !! nonlinear -> linear tolerance factor
        real(wp) :: nrmfac = ONE             !! integrator -> LS norm conversion factor
        class(sun_linear_solver),allocatable,public :: LS  !! the attached linear solver (a copy; public so that solver-specific settings can be changed after attaching)
        class(sun_matrix),allocatable :: A          !! A = I - gamma * df/dy
        class(sun_matrix),allocatable :: savedJ     !! saved Jacobian
        real(wp),dimension(:),allocatable :: ls_ytemp !! temporary (jtimes, psolve)
        real(wp),dimension(:),allocatable :: ls_x     !! solution of the linear solve
        integer(li) :: msbj = CVLS_MSBJ      !! max steps between Jacobian evaluations
        integer(li) :: nje = 0               !! number of Jacobian evaluations
        integer(li) :: nfeDQ = 0             !! f evaluations for DQ Jacobian / J*v
        integer(li) :: nstlj = 0             !! nst at the last Jacobian evaluation
        integer(li) :: npe = 0               !! preconditioner setups
        integer(li) :: nli = 0               !! linear iterations
        integer(li) :: nps = 0               !! preconditioner solves
        integer(li) :: ncfl = 0              !! linear convergence failures
        integer(li) :: njtsetup = 0          !! J*v setup calls
        integer(li) :: njtimes = 0           !! J*v calls
        real(wp) :: tnlj = ZERO              !! tn at the last Jacobian evaluation
        integer :: ls_last_flag = CVLS_SUCCESS !! last flag from a linear solver function
        procedure(cv_psetup_fn),pointer,nopass :: pset_user => null()  !! user preconditioner setup function
        procedure(cv_psolve_fn),pointer,nopass :: psolve_user => null()  !! user preconditioner solve function
        integer :: pre_kind = PRE_NONE  !! which preconditioner is attached: PRE_NONE, PRE_USER, PRE_BANDPRE
        logical :: jtimesDQ = .true.  !! use the internal DQ Jacobian-times-vector product
        procedure(cv_jtsetup_fn),pointer,nopass :: jtsetup_user => null()  !! user J*v setup function
        procedure(cv_jtimes_fn),pointer,nopass :: jtimes_user => null()  !! user J*v function
        procedure(cv_rhs_fn),pointer,nopass :: jt_f => null()  !! f used by the DQ J*v (null: the ODE f)

        !---- band preconditioner (CVBandPrecDataRec, cvode_bandpre_impl.h)
        integer(ip) :: bp_N = 0   !! problem size
        integer(ip) :: bp_mu = 0  !! upper half-bandwidth of the preconditioner
        integer(ip) :: bp_ml = 0  !! lower half-bandwidth of the preconditioner
        type(band_matrix) :: bp_savedJ  !! saved banded Jacobian approximation
        type(band_matrix) :: bp_savedP  !! P = I - gamma*J (factored)
        type(sun_linsol_band) :: bp_LS  !! band LU solver for P
        real(wp),dimension(:),allocatable :: bp_tmp1  !! work vector (f at the perturbed y)
        real(wp),dimension(:),allocatable :: bp_tmp2  !! work vector (perturbed y)
        integer(li) :: bp_nfeBP = 0  !! number of f evaluations for the band Jacobian

        !---- diagonal linear solver (CVDiagMemRec, cvode_diag_impl.h)
        real(wp) :: di_gammasv = ZERO              !! gamma at the last call to setup
        real(wp),dimension(:),allocatable :: di_M  !! M = (I - gamma J)^{-1}
        real(wp),dimension(:),allocatable :: di_bit, di_bitcomp  !! temporary storage
        integer(li) :: di_nfeDI = 0                !! f evaluations for the diagonal J
        integer :: di_last_flag = CVDIAG_SUCCESS   !! last error flag

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

        !---- stability limit detection (STALD)
        logical :: sldeton = .false.  !! is stability limit detection on?
        integer :: nscon = 0         !! counter for STALD method
        integer(li) :: nor = 0       !! number of order reductions
        real(wp),dimension(0:5,0:3) :: ssdat = ZERO  !! scaled data array for STALD (C: ssdat[6][4], used as [1..5][1..3])

        !---- root finding
        procedure(cv_root_fn),pointer,nopass :: gfun => null()  !! function g for the roots sought
        integer :: nrtfn = 0                         !! number of components of g
        integer,dimension(:),allocatable :: iroots   !! root information: +-1 for a root in g(i) (direction), 0 otherwise
        integer,dimension(:),allocatable :: rootdir  !! direction of zero-crossing to look for: +1, -1, or 0 (both)
        real(wp) :: tlo = ZERO    !! nearest endpoint of the interval in the root search
        real(wp) :: thi = ZERO    !! farthest endpoint of the interval in the root search
        real(wp) :: trout = ZERO  !! t value returned by the root search
        real(wp),dimension(:),allocatable :: glo    !! saved g values at t = tlo
        real(wp),dimension(:),allocatable :: ghi    !! saved g values at t = thi
        real(wp),dimension(:),allocatable :: grout  !! g values at t = trout
        real(wp) :: ttol = ZERO   !! tolerance on the root location trout
        integer :: irfnd = 0      !! flag: did the last step find a root? (1: yes)
        integer(li) :: nge = 0    !! number of root function evaluations
        logical,dimension(:),allocatable :: gactive  !! which root functions are active
        integer :: mxgnull = 1    !! number of warnings about root functions that are identically 0

        !---- inequality constraints
        real(wp),dimension(:),allocatable :: constraints  !! constraint flags (allocated: constraints are on)
        integer(li) :: constraint_corrections = 0  !! total constraint corrections
        integer(li) :: constraint_fails = 0        !! total constraint failures
        integer :: max_constraint_fails = MAX_CONSTRAINT_FAILS  !! max constraint failures per step

        !---- error reporting (replaces the SUNContext error handler)
        logical :: print_errors = .true.     !! print error and warning messages
        integer :: errunit = error_unit      !! unit for messages

    contains

        ! user-callable
        procedure,public :: create              => CVodeCreate  !! `CVodeCreate`: set the method (Adams/BDF) and the default options
        procedure,public :: init                => CVodeInit  !! `CVodeInit`: allocate memory and set the problem
        procedure,public :: reinit              => CVodeReInit  !! `CVodeReInit`: restart with a new initial condition
        procedure,public :: ss_tolerances       => CVodeSStolerances  !! `CVodeSStolerances`: scalar relative and absolute tolerances
        procedure,public :: sv_tolerances       => CVodeSVtolerances  !! `CVodeSVtolerances`: scalar relative and vector absolute tolerances
        procedure,public :: wf_tolerances       => CVodeWFtolerances  !! `CVodeWFtolerances`: user error weight function
        procedure,public :: set_nonlinear_solver => CVodeSetNonlinearSolver  !! `CVodeSetNonlinearSolver`: attach a nonlinear solver
        procedure,public :: solve               => CVode  !! `CVode`: integrate to tout, or take one step
        procedure,public :: get_dky             => CVodeGetDky  !! `CVodeGetDky`: interpolate the solution (or a derivative)
        procedure,public :: destroy             => CVodeFree  !! `CVodeFree`: free all memory

        ! optional inputs (cvode_io.c)
        procedure,public :: set_error_output    => CVodeSetErrorOutput  !! `CVodeSetErrorOutput`: control error and warning messages
        procedure,public :: set_max_ord         => CVodeSetMaxOrd  !! `CVodeSetMaxOrd`: maximum method order
        procedure,public :: set_max_num_steps   => CVodeSetMaxNumSteps  !! `CVodeSetMaxNumSteps`: maximum steps per call
        procedure,public :: set_max_hnil_warns  => CVodeSetMaxHnilWarns  !! `CVodeSetMaxHnilWarns`: maximum 't + h = t' warnings
        procedure,public :: set_init_step       => CVodeSetInitStep  !! `CVodeSetInitStep`: initial step size
        procedure,public :: set_min_step        => CVodeSetMinStep  !! `CVodeSetMinStep`: minimum step size
        procedure,public :: set_max_step        => CVodeSetMaxStep  !! `CVodeSetMaxStep`: maximum step size
        procedure,public :: set_stop_time       => CVodeSetStopTime  !! `CVodeSetStopTime`: stop time
        procedure,public :: set_interpolate_stop_time => CVodeSetInterpolateStopTime  !! `CVodeSetInterpolateStopTime`: interpolate at the stop time
        procedure,public :: clear_stop_time     => CVodeClearStopTime  !! `CVodeClearStopTime`: disable the stop time
        procedure,public :: set_max_err_test_fails => CVodeSetMaxErrTestFails  !! `CVodeSetMaxErrTestFails`: maximum error test failures per step
        procedure,public :: set_max_conv_fails  => CVodeSetMaxConvFails  !! `CVodeSetMaxConvFails`: maximum nonlinear convergence failures per step
        procedure,public :: set_nonlin_conv_coef => CVodeSetNonlinConvCoef  !! `CVodeSetNonlinConvCoef`: nonlinear convergence coefficient
        procedure,public :: set_max_nonlin_iters => CVodeSetMaxNonlinIters  !! `CVodeSetMaxNonlinIters`: maximum nonlinear iterations per step
        procedure,public :: set_delta_gamma_max_lsetup => CVodeSetDeltaGammaMaxLSetup  !! `CVodeSetDeltaGammaMaxLSetup`: gamma ratio threshold for the linear solver setup
        procedure,public :: set_lsetup_frequency => CVodeSetLSetupFrequency  !! `CVodeSetLSetupFrequency`: maximum steps between linear solver setups
        procedure,public :: set_eta_fixed_step_bounds => CVodeSetEtaFixedStepBounds  !! `CVodeSetEtaFixedStepBounds`: range of eta in which h is kept
        procedure,public :: set_eta_max_first_step => CVodeSetEtaMaxFirstStep  !! `CVodeSetEtaMaxFirstStep`: maximum eta after the first step
        procedure,public :: set_eta_max_early_step => CVodeSetEtaMaxEarlyStep  !! `CVodeSetEtaMaxEarlyStep`: maximum eta on early steps
        procedure,public :: set_num_steps_eta_max_early_step => CVodeSetNumStepsEtaMaxEarlyStep  !! `CVodeSetNumStepsEtaMaxEarlyStep`: number of early steps
        procedure,public :: set_eta_max         => CVodeSetEtaMax  !! `CVodeSetEtaMax`: maximum eta on a general step
        procedure,public :: set_eta_min         => CVodeSetEtaMin  !! `CVodeSetEtaMin`: minimum eta on a general step
        procedure,public :: set_eta_min_err_fail => CVodeSetEtaMinErrFail  !! `CVodeSetEtaMinErrFail`: minimum eta after an error test failure
        procedure,public :: set_eta_max_err_fail => CVodeSetEtaMaxErrFail  !! `CVodeSetEtaMaxErrFail`: maximum eta after multiple error test failures
        procedure,public :: set_num_fails_eta_max_err_fail => CVodeSetNumFailsEtaMaxErrFail  !! `CVodeSetNumFailsEtaMaxErrFail`: error test failures before eta_max_ef applies
        procedure,public :: set_eta_conv_fail   => CVodeSetEtaConvFail  !! `CVodeSetEtaConvFail`: eta after a nonlinear convergence failure

        ! optional outputs (cvode_io.c)
        procedure,public :: get_num_steps       => CVodeGetNumSteps  !! `CVodeGetNumSteps`: number of steps
        procedure,public :: get_num_rhs_evals   => CVodeGetNumRhsEvals  !! `CVodeGetNumRhsEvals`: number of f evaluations
        procedure,public :: get_num_lin_solv_setups => CVodeGetNumLinSolvSetups  !! `CVodeGetNumLinSolvSetups`: number of linear solver setups
        procedure,public :: get_num_err_test_fails => CVodeGetNumErrTestFails  !! `CVodeGetNumErrTestFails`: number of error test failures
        procedure,public :: get_last_order      => CVodeGetLastOrder  !! `CVodeGetLastOrder`: order used on the last step
        procedure,public :: get_current_order   => CVodeGetCurrentOrder  !! `CVodeGetCurrentOrder`: order to be used on the next step
        procedure,public :: get_actual_init_step => CVodeGetActualInitStep  !! `CVodeGetActualInitStep`: initial step size actually used
        procedure,public :: get_last_step       => CVodeGetLastStep  !! `CVodeGetLastStep`: last step size
        procedure,public :: get_current_step    => CVodeGetCurrentStep  !! `CVodeGetCurrentStep`: step size to be used next
        procedure,public :: get_current_time    => CVodeGetCurrentTime  !! `CVodeGetCurrentTime`: current internal time
        procedure,public :: get_tol_scale_factor => CVodeGetTolScaleFactor  !! `CVodeGetTolScaleFactor`: tolerance scale factor
        procedure,public :: get_num_nonlin_solv_iters => CVodeGetNumNonlinSolvIters  !! `CVodeGetNumNonlinSolvIters`: number of nonlinear iterations
        procedure,public :: get_num_nonlin_solv_conv_fails => CVodeGetNumNonlinSolvConvFails  !! `CVodeGetNumNonlinSolvConvFails`: number of nonlinear convergence failures
        procedure,public :: get_num_step_solve_fails => CVodeGetNumStepSolveFails  !! `CVodeGetNumStepSolveFails`: number of failed steps due to nonlinear solver failures
        procedure,public :: get_current_gamma   => CVodeGetCurrentGamma  !! `CVodeGetCurrentGamma`: current gamma
        procedure,public :: get_num_stab_lim_order_reds => CVodeGetNumStabLimOrderReds  !! `CVodeGetNumStabLimOrderReds`: order reductions by STALD
        procedure,public :: get_num_g_evals     => CVodeGetNumGEvals  !! `CVodeGetNumGEvals`: number of root function evaluations
        procedure,public :: get_current_state   => CVodeGetCurrentState  !! `CVodeGetCurrentState`: current internal state
        procedure,public :: get_err_weights     => CVodeGetErrWeights  !! `CVodeGetErrWeights`: error weights
        procedure,public :: get_est_local_errors => CVodeGetEstLocalErrors  !! `CVodeGetEstLocalErrors`: estimated local errors
        procedure,public :: get_integrator_stats => CVodeGetIntegratorStats  !! `CVodeGetIntegratorStats`: main integrator statistics
        procedure,public :: get_nonlin_solv_stats => CVodeGetNonlinSolvStats  !! `CVodeGetNonlinSolvStats`: nonlinear solver statistics
        procedure,public :: print_all_stats     => CVodePrintAllStats  !! `CVodePrintAllStats`: print all statistics

        ! linear solver interface (cvode_ls.c, submodule fvode_core_ls)
        procedure,public :: set_linear_solver   => CVodeSetLinearSolver  !! `CVodeSetLinearSolver`: attach a linear solver and matrix
        procedure,public :: set_jac_fn          => CVodeSetJacFn  !! `CVodeSetJacFn`: user Jacobian function
        procedure,public :: set_jac_eval_frequency => CVodeSetJacEvalFrequency  !! `CVodeSetJacEvalFrequency`: maximum steps between Jacobian evaluations
        procedure,public :: set_linear_solution_scaling => CVodeSetLinearSolutionScaling  !! `CVodeSetLinearSolutionScaling`: BDF solution scaling on/off
        procedure,public :: set_delta_gamma_max_bad_jac => CVodeSetDeltaGammaMaxBadJac  !! `CVodeSetDeltaGammaMaxBadJac`: gamma ratio threshold for a bad Jacobian
        procedure,public :: get_num_jac_evals   => CVodeGetNumJacEvals  !! `CVodeGetNumJacEvals`: number of Jacobian evaluations
        procedure,public :: get_num_lin_rhs_evals => CVodeGetNumLinRhsEvals  !! `CVodeGetNumLinRhsEvals`: f evaluations in the linear solver interface
        procedure,public :: get_num_lin_iters   => CVodeGetNumLinIters  !! `CVodeGetNumLinIters`: number of linear iterations
        procedure,public :: get_num_lin_conv_fails => CVodeGetNumLinConvFails  !! `CVodeGetNumLinConvFails`: number of linear convergence failures
        procedure,public :: get_num_prec_evals  => CVodeGetNumPrecEvals  !! `CVodeGetNumPrecEvals`: number of preconditioner setups
        procedure,public :: get_num_prec_solves => CVodeGetNumPrecSolves  !! `CVodeGetNumPrecSolves`: number of preconditioner solves
        procedure,public :: get_last_lin_flag   => CVodeGetLastLinFlag  !! `CVodeGetLastLinFlag`: last linear solver flag
        procedure,public :: get_jac_time        => CVodeGetJacTime  !! `CVodeGetJacTime`: time of the last Jacobian evaluation
        procedure,public :: get_jac_num_steps   => CVodeGetJacNumSteps  !! `CVodeGetJacNumSteps`: step number of the last Jacobian evaluation
        procedure,public :: set_preconditioner  => CVodeSetPreconditioner  !! `CVodeSetPreconditioner`: user preconditioner (Krylov solvers)
        procedure,public :: set_jac_times       => CVodeSetJacTimes  !! `CVodeSetJacTimes`: user Jacobian-times-vector product
        procedure,public :: set_jac_times_rhs_fn => CVodeSetJacTimesRhsFn  !! `CVodeSetJacTimesRhsFn`: f used by the DQ J*v
        procedure,public :: set_eps_lin         => CVodeSetEpsLin  !! `CVodeSetEpsLin`: linear tolerance factor
        procedure,public :: set_ls_norm_factor  => CVodeSetLSNormFactor  !! `CVodeSetLSNormFactor`: WRMS -> 2-norm factor
        procedure,public :: get_num_jtsetup_evals => CVodeGetNumJTSetupEvals  !! `CVodeGetNumJTSetupEvals`
        procedure,public :: get_num_jtimes_evals => CVodeGetNumJtimesEvals  !! `CVodeGetNumJtimesEvals`
        procedure,public :: band_prec_init      => CVBandPrecInit  !! `CVBandPrecInit`: band preconditioner (Krylov solvers)
        procedure,public :: get_band_prec_num_rhs_evals => CVBandPrecGetNumRhsEvals  !! `CVBandPrecGetNumRhsEvals`
        ! linear solver callbacks (ls_problem)
        procedure :: ls_atimes => cvLsATimes  !! `cvLsATimes`: z = (I - gamma J) v
        procedure :: ls_psetup => cvLsPSetup  !! `cvLsPSetup`: preconditioner setup
        procedure :: ls_psolve => cvLsPSolve  !! `cvLsPSolve`: preconditioner solve

        ! root finding (cvode.c / cvode_io.c, submodule fvode_core_root)
        procedure,public :: root_init           => CVodeRootInit  !! `CVodeRootInit`: specify root functions g(t,y)
        procedure,public :: set_root_direction  => CVodeSetRootDirection  !! `CVodeSetRootDirection`: directions of zero-crossings to find
        procedure,public :: set_no_inactive_root_warn => CVodeSetNoInactiveRootWarn  !! `CVodeSetNoInactiveRootWarn`: no warning for g(i) = 0 at t0
        procedure,public :: get_root_info       => CVodeGetRootInfo  !! `CVodeGetRootInfo`: which root functions have a root
        ! inequality constraints and stability limit detection (cvode_io.c)
        procedure,public :: set_constraints     => CVodeSetConstraints  !! `CVodeSetConstraints`: inequality constraints on y
        procedure,public :: set_max_num_constraint_fails => CVodeSetMaxNumConstraintFails  !! `CVodeSetMaxNumConstraintFails`
        procedure,public :: get_num_constraint_fails => CVodeGetNumConstraintFails  !! `CVodeGetNumConstraintFails`
        procedure,public :: get_num_constraint_corrections => CVodeGetNumConstraintCorrections  !! `CVodeGetNumConstraintCorrections`
        procedure,public :: set_stab_lim_det    => CVodeSetStabLimDet  !! `CVodeSetStabLimDet`: turn stability limit detection on/off (BDF)

        ! diagonal linear solver (cvode_diag.c, submodule fvode_core_diag)
        procedure,public :: diag                => CVDiag  !! `CVDiag`: attach the diagonal linear solver
        procedure,public :: get_diag_num_rhs_evals => CVDiagGetNumRhsEvals  !! `CVDiagGetNumRhsEvals`: f evaluations for the diagonal Jacobian
        procedure,public :: get_diag_last_flag  => CVDiagGetLastFlag  !! `CVDiagGetLastFlag`: last flag of the diagonal solver

        ! nonlinear solver interface (cvode_nls.c)
        procedure :: nls_sys      => cvNlsSys  !! `cvNlsSys`: nonlinear system function (residual or fixed-point)
        procedure :: nls_lsetup   => cvNlsLSetup  !! `cvNlsLSetup`: linear solver setup for Newton
        procedure :: nls_lsolve   => cvNlsLSolve  !! `cvNlsLSolve`: linear solve for Newton
        procedure :: nls_convtest => cvNlsConvTest  !! `cvNlsConvTest`: CVODE's nonlinear convergence test

    end type cvode_solver

    abstract interface

        function cv_rhs_fn(me, t, y, ydot) result(retval)
            !! Right-hand side function `ydot = f(t, y)` (`CVRhsFn`).
            !! Return 0 on success, > 0 for a recoverable error, < 0 for an
            !! unrecoverable error.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the solver (extend it to carry data)
            real(wp),intent(in) :: t  !! independent variable
            real(wp),dimension(:),intent(in) :: y  !! state vector at t
            real(wp),dimension(:),intent(out) :: ydot  !! right-hand side f(t,y)
            integer :: retval  !! 0 on success, > 0 for a recoverable error, < 0 for an unrecoverable error
        end function cv_rhs_fn

        function cv_ewt_fn(me, y, ewt) result(retval)
            !! User error-weight function (`CVEwtFn`). Return 0 on success,
            !! nonzero on failure.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(in) :: y  !! current state vector
            real(wp),dimension(:),intent(out) :: ewt  !! error weights (must be > 0)
            integer :: retval  !! 0 on success, nonzero on failure
        end function cv_ewt_fn

        function cv_jac_fn(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
            !! User Jacobian function (`CVLsJacFn`): fill `J = df/dy` at (t, y);
            !! `fy = f(t, y)`. J is the matrix given to `set_linear_solver`
            !! (use `select type` to access its data). Return 0 on success,
            !! > 0 for a recoverable error, < 0 for an unrecoverable error.
            import :: cvode_solver, sun_matrix, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),intent(in) :: t  !! independent variable
            real(wp),dimension(:),intent(in) :: y  !! state vector at t
            real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
            class(sun_matrix),intent(inout) :: J  !! Jacobian matrix df/dy to fill
            real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors of length n
            integer :: retval  !! 0 on success, > 0 for a recoverable error, < 0 for an unrecoverable error
        end function cv_jac_fn

        function cv_psetup_fn(me, t, y, fy, jok, jcur, gamma) result(retval)
            !! Preconditioner setup function (`CVLsPrecSetupFn`): prepare the
            !! preconditioner for `P ~ I - gamma*J`. If `jok`, saved Jacobian
            !! data may be reused; set `jcur = .true.` if the Jacobian data was
            !! recomputed. Return 0 on success, > 0 recoverable, < 0 unrecoverable.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver (extend it to carry data)
            real(wp),intent(in) :: t  !! current time
            real(wp),dimension(:),intent(in) :: y  !! current solution
            real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
            logical,intent(in) :: jok  !! may saved Jacobian data be reused?
            logical,intent(out) :: jcur  !! was the Jacobian data recomputed?
            real(wp),intent(in) :: gamma  !! the scalar gamma in P = I - gamma*J
            integer :: retval  !! 0 on success, > 0 recoverable, < 0 unrecoverable
        end function cv_psetup_fn

        function cv_psolve_fn(me, t, y, fy, r, z, gamma, delta, lr) result(retval)
            !! Preconditioner solve function (`CVLsPrecSolveFn`): solve
            !! `P z = r` (`lr` = 1 left, 2 right preconditioner), approximately
            !! to the tolerance `delta` if done iteratively. Return 0 on success,
            !! > 0 recoverable, < 0 unrecoverable.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver (extend it to carry data)
            real(wp),intent(in) :: t  !! current time
            real(wp),dimension(:),intent(in) :: y  !! current solution
            real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
            real(wp),dimension(:),intent(in) :: r  !! the right-hand side
            real(wp),dimension(:),intent(out) :: z  !! the solution
            real(wp),intent(in) :: gamma  !! the scalar gamma in P = I - gamma*J
            real(wp),intent(in) :: delta  !! tolerance for an iterative solve
            integer,intent(in) :: lr  !! 1: left preconditioner, 2: right
            integer :: retval  !! 0 on success, > 0 recoverable, < 0 unrecoverable
        end function cv_psolve_fn

        function cv_jtsetup_fn(me, t, y, fy) result(retval)
            !! Jacobian-times-vector setup function (`CVLsJacTimesSetupFn`):
            !! prepare data for [[cv_jtimes_fn]]. Return 0 on success.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver (extend it to carry data)
            real(wp),intent(in) :: t  !! current time
            real(wp),dimension(:),intent(in) :: y  !! current solution
            real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
            integer :: retval  !! 0 on success, nonzero on failure
        end function cv_jtsetup_fn

        function cv_jtimes_fn(me, v, Jv, t, y, fy, tmp) result(retval)
            !! Jacobian-times-vector function (`CVLsJacTimesVecFn`):
            !! `Jv = J(t,y) v`. Return 0 on success, nonzero on failure.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver (extend it to carry data)
            real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
            real(wp),dimension(:),intent(out) :: Jv  !! the product J*v
            real(wp),intent(in) :: t  !! current time
            real(wp),dimension(:),intent(in) :: y  !! current solution
            real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
            real(wp),dimension(:),intent(inout) :: tmp  !! work vector
            integer :: retval  !! 0 on success, nonzero on failure
        end function cv_jtimes_fn

        function cv_root_fn(me, t, y, gout) result(retval)
            !! Root function (`CVRootFn`): compute `gout(i) = g_i(t, y)`,
            !! i = 1..nrtfn. CVODE finds the roots of these functions.
            !! Return 0 on success, nonzero on failure.
            import :: cvode_solver, wp
            implicit none
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver (extend it to carry data)
            real(wp),intent(in) :: t  !! current value of the independent variable
            real(wp),dimension(:),intent(in) :: y  !! current value of the dependent variable vector
            real(wp),dimension(:),intent(out) :: gout  !! the root function values g_i(t, y)
            integer :: retval  !! 0 on success, nonzero on failure
        end function cv_root_fn

    end interface

    public :: cv_rhs_fn, cv_ewt_fn, cv_jac_fn, cv_root_fn
    public :: cv_psetup_fn, cv_psolve_fn, cv_jtsetup_fn, cv_jtimes_fn

    !-------------------------------------------------------------------------------------
    ! linear solver interface: implemented in submodule fvode_core_ls (cvode_ls.c)
    !-------------------------------------------------------------------------------------
    interface

        module subroutine CVodeSetLinearSolver(me, LS, A, flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            class(sun_linear_solver),intent(in) :: LS  !! the linear solver (copied)
            class(sun_matrix),intent(in),optional :: A  !! the matrix (copied; required for direct solvers)
            integer,intent(out) :: flag  !! return flag: `CVLS_SUCCESS` or a negative error code
        end subroutine CVodeSetLinearSolver

        module subroutine CVodeSetJacFn(me, jac, flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            procedure(cv_jac_fn),optional :: jac  !! user Jacobian function (absent: use the difference-quotient Jacobian)
            integer,intent(out) :: flag  !! return flag: `CVLS_SUCCESS` or a negative error code
        end subroutine CVodeSetJacFn

        module subroutine CVodeSetJacEvalFrequency(me, msbj, flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer(li),intent(in) :: msbj  !! max steps between Jacobian evaluations (0: default)
            integer,intent(out) :: flag  !! return flag: `CVLS_SUCCESS` or a negative error code
        end subroutine CVodeSetJacEvalFrequency

        module subroutine CVodeSetLinearSolutionScaling(me, onoff, flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            logical,intent(in) :: onoff  !! enable (`.true.`) or disable the solution scaling
            integer,intent(out) :: flag  !! return flag: `CVLS_SUCCESS` or a negative error code
        end subroutine CVodeSetLinearSolutionScaling

        module subroutine CVodeSetDeltaGammaMaxBadJac(me, dgmax_jbad, flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),intent(in) :: dgmax_jbad  !! gamma ratio threshold (<= 0: default)
            integer,intent(out) :: flag  !! return flag: `CVLS_SUCCESS` or a negative error code
        end subroutine CVodeSetDeltaGammaMaxBadJac

        module function cvLsInitialize(me) result(flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer :: flag  !! 0 on success, otherwise the linear solver's error flag
        end function cvLsInitialize

        module function cvLsReInitialize(me) result(flag)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer :: flag  !! always `CVLS_SUCCESS`
        end function cvLsReInitialize

        module function cvLsSetup(me, convfail) result(flag)
            !! See the implementation in the submodule (`fvode_core_ls`). Uses
            !! `me%y` (ypred) and `me%ftemp` (fpred), and sets `me%jcur`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer,intent(in) :: convfail  !! convergence failure status (`CV_NO_FAILURES`, `CV_FAIL_BAD_J`, `CV_FAIL_OTHER`)
            integer :: flag  !! 0 on success, > 0 for a recoverable failure, < 0 otherwise
        end function cvLsSetup

        module function cvLsSolve(me, b) result(flag)
            !! See the implementation in the submodule (`fvode_core_ls`). Uses
            !! `me%ewt` (weights), `me%y` (ynow) and `me%ftemp` (fnow).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(inout) :: b  !! right-hand side on input, solution on output
            integer :: flag  !! 0 on success, 1 for a recoverable failure, -1 for an unrecoverable one
        end function cvLsSolve

        module subroutine CVodeSetPreconditioner(me, psetup, psolve, flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            procedure(cv_psetup_fn),optional :: psetup  !! preconditioner setup function (absent: none)
            procedure(cv_psolve_fn),optional :: psolve  !! preconditioner solve function (absent: none)
            integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`
        end subroutine CVodeSetPreconditioner

        module subroutine CVodeSetJacTimes(me, jtsetup, jtimes, flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            procedure(cv_jtsetup_fn),optional :: jtsetup  !! J*v setup function (absent: none)
            procedure(cv_jtimes_fn),optional :: jtimes  !! J*v function (absent: internal DQ)
            integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`
        end subroutine CVodeSetJacTimes

        module subroutine CVodeSetJacTimesRhsFn(me, jtimesRhsFn, flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            procedure(cv_rhs_fn),optional :: jtimesRhsFn  !! f used by the DQ J*v (absent: the ODE f)
            integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`
        end subroutine CVodeSetJacTimesRhsFn

        module subroutine CVodeSetEpsLin(me, eplifac, flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),intent(in) :: eplifac  !! linear tolerance factor (0: default 0.05)
            integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`
        end subroutine CVodeSetEpsLin

        module subroutine CVodeSetLSNormFactor(me, nrmfac, flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),intent(in) :: nrmfac  !! > 0: the factor; < 0: computed as sqrt(dot(1,1)); 0: sqrt(N)
            integer,intent(out) :: flag  !! `CVLS_SUCCESS` or `CVLS_LMEM_NULL`
        end subroutine CVodeSetLSNormFactor

        module function cvLsATimes(me, v, z) result(flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
            real(wp),dimension(:),intent(out) :: z  !! z = (I - gamma J) v
            integer :: flag  !! 0 on success, the jtimes flag otherwise
        end function cvLsATimes

        module function cvLsPSetup(me) result(flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer :: flag  !! 0 on success, > 0 recoverable, < 0 unrecoverable
        end function cvLsPSetup

        module function cvLsPSolve(me, r, z, tol, lr) result(flag)
            !! Implemented in submodule `fvode_core_ls`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(in) :: r  !! the right-hand side
            real(wp),dimension(:),intent(out) :: z  !! the solution of P z = r
            real(wp),intent(in) :: tol  !! tolerance
            integer,intent(in) :: lr  !! `SUN_PREC_LEFT` or `SUN_PREC_RIGHT`
            integer :: flag  !! 0 on success, > 0 recoverable, < 0 unrecoverable
        end function cvLsPSolve

        ! band preconditioner: submodule fvode_core_bandpre (cvode_bandpre.c)

        module subroutine CVBandPrecInit(me, N, mu, ml, flag)
            !! Implemented in submodule `fvode_core_bandpre`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer(ip),intent(in) :: N  !! problem size
            integer(ip),intent(in) :: mu  !! upper half-bandwidth of the preconditioner
            integer(ip),intent(in) :: ml  !! lower half-bandwidth of the preconditioner
            integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL`, `CVLS_SUNLS_FAIL`, ...
        end subroutine CVBandPrecInit

        module function cvBandPrecSetup(me, t, y, fy, jok, jcur, gamma) result(retval)
            !! Implemented in submodule `fvode_core_bandpre`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),intent(in) :: t  !! current time
            real(wp),dimension(:),intent(in) :: y  !! current solution
            real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
            logical,intent(in) :: jok  !! may the saved Jacobian be reused?
            logical,intent(out) :: jcur  !! was the Jacobian recomputed?
            real(wp),intent(in) :: gamma  !! the scalar gamma
            integer :: retval  !! 0 on success, > 0 recoverable, < 0 unrecoverable
        end function cvBandPrecSetup

        module function cvBandPrecSolve(me, r, z) result(retval)
            !! Implemented in submodule `fvode_core_bandpre`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(in) :: r  !! the right-hand side
            real(wp),dimension(:),intent(out) :: z  !! the solution of P z = r
            integer :: retval  !! the band solver flag
        end function cvBandPrecSolve

        module subroutine cvBandPrecFree(me)
            !! Implemented in submodule `fvode_core_bandpre`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        end subroutine cvBandPrecFree

        module subroutine cvLsFree(me)
            !! See the implementation in the submodule (`fvode_core_ls`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        end subroutine cvLsFree

        ! root finding: submodule fvode_core_root (cvode.c)

        module subroutine CVodeRootInit(me, nrtfn, g, flag)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer,intent(in) :: nrtfn  !! number of root functions (<= 0: turn root finding off)
            procedure(cv_root_fn),optional :: g  !! the root function
            integer,intent(out) :: flag  !! `CV_SUCCESS`, `CV_ILL_INPUT` or `CV_MEM_FAIL`
        end subroutine CVodeRootInit

        module subroutine CVodeSetRootDirection(me, rootdir, flag)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer,dimension(:),intent(in) :: rootdir  !! direction for each root function: +1 increasing, -1 decreasing, 0 both
            integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`
        end subroutine CVodeSetRootDirection

        module subroutine CVodeSetNoInactiveRootWarn(me, flag)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer,intent(out) :: flag  !! `CV_SUCCESS`
        end subroutine CVodeSetNoInactiveRootWarn

        module subroutine CVodeGetRootInfo(me, rootsfound)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(in) :: me  !! the CVODE solver
            integer,dimension(:),intent(out) :: rootsfound  !! +1 / -1: root of g(i) found, increasing / decreasing; 0: none
        end subroutine CVodeGetRootInfo

        module function cvRcheck1(me) result(flag)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer :: flag  !! `CV_SUCCESS` or `CV_RTFUNC_FAIL`
        end function cvRcheck1

        module function cvRcheck2(me) result(flag)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer :: flag  !! `CV_SUCCESS`, `RTFOUND`, `CLOSERT` or `CV_RTFUNC_FAIL`
        end function cvRcheck2

        module function cvRcheck3(me, tout, itask) result(flag)
            !! Implemented in submodule `fvode_core_root`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),intent(in) :: tout  !! the user's output time
            integer,intent(in) :: itask  !! `CV_NORMAL` or `CV_ONE_STEP`
            integer :: flag  !! `CV_SUCCESS`, `RTFOUND` or `CV_RTFUNC_FAIL`
        end function cvRcheck3

        ! stability limit detection: submodule fvode_core_stald (cvode.c)

        module subroutine cvBDFStab(me)
            !! Implemented in submodule `fvode_core_stald`.
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        end subroutine cvBDFStab

        ! diagonal linear solver: submodule fvode_core_diag (cvode_diag.c)

        module subroutine CVDiag(me, flag)
            !! See the implementation in the submodule (`fvode_core_diag`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer,intent(out) :: flag  !! return flag: `CVDIAG_SUCCESS` or a negative error code
        end subroutine CVDiag

        module function CVDiagInit(me) result(flag)
            !! See the implementation in the submodule (`fvode_core_diag`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            integer :: flag  !! always 0
        end function CVDiagInit

        module function CVDiagSetup(me, ypred, fpred, jcur) result(flag)
            !! See the implementation in the submodule (`fvode_core_diag`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(in) :: ypred  !! predicted state
            real(wp),dimension(:),intent(in) :: fpred  !! f at the predicted state
            logical,intent(out) :: jcur  !! set to `.true.` on success
            integer :: flag  !! 0 on success, 1 for a recoverable failure, -1 for an unrecoverable one
        end function CVDiagSetup

        module function CVDiagSolve(me, b) result(flag)
            !! See the implementation in the submodule (`fvode_core_diag`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
            real(wp),dimension(:),intent(inout) :: b  !! right-hand side on input, solution on output
            integer :: flag  !! 0 on success, 1 if M is not invertible
        end function CVDiagSolve

        module subroutine CVDiagFree(me)
            !! See the implementation in the submodule (`fvode_core_diag`).
            class(cvode_solver),intent(inout) :: me  !! the CVODE solver
        end subroutine CVDiagFree

    end interface
    public :: CVodeGetReturnFlagName

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: lmm   !! `CV_ADAMS` or `CV_BDF`
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`

    integer :: maxord  !! maximum order of the method

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

    ! initialize constraint variables
    me%constraint_corrections = 0
    me%constraint_fails       = 0
    me%max_constraint_fails   = MAX_CONSTRAINT_FAILS

    ! initialize root finding variables
    me%gfun    => null()
    me%nrtfn   = 0
    me%mxgnull = 1

    ! [not yet ported] projection, resize

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
!  As in C, a Newton nonlinear solver is attached by default; use
!  [[CVodeSetNonlinearSolver]] to attach a different one (e.g. fixed point).
!  A Newton solver also needs a linear solver ([[CVodeSetLinearSolver]] or
!  [[CVDiag]]).

    subroutine CVodeInit(me, f, t0, y0, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_rhs_fn) :: f  !! right-hand side function
    real(wp),intent(in) :: t0  !! initial time
    real(wp),dimension(:),intent(in) :: y0  !! initial state vector
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    integer :: i, k  !! indices into `ssdat`

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

    ! create a Newton nonlinear solver object and attach it to CVode
    call CVodeSetNonlinearSolverInternal(me, sun_nonlinsol_newton(int(me%n)), flag)
    if (flag /= CV_SUCCESS) then
        call cvProcessError(me, flag, 'CVodeInit', 'Setting the nonlinear solver failed')
        call cvFreeVectors(me)
        flag = CV_MEM_FAIL
        return
    end if

    ! initialize all the counters and other optional output values
    me%q      = 1
    me%LL     = 2
    me%qwait  = me%LL
    me%etamax = me%eta_max_fs

    me%qu    = 0
    me%hu    = ZERO
    me%tolsf = ONE

    ! set the linear solver addresses to NULL
    ! (we check != NULL later, in CVode)
    call cvLFree(me)

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
            me%ssdat(i-1,k-1) = ZERO
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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t0  !! initial time
    real(wp),dimension(:),intent(in) :: y0  !! initial state vector (same size as in `init`)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    integer :: i, k  !! indices into `ssdat`

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

    me%constraint_corrections = 0
    me%constraint_fails       = 0

    select case (me%ls_kind)
    case (LS_GENERIC); flag = cvLsReInitialize(me)
    case (LS_DIAG);    flag = CVDiagInit(me)
    end select

    ! initialize other integrator optional outputs
    me%h0u    = ZERO
    me%next_h = ZERO
    me%next_q = 0

    ! initialize stability limit detection data
    me%nor = 0
    do i = 1, 5
        do k = 1, 3
            me%ssdat(i-1,k-1) = ZERO
        end do
    end do

    flag = CV_SUCCESS

    end subroutine CVodeReInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSStolerances`: scalar relative and absolute tolerances.

    subroutine CVodeSStolerances(me, reltol, abstol, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: reltol  !! relative tolerance (>= 0)
    real(wp),intent(in) :: abstol  !! absolute tolerance (>= 0)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: reltol  !! relative tolerance (>= 0)
    real(wp),dimension(:),intent(in) :: abstol  !! absolute tolerances, one per component (>= 0)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    real(wp) :: atolmin  !! smallest absolute tolerance

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_ewt_fn) :: efun  !! user function that computes the error weights
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: tout  !! next time at which a solution is desired
    real(wp),dimension(:),intent(inout) :: yout  !! solution at `tret` (output)
    real(wp),intent(out) :: tret  !! time reached (output)
    integer,intent(in) :: itask  !! `CV_NORMAL` or `CV_ONE_STEP`
    integer,intent(out) :: flag  !! `CV_SUCCESS`, `CV_TSTOP_RETURN`, `CV_ROOT_RETURN`, or a negative error code

    integer(li) :: nstloc  !! number of steps taken in this call
    integer :: retval, hflag, kflag, istate, ier  !! return flags of f, cvHin, cvStep, the call, and CVodeGetDky / cvInitialSetup
    integer :: ewtsetOK  !! return flag of the error weight function (0: success)
    real(wp) :: troundoff, tout_hin, rh, nrm  !! roundoff in t, tout used by cvHin, |h|/hmax, norm of zn(:,0)
    integer :: irfndp  !! value of irfnd on entry (did the previous step find a root?)
    integer :: ir  !! root function index
    logical :: inactive_roots  !! are some root functions inactive after the first step?

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

        ! check for zeros of root function g at and near t0
        if (me%nrtfn > 0) then
            retval = cvRcheck1(me)
            if (retval == CV_RTFUNC_FAIL) then
                call cvProcessError(me, CV_RTFUNC_FAIL, 'cvRcheck1', 'At '//t_str(me%tn)// &
                    ', the rootfinding routine failed in an unrecoverable manner.')
                flag = CV_RTFUNC_FAIL
                return
            end if
        end if

    end if ! end of first call block

    !------------------------------------------------------------------
    ! at following steps, perform stop tests
    !------------------------------------------------------------------

    if (me%nst > 0) then

        ! estimate an infinitesimal time interval to be used as
        ! a roundoff for time quantities (based on current time
        ! and step size)
        troundoff = FUZZ_FACTOR * me%uround * (abs(me%tn) + abs(me%h))

        ! first, check for a root in the last step taken, other than the
        ! last root found, if any. If itask = CV_ONE_STEP and y(tn) was not
        ! returned because of an intervening root, return y(tn) now.
        if (me%nrtfn > 0) then

            irfndp = me%irfnd

            retval = cvRcheck2(me)

            if (retval == CLOSERT) then
                call cvProcessError(me, CV_ILL_INPUT, 'cvRcheck2', &
                    'Root found at and very near '//t_str(me%tlo)//'.')
                flag = CV_ILL_INPUT
                return
            else if (retval == CV_RTFUNC_FAIL) then
                call cvProcessError(me, CV_RTFUNC_FAIL, 'cvRcheck2', 'At '//t_str(me%tlo)// &
                    ', the rootfinding routine failed in an unrecoverable manner.')
                flag = CV_RTFUNC_FAIL
                return
            else if (retval == RTFOUND) then
                me%tretlast = me%tlo
                tret = me%tretlast
                yout = me%y    ! (C: cv_y is yout)
                flag = CV_ROOT_RETURN
                return
            end if

            ! if tn is distinct from tretlast (within roundoff),
            ! check remaining interval for roots
            if (abs(me%tn - me%tretlast) > troundoff) then

                retval = cvRcheck3(me, tout, itask)

                if (retval == CV_SUCCESS) then      ! no root found
                    me%irfnd = 0
                    if ((irfndp == 1) .and. (itask == CV_ONE_STEP)) then
                        me%tretlast = me%tn
                        tret = me%tretlast
                        yout = me%zn(:,0)
                        flag = CV_SUCCESS
                        return
                    end if
                else if (retval == RTFOUND) then    ! a new root was found
                    me%irfnd = 1
                    me%tretlast = me%tlo
                    tret = me%tretlast
                    yout = me%y    ! (C: cv_y is yout)
                    flag = CV_ROOT_RETURN
                    return
                else if (retval == CV_RTFUNC_FAIL) then  ! g failed
                    call cvProcessError(me, CV_RTFUNC_FAIL, 'cvRcheck3', 'At '//t_str(me%tlo)// &
                        ', the rootfinding routine failed in an unrecoverable manner.')
                    flag = CV_RTFUNC_FAIL
                    return
                end if

            end if

        end if ! end of root stop check

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

        ! check for root in last step taken
        if (me%nrtfn > 0) then

            retval = cvRcheck3(me, tout, itask)

            if (retval == RTFOUND) then  ! a new root was found
                me%irfnd = 1
                istate = CV_ROOT_RETURN
                me%tretlast = me%tlo
                tret = me%tretlast
                yout = me%y    ! (C: cv_y is yout)
                exit
            else if (retval == CV_RTFUNC_FAIL) then  ! g failed
                call cvProcessError(me, CV_RTFUNC_FAIL, 'cvRcheck3', 'At '//t_str(me%tlo)// &
                    ', the rootfinding routine failed in an unrecoverable manner.')
                istate = CV_RTFUNC_FAIL
                exit
            end if

            ! if we are at the end of the first step and we still have
            ! some event functions that are inactive, issue a warning
            ! as this may indicate a user error in the implementation
            ! of the root function.
            if (me%nst == 1) then
                inactive_roots = .false.
                do ir = 1, me%nrtfn
                    if (.not. me%gactive(ir)) then
                        inactive_roots = .true.
                        exit
                    end if
                end do
                if ((me%mxgnull > 0) .and. inactive_roots) then
                    call cvProcessError(me, CV_WARNING, 'CVode', &
                        'At the end of the first step, there are still some root functions '// &
                        'identically 0. This warning will not be issued again.')
                end if
            end if

        end if

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! time at which to interpolate (tn - hu <= t <= tn)
    integer,intent(in) :: k  !! derivative order (0 <= k <= q)
    real(wp),dimension(:),intent(inout) :: dky  !! k-th derivative of the interpolating polynomial at t (output)
    integer,intent(out) :: flag  !! `CV_SUCCESS`, `CV_BAD_K`, `CV_BAD_T` or `CV_BAD_DKY`

    real(wp) :: s, r, c  !! scaled time (t - tn)/h, the factor h^(-k), and a polynomial coefficient
    real(wp) :: tfuzz, tp, tn1  !! time slack and the interval limits for t
    integer :: i, j  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    call cvFreeVectors(me)
    if (allocated(me%NLS)) deallocate(me%NLS)
    call cvLFree(me)
    call cvFreeRoots(me)
    me%nrtfn = 0
    me%gfun => null()
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
!  Free the root finding arrays (C: `free(cv_mem->cv_glo)` etc.).

    subroutine cvFreeRoots(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    if (allocated(me%glo))     deallocate(me%glo)
    if (allocated(me%ghi))     deallocate(me%ghi)
    if (allocated(me%grout))   deallocate(me%grout)
    if (allocated(me%iroots))  deallocate(me%iroots)
    if (allocated(me%rootdir)) deallocate(me%rootdir)
    if (allocated(me%gactive)) deallocate(me%gactive)

    end subroutine cvFreeRoots
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cv_mem->cv_lfree(cv_mem)`: free the attached linear solver interface.

    subroutine cvLFree(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    select case (me%ls_kind)
    case (LS_GENERIC); call cvLsFree(me)
    case (LS_DIAG);    call CVDiagFree(me)
    end select
    me%ls_kind       = LS_NONE
    me%lsetup_exists = .false.
    me%lsolve_exists = .false.

    end subroutine cvLFree
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAllocVectors`: allocate the CVODE vectors ewt, acor, tempv, ftemp,
!  vtemp1-3 and zn[0..qmax], all of length `n`.

    subroutine cvAllocVectors(me, n)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(ip),intent(in) :: n  !! problem size

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

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
    if (allocated(me%constraints)) deallocate(me%constraints)
    me%n = 0

    end subroutine cvFreeVectors
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvInitialSetup`: checks done at the first call to `CVode`; computes the
!  initial error weights and initializes the nonlinear solver.

    function cvInitialSetup(me, tout) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: tout  !! first output time
    integer :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    integer :: ier  !! return flag of the called functions
    real(wp) :: tdist, tround  !! distance to tout and roundoff in t
    logical :: conOK  !! does y0 satisfy the constraints?

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

    ! check to see if y0 satisfies constraints
    if (allocated(me%constraints)) then
        conOK = n_vconstrmask(me%constraints, me%zn(:,0), me%tempv)
        if (.not. conOK) then
            call cvProcessError(me, CV_ILL_INPUT, 'cvInitialSetup', &
                                'y0 fails to satisfy constraints.')
            flag = CV_ILL_INPUT
            return
        end if
    end if

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

    ! call linit function (if it exists)
    if (me%ls_kind /= LS_NONE) then
        select case (me%ls_kind)
        case (LS_GENERIC); ier = cvLsInitialize(me)
        case (LS_DIAG);    ier = CVDiagInit(me)
        end select
        if (ier /= 0) then
            call cvProcessError(me, CV_LINIT_FAIL, 'cvInitialSetup', &
                                "The linear solver's init routine failed.")
            flag = CV_LINIT_FAIL
            return
        end if
    end if

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: tout  !! first output time (or tstop)
    integer :: flag  !! `CV_SUCCESS`, `CV_RHSFUNC_FAIL` or `CV_REPTD_RHSFUNC_ERR`

    integer :: retval, sgn, count1, count2  !! return flag of cvYddNorm, direction of integration (+1/-1), and the loop counters
    real(wp) :: tdiff, tdist, tround, hlb, hub  !! tout - t0, its magnitude, roundoff in t, and lower/upper bounds on h0
    real(wp) :: hg, hgs, hs, hnew, hrat, h0, yddnrm  !! trial step, signed trial step, last feasible step, new step, hnew/hg, result, norm of y''
    logical :: hgOK  !! was y'' computed successfully for the trial step?

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: tdist  !! distance |tout - t0|
    real(wp) :: hub  !! upper bound on the initial step size

    real(wp) :: hub_inv  !! inverse of the bound based on |y0|/|y0'|
    integer :: ier  !! return flag of the error weight function (not used, as in C)

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: hg  !! trial step size
    real(wp),intent(out) :: yddnrm  !! WRMS norm of the estimated y'' (output)
    integer :: flag  !! `CV_SUCCESS`, `CV_RHSFUNC_FAIL` or `RHSFUNC_RECVR`

    integer :: retval  !! return flag of f

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

    me%tempv = (ONE / hg) * (me%tempv - me%zn(:,1))

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! `CV_SUCCESS` or a failure flag for cvHandleFailure

    real(wp) :: saved_t  !! time to restore to if a failure occurs
    real(wp) :: dsm  !! local truncation error estimate
    integer :: nflag, kflag  !! nonlinear solver flags
    integer :: eflag  !! error test return flag

    ! are we computing sensitivities with the simultaneous approach?
    integer :: ncf  !! corrector failures in this step attempt
    integer :: nef  !! error test failures in this step attempt
    integer :: step_constraint_fails  !! constraint failures in this step attempt
    integer :: cflag  !! return flag of cvCheckConstraints

    ncf = 0
    nef = 0
    step_constraint_fails = 0

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

        ! check if constraints are satisfied
        if (allocated(me%constraints)) then
            cflag = cvCheckConstraints(me, nflag, saved_t, step_constraint_fails)
            ! go back in loop if we need to predict again (nflag=PREV_CONV_FAIL)
            if (cflag == PREDICT_AGAIN) cycle
            ! return if predict again was not possible
            if (cflag /= CV_SUCCESS) then
                flag = cflag
                return
            end if
        end if

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
    if (me%sldeton) call cvBDFStab(me)

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: deltaq  !! change in order: +1 or -1

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: deltaq  !! change in order: +1 or -1

    integer :: i, j  !! loop indices
    real(wp) :: xi, hsum  !! scaled step history (t_n - t_(n-j))/h and its running sum

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: deltaq  !! change in order: +1 or -1

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    real(wp) :: alpha0, alpha1, prod, xi, xiold, hsum, A1  !! coefficients and step history quantities of the order increase (see cvode.c)
    integer :: i, j  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    real(wp) :: hsum, xi  !! running sum of step sizes and its scaled value
    integer :: i, j  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    integer :: j  !! column index of zn
    real(wp) :: factor  !! eta**j

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    integer :: j, k  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    real(wp),dimension(0:L_MAX-1) :: m  !! coefficients of the product polynomial (cvAdamsStart)
    real(wp),dimension(0:2) :: mm  !! alternating sums M[0..2] (C: M[3], renamed)
    real(wp) :: hsum  !! sum of the last q step sizes

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(0:L_MAX-1),intent(inout) :: m  !! coefficients of the product polynomial (output)
    real(wp) :: hsum  !! sum of the step sizes (result)

    real(wp) :: xi_inv, sum  !! h/hsum and an alternating sum
    integer :: i, j  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(0:L_MAX-1),intent(inout) :: m  !! coefficients of the product polynomial
    real(wp),dimension(0:2),intent(inout) :: mm  !! alternating sums M[0..2]
    real(wp),intent(in) :: hsum  !! sum of the step sizes

    integer :: i  !! loop index
    real(wp) :: M0_inv, xi, xi_inv  !! 1/M[0], hsum/h and its inverse

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

    integer,intent(in) :: iend  !! last index of the sum
    real(wp),dimension(0:),intent(in) :: a  !! terms a(0:iend)
    integer,intent(in) :: k  !! offset in the denominators
    real(wp) :: sum  !! the alternating sum

    integer :: i, sgn  !! loop index and the sign (+1/-1)

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    real(wp) :: alpha0, alpha0_hat, xi_inv, xistar_inv, hsum  !! BDF coefficient quantities (see cvode.c) and the sum of step sizes
    integer :: i, j  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: hsum  !! sum of the step sizes
    real(wp),intent(in) :: alpha0  !! alpha_0 of the BDF formula
    real(wp),intent(in) :: alpha0_hat  !! alpha_0 hat
    real(wp),intent(in) :: xi_inv  !! h / hsum
    real(wp),intent(in) :: xistar_inv  !! 1 / xi*

    real(wp) :: A1, A2, A3, A4, A5, A6  !! intermediate quantities (as in C)
    real(wp) :: C, Cpinv, Cppinv  !! error constants (as in C)
    real(wp) :: hsum_, xi_inv_  !! C modifies its by-value copies

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: nflag  !! `FIRST_CALL`, `PREV_CONV_FAIL` or `PREV_ERR_FAIL`
    integer :: flag  !! `CV_SUCCESS` or a nonlinear solver failure flag

    logical :: callSetup  !! call the linear solver setup in this solve?
    integer(li) :: nni_inc, nnf_inc  !! nonlinear iterations and failures in this solve

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
!  `cvCheckConstraints`: check whether the corrected solution `y` satisfies
!  the inequality constraints. If not, try to correct acor so that the
!  violating components are moved onto the constraint boundary (if the
!  correction is small, compared with the nonlinear tolerance `tq(4)`).
!  Otherwise restore the history, reduce the step size, and return
!  `PREDICT_AGAIN` (or `CV_CONSTR_FAIL` if h cannot be reduced further).

    function cvCheckConstraints(me, nflag, saved_t, step_constraint_fails) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(inout) :: nflag  !! nonlinear solver flag (set to `PREV_CONV_FAIL` on retry)
    real(wp),intent(in) :: saved_t  !! time to restore to on failure
    integer,intent(inout) :: step_constraint_fails  !! constraint failures in this step attempt
    integer :: flag  !! `CV_SUCCESS`, `PREDICT_AGAIN` or `CV_CONSTR_FAIL`

    logical :: constraintsPassed  !! are all constraints satisfied?
    real(wp) :: vnorm  !! weighted norm of the proposed correction

    associate (mm => me%ftemp, tmp => me%tempv)  ! mm: mask of the violated constraints

    ! get mask vector mm, set where constraints failed
    constraintsPassed = n_vconstrmask(me%constraints, me%y, mm)
    if (constraintsPassed) then
        flag = CV_SUCCESS
        return
    end if

    ! constraints not met: compute the proposed correction
    ! tmp = mm * (y - 0.1 * c' / ewt), where c' = c for |c| = 2 (strict
    ! inequality: keep a small margin from zero) and 0 for |c| = 1, and
    ! check if it is small enough to be applied
    tmp = merge(ONE, ZERO, abs(me%constraints) >= ONEPT5)   ! N_VCompare
    tmp = tmp * me%constraints
    tmp = tmp / me%ewt
    me%vtemp1 = -PT1 * tmp
    tmp = me%y - PT1 * tmp
    tmp = tmp * mm

    vnorm = n_vwrmsnorm(tmp, me%ewt)

    ! if vnorm is small enough, correct acor so that y = zn[0] + acor satisfies
    ! the constraints, and continue
    if (vnorm <= me%tq(4)) then

        me%constraint_corrections = me%constraint_corrections + 1

        ! acor = acor - mm*acor - mm*zn[0] - mm*vtemp1, so that the violating
        ! components of y become 0.1*c'/ewt (just inside the constraint)
        tmp = mm * me%acor
        me%acor = me%acor - tmp
        tmp = mm * me%zn(:,0)
        me%acor = me%acor - tmp
        me%vtemp1 = mm * me%vtemp1
        me%acor = me%acor - me%vtemp1

        flag = CV_SUCCESS
        return
    end if

    ! the correction is too large: reduce the step size
    step_constraint_fails = step_constraint_fails + 1
    me%constraint_fails = me%constraint_fails + 1

    ! restore zn
    call cvRestore(me, saved_t)

    ! check for |h| == hmin
    if (abs(me%h) <= me%hmin * ONEPSM) then
        flag = CV_CONSTR_FAIL
        return
    end if

    ! check for max step attempt failures
    if (step_constraint_fails == me%max_constraint_fails) then
        flag = CV_CONSTR_FAIL
        return
    end if

    ! reduce h by computing eta = h'/h: the largest ratio (times 0.9) that keeps
    ! the violating components on the right side, bounded below by 0.1
    tmp = me%zn(:,0) - me%y
    tmp = mm * tmp

    me%eta = PT9 * n_vminquotient(me%zn(:,0), tmp)
    me%eta = max(me%eta, PT1)
    me%eta = max(me%eta, me%hmin / abs(me%h))
    call cvRescale(me)
    nflag = PREV_CONV_FAIL

    flag = PREDICT_AGAIN

    end associate

    end function cvCheckConstraints
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvHandleNFlag`: take action on the return value nflag = cvNls, and
!  return `DO_ERROR_TEST`, `PREDICT_AGAIN` or a failure flag.

    function cvHandleNFlag(me, nflag, saved_t, ncf) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(inout) :: nflag  !! cvNls return flag on input; `PREV_CONV_FAIL` on a retry
    real(wp),intent(in) :: saved_t  !! time at the start of the step
    integer,intent(inout) :: ncf  !! number of convergence failures in this step
    integer :: flag  !! `DO_ERROR_TEST`, `PREDICT_AGAIN`, or a failure flag

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: saved_t  !! time at the start of the step

    integer :: j, k  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(inout) :: nflag  !! set to `PREV_ERR_FAIL` on failure
    real(wp),intent(in) :: saved_t  !! time at the start of the step
    integer,intent(inout) :: nef  !! number of error test failures in this step
    real(wp),intent(out) :: dsm  !! weighted local error norm (output)
    integer :: flag  !! `CV_SUCCESS`, `TRY_AGAIN`, or a failure flag

    integer :: retval  !! return flag of f

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    integer :: i, j  !! loop indices

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: dsm  !! weighted local error norm of the step

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp) :: etaqm1  !! step size ratio for order q-1

    real(wp) :: ddn  !! error norm at order q-1

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp) :: etaqp1  !! step size ratio for order q+1

    real(wp) :: dup, cquot  !! error norm at order q+1 and the correction factor

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    real(wp) :: etam  !! largest of etaqm1, etaq, etaqp1

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: flag  !! the failure flag
    integer :: istate  !! the flag returned to the user

    character(len=*),parameter :: fn = 'CVode'  !! function name in the messages

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ycur  !! current state vector
    real(wp),dimension(:),intent(inout) :: weight  !! error weights (output)
    integer :: flag  !! 0 on success, nonzero if a weight would be <= 0

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ycur  !! current state vector
    real(wp),dimension(:),intent(inout) :: weight  !! error weights (output)
    integer :: flag  !! 0 on success, -1 if a weight would be <= 0

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ycur  !! current state vector
    real(wp),dimension(:),intent(inout) :: weight  !! error weights (output)
    integer :: flag  !! 0 on success, -1 if a weight would be <= 0

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

    class(cvode_solver),intent(in) :: me  !! the CVODE solver (for the message settings)
    integer,intent(in) :: error_code  !! the error (or `CV_WARNING`) code
    character(len=*),intent(in) :: func  !! name of the reporting function
    character(len=*),intent(in) :: msg  !! the message

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    class(sun_nonlinear_solver),intent(in) :: NLS  !! the nonlinear solver (copied)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeSetNonlinearSolver', &
                            'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    call CVodeSetNonlinearSolverInternal(me, NLS, flag)

    end subroutine CVodeSetNonlinearSolver
!*****************************************************************************************

!*****************************************************************************************
!>
!  The body of `CVodeSetNonlinearSolver` (also used by `CVodeInit` to attach
!  the default Newton solver).

    subroutine CVodeSetNonlinearSolverInternal(me, NLS, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    class(sun_nonlinear_solver),intent(in) :: NLS  !! the nonlinear solver (copied)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    integer :: retval, itype  !! return flag and the solver type

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

    end subroutine CVodeSetNonlinearSolverInternal
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsInit`: record whether the linear solver functions exist, then
!  initialize the nonlinear solver.

    function cvNlsInit(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! `CV_SUCCESS` or `CV_NLS_INIT_FAIL`

    integer :: retval  !! return flag of the solver's initialize

    ! set the linear solver setup / solve wrappers (if they exist)
    me%NLS%has_lsetup = me%lsetup_exists
    me%NLS%has_lsolve = me%lsolve_exists

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ycor  !! correction to the predictor
    real(wp),dimension(:),intent(inout) :: res  !! nonlinear system function value (output)
    integer :: flag  !! `CV_SUCCESS`, `RHSFUNC_RECVR` or `CV_RHSFUNC_FAIL`

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ycor  !! correction to the predictor
    real(wp),dimension(:),intent(inout) :: res  !! residual (output)
    integer :: flag  !! `CV_SUCCESS`, `RHSFUNC_RECVR` or `CV_RHSFUNC_FAIL`

    integer :: retval  !! return flag of f

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ycor  !! correction to the predictor
    real(wp),dimension(:),intent(inout) :: res  !! fixed-point function value (output)
    integer :: flag  !! `CV_SUCCESS`, `RHSFUNC_RECVR` or `CV_RHSFUNC_FAIL`

    integer :: retval  !! return flag of f

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
!  `cvNlsLSetup`: set up the linear solver for the Newton iteration.

    function cvNlsLSetup(me, jbad, jcur) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: jbad  !! `.true.` if the nonlinear solver suspects a bad Jacobian
    logical,intent(out) :: jcur  !! `.true.` if the Jacobian is current (output)
    integer :: flag  !! `CV_SUCCESS`, `SUN_NLS_CONV_RECVR` or `CV_LSETUP_FAIL`

    integer :: retval  !! return flag of the linear solver setup

    ! if the nonlinear solver marked the Jacobian as bad update convfail
    if (jbad) me%convfail = CV_FAIL_BAD_J

    ! setup the linear system
    select case (me%ls_kind)
    case (LS_GENERIC); retval = cvLsSetup(me, me%convfail)
    case (LS_DIAG);    retval = CVDiagSetup(me, me%y, me%ftemp, me%jcur)
    case default;      retval = -1
    end select
    me%nsetups = me%nsetups + 1

    ! update Jacobian status
    jcur = me%jcur

    me%gamrat = ONE
    me%gammap = me%gamma
    me%crate  = ONE
    me%delnrm = ZERO
    me%nstlp  = me%nst

    if (retval < 0) then
        flag = CV_LSETUP_FAIL
    else if (retval > 0) then
        flag = SUN_NLS_CONV_RECVR
    else
        flag = CV_SUCCESS
    end if

    end function cvNlsLSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsLSolve`: solve the linear system for the Newton update (in place).

    function cvNlsLSolve(me, b) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(inout) :: b  !! right-hand side on input, Newton update on output
    integer :: flag  !! `CV_SUCCESS`, `SUN_NLS_CONV_RECVR` or `CV_LSOLVE_FAIL`

    integer :: retval  !! return flag of the linear solve

    select case (me%ls_kind)
    case (LS_GENERIC); retval = cvLsSolve(me, b)
    case (LS_DIAG);    retval = CVDiagSolve(me, b)
    case default;      retval = -1
    end select

    if (retval < 0) then
        flag = CV_LSOLVE_FAIL
    else if (retval > 0) then
        flag = SUN_NLS_CONV_RECVR
    else
        flag = CV_SUCCESS
    end if

    end function cvNlsLSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvNlsConvTest`: the CVODE nonlinear convergence test, based on the
!  estimated convergence rate `crate` and the norm of the update `del`.

    function cvNlsConvTest(me, nls, ycor, del, tol, ewt) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    class(sun_nonlinear_solver),intent(in) :: nls  !! the nonlinear solver (for the iteration count)
    real(wp),dimension(:),intent(in) :: ycor  !! current correction
    real(wp),dimension(:),intent(in) :: del  !! last update
    real(wp),intent(in) :: tol  !! convergence tolerance
    real(wp),dimension(:),intent(in) :: ewt  !! error weights
    integer :: flag  !! `CV_SUCCESS` (converged), `SUN_NLS_CONTINUE` or `SUN_NLS_CONV_RECVR`

    integer :: m  !! current iteration number
    real(wp) :: dcon  !! convergence test quantity

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: print_errors  !! print error and warning messages?
    integer,intent(in),optional :: unit  !! output unit (default: `error_unit`)

    me%print_errors = print_errors
    if (present(unit)) me%errunit = unit

    end subroutine CVodeSetErrorOutput
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxOrd`: maximum method order (may only be decreased).

    subroutine CVodeSetMaxOrd(me, maxord, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: maxord  !! maximum order (cannot exceed the value at creation)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    integer :: qmax_alloc  !! order used to allocate zn

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(li),intent(in) :: mxsteps  !! max internal steps per call (0: default, < 0: no limit)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: mxhnil  !! max number of 't + h = t' warnings
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    me%mxhnil = mxhnil

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxHnilWarns
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetInitStep`: initial step size (0 means estimate it).

    subroutine CVodeSetInitStep(me, hin, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: hin  !! initial step size (0: estimate it)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    me%hin = hin

    ! if an initial step size was set, take it into account at the first step
    flag = CV_SUCCESS

    end subroutine CVodeSetInitStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMinStep`: minimum absolute step size (0 means no limit).

    subroutine CVodeSetMinStep(me, hmin, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: hmin  !! minimum absolute step size (0: no limit)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: hmax  !! maximum absolute step size (0: no limit)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    real(wp) :: hmax_inv  !! 1/hmax

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: tstop  !! stop time
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: interp  !! interpolate the solution at tstop?
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    me%tstopinterp = interp

    flag = CV_SUCCESS

    end subroutine CVodeSetInterpolateStopTime
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeClearStopTime`: disable the stop time.

    subroutine CVodeClearStopTime(me, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    me%tstopset = .false.

    flag = CV_SUCCESS

    end subroutine CVodeClearStopTime
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxErrTestFails`: maximum number of error test failures per step
!  (<= 0 restores the default).

    subroutine CVodeSetMaxErrTestFails(me, maxnef, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: maxnef  !! max error test failures per step (<= 0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: maxncf  !! max nonlinear convergence failures per step (<= 0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

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

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: nlscoef  !! nonlinear convergence coefficient (<= 0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (nlscoef <= ZERO) then
        me%nlscoef = CORTES
    else
        me%nlscoef = nlscoef
    end if

    flag = CV_SUCCESS

    end subroutine CVodeSetNonlinConvCoef
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxNonlinIters`: maximum nonlinear iterations per step (default 3).

    subroutine CVodeSetMaxNonlinIters(me, maxcor, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: maxcor  !! max nonlinear iterations per step
    integer,intent(out) :: flag  !! `CV_SUCCESS`, `CV_MEM_FAIL`, or the solver's error flag

    if (.not. allocated(me%NLS)) then
        call cvProcessError(me, CV_MEM_FAIL, 'CVodeSetMaxNonlinIters', 'A memory request failed.')
        flag = CV_MEM_FAIL
        return
    end if

    flag = me%NLS%set_max_iters(maxcor)

    end subroutine CVodeSetMaxNonlinIters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetDeltaGammaMaxLSetup`: gamma ratio threshold for calling the linear
!  solver setup (< 0 restores the default).

    subroutine CVodeSetDeltaGammaMaxLSetup(me, dgmax_lsetup, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: dgmax_lsetup  !! gamma ratio threshold (< 0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (dgmax_lsetup < ZERO) then
        me%dgmax_lsetup = DGMAX_LSETUP_DEFAULT
    else
        me%dgmax_lsetup = dgmax_lsetup
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetDeltaGammaMaxLSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetLSetupFrequency`: maximum steps between linear solver setup calls
!  (0 restores the default).

    subroutine CVodeSetLSetupFrequency(me, msbp, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(li),intent(in) :: msbp  !! max steps between setups (0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (msbp < 0) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetLSetupFrequency', &
                            'A negative setup frequency was provided')
        flag = CV_ILL_INPUT
        return
    end if

    if (msbp == 0) then
        me%msbp = MSBP_DEFAULT
    else
        me%msbp = msbp
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetLSetupFrequency
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaFixedStepBounds`: keep the step size if
!  `eta_min_fx < eta < eta_max_fx` (invalid values restore the defaults).

    subroutine CVodeSetEtaFixedStepBounds(me, eta_min_fx, eta_max_fx, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_min_fx  !! lower bound, in [0, 1]
    real(wp),intent(in) :: eta_max_fx  !! upper bound, >= 1
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_min_fx >= ZERO .and. eta_min_fx <= ONE) then
        me%eta_min_fx = eta_min_fx
    else
        me%eta_min_fx = ETA_MIN_FX_DEFAULT
    end if
    if (eta_max_fx >= ONE) then
        me%eta_max_fx = eta_max_fx
    else
        me%eta_max_fx = ETA_MAX_FX_DEFAULT
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaFixedStepBounds
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaMaxFirstStep`: maximum step size increase after the first step
!  (<= 1 restores the default).

    subroutine CVodeSetEtaMaxFirstStep(me, eta_max_fs, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_max_fs  !! max step size ratio (> 1; else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_max_fs <= ONE) then
        me%eta_max_fs = ETA_MAX_FS_DEFAULT
    else
        me%eta_max_fs = eta_max_fs
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaMaxFirstStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaMaxEarlyStep`: maximum step size increase on early steps
!  (<= 1 restores the default).

    subroutine CVodeSetEtaMaxEarlyStep(me, eta_max_es, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_max_es  !! max step size ratio (> 1; else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_max_es <= ONE) then
        me%eta_max_es = ETA_MAX_ES_DEFAULT
    else
        me%eta_max_es = eta_max_es
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaMaxEarlyStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetNumStepsEtaMaxEarlyStep`: number of early steps that use
!  eta_max_es (< 0 restores the default).

    subroutine CVodeSetNumStepsEtaMaxEarlyStep(me, small_nst, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(li),intent(in) :: small_nst  !! number of early steps (< 0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (small_nst < 0) then
        me%small_nst = SMALL_NST_DEFAULT
    else
        me%small_nst = small_nst
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetNumStepsEtaMaxEarlyStep
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaMax`: maximum step size increase on a general step
!  (<= 1 restores the default).

    subroutine CVodeSetEtaMax(me, eta_max_gs, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_max_gs  !! max step size ratio (> 1; else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_max_gs <= ONE) then
        me%eta_max_gs = ETA_MAX_GS_DEFAULT
    else
        me%eta_max_gs = eta_max_gs
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaMax
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaMin`: minimum step size reduction on a general step
!  (outside (0,1) restores the default).

    subroutine CVodeSetEtaMin(me, eta_min, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_min  !! min step size ratio (in (0,1); else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_min <= ZERO .or. eta_min >= ONE) then
        me%eta_min = ETA_MIN_DEFAULT
    else
        me%eta_min = eta_min
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaMin
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaMinErrFail`: minimum step size reduction after an error test
!  failure (outside (0,1) restores the default).

    subroutine CVodeSetEtaMinErrFail(me, eta_min_ef, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_min_ef  !! min step size ratio (in (0,1); else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_min_ef <= ZERO .or. eta_min_ef >= ONE) then
        me%eta_min_ef = ETA_MIN_EF_DEFAULT
    else
        me%eta_min_ef = eta_min_ef
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaMinErrFail
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaMaxErrFail`: maximum step size ratio after multiple error test
!  failures (outside (0,1) restores the default).

    subroutine CVodeSetEtaMaxErrFail(me, eta_max_ef, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_max_ef  !! max step size ratio (in (0,1); else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_max_ef <= ZERO .or. eta_max_ef >= ONE) then
        me%eta_max_ef = ETA_MAX_EF_DEFAULT
    else
        me%eta_max_ef = eta_max_ef
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaMaxErrFail
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetNumFailsEtaMaxErrFail`: number of error test failures after which
!  eta_max_ef is used (< 0 restores the default).

    subroutine CVodeSetNumFailsEtaMaxErrFail(me, small_nef, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: small_nef  !! number of failures (< 0: default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (small_nef < 0) then
        me%small_nef = SMALL_NEF_DEFAULT
    else
        me%small_nef = small_nef
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetNumFailsEtaMaxErrFail
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEtaConvFail`: step size ratio after a nonlinear convergence
!  failure (outside (0,1) restores the default).

    subroutine CVodeSetEtaConvFail(me, eta_cf, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta_cf  !! step size ratio (in (0,1); else default)
    integer,intent(out) :: flag  !! return flag: `CV_SUCCESS` or a negative error code

    if (eta_cf <= ZERO .or. eta_cf >= ONE) then
        me%eta_cf = ETA_CF_DEFAULT
    else
        me%eta_cf = eta_cf
    end if
    flag = CV_SUCCESS

    end subroutine CVodeSetEtaConvFail
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetConstraints`: inequality constraints on the solution
!  (0: none, 1: y_i >= 0, -1: y_i <= 0, 2: y_i > 0, -2: y_i < 0). Without
!  `constraints`, the constraints are turned off.

    subroutine CVodeSetConstraints(me, constraints, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in),optional :: constraints  !! constraint flag for each component
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`

    real(wp) :: temptest  !! max |constraints(i)|

    if (.not. me%MallocDone) then
        call cvProcessError(me, CV_NO_MALLOC, 'CVodeSetConstraints', 'Attempt to call before CVodeInit.')
        flag = CV_NO_MALLOC
        return
    end if

    ! if there are no constraints, destroy data structures
    if (.not. present(constraints)) then
        if (allocated(me%constraints)) deallocate(me%constraints)
        flag = CV_SUCCESS
        return
    end if

    if (size(constraints, kind=ip) /= me%n) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetConstraints', 'constraints has the wrong size.')
        flag = CV_ILL_INPUT
        return
    end if

    ! check the constraints vector
    temptest = maxval(abs(constraints))
    if ((temptest > TWOPT5) .or. (temptest < HALF)) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetConstraints', &
                            'Illegal values in constraints vector.')
        flag = CV_ILL_INPUT
        return
    end if

    ! load the constraints vector
    me%constraints = constraints

    flag = CV_SUCCESS

    end subroutine CVodeSetConstraints
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxNumConstraintFails`: maximum number of constraint failures in
!  a step attempt (<= 0 restores the default, 10).

    subroutine CVodeSetMaxNumConstraintFails(me, max_fails, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: max_fails  !! maximum number of constraint failures per step
    integer,intent(out) :: flag  !! `CV_SUCCESS`

    if (max_fails <= 0) then
        me%max_constraint_fails = MAX_CONSTRAINT_FAILS
    else
        me%max_constraint_fails = max_fails
    end if

    flag = CV_SUCCESS

    end subroutine CVodeSetMaxNumConstraintFails
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetStabLimDet`: turn the BDF stability limit detection on or off.

    subroutine CVodeSetStabLimDet(me, sldet, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: sldet  !! use stability limit detection?
    integer,intent(out) :: flag  !! `CV_SUCCESS`, or `CV_ILL_INPUT` for the Adams method

    if (sldet .and. (me%lmm /= CV_BDF)) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetStabLimDet', &
            'Attempt to use stability limit detection with the CV_ADAMS method illegal.')
        flag = CV_ILL_INPUT
        return
    end if

    me%sldeton = sldet

    flag = CV_SUCCESS

    end subroutine CVodeSetStabLimDet
!*****************************************************************************************

!=========================================================================================
!  Optional outputs (cvode_io.c)
!=========================================================================================

!>
!  `CVodeGetNumSteps`: number of steps.
    pure integer(li) function CVodeGetNumSteps(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nst
    end function CVodeGetNumSteps

!>
!  `CVodeGetNumRhsEvals`: number of f evaluations.
    pure integer(li) function CVodeGetNumRhsEvals(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nfe
    end function CVodeGetNumRhsEvals

!>
!  `CVodeGetNumLinSolvSetups`: number of linear solver setups.
    pure integer(li) function CVodeGetNumLinSolvSetups(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nsetups
    end function CVodeGetNumLinSolvSetups

!>
!  `CVodeGetNumErrTestFails`: number of error test failures.
    pure integer(li) function CVodeGetNumErrTestFails(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%netf
    end function CVodeGetNumErrTestFails

!>
!  `CVodeGetLastOrder`: order used on the last step.
    pure integer function CVodeGetLastOrder(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%qu
    end function CVodeGetLastOrder

!>
!  `CVodeGetCurrentOrder`: order to be used on the next step.
    pure integer function CVodeGetCurrentOrder(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%next_q
    end function CVodeGetCurrentOrder

!>
!  `CVodeGetActualInitStep`: initial step size actually used.
    pure real(wp) function CVodeGetActualInitStep(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%h0u
    end function CVodeGetActualInitStep

!>
!  `CVodeGetLastStep`: last step size.
    pure real(wp) function CVodeGetLastStep(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%hu
    end function CVodeGetLastStep

!>
!  `CVodeGetCurrentStep`: step size to be used next.
    pure real(wp) function CVodeGetCurrentStep(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%next_h
    end function CVodeGetCurrentStep

!>
!  `CVodeGetCurrentTime`: current internal time.
    pure real(wp) function CVodeGetCurrentTime(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%tn
    end function CVodeGetCurrentTime

!>
!  `CVodeGetTolScaleFactor`: tolerance scale factor.
    pure real(wp) function CVodeGetTolScaleFactor(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%tolsf
    end function CVodeGetTolScaleFactor

!>
!  `CVodeGetNumNonlinSolvIters`: number of nonlinear iterations.
    pure integer(li) function CVodeGetNumNonlinSolvIters(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nni
    end function CVodeGetNumNonlinSolvIters

!>
!  `CVodeGetNumNonlinSolvConvFails`: number of nonlinear convergence failures.
    pure integer(li) function CVodeGetNumNonlinSolvConvFails(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nnf
    end function CVodeGetNumNonlinSolvConvFails

!>
!  `CVodeGetNumStepSolveFails`: number of failed steps due to nonlinear solver failures.
    pure integer(li) function CVodeGetNumStepSolveFails(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%ncfn
    end function CVodeGetNumStepSolveFails

!>
!  `CVodeGetNumConstraintFails`: number of constraint failures.
    pure integer(li) function CVodeGetNumConstraintFails(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%constraint_fails
    end function CVodeGetNumConstraintFails

!>
!  `CVodeGetNumConstraintCorrections`: number of constraint corrections.
    pure integer(li) function CVodeGetNumConstraintCorrections(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%constraint_corrections
    end function CVodeGetNumConstraintCorrections

    pure integer(li) function CVDiagGetNumRhsEvals(me) result(v)
    !! f evaluations for the diagonal Jacobian approximation
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%di_nfeDI
    end function CVDiagGetNumRhsEvals

!>
!  `CVDiagGetLastFlag`: last flag of the diagonal solver.
    pure integer function CVDiagGetLastFlag(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%di_last_flag
    end function CVDiagGetLastFlag

!>
!  `CVodeGetNumJacEvals`: number of Jacobian evaluations.
    pure integer(li) function CVodeGetNumJacEvals(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nje
    end function CVodeGetNumJacEvals

!>
!  `CVodeGetNumLinRhsEvals`: f evaluations in the linear solver interface.
    pure integer(li) function CVodeGetNumLinRhsEvals(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nfeDQ
    end function CVodeGetNumLinRhsEvals

!>
!  `CVodeGetNumLinIters`: number of linear iterations.
    pure integer(li) function CVodeGetNumLinIters(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nli
    end function CVodeGetNumLinIters

!>
!  `CVodeGetNumLinConvFails`: number of linear convergence failures.
    pure integer(li) function CVodeGetNumLinConvFails(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%ncfl
    end function CVodeGetNumLinConvFails

!>
!  `CVodeGetNumPrecEvals`: number of preconditioner setups.
    pure integer(li) function CVodeGetNumPrecEvals(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%npe
    end function CVodeGetNumPrecEvals

!>
!  `CVodeGetNumPrecSolves`: number of preconditioner solves.
    pure integer(li) function CVodeGetNumPrecSolves(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nps
    end function CVodeGetNumPrecSolves

    pure integer(li) function CVodeGetNumJTSetupEvals(me) result(v)
    !! `CVodeGetNumJTSetupEvals`: number of Jacobian-times-vector setup calls.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%njtsetup
    end function CVodeGetNumJTSetupEvals

    pure integer(li) function CVodeGetNumJtimesEvals(me) result(v)
    !! `CVodeGetNumJtimesEvals`: number of Jacobian-times-vector calls.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%njtimes
    end function CVodeGetNumJtimesEvals

    pure integer(li) function CVBandPrecGetNumRhsEvals(me) result(v)
    !! `CVBandPrecGetNumRhsEvals`: f evaluations for the band preconditioner.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%bp_nfeBP
    end function CVBandPrecGetNumRhsEvals

!>
!  `CVodeGetLastLinFlag`: last linear solver flag.
    pure integer function CVodeGetLastLinFlag(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%ls_last_flag
    end function CVodeGetLastLinFlag

    pure real(wp) function CVodeGetJacTime(me) result(v)
    !! time at which the Jacobian was last evaluated
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%tnlj
    end function CVodeGetJacTime

    pure integer(li) function CVodeGetJacNumSteps(me) result(v)
    !! step number at which the Jacobian was last evaluated
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nstlj
    end function CVodeGetJacNumSteps

!>
!  `CVodeGetCurrentGamma`: current gamma.
    pure real(wp) function CVodeGetCurrentGamma(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%gamma
    end function CVodeGetCurrentGamma

!>
!  `CVodeGetNumStabLimOrderReds`: order reductions by STALD.
    pure integer(li) function CVodeGetNumStabLimOrderReds(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nor
    end function CVodeGetNumStabLimOrderReds

!>
!  `CVodeGetNumGEvals`: number of root function evaluations.
    pure integer(li) function CVodeGetNumGEvals(me) result(v)
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    v = me%nge
    end function CVodeGetNumGEvals

    subroutine CVodeGetCurrentState(me, y)
    !! The current internal state vector (C returns a pointer to `cv_y`).
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(out) :: y  !! current internal state vector (output)
    y = me%y
    end subroutine CVodeGetCurrentState

    subroutine CVodeGetErrWeights(me, eweight)
    !! The error weight vector `ewt`.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(out) :: eweight  !! error weight vector (output)
    eweight = me%ewt
    end subroutine CVodeGetErrWeights

    subroutine CVodeGetEstLocalErrors(me, ele)
    !! The estimated local errors of the last step.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(out) :: ele  !! estimated local errors (output)
    ele = me%acor
    end subroutine CVodeGetEstLocalErrors

    subroutine CVodeGetIntegratorStats(me, nsteps, nfevals, nlinsetups, netfails, &
                                       qlast, qcur, hinused, hlast, hcur, tcur)
    !! All main integrator statistics in one call.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    integer(li),intent(out) :: nsteps, nfevals, nlinsetups, netfails  !! steps, f evaluations, linear solver setups, error test failures
    integer,intent(out) :: qlast, qcur  !! order used on the last step and to be used on the next
    real(wp),intent(out) :: hinused, hlast, hcur, tcur  !! initial, last and next step sizes, and the current time
    nsteps     = me%nst
    nfevals    = me%nfe
    nlinsetups = me%nsetups
    netfails   = me%netf
    qlast      = me%qu
    qcur       = me%next_q
    hinused    = me%h0u
    hlast      = me%hu
    hcur       = me%next_h
    tcur       = me%tn
    end subroutine CVodeGetIntegratorStats

    subroutine CVodeGetNonlinSolvStats(me, nniters, nnfails)
    !! Nonlinear solver statistics.
    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    integer(li),intent(out) :: nniters, nnfails  !! nonlinear iterations and convergence failures
    nniters = me%nni
    nnfails = me%nnf
    end subroutine CVodeGetNonlinSolvStats

!*****************************************************************************************
!>
!  `CVodePrintAllStats`: print all integrator statistics, as a table
!  (`SUN_OUTPUTFORMAT_TABLE`) or as one CSV record of name,value pairs
!  (`SUN_OUTPUTFORMAT_CSV`).

    subroutine CVodePrintAllStats(me, iunit, fmt, flag)

    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    integer,intent(in) :: iunit  !! output unit
    integer,intent(in) :: fmt  !! `SUN_OUTPUTFORMAT_TABLE` or `SUN_OUTPUTFORMAT_CSV`
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`

    if (fmt /= SUN_OUTPUTFORMAT_TABLE .and. fmt /= SUN_OUTPUTFORMAT_CSV) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodePrintAllStats', 'Invalid formatting option.')
        flag = CV_ILL_INPUT
        return
    end if

    ! step and method stats
    call put_real('Current time', me%tn, start=.true.)
    call put_long('Steps', me%nst)
    call put_long('Error test fails', me%netf)
    call put_long('NLS step fails', me%ncfn)
    call put_long('Constraint fails', me%constraint_fails)
    call put_long('Constraint corrections', me%constraint_corrections)
    call put_real('Initial step size', me%h0u)
    call put_real('Last step size', me%hu)
    call put_real('Current step size', me%next_h)
    call put_long('Last method order', int(me%qu, li))
    call put_long('Current method order', int(me%next_q, li))
    call put_long('Stab. lim. order reductions', me%nor)

    ! function evaluations
    call put_long('RHS fn evals', me%nfe)

    ! nonlinear solver stats
    call put_long('NLS iters', me%nni)
    call put_long('NLS fails', me%nnf)
    if (me%nst > 0) call put_real('NLS iters per step', real(me%nni, wp) / real(me%nst, wp))

    ! linear solver stats
    call put_long('LS setups', me%nsetups)
    ! (C reads the CVLsMem statistics whenever a linear solver is attached,
    ! which is wrong for CVDiag; FVODE prints them only for the generic interface)
    if (me%ls_kind == LS_GENERIC) then
        call put_long('Jac fn evals', me%nje)
        call put_long('LS RHS fn evals', me%nfeDQ)
        call put_long('Prec setup evals', me%npe)
        call put_long('Prec solves', me%nps)
        call put_long('LS iters', me%nli)
        call put_long('LS fails', me%ncfl)
        call put_long('Jac-times setups', me%njtsetup)
        call put_long('Jac-times evals', me%njtimes)
        if (me%nni > 0) then
            call put_real('LS iters per NLS iter', real(me%nli, wp) / real(me%nni, wp))
            call put_real('Jac evals per NLS iter', real(me%nje, wp) / real(me%nni, wp))
            call put_real('Prec evals per NLS iter', real(me%npe, wp) / real(me%nni, wp))
        end if
    end if

    ! rootfinding stats
    call put_long('Root fn evals', me%nge)

    ! [not yet ported] projection stats

    if (fmt == SUN_OUTPUTFORMAT_CSV) write(iunit, '(a)') ''   ! end the CSV record
    flag = CV_SUCCESS

    contains

        subroutine put_real(name, value, start)
        !! write a real statistic (table row or CSV item)
        character(len=*),intent(in) :: name  !! statistic name
        real(wp),intent(in) :: value  !! its value
        logical,intent(in),optional :: start  !! first item of a CSV record?
        if (fmt == SUN_OUTPUTFORMAT_TABLE) then
            write(iunit, '(a," = ",a)') pad(name), g_str(value)
        else
            if (.not. present(start)) write(iunit, '(a)', advance='no') ','
            write(iunit, '(a,",",a)', advance='no') name, e_str(value)
        end if
        end subroutine put_real

        subroutine put_long(name, value)
        !! write an integer statistic (table row or CSV item)
        character(len=*),intent(in) :: name  !! statistic name
        integer(li),intent(in) :: value  !! its value
        if (fmt == SUN_OUTPUTFORMAT_TABLE) then
            write(iunit, '(a," = ",i0)') pad(name), value
        else
            write(iunit, '(",",a,",",i0)', advance='no') name, value
        end if
        end subroutine put_long

        function pad(name) result(s)
        !! left-justified in SUN_TABLE_WIDTH = 29 characters
        character(len=*),intent(in) :: name  !! text to pad
        character(len=max(29, len(name))) :: s  !! the text, left-justified in 29 characters
        s = name
        end function pad

    end subroutine CVodePrintAllStats
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeGetReturnFlagName`: the name of a CVODE return flag.
!
!  (Upstream returns the misspelled "CV_NLS_SETUPT_FAIL" for
!  `CV_NLS_SETUP_FAIL`; the correct name is returned here.)

    function CVodeGetReturnFlagName(flag) result(name)

    integer,intent(in) :: flag  !! a CVODE return flag
    character(len=:),allocatable :: name  !! its name, or 'NONE'

    select case (flag)
    case (CV_SUCCESS);            name = 'CV_SUCCESS'
    case (CV_TSTOP_RETURN);       name = 'CV_TSTOP_RETURN'
    case (CV_ROOT_RETURN);        name = 'CV_ROOT_RETURN'
    case (CV_TOO_MUCH_WORK);      name = 'CV_TOO_MUCH_WORK'
    case (CV_TOO_MUCH_ACC);       name = 'CV_TOO_MUCH_ACC'
    case (CV_ERR_FAILURE);        name = 'CV_ERR_FAILURE'
    case (CV_CONV_FAILURE);       name = 'CV_CONV_FAILURE'
    case (CV_LINIT_FAIL);         name = 'CV_LINIT_FAIL'
    case (CV_LSETUP_FAIL);        name = 'CV_LSETUP_FAIL'
    case (CV_LSOLVE_FAIL);        name = 'CV_LSOLVE_FAIL'
    case (CV_RHSFUNC_FAIL);       name = 'CV_RHSFUNC_FAIL'
    case (CV_FIRST_RHSFUNC_ERR);  name = 'CV_FIRST_RHSFUNC_ERR'
    case (CV_REPTD_RHSFUNC_ERR);  name = 'CV_REPTD_RHSFUNC_ERR'
    case (CV_UNREC_RHSFUNC_ERR);  name = 'CV_UNREC_RHSFUNC_ERR'
    case (CV_RTFUNC_FAIL);        name = 'CV_RTFUNC_FAIL'
    case (CV_MEM_FAIL);           name = 'CV_MEM_FAIL'
    case (CV_MEM_NULL);           name = 'CV_MEM_NULL'
    case (CV_ILL_INPUT);          name = 'CV_ILL_INPUT'
    case (CV_NO_MALLOC);          name = 'CV_NO_MALLOC'
    case (CV_BAD_K);              name = 'CV_BAD_K'
    case (CV_BAD_T);              name = 'CV_BAD_T'
    case (CV_BAD_DKY);            name = 'CV_BAD_DKY'
    case (CV_TOO_CLOSE);          name = 'CV_TOO_CLOSE'
    case (CV_NLS_INIT_FAIL);      name = 'CV_NLS_INIT_FAIL'
    case (CV_NLS_SETUP_FAIL);     name = 'CV_NLS_SETUP_FAIL'
    case (CV_NLS_FAIL);           name = 'CV_NLS_FAIL'
    case (CV_PROJ_MEM_NULL);      name = 'CV_PROJ_MEM_NULL'
    case (CV_PROJFUNC_FAIL);      name = 'CV_PROJFUNC_FAIL'
    case (CV_REPTD_PROJFUNC_ERR); name = 'CV_REPTD_PROJFUNC_ERR'
    case default;                 name = 'NONE'
    end select

    end function CVodeGetReturnFlagName
!*****************************************************************************************

!=========================================================================================
!  Message formatting helpers (SUN_FORMAT_G is "%.15g")
!=========================================================================================

    function g_str(x) result(s)
    !! A real formatted like C's `%.15g` (`SUN_FORMAT_G` in double precision):
    !! 15 significant digits, fixed notation if the decimal exponent is in
    !! [-4, 15), otherwise exponential notation, trailing zeros removed.
    real(wp),intent(in) :: x  !! value to format
    character(len=:),allocatable :: s  !! formatted value
    integer,parameter :: p = 15  !! number of significant digits
    character(len=64) :: buf  !! write buffer
    character(len=16) :: fmt  !! run-time format for the fixed notation
    character(len=:),allocatable :: mant  !! mantissa digits (exponential notation)
    integer :: ie  !! position of the exponent letter in buf
    integer :: iexp  !! decimal exponent (after rounding to p digits)
    if (x /= x) then
        s = 'nan'
        return
    else if (abs(x) > huge(x)) then
        s = merge('-inf', 'inf ', x < 0.0_wp)
        s = trim(s)
        return
    else if (x == 0.0_wp) then
        s = merge('-0', '0 ', sign(1.0_wp, x) < 0.0_wp)
        s = trim(s)
        return
    end if
    ! exponent after rounding to p significant digits
    write(buf, '(es30.14e3)') x
    buf = adjustl(buf)
    ie = index(buf, 'E')
    read(buf(ie+1:), *) iexp
    if (iexp < -4 .or. iexp >= p) then
        mant = strip_zeros(buf(1:ie-1))
        write(buf, '(a,"e",sp,i0.2)') mant, iexp
        s = trim(buf)
    else
        write(fmt, '(a,i0,a)') '(f0.', p - 1 - iexp, ')'
        write(buf, fmt) x
        buf = adjustl(buf)
        if (buf(1:1) == '.') buf = '0'//buf
        if (buf(1:2) == '-.') buf = '-0'//buf(2:)
        s = strip_zeros(trim(buf))
    end if
    end function g_str

    function e_str(x) result(s)
    !! A real formatted like C's `% .15e` (`SUN_FORMAT_E` in double precision):
    !! a space or minus sign, 16 significant digits, a lowercase `e` and an
    !! exponent of at least two digits.
    real(wp),intent(in) :: x  !! value to format
    character(len=:),allocatable :: s  !! formatted value
    character(len=64) :: buf  !! write buffer
    integer :: ie  !! position of the exponent letter in buf
    integer :: iexp  !! decimal exponent
    if (x /= x .or. abs(x) > huge(x)) then
        s = ' '//g_str(x)
        return
    end if
    write(buf, '(es30.15e3)') x
    buf = adjustl(buf)
    ie = index(buf, 'E')
    read(buf(ie+1:), *) iexp
    write(buf, '(a,"e",sp,i0.2)') buf(1:ie-1), iexp
    s = trim(buf)
    if (s(1:1) /= '-') s = ' '//s
    end function e_str

    pure function strip_zeros(str) result(s)
    !! Remove trailing zeros after the decimal point, and the point itself if
    !! nothing follows it (C `%g` behavior).
    character(len=*),intent(in) :: str  !! number with a decimal point
    character(len=:),allocatable :: s  !! the stripped number
    integer :: n  !! length of the result
    s = trim(str)
    if (index(s, '.') == 0) return
    n = len(s)
    do while (n > 0)
        if (s(n:n) /= '0') exit
        n = n - 1
    end do
    if (s(n:n) == '.') n = n - 1
    s = s(1:n)
    end function strip_zeros

    function t_str(t) result(s)
    !! `MSG_TIME`: "t = <t>"
    real(wp),intent(in) :: t  !! time value
    character(len=:),allocatable :: s  !! "t = <t>"
    s = 't = '//g_str(t)
    end function t_str

    function th_str(t, h) result(s)
    !! `MSG_TIME_H`: "t = <t> and h = <h>"
    real(wp),intent(in) :: t  !! time value
    real(wp),intent(in) :: h  !! step size
    character(len=:),allocatable :: s  !! "t = <t> and h = <h>"
    s = 't = '//g_str(t)//' and h = '//g_str(h)
    end function th_str

!*****************************************************************************************
    end module fvode_core
!*****************************************************************************************

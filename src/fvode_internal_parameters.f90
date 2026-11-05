!*****************************************************************************************
!>
!  Internal constants of CVODE (not part of its public API).
!
!  Copied from SUNDIALS v7.9.0 with the C names and values (see
!  [[fvode_parameters]] for how they are verified). These are not re-exported
!  by the [[fvode]] umbrella module, so names like `ONE` and `HALF` do not
!  leak into user code.
!
!  Renamed: `TINY` (cvode.c) is `CV_TINY`, because `tiny` is a Fortran intrinsic.

    module fvode_internal_parameters

    use fvode_kinds, only: wp

    implicit none

    public


    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : method order limits
    !-------------------------------------------------------------------------------------
    integer,parameter :: ADAMS_Q_MAX = 12          !! max value of q for lmm == ADAMS
    integer,parameter :: BDF_Q_MAX   = 5           !! max value of q for lmm == BDF
    integer,parameter :: Q_MAX       = ADAMS_Q_MAX !! max value of q for either lmm
    integer,parameter :: L_MAX       = Q_MAX + 1   !! max value of L for either lmm
    integer,parameter :: NUM_TESTS   = 5           !! number of error test quantities

    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : defaults
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: HMIN_DEFAULT     = 0.0_wp !! hmin default value
    real(wp),parameter :: HMAX_INV_DEFAULT = 0.0_wp !! hmax_inv default value
    integer,parameter  :: MXHNIL_DEFAULT   = 10     !! mxhnil default value
    integer,parameter  :: MXSTEP_DEFAULT   = 500    !! mxstep default value

    integer,parameter  :: MSBP_DEFAULT         = 20     !! max steps between lsetup calls
    real(wp),parameter :: DGMAX_LSETUP_DEFAULT = 0.3_wp !! gamma threshold to call lsetup

    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : step size change constants
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: ETA_MIN_FX_DEFAULT = 0.0_wp
    real(wp),parameter :: ETA_MAX_FX_DEFAULT = 1.5_wp
    real(wp),parameter :: ETA_MAX_FS_DEFAULT = 10000.0_wp
    real(wp),parameter :: ETA_MAX_ES_DEFAULT = 10.0_wp
    real(wp),parameter :: ETA_MAX_GS_DEFAULT = 10.0_wp
    real(wp),parameter :: ETA_MIN_DEFAULT    = 0.1_wp
    real(wp),parameter :: ETA_MAX_EF_DEFAULT = 0.2_wp
    real(wp),parameter :: ETA_MIN_EF_DEFAULT = 0.1_wp
    real(wp),parameter :: ETA_CF_DEFAULT     = 0.25_wp
    integer,parameter  :: SMALL_NST_DEFAULT  = 10
    integer,parameter  :: SMALL_NEF_DEFAULT  = 2
    real(wp),parameter :: ONEPSM             = 1.000001_wp

    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : step size controller constants
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: ADDON = 0.000001_wp !! safety factor in computing eta
    real(wp),parameter :: BIAS1 = 6.0_wp      !! bias factors in eta selection
    real(wp),parameter :: BIAS2 = 6.0_wp
    real(wp),parameter :: BIAS3 = 10.0_wp

    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : order selection and failure limits
    !-------------------------------------------------------------------------------------
    integer,parameter :: LONG_WAIT            = 10
    integer,parameter :: MXNCF                = 10 !! max convergence failures per step try
    integer,parameter :: MXNEF                = 7  !! max error test failures per step try
    integer,parameter :: MXNEF1               = 3  !! error test failures before order reduction
    integer,parameter :: MAX_CONSTRAINT_FAILS = 10

    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : control constants for lower-level functions used by cvStep
    !-------------------------------------------------------------------------------------
    integer,parameter :: DO_ERROR_TEST = 2
    integer,parameter :: PREDICT_AGAIN = 3

    integer,parameter :: TRY_AGAIN      = 5
    integer,parameter :: FIRST_CALL     = 6
    integer,parameter :: PREV_CONV_FAIL = 7
    integer,parameter :: PREV_PROJ_FAIL = 8
    integer,parameter :: PREV_ERR_FAIL  = 9

    integer,parameter :: RHSFUNC_RECVR    = 10
    integer,parameter :: CONSTRFUNC_RECVR = 11
    integer,parameter :: PROJFUNC_RECVR   = 12

    !-------------------------------------------------------------------------------------
    ! cvode_impl.h : convfail values passed to lsetup
    !-------------------------------------------------------------------------------------
    integer,parameter :: CV_NO_FAILURES = 0
    integer,parameter :: CV_FAIL_BAD_J  = 1
    integer,parameter :: CV_FAIL_OTHER  = 2

    !-------------------------------------------------------------------------------------
    ! cvode.c : real constants
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: ZERO    = 0.0_wp
    real(wp),parameter :: CV_TINY = 1.0e-10_wp !! `TINY` in cvode.c (renamed: clashes with intrinsic)
    real(wp),parameter :: PT1     = 0.1_wp
    real(wp),parameter :: POINT2  = 0.2_wp
    real(wp),parameter :: FOURTH  = 0.25_wp
    real(wp),parameter :: HALF    = 0.5_wp
    real(wp),parameter :: PT9     = 0.9_wp
    real(wp),parameter :: ONE     = 1.0_wp
    real(wp),parameter :: ONEPT5  = 1.5_wp
    real(wp),parameter :: TWO     = 2.0_wp
    real(wp),parameter :: THREE   = 3.0_wp
    real(wp),parameter :: FOUR    = 4.0_wp
    real(wp),parameter :: FIVE    = 5.0_wp
    real(wp),parameter :: TWELVE  = 12.0_wp
    real(wp),parameter :: HUNDRED = 100.0_wp

    !-------------------------------------------------------------------------------------
    ! cvode.c : root finding, tolerance types, initial step
    !-------------------------------------------------------------------------------------
    integer,parameter :: RTFOUND = 1
    integer,parameter :: CLOSERT = 3

    integer,parameter :: CV_NN = 0 !! no tolerances specified
    integer,parameter :: CV_SS = 1 !! scalar reltol, scalar abstol
    integer,parameter :: CV_SV = 2 !! scalar reltol, vector abstol
    integer,parameter :: CV_WF = 3 !! user-supplied ewt function

    real(wp),parameter :: FUZZ_FACTOR = 100.0_wp
    real(wp),parameter :: HLB_FACTOR  = 100.0_wp
    real(wp),parameter :: HUB_FACTOR  = 0.1_wp
    real(wp),parameter :: H_BIAS      = HALF
    integer,parameter  :: MAX_ITERS   = 4

    real(wp),parameter :: CORTES = 0.1_wp !! default nonlinear convergence coefficient

    !-------------------------------------------------------------------------------------
    ! cvode_nls.c
    !-------------------------------------------------------------------------------------
    integer,parameter  :: NLS_MAXCOR = 3      !! max corrector iterations
    real(wp),parameter :: CRDOWN     = 0.3_wp !! constant used in convergence rate estimate
    real(wp),parameter :: RDIV       = 2.0_wp !! declare divergence if ratio del/delp > RDIV

    !-------------------------------------------------------------------------------------
    ! cvode_ls_impl.h, cvode_ls.c, cvode_bandpre.c
    !-------------------------------------------------------------------------------------
    integer,parameter  :: CVLS_MSBJ    = 51       !! max steps between Jacobian evaluations
    real(wp),parameter :: CVLS_DGMAX   = 0.2_wp   !! max |gamma/gammap - 1| for Jacobian reuse
    real(wp),parameter :: CVLS_EPLIN   = 0.05_wp  !! default eplifac
    real(wp),parameter :: MIN_INC_MULT = 1000.0_wp
    integer,parameter  :: MAX_DQITERS  = 3        !! max attempts to recover in DQ J*v
    real(wp),parameter :: PT25         = 0.25_wp

    !-------------------------------------------------------------------------------------
    ! cvode_diag.c
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: FRACT = 0.1_wp

    !-------------------------------------------------------------------------------------
    ! cvode_proj_impl.h
    !-------------------------------------------------------------------------------------
    integer,parameter  :: PROJ_MAX_FAILS = 10
    real(wp),parameter :: PROJ_EPS       = 0.1_wp
    real(wp),parameter :: PROJ_FAIL_ETA  = 0.25_wp

    end module fvode_internal_parameters
!*****************************************************************************************

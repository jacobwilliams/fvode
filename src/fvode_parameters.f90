!*****************************************************************************************
!>
!  Constants used by FVODE.
!
!  Public CVODE/SUNDIALS constants (method choices, task modes, return codes).
!  Internal constants are in [[fvode_internal_parameters]].
!
!  Every value is copied from the SUNDIALS v7.9.0 headers/sources, keeping
!  the C names and integer values. CVODE's control flow branches on these codes,
!  so they must not change. `test/test_parameters.f90` checks them against a
!  table generated from the C headers (`reference/constants.c`).
!
!  Sources:
!
!  * `include/cvode/cvode.h`, `cvode_ls.h`, `cvode_diag.h`
!  * `include/sundials/sundials_errors.h`, `sundials_linearsolver.h`,
!    `sundials_nonlinearsolver.h`, `sundials_iterative.h`, `sundials_types.h`,
!    `include/sunlinsol/sunlinsol_sp*.h`, `sunlinsol_pcg.h`

    module fvode_parameters

    use fvode_kinds, only: wp

    implicit none

    public

    !-------------------------------------------------------------------------------------
    ! cvode.h : lmm, itask
    !-------------------------------------------------------------------------------------
    integer,parameter :: CV_ADAMS = 1
    integer,parameter :: CV_BDF   = 2

    integer,parameter :: CV_NORMAL   = 1
    integer,parameter :: CV_ONE_STEP = 2

    !-------------------------------------------------------------------------------------
    ! cvode.h : return values
    !-------------------------------------------------------------------------------------
    integer,parameter :: CV_SUCCESS      = 0
    integer,parameter :: CV_TSTOP_RETURN = 1
    integer,parameter :: CV_ROOT_RETURN  = 2

    integer,parameter :: CV_WARNING = 99

    integer,parameter :: CV_TOO_MUCH_WORK = -1
    integer,parameter :: CV_TOO_MUCH_ACC  = -2
    integer,parameter :: CV_ERR_FAILURE   = -3
    integer,parameter :: CV_CONV_FAILURE  = -4

    integer,parameter :: CV_LINIT_FAIL        = -5
    integer,parameter :: CV_LSETUP_FAIL       = -6
    integer,parameter :: CV_LSOLVE_FAIL       = -7
    integer,parameter :: CV_RHSFUNC_FAIL      = -8
    integer,parameter :: CV_FIRST_RHSFUNC_ERR = -9
    integer,parameter :: CV_REPTD_RHSFUNC_ERR = -10
    integer,parameter :: CV_UNREC_RHSFUNC_ERR = -11
    integer,parameter :: CV_RTFUNC_FAIL       = -12
    integer,parameter :: CV_NLS_INIT_FAIL     = -13
    integer,parameter :: CV_NLS_SETUP_FAIL    = -14
    integer,parameter :: CV_CONSTR_FAIL       = -15
    integer,parameter :: CV_NLS_FAIL          = -16

    integer,parameter :: CV_MEM_FAIL     = -20
    integer,parameter :: CV_MEM_NULL     = -21
    integer,parameter :: CV_ILL_INPUT    = -22
    integer,parameter :: CV_NO_MALLOC    = -23
    integer,parameter :: CV_BAD_K        = -24
    integer,parameter :: CV_BAD_T        = -25
    integer,parameter :: CV_BAD_DKY      = -26
    integer,parameter :: CV_TOO_CLOSE    = -27
    integer,parameter :: CV_VECTOROP_ERR = -28

    integer,parameter :: CV_PROJ_MEM_NULL      = -29
    integer,parameter :: CV_PROJFUNC_FAIL      = -30
    integer,parameter :: CV_REPTD_PROJFUNC_ERR = -31

    integer,parameter :: CV_CONTEXT_ERR = -32

    integer,parameter :: CV_UNRECOGNIZED_ERR = -99

    !-------------------------------------------------------------------------------------
    ! cvode_ls.h : linear solver interface return values
    !-------------------------------------------------------------------------------------
    integer,parameter :: CVLS_SUCCESS         = 0
    integer,parameter :: CVLS_MEM_NULL        = -1
    integer,parameter :: CVLS_LMEM_NULL       = -2
    integer,parameter :: CVLS_ILL_INPUT       = -3
    integer,parameter :: CVLS_MEM_FAIL        = -4
    integer,parameter :: CVLS_PMEM_NULL       = -5
    integer,parameter :: CVLS_JACFUNC_UNRECVR = -6
    integer,parameter :: CVLS_JACFUNC_RECVR   = -7
    integer,parameter :: CVLS_SUNMAT_FAIL     = -8
    integer,parameter :: CVLS_SUNLS_FAIL      = -9

    !-------------------------------------------------------------------------------------
    ! cvode_diag.h : diagonal linear solver return values
    !-------------------------------------------------------------------------------------
    integer,parameter :: CVDIAG_SUCCESS   = 0
    integer,parameter :: CVDIAG_MEM_NULL  = -1
    integer,parameter :: CVDIAG_LMEM_NULL = -2
    integer,parameter :: CVDIAG_ILL_INPUT = -3
    integer,parameter :: CVDIAG_MEM_FAIL  = -4

    integer,parameter :: CVDIAG_INV_FAIL        = -5
    integer,parameter :: CVDIAG_RHSFUNC_UNRECVR = -6
    integer,parameter :: CVDIAG_RHSFUNC_RECVR   = -7

    !-------------------------------------------------------------------------------------
    ! sundials_errors.h : generic SUNErrCode values (enum starting at -10000)
    !-------------------------------------------------------------------------------------
    integer,parameter :: SUN_SUCCESS              = 0
    integer,parameter :: SUN_ERR_MINIMUM          = -10000
    integer,parameter :: SUN_ERR_ARG_CORRUPT      = -9999
    integer,parameter :: SUN_ERR_ARG_INCOMPATIBLE = -9998
    integer,parameter :: SUN_ERR_ARG_OUTOFRANGE   = -9997
    integer,parameter :: SUN_ERR_ARG_WRONGTYPE    = -9996
    integer,parameter :: SUN_ERR_ARG_DIMSMISMATCH = -9995
    integer,parameter :: SUN_ERR_GENERIC          = -9994
    integer,parameter :: SUN_ERR_CORRUPT          = -9993
    integer,parameter :: SUN_ERR_OUTOFRANGE       = -9992
    integer,parameter :: SUN_ERR_FILE_OPEN        = -9991
    integer,parameter :: SUN_ERR_OP_FAIL          = -9990
    integer,parameter :: SUN_ERR_MEM_FAIL         = -9989
    integer,parameter :: SUN_ERR_MALLOC_FAIL      = -9988
    integer,parameter :: SUN_ERR_EXT_FAIL         = -9987
    integer,parameter :: SUN_ERR_DESTROY_FAIL     = -9986
    integer,parameter :: SUN_ERR_NOT_IMPLEMENTED  = -9985
    integer,parameter :: SUN_ERR_USER_FCN_FAIL    = -9984
    integer,parameter :: SUN_ERR_UNKNOWN          = -9970
    integer,parameter :: SUN_ERR_MAXIMUM          = -1000

    !-------------------------------------------------------------------------------------
    ! sundials_linearsolver.h
    !-------------------------------------------------------------------------------------
    ! SUNLinearSolver_Type
    integer,parameter :: SUNLINEARSOLVER_DIRECT           = 0
    integer,parameter :: SUNLINEARSOLVER_ITERATIVE        = 1
    integer,parameter :: SUNLINEARSOLVER_MATRIX_ITERATIVE = 2
    integer,parameter :: SUNLINEARSOLVER_MATRIX_EMBEDDED  = 3

    ! unrecoverable failure codes
    integer,parameter :: SUNLS_ATIMES_NULL       = -804
    integer,parameter :: SUNLS_ATIMES_FAIL_UNREC = -805
    integer,parameter :: SUNLS_PSET_FAIL_UNREC   = -806
    integer,parameter :: SUNLS_PSOLVE_NULL       = -807
    integer,parameter :: SUNLS_PSOLVE_FAIL_UNREC = -808
    integer,parameter :: SUNLS_GS_FAIL           = -810
    integer,parameter :: SUNLS_QRSOL_FAIL        = -811

    ! recoverable failure codes
    integer,parameter :: SUNLS_RECOV_FAILURE    = 800
    integer,parameter :: SUNLS_RES_REDUCED      = 801
    integer,parameter :: SUNLS_CONV_FAIL        = 802
    integer,parameter :: SUNLS_ATIMES_FAIL_REC  = 803
    integer,parameter :: SUNLS_PSET_FAIL_REC    = 804
    integer,parameter :: SUNLS_PSOLVE_FAIL_REC  = 805
    integer,parameter :: SUNLS_PACKAGE_FAIL_REC = 806
    integer,parameter :: SUNLS_QRFACT_FAIL      = 807
    integer,parameter :: SUNLS_LUFACT_FAIL      = 808

    !-------------------------------------------------------------------------------------
    ! sundials_nonlinearsolver.h
    !-------------------------------------------------------------------------------------
    ! SUNNonlinearSolver_Type
    integer,parameter :: SUNNONLINEARSOLVER_ROOTFIND   = 0
    integer,parameter :: SUNNONLINEARSOLVER_FIXEDPOINT = 1
    integer,parameter :: SUNNONLINEARSOLVER_HYBRID     = 2

    integer,parameter :: SUN_NLS_CONTINUE   = 901
    integer,parameter :: SUN_NLS_CONV_RECVR = 902
    integer,parameter :: SUN_NLS_SWITCH     = 903

    !-------------------------------------------------------------------------------------
    ! sundials_iterative.h
    !-------------------------------------------------------------------------------------
    ! SUNPrecType
    integer,parameter :: SUN_PREC_NONE  = 0
    integer,parameter :: SUN_PREC_LEFT  = 1
    integer,parameter :: SUN_PREC_RIGHT = 2
    integer,parameter :: SUN_PREC_BOTH  = 3

    ! SUNGramSchmidtType
    integer,parameter :: SUN_MODIFIED_GS  = 1
    integer,parameter :: SUN_CLASSICAL_GS = 2

    !-------------------------------------------------------------------------------------
    ! sunlinsol_sp*.h, sunlinsol_pcg.h : Krylov defaults
    !-------------------------------------------------------------------------------------
    integer,parameter :: SUNSPGMR_MAXL_DEFAULT    = 5
    integer,parameter :: SUNSPGMR_MAXRS_DEFAULT   = 0
    integer,parameter :: SUNSPGMR_GSTYPE_DEFAULT  = SUN_MODIFIED_GS
    integer,parameter :: SUNSPFGMR_MAXL_DEFAULT   = 5
    integer,parameter :: SUNSPFGMR_MAXRS_DEFAULT  = 0
    integer,parameter :: SUNSPFGMR_GSTYPE_DEFAULT = SUN_MODIFIED_GS
    integer,parameter :: SUNSPBCGS_MAXL_DEFAULT   = 5
    integer,parameter :: SUNSPTFQMR_MAXL_DEFAULT  = 5
    integer,parameter :: SUNPCG_MAXL_DEFAULT      = 5

    !-------------------------------------------------------------------------------------
    ! sundials_types.h : SUNOutputFormat
    !-------------------------------------------------------------------------------------
    integer,parameter :: SUN_OUTPUTFORMAT_TABLE = 0
    integer,parameter :: SUN_OUTPUTFORMAT_CSV   = 1

    !-------------------------------------------------------------------------------------
    ! sundials_types.h : floating point limits
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: SUN_UNIT_ROUNDOFF = epsilon(1.0_wp) !! `DBL_EPSILON`
    real(wp),parameter :: SUN_BIG_REAL      = huge(1.0_wp)    !! `DBL_MAX`
    real(wp),parameter :: SUN_SMALL_REAL    = tiny(1.0_wp)    !! `DBL_MIN`

    end module fvode_parameters
!*****************************************************************************************

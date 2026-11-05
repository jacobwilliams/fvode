!*****************************************************************************************
!>
!  Checks every constant in [[fvode_parameters]] and
!  [[fvode_internal_parameters]] against the values in the
!  SUNDIALS v7.9.0 headers and sources.
!
!  `test/data/sundials_constants.txt` is generated from the C code by
!  `pixi run -e reference ref-run` (see `reference/constants.c` and
!  `reference/parse_defines.py`). The calls below were generated from that file,
!  so every upstream constant is checked.

    module test_parameters

    use fvode_kinds, only: wp
    use fvode_parameters
    use fvode_internal_parameters
    use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

    implicit none

    private

    public :: collect_parameters

    character(len=*),parameter :: table_file = 'test/data/sundials_constants.txt'

    type :: entry_type
        character(len=:),allocatable :: name
        character(len=:),allocatable :: value
        logical :: used = .false.
    end type entry_type

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Collect all tests in this suite.

    subroutine collect_parameters(testsuite)

    type(unittest_type),allocatable,intent(out) :: testsuite(:)

    testsuite = [ new_unittest("constants_match_sundials", test_constants) ]

    end subroutine collect_parameters
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the reference table and compare each constant.

    subroutine test_constants(error)

    type(error_type),allocatable,intent(out) :: error

    type(entry_type),allocatable :: table(:)
    integer :: i

    call read_table(error, table); if (allocated(error)) return

    ! integer constants
    call check_int(error, table, 'ADAMS_Q_MAX', ADAMS_Q_MAX); if (allocated(error)) return
    call check_int(error, table, 'BDF_Q_MAX', BDF_Q_MAX); if (allocated(error)) return
    call check_int(error, table, 'CLOSERT', CLOSERT); if (allocated(error)) return
    call check_int(error, table, 'CONSTRFUNC_RECVR', CONSTRFUNC_RECVR); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_ILL_INPUT', CVDIAG_ILL_INPUT); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_INV_FAIL', CVDIAG_INV_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_LMEM_NULL', CVDIAG_LMEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_MEM_FAIL', CVDIAG_MEM_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_MEM_NULL', CVDIAG_MEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_RHSFUNC_RECVR', CVDIAG_RHSFUNC_RECVR); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_RHSFUNC_UNRECVR', CVDIAG_RHSFUNC_UNRECVR); if (allocated(error)) return
    call check_int(error, table, 'CVDIAG_SUCCESS', CVDIAG_SUCCESS); if (allocated(error)) return
    call check_int(error, table, 'CVLS_ILL_INPUT', CVLS_ILL_INPUT); if (allocated(error)) return
    call check_int(error, table, 'CVLS_JACFUNC_RECVR', CVLS_JACFUNC_RECVR); if (allocated(error)) return
    call check_int(error, table, 'CVLS_JACFUNC_UNRECVR', CVLS_JACFUNC_UNRECVR); if (allocated(error)) return
    call check_int(error, table, 'CVLS_LMEM_NULL', CVLS_LMEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CVLS_MEM_FAIL', CVLS_MEM_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CVLS_MEM_NULL', CVLS_MEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CVLS_MSBJ', CVLS_MSBJ); if (allocated(error)) return
    call check_int(error, table, 'CVLS_PMEM_NULL', CVLS_PMEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CVLS_SUCCESS', CVLS_SUCCESS); if (allocated(error)) return
    call check_int(error, table, 'CVLS_SUNLS_FAIL', CVLS_SUNLS_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CVLS_SUNMAT_FAIL', CVLS_SUNMAT_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_ADAMS', CV_ADAMS); if (allocated(error)) return
    call check_int(error, table, 'CV_BAD_DKY', CV_BAD_DKY); if (allocated(error)) return
    call check_int(error, table, 'CV_BAD_K', CV_BAD_K); if (allocated(error)) return
    call check_int(error, table, 'CV_BAD_T', CV_BAD_T); if (allocated(error)) return
    call check_int(error, table, 'CV_BDF', CV_BDF); if (allocated(error)) return
    call check_int(error, table, 'CV_CONSTR_FAIL', CV_CONSTR_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_CONTEXT_ERR', CV_CONTEXT_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_CONV_FAILURE', CV_CONV_FAILURE); if (allocated(error)) return
    call check_int(error, table, 'CV_ERR_FAILURE', CV_ERR_FAILURE); if (allocated(error)) return
    call check_int(error, table, 'CV_FAIL_BAD_J', CV_FAIL_BAD_J); if (allocated(error)) return
    call check_int(error, table, 'CV_FAIL_OTHER', CV_FAIL_OTHER); if (allocated(error)) return
    call check_int(error, table, 'CV_FIRST_RHSFUNC_ERR', CV_FIRST_RHSFUNC_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_ILL_INPUT', CV_ILL_INPUT); if (allocated(error)) return
    call check_int(error, table, 'CV_LINIT_FAIL', CV_LINIT_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_LSETUP_FAIL', CV_LSETUP_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_LSOLVE_FAIL', CV_LSOLVE_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_MEM_FAIL', CV_MEM_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_MEM_NULL', CV_MEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CV_NLS_FAIL', CV_NLS_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_NLS_INIT_FAIL', CV_NLS_INIT_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_NLS_SETUP_FAIL', CV_NLS_SETUP_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_NN', CV_NN); if (allocated(error)) return
    call check_int(error, table, 'CV_NORMAL', CV_NORMAL); if (allocated(error)) return
    call check_int(error, table, 'CV_NO_FAILURES', CV_NO_FAILURES); if (allocated(error)) return
    call check_int(error, table, 'CV_NO_MALLOC', CV_NO_MALLOC); if (allocated(error)) return
    call check_int(error, table, 'CV_ONE_STEP', CV_ONE_STEP); if (allocated(error)) return
    call check_int(error, table, 'CV_PROJFUNC_FAIL', CV_PROJFUNC_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_PROJ_MEM_NULL', CV_PROJ_MEM_NULL); if (allocated(error)) return
    call check_int(error, table, 'CV_REPTD_PROJFUNC_ERR', CV_REPTD_PROJFUNC_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_REPTD_RHSFUNC_ERR', CV_REPTD_RHSFUNC_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_RHSFUNC_FAIL', CV_RHSFUNC_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_ROOT_RETURN', CV_ROOT_RETURN); if (allocated(error)) return
    call check_int(error, table, 'CV_RTFUNC_FAIL', CV_RTFUNC_FAIL); if (allocated(error)) return
    call check_int(error, table, 'CV_SS', CV_SS); if (allocated(error)) return
    call check_int(error, table, 'CV_SUCCESS', CV_SUCCESS); if (allocated(error)) return
    call check_int(error, table, 'CV_SV', CV_SV); if (allocated(error)) return
    call check_int(error, table, 'CV_TOO_CLOSE', CV_TOO_CLOSE); if (allocated(error)) return
    call check_int(error, table, 'CV_TOO_MUCH_ACC', CV_TOO_MUCH_ACC); if (allocated(error)) return
    call check_int(error, table, 'CV_TOO_MUCH_WORK', CV_TOO_MUCH_WORK); if (allocated(error)) return
    call check_int(error, table, 'CV_TSTOP_RETURN', CV_TSTOP_RETURN); if (allocated(error)) return
    call check_int(error, table, 'CV_UNRECOGNIZED_ERR', CV_UNRECOGNIZED_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_UNREC_RHSFUNC_ERR', CV_UNREC_RHSFUNC_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_VECTOROP_ERR', CV_VECTOROP_ERR); if (allocated(error)) return
    call check_int(error, table, 'CV_WARNING', CV_WARNING); if (allocated(error)) return
    call check_int(error, table, 'CV_WF', CV_WF); if (allocated(error)) return
    call check_int(error, table, 'DO_ERROR_TEST', DO_ERROR_TEST); if (allocated(error)) return
    call check_int(error, table, 'FIRST_CALL', FIRST_CALL); if (allocated(error)) return
    call check_int(error, table, 'LONG_WAIT', LONG_WAIT); if (allocated(error)) return
    call check_int(error, table, 'L_MAX', L_MAX); if (allocated(error)) return
    call check_int(error, table, 'MAX_CONSTRAINT_FAILS', MAX_CONSTRAINT_FAILS); if (allocated(error)) return
    call check_int(error, table, 'MAX_DQITERS', MAX_DQITERS); if (allocated(error)) return
    call check_int(error, table, 'MAX_ITERS', MAX_ITERS); if (allocated(error)) return
    call check_int(error, table, 'MSBP_DEFAULT', MSBP_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'MXHNIL_DEFAULT', MXHNIL_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'MXNCF', MXNCF); if (allocated(error)) return
    call check_int(error, table, 'MXNEF', MXNEF); if (allocated(error)) return
    call check_int(error, table, 'MXNEF1', MXNEF1); if (allocated(error)) return
    call check_int(error, table, 'MXSTEP_DEFAULT', MXSTEP_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'NLS_MAXCOR', NLS_MAXCOR); if (allocated(error)) return
    call check_int(error, table, 'NUM_TESTS', NUM_TESTS); if (allocated(error)) return
    call check_int(error, table, 'PREDICT_AGAIN', PREDICT_AGAIN); if (allocated(error)) return
    call check_int(error, table, 'PREV_CONV_FAIL', PREV_CONV_FAIL); if (allocated(error)) return
    call check_int(error, table, 'PREV_ERR_FAIL', PREV_ERR_FAIL); if (allocated(error)) return
    call check_int(error, table, 'PREV_PROJ_FAIL', PREV_PROJ_FAIL); if (allocated(error)) return
    call check_int(error, table, 'PROJFUNC_RECVR', PROJFUNC_RECVR); if (allocated(error)) return
    call check_int(error, table, 'PROJ_MAX_FAILS', PROJ_MAX_FAILS); if (allocated(error)) return
    call check_int(error, table, 'Q_MAX', Q_MAX); if (allocated(error)) return
    call check_int(error, table, 'RHSFUNC_RECVR', RHSFUNC_RECVR); if (allocated(error)) return
    call check_int(error, table, 'RTFOUND', RTFOUND); if (allocated(error)) return
    call check_int(error, table, 'SMALL_NEF_DEFAULT', SMALL_NEF_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SMALL_NST_DEFAULT', SMALL_NST_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNLINEARSOLVER_DIRECT', SUNLINEARSOLVER_DIRECT); if (allocated(error)) return
    call check_int(error, table, 'SUNLINEARSOLVER_ITERATIVE', SUNLINEARSOLVER_ITERATIVE); if (allocated(error)) return
    call check_int(error, table, 'SUNLINEARSOLVER_MATRIX_EMBEDDED', SUNLINEARSOLVER_MATRIX_EMBEDDED); if (allocated(error)) return
    call check_int(error, table, 'SUNLINEARSOLVER_MATRIX_ITERATIVE', SUNLINEARSOLVER_MATRIX_ITERATIVE); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_ATIMES_FAIL_REC', SUNLS_ATIMES_FAIL_REC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_ATIMES_FAIL_UNREC', SUNLS_ATIMES_FAIL_UNREC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_ATIMES_NULL', SUNLS_ATIMES_NULL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_CONV_FAIL', SUNLS_CONV_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_GS_FAIL', SUNLS_GS_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_LUFACT_FAIL', SUNLS_LUFACT_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_PACKAGE_FAIL_REC', SUNLS_PACKAGE_FAIL_REC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_PSET_FAIL_REC', SUNLS_PSET_FAIL_REC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_PSET_FAIL_UNREC', SUNLS_PSET_FAIL_UNREC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_PSOLVE_FAIL_REC', SUNLS_PSOLVE_FAIL_REC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_PSOLVE_FAIL_UNREC', SUNLS_PSOLVE_FAIL_UNREC); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_PSOLVE_NULL', SUNLS_PSOLVE_NULL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_QRFACT_FAIL', SUNLS_QRFACT_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_QRSOL_FAIL', SUNLS_QRSOL_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_RECOV_FAILURE', SUNLS_RECOV_FAILURE); if (allocated(error)) return
    call check_int(error, table, 'SUNLS_RES_REDUCED', SUNLS_RES_REDUCED); if (allocated(error)) return
    call check_int(error, table, 'SUNNONLINEARSOLVER_FIXEDPOINT', SUNNONLINEARSOLVER_FIXEDPOINT); if (allocated(error)) return
    call check_int(error, table, 'SUNNONLINEARSOLVER_HYBRID', SUNNONLINEARSOLVER_HYBRID); if (allocated(error)) return
    call check_int(error, table, 'SUNNONLINEARSOLVER_ROOTFIND', SUNNONLINEARSOLVER_ROOTFIND); if (allocated(error)) return
    call check_int(error, table, 'SUNPCG_MAXL_DEFAULT', SUNPCG_MAXL_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPBCGS_MAXL_DEFAULT', SUNSPBCGS_MAXL_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPFGMR_GSTYPE_DEFAULT', SUNSPFGMR_GSTYPE_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPFGMR_MAXL_DEFAULT', SUNSPFGMR_MAXL_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPFGMR_MAXRS_DEFAULT', SUNSPFGMR_MAXRS_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPGMR_GSTYPE_DEFAULT', SUNSPGMR_GSTYPE_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPGMR_MAXL_DEFAULT', SUNSPGMR_MAXL_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPGMR_MAXRS_DEFAULT', SUNSPGMR_MAXRS_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUNSPTFQMR_MAXL_DEFAULT', SUNSPTFQMR_MAXL_DEFAULT); if (allocated(error)) return
    call check_int(error, table, 'SUN_CLASSICAL_GS', SUN_CLASSICAL_GS); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_ARG_CORRUPT', SUN_ERR_ARG_CORRUPT); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_ARG_DIMSMISMATCH', SUN_ERR_ARG_DIMSMISMATCH); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_ARG_INCOMPATIBLE', SUN_ERR_ARG_INCOMPATIBLE); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_ARG_OUTOFRANGE', SUN_ERR_ARG_OUTOFRANGE); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_ARG_WRONGTYPE', SUN_ERR_ARG_WRONGTYPE); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_CORRUPT', SUN_ERR_CORRUPT); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_DESTROY_FAIL', SUN_ERR_DESTROY_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_EXT_FAIL', SUN_ERR_EXT_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_FILE_OPEN', SUN_ERR_FILE_OPEN); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_GENERIC', SUN_ERR_GENERIC); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_MALLOC_FAIL', SUN_ERR_MALLOC_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_MAXIMUM', SUN_ERR_MAXIMUM); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_MEM_FAIL', SUN_ERR_MEM_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_MINIMUM', SUN_ERR_MINIMUM); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_NOT_IMPLEMENTED', SUN_ERR_NOT_IMPLEMENTED); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_OP_FAIL', SUN_ERR_OP_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_OUTOFRANGE', SUN_ERR_OUTOFRANGE); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_UNKNOWN', SUN_ERR_UNKNOWN); if (allocated(error)) return
    call check_int(error, table, 'SUN_ERR_USER_FCN_FAIL', SUN_ERR_USER_FCN_FAIL); if (allocated(error)) return
    call check_int(error, table, 'SUN_MODIFIED_GS', SUN_MODIFIED_GS); if (allocated(error)) return
    call check_int(error, table, 'SUN_NLS_CONTINUE', SUN_NLS_CONTINUE); if (allocated(error)) return
    call check_int(error, table, 'SUN_NLS_CONV_RECVR', SUN_NLS_CONV_RECVR); if (allocated(error)) return
    call check_int(error, table, 'SUN_NLS_SWITCH', SUN_NLS_SWITCH); if (allocated(error)) return
    call check_int(error, table, 'SUN_OUTPUTFORMAT_CSV', SUN_OUTPUTFORMAT_CSV); if (allocated(error)) return
    call check_int(error, table, 'SUN_OUTPUTFORMAT_TABLE', SUN_OUTPUTFORMAT_TABLE); if (allocated(error)) return
    call check_int(error, table, 'SUN_PREC_BOTH', SUN_PREC_BOTH); if (allocated(error)) return
    call check_int(error, table, 'SUN_PREC_LEFT', SUN_PREC_LEFT); if (allocated(error)) return
    call check_int(error, table, 'SUN_PREC_NONE', SUN_PREC_NONE); if (allocated(error)) return
    call check_int(error, table, 'SUN_PREC_RIGHT', SUN_PREC_RIGHT); if (allocated(error)) return
    call check_int(error, table, 'SUN_SUCCESS', SUN_SUCCESS); if (allocated(error)) return
    call check_int(error, table, 'TRY_AGAIN', TRY_AGAIN); if (allocated(error)) return

    ! real constants
    call check_real(error, table, 'ADDON', ADDON); if (allocated(error)) return
    call check_real(error, table, 'BIAS1', BIAS1); if (allocated(error)) return
    call check_real(error, table, 'BIAS2', BIAS2); if (allocated(error)) return
    call check_real(error, table, 'BIAS3', BIAS3); if (allocated(error)) return
    call check_real(error, table, 'CORTES', CORTES); if (allocated(error)) return
    call check_real(error, table, 'CRDOWN', CRDOWN); if (allocated(error)) return
    call check_real(error, table, 'CVLS_DGMAX', CVLS_DGMAX); if (allocated(error)) return
    call check_real(error, table, 'CVLS_EPLIN', CVLS_EPLIN); if (allocated(error)) return
    call check_real(error, table, 'DGMAX_LSETUP_DEFAULT', DGMAX_LSETUP_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_CF_DEFAULT', ETA_CF_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MAX_EF_DEFAULT', ETA_MAX_EF_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MAX_ES_DEFAULT', ETA_MAX_ES_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MAX_FS_DEFAULT', ETA_MAX_FS_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MAX_FX_DEFAULT', ETA_MAX_FX_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MAX_GS_DEFAULT', ETA_MAX_GS_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MIN_DEFAULT', ETA_MIN_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MIN_EF_DEFAULT', ETA_MIN_EF_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'ETA_MIN_FX_DEFAULT', ETA_MIN_FX_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'FIVE', FIVE); if (allocated(error)) return
    call check_real(error, table, 'FOUR', FOUR); if (allocated(error)) return
    call check_real(error, table, 'FOURTH', FOURTH); if (allocated(error)) return
    call check_real(error, table, 'FRACT', FRACT); if (allocated(error)) return
    call check_real(error, table, 'FUZZ_FACTOR', FUZZ_FACTOR); if (allocated(error)) return
    call check_real(error, table, 'HALF', HALF); if (allocated(error)) return
    call check_real(error, table, 'HLB_FACTOR', HLB_FACTOR); if (allocated(error)) return
    call check_real(error, table, 'HMAX_INV_DEFAULT', HMAX_INV_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'HMIN_DEFAULT', HMIN_DEFAULT); if (allocated(error)) return
    call check_real(error, table, 'HUB_FACTOR', HUB_FACTOR); if (allocated(error)) return
    call check_real(error, table, 'HUNDRED', HUNDRED); if (allocated(error)) return
    call check_real(error, table, 'H_BIAS', H_BIAS); if (allocated(error)) return
    call check_real(error, table, 'MIN_INC_MULT', MIN_INC_MULT); if (allocated(error)) return
    call check_real(error, table, 'ONE', ONE); if (allocated(error)) return
    call check_real(error, table, 'ONEPSM', ONEPSM); if (allocated(error)) return
    call check_real(error, table, 'ONEPT5', ONEPT5); if (allocated(error)) return
    call check_real(error, table, 'POINT2', POINT2); if (allocated(error)) return
    call check_real(error, table, 'PROJ_EPS', PROJ_EPS); if (allocated(error)) return
    call check_real(error, table, 'PROJ_FAIL_ETA', PROJ_FAIL_ETA); if (allocated(error)) return
    call check_real(error, table, 'PT1', PT1); if (allocated(error)) return
    call check_real(error, table, 'PT25', PT25); if (allocated(error)) return
    call check_real(error, table, 'PT9', PT9); if (allocated(error)) return
    call check_real(error, table, 'RDIV', RDIV); if (allocated(error)) return
    call check_real(error, table, 'SUN_BIG_REAL', SUN_BIG_REAL); if (allocated(error)) return
    call check_real(error, table, 'SUN_SMALL_REAL', SUN_SMALL_REAL); if (allocated(error)) return
    call check_real(error, table, 'SUN_UNIT_ROUNDOFF', SUN_UNIT_ROUNDOFF); if (allocated(error)) return
    call check_real(error, table, 'THREE', THREE); if (allocated(error)) return
    call check_real(error, table, 'TINY', CV_TINY); if (allocated(error)) return
    call check_real(error, table, 'TWELVE', TWELVE); if (allocated(error)) return
    call check_real(error, table, 'TWO', TWO); if (allocated(error)) return
    call check_real(error, table, 'ZERO', ZERO); if (allocated(error)) return

    ! every upstream constant must have been checked
    do i = 1, size(table)
        if (.not. table(i)%used) then
            call test_failed(error, 'constant not checked: '//table(i)%name)
            return
        end if
    end do

    end subroutine test_constants
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read `NAME VALUE` lines from the reference table.

    subroutine read_table(error, table)

    type(error_type),allocatable,intent(out) :: error
    type(entry_type),allocatable,intent(out) :: table(:)

    integer :: iunit, istat, n, i
    character(len=256) :: line, name, value

    open(newunit=iunit, file=table_file, status='old', action='read', iostat=istat)
    if (istat /= 0) then
        call test_failed(error, 'cannot open '//table_file)
        return
    end if

    allocate(table(0))
    do
        read(iunit, '(a)', iostat=istat) line
        if (istat /= 0) exit
        if (len_trim(line) == 0) cycle
        read(line, *) name, value
        table = [table, entry_type(trim(name), trim(value))]
    end do
    close(iunit)

    n = size(table)
    do i = 1, n
        if (len(table(i)%name) == 0) then
            call test_failed(error, 'bad line in '//table_file)
            return
        end if
    end do

    end subroutine read_table
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find a name in the table (marks it as used). Returns 0 if not found.

    function lookup(table, name) result(idx)

    type(entry_type),intent(inout) :: table(:)
    character(len=*),intent(in) :: name
    integer :: idx

    do idx = 1, size(table)
        if (table(idx)%name == name) then
            table(idx)%used = .true.
            return
        end if
    end do
    idx = 0

    end function lookup
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check an integer constant (must match exactly).

    subroutine check_int(error, table, name, val)

    type(error_type),allocatable,intent(out) :: error
    type(entry_type),intent(inout) :: table(:)
    character(len=*),intent(in) :: name
    integer,intent(in) :: val

    integer :: idx, ref, istat

    idx = lookup(table, name)
    if (idx == 0) then
        call test_failed(error, name//' missing from '//table_file)
        return
    end if
    read(table(idx)%value, *, iostat=istat) ref
    if (istat /= 0) then
        call test_failed(error, name//': reference value is not an integer: '//table(idx)%value)
        return
    end if
    call check(error, val, ref, name)

    end subroutine check_int
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check a real constant. Values are printed with 17 significant digits,
!  so they read back exactly: the comparison is bitwise.

    subroutine check_real(error, table, name, val)

    type(error_type),allocatable,intent(out) :: error
    type(entry_type),intent(inout) :: table(:)
    character(len=*),intent(in) :: name
    real(wp),intent(in) :: val

    integer :: idx, istat
    real(wp) :: ref

    idx = lookup(table, name)
    if (idx == 0) then
        call test_failed(error, name//' missing from '//table_file)
        return
    end if
    read(table(idx)%value, *, iostat=istat) ref
    if (istat /= 0) then
        call test_failed(error, name//': reference value is not a real: '//table(idx)%value)
        return
    end if
    if (val /= ref) call test_failed(error, name//' differs from SUNDIALS value '//table(idx)%value)

    end subroutine check_real
!*****************************************************************************************

    end module test_parameters
!*****************************************************************************************

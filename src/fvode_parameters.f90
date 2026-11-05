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
    integer,parameter :: CV_ADAMS = 1 !! Adams-Moulton methods (nonstiff problems)
    integer,parameter :: CV_BDF   = 2 !! backward differentiation formulas (stiff problems)

    integer,parameter :: CV_NORMAL   = 1 !! itask: integrate to tout and interpolate the solution there
    integer,parameter :: CV_ONE_STEP = 2 !! `itask` value: take one internal step and return

    !-------------------------------------------------------------------------------------
    ! cvode.h : return values
    !-------------------------------------------------------------------------------------
    integer,parameter :: CV_SUCCESS      = 0 !! successful return
    integer,parameter :: CV_TSTOP_RETURN = 1 !! the stop time tstop was reached
    integer,parameter :: CV_ROOT_RETURN  = 2 !! a root of the root function was found

    integer,parameter :: CV_WARNING = 99 !! warning (not an error)

    integer,parameter :: CV_TOO_MUCH_WORK = -1 !! mxstep steps taken before reaching tout
    integer,parameter :: CV_TOO_MUCH_ACC  = -2 !! the requested accuracy is too high for the machine precision
    integer,parameter :: CV_ERR_FAILURE   = -3 !! the error test failed repeatedly or with |h| = hmin
    integer,parameter :: CV_CONV_FAILURE  = -4 !! the nonlinear solver failed to converge repeatedly or with |h| = hmin

    integer,parameter :: CV_LINIT_FAIL        = -5 !! the linear solver's initialization function failed
    integer,parameter :: CV_LSETUP_FAIL       = -6 !! the linear solver's setup function failed unrecoverably
    integer,parameter :: CV_LSOLVE_FAIL       = -7 !! the linear solver's solve function failed unrecoverably
    integer,parameter :: CV_RHSFUNC_FAIL      = -8 !! the right-hand side function failed unrecoverably
    integer,parameter :: CV_FIRST_RHSFUNC_ERR = -9 !! the right-hand side function failed at the first call
    integer,parameter :: CV_REPTD_RHSFUNC_ERR = -10 !! repeated recoverable right-hand side function errors
    integer,parameter :: CV_UNREC_RHSFUNC_ERR = -11 !! recoverable right-hand side failure, but no recovery possible
    integer,parameter :: CV_RTFUNC_FAIL       = -12 !! the root function failed unrecoverably
    integer,parameter :: CV_NLS_INIT_FAIL     = -13 !! the nonlinear solver's initialization function failed
    integer,parameter :: CV_NLS_SETUP_FAIL    = -14 !! the nonlinear solver's setup function failed unrecoverably
    integer,parameter :: CV_CONSTR_FAIL       = -15 !! the inequality constraints could not be satisfied
    integer,parameter :: CV_NLS_FAIL          = -16 !! the nonlinear solver failed unrecoverably

    integer,parameter :: CV_MEM_FAIL     = -20 !! a memory allocation failed
    integer,parameter :: CV_MEM_NULL     = -21 !! the solver object was not created
    integer,parameter :: CV_ILL_INPUT    = -22 !! an input argument is illegal
    integer,parameter :: CV_NO_MALLOC    = -23 !! the solver was not initialized (init not called)
    integer,parameter :: CV_BAD_K        = -24 !! illegal derivative order k in get_dky
    integer,parameter :: CV_BAD_T        = -25 !! illegal time t in get_dky (outside the last step)
    integer,parameter :: CV_BAD_DKY      = -26 !! illegal output vector in get_dky
    integer,parameter :: CV_TOO_CLOSE    = -27 !! tout too close to t0 to start the integration
    integer,parameter :: CV_VECTOROP_ERR = -28 !! a vector operation failed

    integer,parameter :: CV_PROJ_MEM_NULL      = -29 !! projection memory is not allocated
    integer,parameter :: CV_PROJFUNC_FAIL      = -30 !! the projection function failed unrecoverably
    integer,parameter :: CV_REPTD_PROJFUNC_ERR = -31 !! repeated recoverable projection function errors

    integer,parameter :: CV_CONTEXT_ERR = -32 !! SUNContext error (not used by FVODE)

    integer,parameter :: CV_UNRECOGNIZED_ERR = -99 !! an unrecognized error occurred

    !-------------------------------------------------------------------------------------
    ! cvode_ls.h : linear solver interface return values
    !-------------------------------------------------------------------------------------
    integer,parameter :: CVLS_SUCCESS         = 0 !! successful return
    integer,parameter :: CVLS_MEM_NULL        = -1 !! the solver object was not created
    integer,parameter :: CVLS_LMEM_NULL       = -2 !! no linear solver is attached
    integer,parameter :: CVLS_ILL_INPUT       = -3 !! an input argument is illegal or incompatible
    integer,parameter :: CVLS_MEM_FAIL        = -4 !! a memory allocation failed
    integer,parameter :: CVLS_PMEM_NULL       = -5 !! the preconditioner module was not initialized
    integer,parameter :: CVLS_JACFUNC_UNRECVR = -6 !! the Jacobian function failed unrecoverably
    integer,parameter :: CVLS_JACFUNC_RECVR   = -7 !! the Jacobian function failed recoverably
    integer,parameter :: CVLS_SUNMAT_FAIL     = -8 !! a matrix operation failed
    integer,parameter :: CVLS_SUNLS_FAIL      = -9 !! a linear solver operation failed

    !-------------------------------------------------------------------------------------
    ! cvode_diag.h : diagonal linear solver return values
    !-------------------------------------------------------------------------------------
    integer,parameter :: CVDIAG_SUCCESS   = 0 !! successful return
    integer,parameter :: CVDIAG_MEM_NULL  = -1 !! the solver object was not created
    integer,parameter :: CVDIAG_LMEM_NULL = -2 !! the diagonal linear solver is not attached
    integer,parameter :: CVDIAG_ILL_INPUT = -3 !! an input argument is illegal
    integer,parameter :: CVDIAG_MEM_FAIL  = -4 !! a memory allocation failed

    integer,parameter :: CVDIAG_INV_FAIL        = -5 !! the diagonal matrix M has a zero element (not invertible)
    integer,parameter :: CVDIAG_RHSFUNC_UNRECVR = -6 !! the right-hand side function failed unrecoverably
    integer,parameter :: CVDIAG_RHSFUNC_RECVR   = -7 !! the right-hand side function failed recoverably

    !-------------------------------------------------------------------------------------
    ! sundials_errors.h : generic SUNErrCode values (enum starting at -10000)
    !-------------------------------------------------------------------------------------
    integer,parameter :: SUN_SUCCESS              = 0 !! successful return
    integer,parameter :: SUN_ERR_MINIMUM          = -10000 !! lower bound of the SUNErrCode range (not an error code)
    integer,parameter :: SUN_ERR_ARG_CORRUPT      = -9999 !! argument provided is NULL or corrupted
    integer,parameter :: SUN_ERR_ARG_INCOMPATIBLE = -9998 !! argument provided is not compatible
    integer,parameter :: SUN_ERR_ARG_OUTOFRANGE   = -9997 !! argument is out of the valid range
    integer,parameter :: SUN_ERR_ARG_WRONGTYPE    = -9996 !! argument provided is not the right type
    integer,parameter :: SUN_ERR_ARG_DIMSMISMATCH = -9995 !! argument dimensions do not agree
    integer,parameter :: SUN_ERR_GENERIC          = -9994 !! an error occurred
    integer,parameter :: SUN_ERR_CORRUPT          = -9993 !! value is NULL or corrupt
    integer,parameter :: SUN_ERR_OUTOFRANGE       = -9992 !! value is out of the expected range
    integer,parameter :: SUN_ERR_FILE_OPEN        = -9991 !! unable to open file
    integer,parameter :: SUN_ERR_OP_FAIL          = -9990 !! an operation failed
    integer,parameter :: SUN_ERR_MEM_FAIL         = -9989 !! a memory operation failed
    integer,parameter :: SUN_ERR_MALLOC_FAIL      = -9988 !! memory allocation failed
    integer,parameter :: SUN_ERR_EXT_FAIL         = -9987 !! a failure occurred in an external library
    integer,parameter :: SUN_ERR_DESTROY_FAIL     = -9986 !! a destroy function returned an error
    integer,parameter :: SUN_ERR_NOT_IMPLEMENTED  = -9985 !! operation is not implemented
    integer,parameter :: SUN_ERR_USER_FCN_FAIL    = -9984 !! the user provided callback function failed
    integer,parameter :: SUN_ERR_UNKNOWN          = -9970 !! unknown error occurred
    integer,parameter :: SUN_ERR_MAXIMUM          = -1000 !! upper bound of the SUNErrCode range (not an error code)

    !-------------------------------------------------------------------------------------
    ! sundials_linearsolver.h
    !-------------------------------------------------------------------------------------
    ! SUNLinearSolver_Type
    integer,parameter :: SUNLINEARSOLVER_DIRECT           = 0 !! direct solver (needs a matrix)
    integer,parameter :: SUNLINEARSOLVER_ITERATIVE        = 1 !! iterative solver (matrix-free)
    integer,parameter :: SUNLINEARSOLVER_MATRIX_ITERATIVE = 2 !! iterative solver that uses a matrix
    integer,parameter :: SUNLINEARSOLVER_MATRIX_EMBEDDED  = 3 !! solver with an embedded (internal) matrix

    ! unrecoverable failure codes
    integer,parameter :: SUNLS_ATIMES_NULL       = -804 !! the atimes (A*v) function is NULL
    integer,parameter :: SUNLS_ATIMES_FAIL_UNREC = -805 !! the atimes function failed unrecoverably
    integer,parameter :: SUNLS_PSET_FAIL_UNREC   = -806 !! the preconditioner setup failed unrecoverably
    integer,parameter :: SUNLS_PSOLVE_NULL       = -807 !! the preconditioner solve function is NULL
    integer,parameter :: SUNLS_PSOLVE_FAIL_UNREC = -808 !! the preconditioner solve failed unrecoverably
    integer,parameter :: SUNLS_GS_FAIL           = -810 !! Gram-Schmidt orthogonalization failed
    integer,parameter :: SUNLS_QRSOL_FAIL        = -811 !! QR solve found a singular R

    ! recoverable failure codes
    integer,parameter :: SUNLS_RECOV_FAILURE    = 800 !! generic recoverable failure
    integer,parameter :: SUNLS_RES_REDUCED      = 801 !! nonconvergent solve, but the residual was reduced
    integer,parameter :: SUNLS_CONV_FAIL        = 802 !! nonconvergent solve
    integer,parameter :: SUNLS_ATIMES_FAIL_REC  = 803 !! the atimes function failed recoverably
    integer,parameter :: SUNLS_PSET_FAIL_REC    = 804 !! the preconditioner setup failed recoverably
    integer,parameter :: SUNLS_PSOLVE_FAIL_REC  = 805 !! the preconditioner solve failed recoverably
    integer,parameter :: SUNLS_PACKAGE_FAIL_REC = 806 !! recoverable failure in an external package
    integer,parameter :: SUNLS_QRFACT_FAIL      = 807 !! QR factorization found a singular matrix
    integer,parameter :: SUNLS_LUFACT_FAIL      = 808 !! LU factorization found a singular matrix (zero pivot)

    !-------------------------------------------------------------------------------------
    ! sundials_nonlinearsolver.h
    !-------------------------------------------------------------------------------------
    ! SUNNonlinearSolver_Type
    integer,parameter :: SUNNONLINEARSOLVER_ROOTFIND   = 0 !! solves F(y) = 0 (e.g. Newton)
    integer,parameter :: SUNNONLINEARSOLVER_FIXEDPOINT = 1 !! solves G(y) = y (fixed-point iteration)
    integer,parameter :: SUNNONLINEARSOLVER_HYBRID     = 2 !! supports both formulations

    integer,parameter :: SUN_NLS_CONTINUE   = 901 !! not converged, keep iterating
    integer,parameter :: SUN_NLS_CONV_RECVR = 902 !! convergence failure, try to recover
    integer,parameter :: SUN_NLS_SWITCH     = 903 !! switch nonlinear solver algorithms (auto solver)

    !-------------------------------------------------------------------------------------
    ! sundials_iterative.h
    !-------------------------------------------------------------------------------------
    ! SUNPrecType
    integer,parameter :: SUN_PREC_NONE  = 0 !! no preconditioning
    integer,parameter :: SUN_PREC_LEFT  = 1 !! preconditioning on the left only
    integer,parameter :: SUN_PREC_RIGHT = 2 !! preconditioning on the right only
    integer,parameter :: SUN_PREC_BOTH  = 3 !! preconditioning on both sides

    ! SUNGramSchmidtType
    integer,parameter :: SUN_MODIFIED_GS  = 1 !! modified Gram-Schmidt orthogonalization
    integer,parameter :: SUN_CLASSICAL_GS = 2 !! classical Gram-Schmidt orthogonalization

    !-------------------------------------------------------------------------------------
    ! sunlinsol_sp*.h, sunlinsol_pcg.h : Krylov defaults
    !-------------------------------------------------------------------------------------
    integer,parameter :: SUNSPGMR_MAXL_DEFAULT    = 5 !! SPGMR: default maximum Krylov dimension
    integer,parameter :: SUNSPGMR_MAXRS_DEFAULT   = 0 !! SPGMR: default maximum number of restarts
    integer,parameter :: SUNSPGMR_GSTYPE_DEFAULT  = SUN_MODIFIED_GS !! SPGMR: default Gram-Schmidt type
    integer,parameter :: SUNSPFGMR_MAXL_DEFAULT   = 5 !! SPFGMR: default maximum Krylov dimension
    integer,parameter :: SUNSPFGMR_MAXRS_DEFAULT  = 0 !! SPFGMR: default maximum number of restarts
    integer,parameter :: SUNSPFGMR_GSTYPE_DEFAULT = SUN_MODIFIED_GS !! SPFGMR: default Gram-Schmidt type
    integer,parameter :: SUNSPBCGS_MAXL_DEFAULT   = 5 !! SPBCGS: default maximum number of iterations
    integer,parameter :: SUNSPTFQMR_MAXL_DEFAULT  = 5 !! SPTFQMR: default maximum number of iterations
    integer,parameter :: SUNPCG_MAXL_DEFAULT      = 5 !! default maximum number of PCG iterations

    !-------------------------------------------------------------------------------------
    ! sundials_types.h : SUNOutputFormat
    !-------------------------------------------------------------------------------------
    integer,parameter :: SUN_OUTPUTFORMAT_TABLE = 0 !! print statistics as a table
    integer,parameter :: SUN_OUTPUTFORMAT_CSV   = 1 !! print statistics as comma-separated values

    !-------------------------------------------------------------------------------------
    ! sundials_types.h : floating point limits
    !-------------------------------------------------------------------------------------
    real(wp),parameter :: SUN_UNIT_ROUNDOFF = epsilon(1.0_wp) !! `DBL_EPSILON`
    real(wp),parameter :: SUN_BIG_REAL      = huge(1.0_wp)    !! `DBL_MAX`
    real(wp),parameter :: SUN_SMALL_REAL    = tiny(1.0_wp)    !! `DBL_MIN`

    end module fvode_parameters
!*****************************************************************************************

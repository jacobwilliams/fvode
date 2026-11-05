!*****************************************************************************************
!>
!  CVODE linear solver interface (SUNDIALS v7.9.0 `cvode_ls.c`).
!
!  Connects a [[sun_linear_solver]] (and, for matrix-based solvers, a
!  [[sun_matrix]]) to the Newton iteration: decides when to re-evaluate the
!  Jacobian, builds `A = I - gamma*J`, and solves the linear systems. The
!  Jacobian is either user-supplied ([[cv_jac_fn]]) or approximated by
!  difference quotients (dense and band matrices).
!
!  Only direct (matrix-based) solvers are supported so far; the paths for
!  iterative solvers and preconditioning are marked `[not yet ported]`.

    submodule (fvode_core) fvode_core_ls

    use fvode_parameters, only: SUNLINEARSOLVER_DIRECT, SUNLINEARSOLVER_ITERATIVE, &
                                SUNLINEARSOLVER_MATRIX_ITERATIVE, SUNLINEARSOLVER_MATRIX_EMBEDDED

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetLinearSolver`: attach a linear solver `LS` and, for matrix-based
!  solvers, a matrix `A`. Both are copied into the CVODE object. For direct
!  solvers the internal difference-quotient Jacobian is used unless
!  [[CVodeSetJacFn]] is called; for iterative (Krylov) solvers the internal
!  difference-quotient Jacobian-times-vector product is used unless
!  [[CVodeSetJacTimes]] is called, and a preconditioner can be attached with
!  [[CVodeSetPreconditioner]] or [[CVBandPrecInit]].

    module subroutine CVodeSetLinearSolver(me, LS, A, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    class(sun_linear_solver),intent(in) :: LS  !! the linear solver (copied)
    class(sun_matrix),intent(in),optional :: A  !! the matrix (copied; direct solvers only)
    integer,intent(out) :: flag  !! `CVLS_SUCCESS` or a negative `CVLS_*` error flag

    integer :: LSType  !! the linear solver type
    logical :: iterative    !! is the solver iterative?
    logical :: matrixbased  !! is a matrix structure used?

    if (.not. me%MallocDone) then
        call cvProcessError(me, CVLS_MEM_NULL, 'CVodeSetLinearSolver', 'Integrator memory is NULL.')
        flag = CVLS_MEM_NULL
        return
    end if

    ! retrieve the LS type
    LSType = LS%get_type()

    ! set flags based on LS type
    iterative   = (LSType /= SUNLINEARSOLVER_DIRECT)
    matrixbased = ((LSType /= SUNLINEARSOLVER_ITERATIVE) .and. &
                   (LSType /= SUNLINEARSOLVER_MATRIX_EMBEDDED))

    ! check for compatible LS type, matrix and "atimes" support
    if ((LSType == SUNLINEARSOLVER_MATRIX_EMBEDDED) .and. present(A)) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetLinearSolver', &
                            'Incompatible inputs: matrix-embedded LS requires NULL matrix')
        flag = CVLS_ILL_INPUT
        return
    end if

    if (iterative) then
        if (.not. matrixbased .and. (LSType /= SUNLINEARSOLVER_MATRIX_EMBEDDED) .and. &
            .not. LS%supports_atimes()) then
            call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetLinearSolver', &
                                'Incompatible inputs: iterative LS must support ATimes routine')
            flag = CVLS_ILL_INPUT
            return
        end if
        if (matrixbased .and. .not. present(A)) then
            call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetLinearSolver', &
                                'Incompatible inputs: matrix-iterative LS requires non-NULL matrix')
            flag = CVLS_ILL_INPUT
            return
        end if
    else if (.not. present(A)) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetLinearSolver', &
                            'Incompatible inputs: direct LS requires non-NULL matrix')
        flag = CVLS_ILL_INPUT
        return
    end if

    ! free any existing system solver attached to CVode
    call cvLFree(me)

    ! set four main system linear solver function fields in cv_mem
    me%ls_kind       = LS_GENERIC
    me%lsetup_exists = .true.
    me%lsolve_exists = .true.

    ! set SUNLinearSolver pointer
    allocate(me%LS, source=LS)

    ! set "iterative" and "matrixbased" flags
    me%ls_iterative   = iterative
    me%ls_matrixbased = matrixbased

    ! initialize Jacobian-related data
    if (present(A)) then
        me%jacDQ = .true.
    else
        me%jacDQ = .false.
    end if
    me%jac_user => null()

    ! initialize the linear system function (internal)
    me%user_linsys = .false.
    me%linsys_user => null()

    ! initialize Jacobian-times-vector data (internal DQ)
    me%jtimesDQ = .true.
    me%jtsetup_user => null()
    me%jtimes_user => null()
    me%jt_f => null()

    ! initialize the preconditioner (none)
    me%pset_user => null()
    me%psolve_user => null()
    me%pre_kind = PRE_NONE

    ! initialize counters
    call cvLsInitializeCounters(me)

    ! set default values for the rest of the LS parameters
    me%msbj         = CVLS_MSBJ
    me%jbad         = .true.
    me%dgmax_jbad   = CVLS_DGMAX
    me%eplifac      = CVLS_EPLIN
    me%ls_last_flag = CVLS_SUCCESS

    ! if LS supports ATimes, attach CVLs routine
    if (me%LS%supports_atimes()) me%LS%has_atimes = .true.

    ! if LS supports preconditioning, initialize pset/psol to NULL
    if (me%LS%supports_preconditioner()) then
        me%LS%has_psetup = .false.
        me%LS%has_psolve = .false.
    end if

    ! when using a SUNMatrix object, store pointer to A (the saved Jacobian
    ! is allocated in cvLsInitialize)
    if (present(A)) allocate(me%A, source=A)

    ! allocate memory for ytemp and x
    allocate(me%ls_ytemp(me%n), me%ls_x(me%n))

    ! compute sqrtN from a dot product
    if (iterative) me%nrmfac = sqrt(real(me%n, wp))

    ! for matrix-based LS, enable solution scaling with BDF
    me%scalesol = (matrixbased .and. me%lmm == CV_BDF)

    flag = CVLS_SUCCESS

    end subroutine CVodeSetLinearSolver
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetPreconditioner`: set the preconditioner setup and solve
!  functions for an iterative linear solver (absent: none).

    module subroutine CVodeSetPreconditioner(me, psetup, psolve, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_psetup_fn),optional :: psetup  !! preconditioner setup function
    procedure(cv_psolve_fn),optional :: psolve  !! preconditioner solve function
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetPreconditioner', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! store function pointers for user-supplied routines in CVLs interface
    me%pset_user => null()
    me%psolve_user => null()
    if (present(psetup)) me%pset_user => psetup
    if (present(psolve)) me%psolve_user => psolve
    if (present(psetup) .or. present(psolve)) then
        me%pre_kind = PRE_USER
    else
        me%pre_kind = PRE_NONE
    end if

    ! issue error if LS object does not allow user-supplied preconditioning
    if (.not. me%LS%supports_preconditioner()) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetPreconditioner', &
            'SUNLinearSolver object does not support user-supplied preconditioning')
        flag = CVLS_ILL_INPUT
        return
    end if

    ! notify iterative linear solver to call CVLs interface routines
    call set_ls_preconditioner(me, present(psetup), present(psolve))

    flag = CVLS_SUCCESS

    end subroutine CVodeSetPreconditioner
!*****************************************************************************************

!*****************************************************************************************
!>
!  Attach the CVLs preconditioner routines to the linear solver (C:
!  `SUNLinSolSetPreconditioner(LS, cv_mem, cvLsPSetup, cvLsPSolve)`).

    subroutine set_ls_preconditioner(me, has_psetup, has_psolve)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: has_psetup  !! is there a preconditioner setup function?
    logical,intent(in) :: has_psolve  !! is there a preconditioner solve function?

    me%LS%has_psetup = has_psetup
    me%LS%has_psolve = has_psolve

    end subroutine set_ls_preconditioner
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetJacTimes`: set the Jacobian-times-vector setup and product
!  functions. Without `jtimes`, the internal difference-quotient
!  approximation is used.

    module subroutine CVodeSetJacTimes(me, jtsetup, jtimes, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_jtsetup_fn),optional :: jtsetup  !! J*v setup function
    procedure(cv_jtimes_fn),optional :: jtimes  !! J*v function
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetJacTimes', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! issue error if LS object does not allow user-supplied ATimes
    if (.not. me%LS%supports_atimes()) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetJacTimes', &
            'SUNLinearSolver object does not support user-supplied ATimes routine')
        flag = CVLS_ILL_INPUT
        return
    end if

    ! store function pointers for user-supplied routines in CVLs
    ! interface (NULL jtimes implies use of DQ default)
    if (present(jtimes)) then
        me%jtimesDQ = .false.
        me%jtsetup_user => null()
        if (present(jtsetup)) me%jtsetup_user => jtsetup
        me%jtimes_user => jtimes
    else
        me%jtimesDQ = .true.
        me%jtsetup_user => null()
        me%jtimes_user => null()
        me%jt_f => null()
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetJacTimes
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetJacTimesRhsFn`: the right-hand side function used in the
!  difference-quotient Jacobian-times-vector approximation (absent: the ODE
!  right-hand side).

    module subroutine CVodeSetJacTimesRhsFn(me, jtimesRhsFn, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_rhs_fn),optional :: jtimesRhsFn  !! f used by the DQ J*v
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetJacTimesRhsFn', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! check if using internal finite difference approximation
    if (.not. me%jtimesDQ) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetJacTimesRhsFn', &
            'Internal finite-difference Jacobian-vector product is disabled.')
        flag = CVLS_ILL_INPUT
        return
    end if

    ! store function pointers for RHS function (NULL implies use ODE RHS)
    me%jt_f => null()
    if (present(jtimesRhsFn)) me%jt_f => jtimesRhsFn

    flag = CVLS_SUCCESS

    end subroutine CVodeSetJacTimesRhsFn
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEpsLin`: factor in the linear solver tolerance
!  `delta = eplifac * tq(4)` (0 restores the default, 0.05).

    module subroutine CVodeSetEpsLin(me, eplifac, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eplifac  !! linear tolerance factor (>= 0)
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetEpsLin', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! check for legal eplifac
    if (eplifac < ZERO) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetEpsLin', 'eplifac < 0 illegal.')
        flag = CVLS_ILL_INPUT
        return
    end if

    if (eplifac == ZERO) then
        me%eplifac = CVLS_EPLIN
    else
        me%eplifac = eplifac
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetEpsLin
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetLSNormFactor`: factor converting the WRMS norm to the 2-norm
!  used by the iterative solver (> 0: given; < 0: `sqrt(dot(1,1))`;
!  0: `sqrt(N)`, the default).

    module subroutine CVodeSetLSNormFactor(me, nrmfac, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: nrmfac  !! the norm factor
    integer,intent(out) :: flag  !! `CVLS_SUCCESS` or `CVLS_LMEM_NULL`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetLSNormFactor', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    if (nrmfac > ZERO) then
        ! user-provided factor
        me%nrmfac = nrmfac
    else if (nrmfac < ZERO) then
        ! compute factor for WRMS norm with dot product
        me%ls_ytemp = ONE
        me%nrmfac = sqrt(dot_product(me%ls_ytemp, me%ls_ytemp))
    else
        ! compute default factor for WRMS norm from vector length
        me%nrmfac = sqrt(real(me%n, wp))
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetLSNormFactor
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetLinSysFn`: set a user function that computes the linear system
!  matrix `A = I - gamma*J` directly (absent: the internal one, from the
!  Jacobian). Requires a matrix-based linear solver.

    module subroutine CVodeSetLinSysFn(me, linsys, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_linsys_fn),optional :: linsys  !! user linear system function
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetLinSysFn', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! return with failure if linsys cannot be used
    if (present(linsys) .and. .not. allocated(me%A)) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetLinSysFn', &
            'Linear system setup routine cannot be supplied for NULL SUNMatrix')
        flag = CVLS_ILL_INPUT
        return
    end if

    ! set the linear system routine pointer, and update relevant flags
    if (present(linsys)) then
        me%user_linsys = .true.
        me%linsys_user => linsys
    else
        me%user_linsys = .false.
        me%linsys_user => null()
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetLinSysFn
!*****************************************************************************************
!>
!  `CVodeSetJacFn`: set a user Jacobian function, or (if `jac` is not
!  present) go back to the internal difference-quotient Jacobian.

    module subroutine CVodeSetJacFn(me, jac, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_jac_fn),optional :: jac  !! user Jacobian (absent: difference quotients)
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetJacFn', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! return with failure if jac cannot be used
    if (present(jac) .and. .not. allocated(me%A)) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetJacFn', &
                            'Jacobian routine cannot be supplied for NULL SUNMatrix')
        flag = CVLS_ILL_INPUT
        return
    end if

    ! set Jacobian routine pointer, and update relevant flags
    if (present(jac)) then
        me%jacDQ = .false.
        me%jac_user => jac
    else
        me%jacDQ = .true.
        me%jac_user => null()
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetJacFn
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetJacEvalFrequency`: maximum number of steps between Jacobian
!  evaluations (0 restores the default, 51).

    module subroutine CVodeSetJacEvalFrequency(me, msbj, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(li),intent(in) :: msbj  !! max steps between Jacobian evaluations (0: default)
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetJacEvalFrequency', &
                            'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    if (msbj < 0) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetJacEvalFrequency', &
                            'A negative evaluation frequency was provided.')
        flag = CVLS_ILL_INPUT
        return
    end if

    if (msbj == 0) then
        me%msbj = CVLS_MSBJ
    else
        me%msbj = msbj
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetJacEvalFrequency
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetLinearSolutionScaling`: turn the BDF solution scaling for a change
!  in gamma on/off (matrix-based solvers with BDF only).

    module subroutine CVodeSetLinearSolutionScaling(me, onoff, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: onoff  !! turn the solution scaling on or off
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL` or `CVLS_ILL_INPUT`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetLinearSolutionScaling', &
                            'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! check for valid solver and method type
    if (.not. me%ls_matrixbased .or. me%lmm /= CV_BDF) then
        flag = CVLS_ILL_INPUT
        return
    end if

    me%scalesol = onoff

    flag = CVLS_SUCCESS

    end subroutine CVodeSetLinearSolutionScaling
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetDeltaGammaMaxBadJac`: gamma ratio threshold below which a Jacobian
!  is considered bad after a convergence failure (<= 0 restores the default).

    module subroutine CVodeSetDeltaGammaMaxBadJac(me, dgmax_jbad, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: dgmax_jbad  !! gamma ratio threshold (<= 0: default, 0.2)
    integer,intent(out) :: flag  !! `CVLS_SUCCESS` or `CVLS_LMEM_NULL`

    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVodeSetDeltaGammaMaxBadJac', &
                            'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    if (dgmax_jbad <= ZERO) then
        me%dgmax_jbad = CVLS_DGMAX
    else
        me%dgmax_jbad = dgmax_jbad
    end if

    flag = CVLS_SUCCESS

    end subroutine CVodeSetDeltaGammaMaxBadJac
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsInitialize`: performs remaining initializations specific to the
!  linear solver interface (and solver itself).

    module function cvLsInitialize(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! the flag of the linear solver's initialize, or `CVLS_ILL_INPUT`

    ! test for valid combinations of matrix & Jacobian routines
    if (allocated(me%A) .and. me%user_linsys) then
        ! matrix-based case with a user-supplied linear system function
        continue
    else if (allocated(me%A)) then
        ! matrix-based case with the internal linear system function:
        ! the internal DQ Jacobian needs a dense or band matrix
        if (me%jacDQ) then
            if (me%A%get_id() /= SUNMATRIX_DENSE .and. me%A%get_id() /= SUNMATRIX_BAND) then
                call cvProcessError(me, CVLS_ILL_INPUT, 'cvLsInitialize', &
                                    'No Jacobian constructor available for SUNMatrix type')
                me%ls_last_flag = CVLS_ILL_INPUT
                flag = CVLS_ILL_INPUT
                return
            end if
        end if

        ! allocate internally saved Jacobian if not already done
        if (.not. allocated(me%savedJ)) allocate(me%savedJ, source=me%A)
    else
        ! matrix-free: no Jacobian or linear system function
        me%jacDQ = .false.
        me%jac_user => null()
        me%user_linsys = .false.
        me%linsys_user => null()
    end if

    ! reset counters
    call cvLsInitializeCounters(me)

    ! set Jacobian-vector product fields, based on jtimesDQ
    if (me%jtimesDQ) me%jtsetup_user => null()

    ! if A is NULL and psetup is not present, then cvLsSetup does
    ! not need to be called, so set the lsetup function to NULL
    if (.not. allocated(me%A) .and. .not. me%LS%has_psetup) me%lsetup_exists = .false.

    ! call LS initialize routine, and return result
    me%ls_last_flag = me%LS%initialize()
    flag = me%ls_last_flag

    end function cvLsInitialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsReInitialize`: reset the counters.

    module function cvLsReInitialize(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! always `CVLS_SUCCESS`

    call cvLsInitializeCounters(me)
    flag = CVLS_SUCCESS

    end function cvLsReInitialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsSetup`: decide whether the Jacobian is out of date (`jbad`), update
!  the matrix `A = I - gamma*J` (re-evaluating J if needed) or the
!  preconditioner, and set up the linear solver.
!
!  The Jacobian is re-evaluated on the first step, after `msbj` steps, after
!  a convergence failure with a bad Jacobian and little change in gamma, or
!  after any other convergence failure. Uses `me%y` and `me%ftemp` (C:
!  ypred, fpred) and sets `me%jcur` (C: `*jcurPtr`, which is `&cv_jcur`).

    module function cvLsSetup(me, convfail) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: convfail  !! convergence failure flag (`CV_NO_FAILURES`, `CV_FAIL_BAD_J`, `CV_FAIL_OTHER`)
    integer :: flag  !! 0 on success, > 0 for a recoverable failure, < 0 otherwise

    real(wp) :: dgamma  !! |gamma/gammap - 1|
    integer :: retval  !! return flag of the linear system function
    logical :: jc  !! was the Jacobian re-evaluated? (user linear system function)

    ! use nst, gamma/gammap, and convfail to set J/P eval. flag jok
    dgamma  = abs((me%gamma / me%gammap) - ONE)
    me%jbad = (me%nst == 0) .or. me%first_step_after_resize .or. &
              (me%nst >= me%nstlj + me%msbj) .or. &
              ((convfail == CV_FAIL_BAD_J) .and. (dgamma < me%dgmax_jbad)) .or. &
              (convfail == CV_FAIL_OTHER)

    ! setup the linear system if necessary
    if (allocated(me%A)) then

        ! update J if appropriate and evaluate A = I - gamma J
        if (me%user_linsys) then
            jc = me%jcur
            retval = me%linsys_user(me, me%tn, me%y, me%ftemp, me%A, .not. me%jbad, jc, &
                                    me%gamma, me%vtemp1, me%vtemp2, me%vtemp3)
            me%jcur = jc
        else
            retval = cvLsLinSys(me, me%tn, me%y, me%ftemp, .not. me%jbad, me%gamma)
        end if

        ! update counters
        if (me%jcur) then
            me%nje   = me%nje + 1
            me%nstlj = me%nst
            me%tnlj  = me%tn
        end if

        ! check linsys return flag
        if (retval /= CVLS_SUCCESS) then
            if (me%user_linsys) then
                if (retval < 0) then
                    call cvProcessError(me, CVLS_JACFUNC_UNRECVR, 'cvLsSetup', &
                        'The Jacobian routine failed in an unrecoverable manner.')
                    me%ls_last_flag = CVLS_JACFUNC_UNRECVR
                    flag = -1
                else
                    me%ls_last_flag = CVLS_JACFUNC_RECVR
                    flag = 1
                end if
            else
                flag = retval
            end if
            return
        end if

    else
        ! matrix-free: assume the preconditioner is recomputed if jbad
        me%jcur = me%jbad
    end if

    ! call LS setup routine -- the LS may call cvLsPSetup, who will
    ! pass the heuristic suggestions above to the user code(s)
    if (allocated(me%A)) then
        me%ls_last_flag = me%LS%setup(me%A, prob=me)
    else
        me%ls_last_flag = me%LS%setup(prob=me)
    end if

    ! if the SUNMatrix was NULL, update heuristics flags
    if (.not. allocated(me%A)) then

        ! if user set jcur to SUNTRUE, increment npe and save nst value
        if (me%jcur) then
            me%npe   = me%npe + 1
            me%nstlj = me%nst
            me%tnlj  = me%tn
        end if

        ! update jcur flag if we suggested an update
        if (me%jbad) me%jcur = .true.

    end if

    flag = me%ls_last_flag

    end function cvLsSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsSolve`: solve the linear system `A x = b`; b is overwritten by x.
!
!  For iterative solvers the tolerance is `eplifac * tq(4)` (in the WRMS
!  norm, converted to the solver's 2-norm with `nrmfac`), and the weights
!  are passed as scaling vectors (or folded into the tolerance). For BDF
!  with a matrix-based solver the solution is scaled by `2/(1 + gamrat)`
!  to account for a change in gamma since the last setup. Uses `me%ewt`,
!  `me%y` and `me%ftemp` (C: weight, ynow, fnow).
!
!  Returns 0 on success, 1 for a recoverable failure, -1 otherwise.

    module function cvLsSolve(me, b) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(inout) :: b  !! right-hand side on input, solution on output
    integer :: flag  !! 0 on success, 1 for a recoverable failure, -1 for an unrecoverable one

    integer :: retval  !! return flag of the linear solver
    integer :: curiter  !! current nonlinear iteration
    integer :: nli_inc  !! number of linear iterations in this solve
    real(wp) :: bnorm  !! WRMS norm of b
    real(wp) :: deltar  !! tolerance in the WRMS norm
    real(wp) :: delta  !! tolerance passed to the linear solver
    real(wp) :: w_mean  !! RMS of the weights

    ! retrieve the current nonlinear solver iteration count
    curiter = me%NLS%get_cur_iter()

    ! if the linear solver is iterative:
    ! test norm(b), if small, return x = 0 or x = b;
    ! set linear solver tolerance (in left/right scaled 2-norm)
    if (me%ls_iterative) then
        deltar = me%eplifac * me%tq(4)
        bnorm = n_vwrmsnorm(b, me%ewt)
        if (bnorm <= deltar) then
            if (curiter > 0) b = ZERO
            me%ls_last_flag = CVLS_SUCCESS
            flag = me%ls_last_flag
            return
        end if
        ! adjust tolerance for 2-norm
        delta = deltar * me%nrmfac
    else
        delta = ZERO
    end if

    ! set scaling vectors for LS to use (if applicable)
    if (me%LS%supports_scaling_vectors()) then
        retval = me%LS%set_scaling_vectors(me%ewt, me%ewt)
        if (retval /= SUN_SUCCESS) then
            call cvProcessError(me, CVLS_SUNLS_FAIL, 'cvLsSolve', &
                                'Error in calling SUNLinSolSetScalingVectors')
            me%ls_last_flag = CVLS_SUNLS_FAIL
            flag = me%ls_last_flag
            return
        end if

    ! if solver is iterative and does not support scaling vectors, update the
    ! tolerance in an attempt to account for weight vector (assuming
    ! homogeneous weights and a basic 2-norm in the solver)
    else if (me%ls_iterative) then
        me%ls_x = ONE
        w_mean = n_vwrmsnorm(me%ewt, me%ls_x)
        delta = delta / w_mean
    end if

    ! set initial guess x = 0 to LS
    me%ls_x = ZERO

    ! set zero initial guess flag
    retval = me%LS%set_zero_guess(.true.)
    if (retval /= SUN_SUCCESS) then
        flag = -1
        return
    end if

    ! if a user-provided jtsetup routine is supplied, call that here
    if (associated(me%jtsetup_user)) then
        me%ls_last_flag = me%jtsetup_user(me, me%tn, me%y, me%ftemp)
        me%njtsetup = me%njtsetup + 1
        if (me%ls_last_flag /= 0) then
            call cvProcessError(me, me%ls_last_flag, 'cvLsSolve', &
                'The Jacobian x vector setup routine failed in an unrecoverable manner.')
            flag = me%ls_last_flag
            return
        end if
    end if

    ! call solver, and copy x to b
    if (allocated(me%A)) then
        retval = me%LS%solve(me%A, me%ls_x, b, delta, prob=me)
    else
        retval = me%LS%solve(x=me%ls_x, b=b, tol=delta, prob=me)
    end if
    b = me%ls_x

    ! if using a direct or matrix-iterative solver, BDF method, and gamma has
    ! changed, scale the correction to account for change in gamma
    if (me%scalesol .and. me%gamrat /= ONE) b = (TWO / (ONE + me%gamrat)) * b

    ! retrieve statistics from iterative linear solvers
    nli_inc = 0
    if (me%ls_iterative) nli_inc = me%LS%num_iters()

    ! increment counters nli and ncfl
    me%nli = me%nli + nli_inc
    if (retval /= SUN_SUCCESS) me%ncfl = me%ncfl + 1

    ! interpret solver return value
    me%ls_last_flag = retval

    select case (retval)
    case (SUN_SUCCESS)
        flag = 0
    case (SUNLS_RES_REDUCED)
        ! allow reduction but not solution on first nonlinear iteration,
        ! otherwise return with a recoverable failure
        if (curiter == 0) then
            flag = 0
        else
            flag = 1
        end if
    case (SUNLS_CONV_FAIL, SUNLS_ATIMES_FAIL_REC, SUNLS_PSOLVE_FAIL_REC, &
          SUNLS_PACKAGE_FAIL_REC, SUNLS_QRFACT_FAIL, SUNLS_LUFACT_FAIL)
        flag = 1
    case (SUN_ERR_ARG_CORRUPT, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_MEM_FAIL, &
          SUNLS_GS_FAIL, SUNLS_QRSOL_FAIL)
        flag = -1
    case (SUN_ERR_EXT_FAIL)
        call cvProcessError(me, SUN_ERR_EXT_FAIL, 'cvLsSolve', &
                            'Failure in SUNLinSol external package')
        flag = -1
    case (SUNLS_ATIMES_FAIL_UNREC)
        call cvProcessError(me, SUNLS_ATIMES_FAIL_UNREC, 'cvLsSolve', &
            'The Jacobian x vector routine failed in an unrecoverable manner.')
        flag = -1
    case (SUNLS_PSOLVE_FAIL_UNREC)
        call cvProcessError(me, SUNLS_PSOLVE_FAIL_UNREC, 'cvLsSolve', &
            'The preconditioner solve routine failed in an unrecoverable manner.')
        flag = -1
    case default
        call cvProcessError(me, retval, 'cvLsSolve', &
                            'Unrecognized error return value from SUNLinSolSolve')
        flag = -1
    end select

    end function cvLsSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsATimes`: the matrix-vector product for iterative solvers,
!  `z = (I - gamma*J) v`, with J*v from the user's jtimes or the
!  difference-quotient approximation [[cvLsDQJtimes]].

    module function cvLsATimes(me, v, z) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
    real(wp),dimension(:),intent(out) :: z  !! z = (I - gamma J) v
    integer :: flag  !! 0 on success, the jtimes flag otherwise

    integer :: retval  !! return flag of jtimes

    ! call Jacobian-times-vector product routine
    ! (either user-supplied or internal DQ)
    if (me%jtimesDQ) then
        retval = cvLsDQJtimes(me, v, z, me%tn, me%y, me%ftemp)
    else
        retval = me%jtimes_user(me, v, z, me%tn, me%y, me%ftemp, me%ls_ytemp)
    end if
    me%njtimes = me%njtimes + 1
    if (retval /= 0) then
        flag = retval
        return
    end if

    ! add contribution from identity matrix
    z = v - me%gamma * z

    flag = 0

    end function cvLsATimes
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsPSetup`: call the preconditioner setup function (user or band
!  preconditioner) with `jok = .not. jbad`; sets `me%jcur`.

    module function cvLsPSetup(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! 0 on success, > 0 recoverable, < 0 unrecoverable

    logical :: jcur  !! was the Jacobian data recomputed?

    ! call user-supplied pset routine
    jcur = me%jcur
    select case (me%pre_kind)
    case (PRE_BANDPRE)
        flag = cvBandPrecSetup(me, me%tn, me%y, me%ftemp, .not. me%jbad, jcur, me%gamma)
    case default
        flag = me%pset_user(me, me%tn, me%y, me%ftemp, .not. me%jbad, jcur, me%gamma)
    end select
    me%jcur = jcur

    end function cvLsPSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsPSolve`: call the preconditioner solve function (user or band
!  preconditioner) and count it.

    module function cvLsPSolve(me, r, z, tol, lr) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: r  !! the right-hand side
    real(wp),dimension(:),intent(out) :: z  !! the solution of P z = r
    real(wp),intent(in) :: tol  !! tolerance
    integer,intent(in) :: lr  !! `SUN_PREC_LEFT` or `SUN_PREC_RIGHT`
    integer :: flag  !! 0 on success, > 0 recoverable, < 0 unrecoverable

    ! call the user-supplied psolve routine, and accumulate count
    select case (me%pre_kind)
    case (PRE_BANDPRE)
        flag = cvBandPrecSolve(me, r, z)
    case default
        flag = me%psolve_user(me, me%tn, me%y, me%ftemp, r, z, me%gamma, tol, lr)
    end select
    me%nps = me%nps + 1

    end function cvLsPSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsDQJtimes`: difference-quotient approximation of `Jv = J(t,y) v`:
!  `Jv = [f(t, y + sig*v) - fy] / sig` with `sig = 1/||v||_WRMS`; on a
!  recoverable f failure, sig is reduced by 0.25 (up to 3 tries).

    function cvLsDQJtimes(me, v, Jv, t, y, fy) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: v  !! the vector to multiply
    real(wp),dimension(:),intent(out) :: Jv  !! the product J*v
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    integer :: flag  !! 0 on success, 1 recoverable, -1 unrecoverable

    real(wp) :: sig  !! the increment
    real(wp) :: siginv  !! 1/sig
    integer :: iter  !! retry counter
    integer :: retval  !! return flag of f

    ! initialize perturbation to 1/||v||
    sig = ONE / n_vwrmsnorm(v, me%ewt)

    retval = 0
    do iter = 0, MAX_DQITERS - 1

        ! set work = y + sig*v
        me%ls_ytemp = sig * v + y

        ! set Jv = f(tn, y+sig*v)
        if (associated(me%jt_f)) then
            retval = me%jt_f(me, t, me%ls_ytemp, Jv)
        else
            retval = me%f(me, t, me%ls_ytemp, Jv)
        end if
        me%nfeDQ = me%nfeDQ + 1
        if (retval == 0) exit
        if (retval < 0) then
            flag = -1
            return
        end if

        ! if f failed recoverably, shrink sig and retry
        sig = sig * PT25
    end do

    ! if retval still isn't 0, return with a recoverable failure
    if (retval > 0) then
        flag = 1
        return
    end if

    ! replace Jv by (Jv - fy)/sig
    siginv = ONE / sig
    Jv = siginv * (Jv - fy)

    flag = 0

    end function cvLsDQJtimes
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsFree`: free the linear solver interface memory.

    module subroutine cvLsFree(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    if (allocated(me%ls_ytemp)) deallocate(me%ls_ytemp)
    if (allocated(me%ls_x))     deallocate(me%ls_x)
    if (allocated(me%savedJ))   deallocate(me%savedJ)
    if (allocated(me%A))        deallocate(me%A)
    if (allocated(me%LS))       deallocate(me%LS)
    me%jac_user => null()
    me%user_linsys = .false.
    me%linsys_user => null()
    me%pset_user => null()
    me%psolve_user => null()
    me%jtsetup_user => null()
    me%jtimes_user => null()
    me%jt_f => null()
    ! free preconditioner memory (if applicable)
    call cvBandPrecFree(me)
    me%pre_kind = PRE_NONE

    end subroutine cvLsFree
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsInitializeCounters`

    subroutine cvLsInitializeCounters(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    me%nje      = 0
    me%nfeDQ    = 0
    me%nstlj    = 0
    me%npe      = 0
    me%nli      = 0
    me%nps      = 0
    me%ncfl     = 0
    me%njtsetup = 0
    me%njtimes  = 0

    end subroutine cvLsInitializeCounters
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsLinSys`: set up the linear system `A = I - gamma*J`. If `jok`, the
!  saved Jacobian is reused; otherwise J is re-evaluated (and saved).
!  Sets `me%jcur`.

    function cvLsLinSys(me, t, y, fy, jok, gamma) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    logical,intent(in) :: jok  !! the saved Jacobian may be reused
    real(wp),intent(in) :: gamma  !! gamma = h*rl1 in A = I - gamma*J
    integer :: flag  !! `CVLS_SUCCESS`, 1 recoverable, -1 or `CVLS_SUNMAT_FAIL` unrecoverable

    integer :: retval  !! return flag of the Jacobian or matrix routines

    ! if jok = .true., use saved copy of J
    if (jok) then

        ! indicate J is not current
        me%jcur = .false.

        ! copy saved J into A
        call me%savedJ%copy_to(me%A, retval)
        if (retval /= 0) then
            call cvProcessError(me, CVLS_SUNMAT_FAIL, 'cvLsLinSys', &
                                'A SUNMatrix routine failed in an unrecoverable manner.')
            me%ls_last_flag = CVLS_SUNMAT_FAIL
            flag = me%ls_last_flag
            return
        end if

    else

        ! call jac() to update the function
        me%jcur = .true.

        ! clear the linear system matrix if necessary
        if (me%LS%get_type() == SUNLINEARSOLVER_DIRECT) call me%A%zero()

        ! compute new Jacobian matrix
        if (me%jacDQ) then
            retval = cvLsDQJac(me, t, y, fy, me%A)
        else
            retval = me%jac_user(me, t, y, fy, me%A, me%vtemp1, me%vtemp2, me%vtemp3)
        end if
        if (retval < 0) then
            call cvProcessError(me, CVLS_JACFUNC_UNRECVR, 'cvLsLinSys', &
                                'The Jacobian routine failed in an unrecoverable manner.')
            me%ls_last_flag = CVLS_JACFUNC_UNRECVR
            flag = -1
            return
        end if
        if (retval > 0) then
            me%ls_last_flag = CVLS_JACFUNC_RECVR
            flag = 1
            return
        end if

        ! update saved copy of the Jacobian
        call me%A%copy_to(me%savedJ, retval)
        if (retval /= 0) then
            call cvProcessError(me, CVLS_SUNMAT_FAIL, 'cvLsLinSys', &
                                'A SUNMatrix routine failed in an unrecoverable manner.')
            me%ls_last_flag = CVLS_SUNMAT_FAIL
            flag = me%ls_last_flag
            return
        end if

    end if

    ! perform linear combination A = I - gamma*J
    call me%A%scale_add_i(-gamma)

    flag = CVLS_SUCCESS

    end function cvLsLinSys
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsDQJac`: difference-quotient Jacobian (dense or band matrix).

    function cvLsDQJac(me, t, y, fy, Jac) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    class(sun_matrix),intent(inout) :: Jac  !! the Jacobian (a dense or band matrix)
    integer :: flag  !! 0 on success, f's flag on failure, `CVLS_ILL_INPUT` for other matrices

    select type (Jac)
    type is (dense_matrix)
        flag = cvLsDenseDQJac(me, t, y, fy, Jac)
    type is (band_matrix)
        flag = cvLsBandDQJac(me, t, y, fy, Jac)
    class default
        call cvProcessError(me, CVLS_ILL_INPUT, 'cvLsDQJac', &
                            'unrecognized matrix type for cvLsDQJac')
        flag = CVLS_ILL_INPUT
    end select

    end function cvLsDQJac
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsDenseDQJac`: dense difference-quotient Jacobian, one column at a time:
!  `J(:,j) = [f(t, y + inc*e_j) - fy] / inc`, with
!  `inc = max(sqrt(uround)*|y_j|, minInc/ewt_j)`.
!
!  (C perturbs the y passed in and restores it; here a copy in `vtemp1`
!  is perturbed, and `vtemp2` holds f.)

    function cvLsDenseDQJac(me, t, y, fy, Jac) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    type(dense_matrix),intent(inout) :: Jac  !! the Jacobian
    integer :: retval  !! 0 on success, or the flag of a failed f evaluation

    real(wp) :: fnorm    !! WRMS norm of fy
    real(wp) :: minInc   !! minimum increment
    real(wp) :: inc      !! increment in y_j
    real(wp) :: inc_inv  !! 1/inc
    real(wp) :: yjsaved  !! saved y_j
    real(wp) :: conj     !! constraint flag of y_j
    real(wp) :: srur     !! sqrt(uround)
    integer(ip) :: j, N  !! column index, number of equations

    retval = 0

    ! obtain the number of columns
    N = Jac%n

    ! work vectors: ytmp (perturbed y), ftemp (f at the perturbed y)
    associate (ytmp => me%vtemp1, ftemp => me%vtemp2)

    ytmp = y

    ! set minimum increment based on uround and norm of f
    srur   = sqrt(me%uround)
    fnorm  = n_vwrmsnorm(fy, me%ewt)
    if (fnorm /= ZERO) then
        minInc = MIN_INC_MULT * abs(me%h) * me%uround * N * fnorm
    else
        minInc = ONE
    end if

    do j = 1, N

        ! generate the jth col of J(tn,y)
        yjsaved = ytmp(j)
        inc = max(srur * abs(yjsaved), minInc / me%ewt(j))

        ! adjust sign(inc) if y_j has an inequality constraint
        if (allocated(me%constraints)) then
            conj = me%constraints(j)
            if (abs(conj) == ONE) then
                if ((yjsaved + inc) * conj < ZERO) inc = -inc
            else if (abs(conj) == TWO) then
                if ((yjsaved + inc) * conj <= ZERO) inc = -inc
            end if
        end if

        ytmp(j) = ytmp(j) + inc

        retval = me%f(me, t, ytmp, ftemp)
        me%nfeDQ = me%nfeDQ + 1
        if (retval /= 0) exit

        ytmp(j) = yjsaved

        inc_inv = ONE / inc
        Jac%data(:,j) = inc_inv * (ftemp - fy)

    end do

    end associate

    end function cvLsDenseDQJac
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsBandDQJac`: band difference-quotient Jacobian. Columns that are more
!  than `ml + mu + 1` apart do not overlap in the band, so they are perturbed
!  together: `min(ml + mu + 1, N)` f evaluations in total.

    function cvLsBandDQJac(me, t, y, fy, Jac) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    type(band_matrix),intent(inout) :: Jac  !! the Jacobian
    integer :: retval  !! 0 on success, or the flag of a failed f evaluation

    real(wp) :: fnorm    !! WRMS norm of fy
    real(wp) :: minInc   !! minimum increment
    real(wp) :: inc      !! increment in y_j
    real(wp) :: inc_inv  !! 1/inc
    real(wp) :: srur     !! sqrt(uround)
    integer(ip) :: group    !! column group (1..ngroups)
    integer(ip) :: i, j     !! row and column indices (0-based)
    integer(ip) :: width    !! spacing of the columns in a group, ml + mu + 1
    integer(ip) :: ngroups  !! number of column groups
    integer(ip) :: i1, i2   !! first and last row in the band of column j
    integer(ip) :: N, mupper, mlower  !! number of equations, upper and lower bandwidths

    retval = 0

    ! obtain matrix dimensions
    N      = Jac%n
    mupper = Jac%mu
    mlower = Jac%ml

    ! C: tmp1 = ftemp, tmp2 = ytemp
    associate (ftemp => me%vtemp1, ytemp => me%vtemp2)

    ! load ytemp with y = predicted y vector
    ytemp = y

    ! set minimum increment based on uround and norm of f
    srur   = sqrt(me%uround)
    fnorm  = n_vwrmsnorm(fy, me%ewt)
    if (fnorm /= ZERO) then
        minInc = MIN_INC_MULT * abs(me%h) * me%uround * N * fnorm
    else
        minInc = ONE
    end if

    ! set bandwidth and number of column groups for band differencing
    width   = mlower + mupper + 1
    ngroups = min(width, N)

    ! loop over column groups
    do group = 1, ngroups

        ! increment all y_j in group (0-based j as in C)
        do j = group - 1, N - 1, width
            inc = max(srur * abs(y(j+1)), minInc / me%ewt(j+1))
            ! adjust sign(inc) if y_j has an inequality constraint
            if (allocated(me%constraints)) call adjust_inc(j)
            ytemp(j+1) = ytemp(j+1) + inc
        end do

        ! evaluate f with incremented y
        retval = me%f(me, t, ytemp, ftemp)
        me%nfeDQ = me%nfeDQ + 1
        if (retval /= 0) exit

        ! restore ytemp, then form and load difference quotients
        do j = group - 1, N - 1, width
            ytemp(j+1) = y(j+1)
            inc = max(srur * abs(y(j+1)), minInc / me%ewt(j+1))
            ! adjust sign(inc) as before
            if (allocated(me%constraints)) call adjust_inc(j)
            inc_inv = ONE / inc
            i1 = max(0_ip, j - mupper)
            i2 = min(j + mlower, N - 1)
            do i = i1, i2
                Jac%data(i - j + Jac%smu, j) = inc_inv * (ftemp(i+1) - fy(i+1))
            end do
        end do

    end do

    end associate

    contains

        subroutine adjust_inc(j)
        !! Reverse the increment of column j (0-based) if `ytemp_j + inc` would
        !! violate the constraint on y_j (C: the `conj` tests in cvLsBandDQJac).
        integer(ip),intent(in) :: j  !! column index (0-based)
        real(wp) :: conj  !! constraint flag of y_j
        conj = me%constraints(j+1)
        if (abs(conj) == ONE) then
            if ((me%vtemp2(j+1) + inc) * conj < ZERO) inc = -inc
        else if (abs(conj) == TWO) then
            if ((me%vtemp2(j+1) + inc) * conj <= ZERO) inc = -inc
        end if
        end subroutine adjust_inc

    end function cvLsBandDQJac
!*****************************************************************************************

    end submodule fvode_core_ls
!*****************************************************************************************

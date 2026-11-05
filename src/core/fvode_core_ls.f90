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
!  `CVodeSetLinearSolver`: attach a linear solver `LS` and (for matrix-based
!  solvers) a matrix `A`. Both are copied into the CVODE object. The internal
!  difference-quotient Jacobian is used unless [[CVodeSetJacFn]] is called.

    module subroutine CVodeSetLinearSolver(me, LS, A, flag)

    class(cvode_solver),intent(inout) :: me
    class(sun_linear_solver),intent(in) :: LS
    class(sun_matrix),intent(in),optional :: A
    integer,intent(out) :: flag

    integer :: LSType
    logical :: iterative    ! is the solver iterative?
    logical :: matrixbased  ! is a matrix structure used?

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

    if (iterative) then
        ! [not yet ported] iterative linear solvers
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetLinearSolver', &
                            'Iterative linear solvers are not supported yet.')
        flag = CVLS_ILL_INPUT
        return
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

    ! [not yet ported] jtimes, linsys and preconditioner functions

    ! initialize counters
    call cvLsInitializeCounters(me)

    ! set default values for the rest of the LS parameters
    me%msbj         = CVLS_MSBJ
    me%jbad         = .true.
    me%dgmax_jbad   = CVLS_DGMAX
    me%eplifac      = CVLS_EPLIN
    me%ls_last_flag = CVLS_SUCCESS

    ! initialize the matrix (the saved Jacobian is allocated in cvLsInitialize)
    if (present(A)) allocate(me%A, source=A)

    ! allocate memory for ytemp and x
    allocate(me%ls_ytemp(me%n), me%ls_x(me%n))

    ! for matrix-based LS, enable solution scaling with BDF
    me%scalesol = (matrixbased .and. me%lmm == CV_BDF)

    flag = CVLS_SUCCESS

    end subroutine CVodeSetLinearSolver
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetJacFn`: set a user Jacobian function, or (if `jac` is not
!  present) go back to the internal difference-quotient Jacobian.

    module subroutine CVodeSetJacFn(me, jac, flag)

    class(cvode_solver),intent(inout) :: me
    procedure(cv_jac_fn),optional :: jac
    integer,intent(out) :: flag

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

    class(cvode_solver),intent(inout) :: me
    integer(li),intent(in) :: msbj
    integer,intent(out) :: flag

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

    class(cvode_solver),intent(inout) :: me
    logical,intent(in) :: onoff
    integer,intent(out) :: flag

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

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: dgmax_jbad
    integer,intent(out) :: flag

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

    class(cvode_solver),intent(inout) :: me
    integer :: flag

    ! set up the Jacobian (matrix-based case)
    if (allocated(me%A)) then
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
    end if

    ! reset counters
    call cvLsInitializeCounters(me)

    ! [not yet ported] jtimes setup; no lsetup without a matrix or preconditioner

    ! call LS initialize routine, and return result
    me%ls_last_flag = me%LS%initialize()
    flag = me%ls_last_flag

    end function cvLsInitialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsReInitialize`: reset the counters.

    module function cvLsReInitialize(me) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer :: flag

    call cvLsInitializeCounters(me)
    flag = CVLS_SUCCESS

    end function cvLsReInitialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsSetup`: decide whether the Jacobian is out of date (`jbad`), update
!  the matrix `A = I - gamma*J` (re-evaluating J if needed), and set up the
!  linear solver (e.g. LU factorization).
!
!  The Jacobian is re-evaluated on the first step, after `msbj` steps, after
!  a convergence failure with a bad Jacobian and little change in gamma, or
!  after any other convergence failure.

    module function cvLsSetup(me, convfail, ypred, fpred, jcur) result(flag)

    class(cvode_solver),intent(inout) :: me
    integer,intent(in) :: convfail
    real(wp),dimension(:),intent(in) :: ypred
    real(wp),dimension(:),intent(in) :: fpred
    logical,intent(out) :: jcur
    integer :: flag

    real(wp) :: dgamma
    integer :: retval

    ! use nst, gamma/gammap, and convfail to set J/P eval. flag jok
    dgamma  = abs((me%gamma / me%gammap) - ONE)
    me%jbad = (me%nst == 0) .or. &
              (me%nst >= me%nstlj + me%msbj) .or. &
              ((convfail == CV_FAIL_BAD_J) .and. (dgamma < me%dgmax_jbad)) .or. &
              (convfail == CV_FAIL_OTHER)

    ! setup the linear system if necessary
    if (allocated(me%A)) then

        ! solver calls linsys to set up the linear system
        retval = cvLsLinSys(me, me%tn, ypred, fpred, .not. me%jbad, jcur, me%gamma)

        ! update counters
        if (jcur) then
            me%nje   = me%nje + 1
            me%nstlj = me%nst
            me%tnlj  = me%tn
        end if

        ! check linsys return flag
        if (retval /= CVLS_SUCCESS) then
            flag = retval
            return
        end if

    else
        ! [not yet ported] matrix-free (preconditioner) case
        jcur = me%jbad
    end if

    ! call LS setup routine -- the LS may call cvLsPSetup, who will
    ! pass the heuristic suggestions above to the user code(s)
    me%ls_last_flag = me%LS%setup(me%A)

    flag = me%ls_last_flag

    end function cvLsSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsSolve`: solve the linear system `A x = b`; b is overwritten by x.
!  For BDF with a matrix-based solver the solution is scaled by
!  `2/(1 + gamrat)` to account for a change in gamma since the last setup.
!
!  Returns 0 on success, 1 for a recoverable failure, -1 otherwise.

    module function cvLsSolve(me, b, weight, ynow, fnow) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),dimension(:),intent(inout) :: b
    real(wp),dimension(:),intent(in) :: weight
    real(wp),dimension(:),intent(in) :: ynow
    real(wp),dimension(:),intent(in) :: fnow
    integer :: flag

    integer :: retval
    real(wp) :: delta

    ! [not yet ported] iterative solvers: tolerance and quick-return test
    delta = ZERO

    ! set initial guess x = 0 to LS
    me%ls_x = ZERO

    ! call solver, and copy x to b
    retval = me%LS%solve(me%A, me%ls_x, b, delta)
    b = me%ls_x

    ! if using a direct or matrix-iterative solver, BDF method, and gamma has
    ! changed, scale the correction to account for change in gamma
    if (me%scalesol .and. me%gamrat /= ONE) b = (TWO / (ONE + me%gamrat)) * b

    ! increment counters nli and ncfl
    ! ([not yet ported] nli for iterative solvers)
    if (retval /= SUN_SUCCESS) me%ncfl = me%ncfl + 1

    ! interpret solver return value
    me%ls_last_flag = retval

    select case (retval)
    case (SUN_SUCCESS)
        flag = 0
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
    case default
        call cvProcessError(me, retval, 'cvLsSolve', &
                            'Unrecognized error return value from SUNLinSolSolve')
        flag = -1
    end select

    end function cvLsSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsFree`: free the linear solver interface memory.

    module subroutine cvLsFree(me)

    class(cvode_solver),intent(inout) :: me

    if (allocated(me%ls_ytemp)) deallocate(me%ls_ytemp)
    if (allocated(me%ls_x))     deallocate(me%ls_x)
    if (allocated(me%savedJ))   deallocate(me%savedJ)
    if (allocated(me%A))        deallocate(me%A)
    if (allocated(me%LS))       deallocate(me%LS)
    me%jac_user => null()

    end subroutine cvLsFree
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvLsInitializeCounters`

    subroutine cvLsInitializeCounters(me)

    class(cvode_solver),intent(inout) :: me

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

    function cvLsLinSys(me, t, y, fy, jok, jcur, gamma) result(flag)

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    real(wp),dimension(:),intent(in) :: y
    real(wp),dimension(:),intent(in) :: fy
    logical,intent(in) :: jok
    logical,intent(out) :: jcur
    real(wp),intent(in) :: gamma
    integer :: flag

    integer :: retval

    ! if jok = .true., use saved copy of J
    if (jok) then

        ! indicate J is not current
        jcur = .false.

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
        jcur = .true.

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

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    real(wp),dimension(:),intent(in) :: y
    real(wp),dimension(:),intent(in) :: fy
    class(sun_matrix),intent(inout) :: Jac
    integer :: flag

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

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    real(wp),dimension(:),intent(in) :: y
    real(wp),dimension(:),intent(in) :: fy
    type(dense_matrix),intent(inout) :: Jac
    integer :: retval

    real(wp) :: fnorm, minInc, inc, inc_inv, yjsaved, srur
    integer(ip) :: j, N

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

        ! [not yet ported] adjust sign(inc) if y_j has an inequality constraint

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

    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    real(wp),dimension(:),intent(in) :: y
    real(wp),dimension(:),intent(in) :: fy
    type(band_matrix),intent(inout) :: Jac
    integer :: retval

    real(wp) :: fnorm, minInc, inc, inc_inv, srur
    integer(ip) :: group, i, j, width, ngroups, i1, i2
    integer(ip) :: N, mupper, mlower

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
            ! [not yet ported] adjust sign(inc) if y_j has an inequality constraint
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
            ! [not yet ported] adjust sign(inc) if y_j has an inequality constraint
            inc_inv = ONE / inc
            i1 = max(0_ip, j - mupper)
            i2 = min(j + mlower, N - 1)
            do i = i1, i2
                Jac%data(i - j + Jac%smu, j) = inc_inv * (ftemp(i+1) - fy(i+1))
            end do
        end do

    end do

    end associate

    end function cvLsBandDQJac
!*****************************************************************************************

    end submodule fvode_core_ls
!*****************************************************************************************

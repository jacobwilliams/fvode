!*****************************************************************************************
!>
!  CVODE diagonal linear solver (SUNDIALS v7.9.0 `cvode_diag.c`).
!
!  Approximates the Jacobian by a diagonal matrix built from one extra
!  f evaluation (a difference quotient along the current Newton direction),
!  and stores `M = (I - gamma*D)^{-1}`. Cheap, and useful when the Jacobian
!  is roughly diagonal.

    submodule (fvode_core) fvode_core_diag

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVDiag`: attach the diagonal linear solver.

    module subroutine CVDiag(me, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(out) :: flag  !! `CVDIAG_SUCCESS` or `CVDIAG_MEM_NULL` (init not called)

    if (.not. me%MallocDone) then
        call cvProcessError(me, CVDIAG_MEM_NULL, 'CVDiag', 'Integrator memory is NULL.')
        flag = CVDIAG_MEM_NULL
        return
    end if

    ! free any existing linear solver attached to CVode
    call cvLFree(me)

    ! set four main function fields in cv_mem
    me%ls_kind       = LS_DIAG
    me%lsetup_exists = .true.
    me%lsolve_exists = .true.

    ! set default values for the rest of the diag parameters
    me%di_last_flag = CVDIAG_SUCCESS

    ! allocate memory for M, bit, and bitcomp
    allocate(me%di_M(me%n), me%di_bit(me%n), me%di_bitcomp(me%n))

    flag = CVDIAG_SUCCESS

    end subroutine CVDiag
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVDiagInit`: reset the counter (also the `lreinit` function).

    module function CVDiagInit(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! always 0

    me%di_nfeDI     = 0
    me%di_last_flag = CVDIAG_SUCCESS
    flag = 0

    end function CVDiagInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVDiagSetup`: form `M = (I - gamma*D)^{-1}`, where D is the diagonal
!  difference-quotient approximation of the Jacobian along the direction
!  `ftemp = h*fpred - zn[1]`:
!
!  `D(i) = [f(ypred + 0.1*rl1*ftemp)(i) - fpred(i)] / (0.1*rl1*ftemp(i))`
!
!  Components with `|ftemp(i)*ewt(i)| < uround` get `D(i) = 0`.
!  Returns 0 on success, 1 for a recoverable failure (M not invertible or
!  f failed recoverably), -1 for an unrecoverable f failure.

    module function CVDiagSetup(me, ypred, fpred, jcur) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: ypred  !! predicted solution y at the new time
    real(wp),dimension(:),intent(in) :: fpred  !! f(tn, ypred)
    logical,intent(out) :: jcur  !! set to .true.: the Jacobian approximation is now current
    integer :: flag  !! 0 on success, 1 recoverable failure, -1 unrecoverable failure

    real(wp) :: r  !! relative size of the perturbation, 0.1*rl1
    integer :: retval  !! return flag of the right-hand side function
    logical :: invOK  !! M could be inverted (no zero diagonal entries)

    jcur = .false.

    ! use vtemp1 and vtemp2 for temporary storage
    associate (ftemp => me%vtemp1, y => me%vtemp2, M => me%di_M, &
               bit => me%di_bit, bitcomp => me%di_bitcomp)

    ! form y with perturbation = FRACT*(func. iter. correction)
    r = FRACT * me%rl1
    ftemp = me%h * fpred - me%zn(:,1)
    y = r * ftemp + ypred

    ! evaluate f at perturbed y
    retval = me%f(me, me%tn, y, M)
    me%di_nfeDI = me%di_nfeDI + 1
    if (retval < 0) then
        call cvProcessError(me, CVDIAG_RHSFUNC_UNRECVR, 'CVDiagSetup', &
            'The right-hand side routine failed in an unrecoverable manner.')
        me%di_last_flag = CVDIAG_RHSFUNC_UNRECVR
        flag = -1
        return
    end if
    if (retval > 0) then
        me%di_last_flag = CVDIAG_RHSFUNC_RECVR
        flag = 1
        return
    end if

    ! construct M = I - gamma*J with J = diag(deltaf_i/deltay_i);
    ! protect against deltay_i being at roundoff level
    M = M - fpred
    M = FRACT * ftemp - me%h * M
    y = ftemp * me%ewt
    bit = merge(ONE, ZERO, abs(y) >= me%uround)     ! N_VCompare
    bitcomp = bit - ONE
    y = ftemp * bit
    y = FRACT * y - bitcomp
    M = M / y
    M = M * bit
    M = M - bitcomp

    ! invert M with test for zero components
    invOK = n_vinvtest(M, M)
    if (.not. invOK) then
        me%di_last_flag = CVDIAG_INV_FAIL
        flag = 1
        return
    end if

    end associate

    ! set jcur = .true., save gamma in gammasv, and return
    jcur = .true.
    me%di_gammasv   = me%gamma
    me%di_last_flag = CVDIAG_SUCCESS
    flag = 0

    end function CVDiagSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVDiagSolve`: `b = M*b`. If gamma has changed since the setup, M is
!  updated first: `M = 1 / (1 + (gamma/gammasv)*(1/M - 1))`.

    module function CVDiagSolve(me, b) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(inout) :: b  !! in: right-hand side; out: solution M*b
    integer :: flag  !! 0 on success, 1 if the updated M is not invertible

    logical :: invOK  !! the updated M could be inverted
    real(wp) :: r  !! gamma ratio, gamma/gammasv

    ! if gamma has changed, update factor in M, and save gamma value
    if (me%di_gammasv /= me%gamma) then
        r = me%gamma / me%di_gammasv
        me%di_M = ONE / me%di_M
        me%di_M = me%di_M - ONE
        me%di_M = r * me%di_M
        me%di_M = me%di_M + ONE
        invOK = n_vinvtest(me%di_M, me%di_M)
        if (.not. invOK) then
            me%di_last_flag = CVDIAG_INV_FAIL
            flag = 1
            return
        end if
        me%di_gammasv = me%gamma
    end if

    ! apply M-inverse to b
    b = b * me%di_M

    me%di_last_flag = CVDIAG_SUCCESS
    flag = 0

    end function CVDiagSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVDiagFree`

    module subroutine CVDiagFree(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    if (allocated(me%di_M))       deallocate(me%di_M)
    if (allocated(me%di_bit))     deallocate(me%di_bit)
    if (allocated(me%di_bitcomp)) deallocate(me%di_bitcomp)

    end subroutine CVDiagFree
!*****************************************************************************************

    end submodule fvode_core_diag
!*****************************************************************************************

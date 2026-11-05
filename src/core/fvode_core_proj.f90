!*****************************************************************************************
!>
!  CVODE projection (SUNDIALS v7.9.0 `cvode_proj.c`).
!
!  For ODEs with invariants (e.g. a DAE written as an ODE), a user function
!  projects the corrected solution of each step onto the constraint
!  manifold, `y_n -> y_n + corr`. The correction is added to the Nordsieck
!  history with the coefficients `proj_p` (set in `cvSetBDF`), and the
!  error estimate can be projected too. Projection is only supported with
!  BDF. It is turned on with [[CVodeSetProjFn]].

    submodule (fvode_core) fvode_core_proj

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetProjFn`: set the user projection function and enable
!  projection (BDF only). The projection memory is created with default
!  settings if needed.

    module subroutine CVodeSetProjFn(me, pfun, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    procedure(cv_proj_fn) :: pfun  !! the projection function
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`

    ! check for compatible method
    if (me%lmm /= CV_BDF) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetProjFn', &
                            'Projection is only supported with BDF methods.')
        flag = CV_ILL_INPUT
        return
    end if

    ! create the projection memory (if necessary)
    call cvProjCreate(me)

    ! user-defined projection
    me%internal_proj = .false.

    ! set the projection function
    me%pfun => pfun

    ! enable projection
    me%proj_enabled = .true.

    flag = CV_SUCCESS

    end subroutine CVodeSetProjFn
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetProjErrEst`: turn the projection of the error estimate on or
!  off (default on).

    module subroutine CVodeSetProjErrEst(me, onoff, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    logical,intent(in) :: onoff  !! project the error estimate?
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`

    if (.not. access_proj_mem(me, 'CVodeSetProjErrEst', flag)) return

    ! set projection error flag
    me%err_proj = onoff

    end subroutine CVodeSetProjErrEst
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetProjFrequency`: project every `freq` steps (< 0: the default,
!  every step; 0: turn projection off).

    module subroutine CVodeSetProjFrequency(me, freq, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(li),intent(in) :: freq  !! projection frequency
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`

    if (.not. access_proj_mem(me, 'CVodeSetProjFrequency', flag)) return

    if (freq < 0) then
        ! restore default
        me%proj_freq    = 1
        me%proj_enabled = .true.
    else if (freq == 0) then
        ! disable projection
        me%proj_freq    = 0
        me%proj_enabled = .false.
    else
        ! enable projection at given frequency
        me%proj_freq    = freq
        me%proj_enabled = .true.
    end if

    end subroutine CVodeSetProjFrequency
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetMaxNumProjFails`: maximum number of projection failures in one
!  step attempt (< 1: the default, 10).

    module subroutine CVodeSetMaxNumProjFails(me, max_fails, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: max_fails  !! maximum failures in one step attempt
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`

    if (.not. access_proj_mem(me, 'CVodeSetMaxNumProjFails', flag)) return

    if (max_fails < 1) then
        me%proj_max_fails = PROJ_MAX_FAILS
    else
        me%proj_max_fails = max_fails
    end if

    end subroutine CVodeSetMaxNumProjFails
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetEpsProj`: tolerance for the projection solve (<= 0: the
!  default, 0.1).

    module subroutine CVodeSetEpsProj(me, eps, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eps  !! projection tolerance
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`

    if (.not. access_proj_mem(me, 'CVodeSetEpsProj', flag)) return

    if (eps <= ZERO) then
        me%eps_proj = PROJ_EPS
    else
        me%eps_proj = eps
    end if

    end subroutine CVodeSetEpsProj
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetProjFailEta`: step size reduction factor after a projection
!  failure (outside (0,1]: the default, 0.25).

    module subroutine CVodeSetProjFailEta(me, eta, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: eta  !! step size reduction factor
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`

    if (.not. access_proj_mem(me, 'CVodeSetProjFailEta', flag)) return

    if ((eta <= ZERO) .or. (eta > ONE)) then
        me%eta_pfail = PROJ_FAIL_ETA
    else
        me%eta_pfail = eta
    end if

    end subroutine CVodeSetProjFailEta
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvDoProjection`: call the user projection function on the corrected
!  solution `y`; `tempv` receives the correction acorP, and, with error
!  projection, `ftemp` holds the projected error estimate (from which
!  `acnrm` is recomputed). On a recoverable failure, the step is
!  restored and retried with a smaller step (`PREDICT_AGAIN`).

    module function cvDoProjection(me, nflag, saved_t, npfail) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(inout) :: nflag  !! set to `PREV_PROJ_FAIL` on a recoverable failure
    real(wp),intent(in) :: saved_t  !! time to restore to on failure
    integer,intent(inout) :: npfail  !! projection failures in this step attempt
    integer :: flag  !! `CV_SUCCESS`, `PREDICT_AGAIN` or a failure flag

    integer :: retval  !! return flag of the projection function

    ! access projection memory
    if (.not. me%proj_mem) then
        call cvProcessError(me, CV_PROJ_MEM_NULL, 'cvDoProjection', 'proj_mem = NULL illegal.')
        flag = CV_PROJ_MEM_NULL
        return
    end if

    ! use tempv to store acorP and, if projecting the error, ftemp to store
    ! errP; copy acor into errP (if projecting the error), and call the
    ! user projection function
    if (me%err_proj) then
        me%ftemp = me%acor
        retval = me%pfun(me, me%tn, me%y, me%tempv, me%eps_proj, me%ftemp)
    else
        retval = me%pfun(me, me%tn, me%y, me%tempv, me%eps_proj)
    end if
    me%nproj = me%nproj + 1

    ! this is not the first projection anymore
    me%first_proj = .false.

    ! check the return value
    if (retval == CV_SUCCESS) then
        ! recompute acnrm to be used in error test (if projecting the error)
        if (me%err_proj) me%acnrm = n_vwrmsnorm(me%ftemp, me%ewt)
        ! the projection was successful, return now
        me%proj_applied = .true.
        flag = CV_SUCCESS
        return
    end if

    ! the projection failed, update the return value
    if (retval < 0) retval = CV_PROJFUNC_FAIL
    if (retval > 0) retval = PROJFUNC_RECVR

    ! increment cumulative failure count and restore zn
    me%npfails = me%npfails + 1
    call cvRestore(me, saved_t)

    ! return if failed unrecoverably
    if (retval == CV_PROJFUNC_FAIL) then
        flag = CV_PROJFUNC_FAIL
        return
    end if

    ! recoverable failure, increment failure count for this step attempt
    npfail = npfail + 1
    me%etamax = ONE

    ! check for maximum number of failures or |h| = hmin
    if ((abs(me%h) <= me%hmin * ONEPSM) .or. (npfail == me%proj_max_fails)) then
        if (retval == PROJFUNC_RECVR) then
            flag = CV_REPTD_PROJFUNC_ERR
            return
        end if
    end if

    ! reduce step size; return to reattempt the step
    me%eta = max(me%eta_pfail, me%hmin / abs(me%h))
    nflag  = PREV_PROJ_FAIL
    call cvRescale(me)
    flag = PREDICT_AGAIN

    end function cvDoProjection
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvProjInit`: reset the projection flags and counters (at the first
!  call to `CVode` after `CVodeInit`/`CVodeReInit`).

    module function cvProjInit(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`

    ! check if projection memory exists
    if (.not. me%proj_mem) then
        flag = CV_PROJ_MEM_NULL
        return
    end if

    ! reset flags and counters
    me%first_proj = .true.
    me%nstlprj    = 0
    me%nproj      = 0
    me%npfails    = 0

    flag = CV_SUCCESS

    end function cvProjInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvProjFree`: free the projection memory (and disable projection).

    module subroutine cvProjFree(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    me%proj_mem = .false.
    me%proj_enabled = .false.
    me%pfun => null()

    end subroutine cvProjFree
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvProjCreate`: create the projection memory with the default settings
!  (`cvProjSetDefaults`), unless it already exists.

    subroutine cvProjCreate(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    if (me%proj_mem) return

    me%proj_mem = .true.
    me%internal_proj  = .true.
    me%err_proj       = .true.
    me%first_proj     = .true.
    me%proj_freq      = 1
    me%nstlprj        = 0
    me%proj_max_fails = PROJ_MAX_FAILS
    me%pfun => null()
    me%eps_proj       = PROJ_EPS
    me%eta_pfail      = PROJ_FAIL_ETA
    me%nproj          = 0
    me%npfails        = 0

    end subroutine cvProjCreate
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvAccessProjMem`: check that the projection memory exists (error
!  `CV_PROJ_MEM_NULL` otherwise).

    function access_proj_mem(me, fname, flag) result(ok)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    character(len=*),intent(in) :: fname  !! name of the calling function
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_PROJ_MEM_NULL`
    logical :: ok  !! does the projection memory exist?

    ok = me%proj_mem
    if (ok) then
        flag = CV_SUCCESS
    else
        call cvProcessError(me, CV_PROJ_MEM_NULL, fname, 'proj_mem = NULL illegal.')
        flag = CV_PROJ_MEM_NULL
    end if

    end function access_proj_mem
!*****************************************************************************************

    end submodule fvode_core_proj
!*****************************************************************************************

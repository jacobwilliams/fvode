!*****************************************************************************************
!>
!  CVODE banded preconditioner (SUNDIALS v7.9.0 `cvode_bandpre.c`).
!
!  For use with a Krylov linear solver: the preconditioner is
!  `P = I - gamma*J`, where J is a banded difference-quotient approximation
!  of the Jacobian with (possibly reduced) half-bandwidths `mu` and `ml`.
!  P is factored with the band LU solver. It is attached with
!  [[CVBandPrecInit]] after [[CVodeSetLinearSolver]].

    submodule (fvode_core) fvode_core_bandpre

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVBandPrecInit`: allocate and initialize the band preconditioner and
!  attach it to the (iterative) linear solver. The half-bandwidths are
!  clipped to `0 <= mu, ml <= N-1`.

    module subroutine CVBandPrecInit(me, N, mu, ml, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer(ip),intent(in) :: N  !! problem size
    integer(ip),intent(in) :: mu  !! upper half-bandwidth of the preconditioner
    integer(ip),intent(in) :: ml  !! lower half-bandwidth of the preconditioner
    integer,intent(out) :: flag  !! `CVLS_SUCCESS`, `CVLS_LMEM_NULL`, `CVLS_ILL_INPUT` or `CVLS_SUNLS_FAIL`

    integer(ip) :: mup  !! clipped upper half-bandwidth
    integer(ip) :: mlp  !! clipped lower half-bandwidth
    integer(ip) :: storagemu  !! storage upper bandwidth of the preconditioner matrix

    ! test if the CVLS linear solver interface has been attached
    if (me%ls_kind /= LS_GENERIC) then
        call cvProcessError(me, CVLS_LMEM_NULL, 'CVBandPrecInit', 'Linear solver memory is NULL.')
        flag = CVLS_LMEM_NULL
        return
    end if

    ! make sure the preconditioner data is free from any previous allocations
    call cvBandPrecFree(me)

    ! load bandwidths into the preconditioner data
    me%bp_N  = N
    mup = min(N - 1, max(0_ip, mu))
    mlp = min(N - 1, max(0_ip, ml))
    me%bp_mu = mup
    me%bp_ml = mlp

    ! initialize nfeBP counter
    me%bp_nfeBP = 0

    ! allocate memory for saved banded Jacobian approximation
    me%bp_savedJ = band_matrix(N, mup, mlp, smu=mup)

    ! allocate memory for banded preconditioner
    storagemu = min(N - 1, mup + mlp)
    me%bp_savedP = band_matrix(N, mup, mlp, smu=storagemu)

    ! allocate memory for banded linear solver
    me%bp_LS = sun_linsol_band(N)

    ! allocate memory for temporary vectors
    allocate(me%bp_tmp1(N), me%bp_tmp2(N))

    ! initialize band linear solver object
    if (me%bp_LS%initialize() /= SUN_SUCCESS) then
        call cvBandPrecFree(me)
        call cvProcessError(me, CVLS_SUNLS_FAIL, 'CVBandPrecInit', &
                            'An error arose from a SUNBandLinearSolver routine.')
        flag = CVLS_SUNLS_FAIL
        return
    end if

    ! attach preconditioner solve and setup functions
    ! (C: CVodeSetPreconditioner(cvode_mem, cvBandPrecSetup, cvBandPrecSolve))
    me%pset_user => null()
    me%psolve_user => null()
    me%pre_kind = PRE_BANDPRE
    if (.not. me%LS%supports_preconditioner()) then
        call cvProcessError(me, CVLS_ILL_INPUT, 'CVodeSetPreconditioner', &
            'SUNLinearSolver object does not support user-supplied preconditioning')
        flag = CVLS_ILL_INPUT
        return
    end if
    me%LS%has_psetup = .true.
    me%LS%has_psolve = .true.

    flag = CVLS_SUCCESS

    end subroutine CVBandPrecInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvBandPrecSetup`: form `P = I - gamma*J` with a banded
!  difference-quotient J (reused from the saved copy if `jok`), and
!  LU-factor it.

    module function cvBandPrecSetup(me, t, y, fy, jok, jcur, gamma) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    logical,intent(in) :: jok  !! may the saved Jacobian be reused?
    logical,intent(out) :: jcur  !! was the Jacobian recomputed?
    real(wp),intent(in) :: gamma  !! the scalar gamma
    integer :: retval  !! 0 on success, > 0 recoverable, < 0 unrecoverable

    integer :: iflag  !! return flag of the matrix copy

    if (jok) then

        ! if jok = SUNTRUE, use saved copy of J
        jcur = .false.
        call me%bp_savedJ%copy_to(me%bp_savedP, iflag)
        if (iflag < 0) then
            call cvProcessError(me, -1, 'cvBandPrecSetup', 'An error arose from a SUNBandMatrix routine.')
            retval = -1
            return
        end if
        if (iflag > 0) then
            retval = 1
            return
        end if

    else

        ! if jok = SUNFALSE, call cvBandPrecDQJac for new J value
        jcur = .true.
        call me%bp_savedJ%zero()

        retval = cvBandPrecDQJac(me, t, y, fy)
        if (retval < 0) then
            call cvProcessError(me, -1, 'cvBandPrecSetup', &
                                'The right-hand side routine failed in an unrecoverable manner.')
            retval = -1
            return
        end if
        if (retval > 0) then
            retval = 1
            return
        end if

        call me%bp_savedJ%copy_to(me%bp_savedP, iflag)
        if (iflag < 0) then
            call cvProcessError(me, -1, 'cvBandPrecSetup', 'An error arose from a SUNBandMatrix routine.')
            retval = -1
            return
        end if
        if (iflag > 0) then
            retval = 1
            return
        end if

    end if

    ! scale and add identity to get savedP = I - gamma*J
    call me%bp_savedP%scale_add_i(-gamma)

    ! do LU factorization of matrix and return error flag
    retval = me%bp_LS%setup(me%bp_savedP)

    end function cvBandPrecSetup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvBandPrecSolve`: solve `P z = r` with the band LU factors of P.

    module function cvBandPrecSolve(me, r, z) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),dimension(:),intent(in) :: r  !! the right-hand side
    real(wp),dimension(:),intent(out) :: z  !! the solution of P z = r
    integer :: retval  !! the band solver flag

    ! call banded solver object to do the work
    retval = me%bp_LS%solve(me%bp_savedP, z, r, ZERO)

    end function cvBandPrecSolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvBandPrecFree`: release the band preconditioner data.

    module subroutine cvBandPrecFree(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    if (allocated(me%bp_savedJ%data)) deallocate(me%bp_savedJ%data)
    if (allocated(me%bp_savedP%data)) deallocate(me%bp_savedP%data)
    if (allocated(me%bp_tmp1)) deallocate(me%bp_tmp1)
    if (allocated(me%bp_tmp2)) deallocate(me%bp_tmp2)
    if (me%pre_kind == PRE_BANDPRE) me%pre_kind = PRE_NONE

    end subroutine cvBandPrecFree
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvBandPrecDQJac`: banded difference-quotient approximation of the
!  Jacobian into `bp_savedJ`. Columns more than `ml + mu + 1` apart are
!  perturbed together, so `min(ml + mu + 1, N)` f evaluations are used
!  (counted in `bp_nfeBP`).

    function cvBandPrecDQJac(me, t, y, fy) result(retval)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    integer :: retval  !! 0 on success, or the flag of a failed f evaluation

    real(wp) :: fnorm    !! WRMS norm of fy
    real(wp) :: minInc   !! minimum increment
    real(wp) :: inc      !! increment in y_j
    real(wp) :: inc_inv  !! 1/inc
    real(wp) :: srur     !! sqrt(uround)
    real(wp) :: yj       !! y_j
    integer(ip) :: group    !! column group (1..ngroups)
    integer(ip) :: i, j     !! row and column indices (0-based)
    integer(ip) :: width    !! spacing of the columns in a group, ml + mu + 1
    integer(ip) :: ngroups  !! number of column groups
    integer(ip) :: i1, i2   !! first and last row in the band of column j

    retval = 0

    ! C: ftemp = tmp1, ytemp = tmp2
    associate (ftemp => me%bp_tmp1, ytemp => me%bp_tmp2, Jac => me%bp_savedJ)

    ! load ytemp with y = predicted y vector
    ytemp = y

    ! set minimum increment based on uround and norm of f
    srur  = sqrt(me%uround)
    fnorm = n_vwrmsnorm(fy, me%ewt)
    if (fnorm /= ZERO) then
        minInc = MIN_INC_MULT * abs(me%h) * me%uround * me%bp_N * fnorm
    else
        minInc = ONE
    end if

    ! set bandwidth and number of column groups for band differencing
    width   = me%bp_ml + me%bp_mu + 1
    ngroups = min(width, me%bp_N)

    do group = 1, ngroups

        ! increment all y_j in group
        do j = group - 1, me%bp_N - 1, width
            inc = max(srur * abs(y(j+1)), minInc / me%ewt(j+1))
            yj  = y(j+1)
            ! adjust sign(inc) again if yj has an inequality constraint
            if (allocated(me%constraints)) call adjust_inc(j)
            ytemp(j+1) = ytemp(j+1) + inc
        end do

        ! evaluate f with incremented y
        retval = me%f(me, t, ytemp, ftemp)
        me%bp_nfeBP = me%bp_nfeBP + 1
        if (retval /= 0) exit

        ! restore ytemp, then form and load difference quotients
        do j = group - 1, me%bp_N - 1, width
            ytemp(j+1) = y(j+1)
            yj  = y(j+1)
            inc = max(srur * abs(y(j+1)), minInc / me%ewt(j+1))
            ! adjust sign(inc) as before
            if (allocated(me%constraints)) call adjust_inc(j)
            inc_inv = ONE / inc
            i1 = max(0_ip, j - me%bp_mu)
            i2 = min(j + me%bp_ml, me%bp_N - 1)
            do i = i1, i2
                Jac%data(i - j + Jac%smu, j) = inc_inv * (ftemp(i+1) - fy(i+1))
            end do
        end do

    end do

    end associate

    contains

        subroutine adjust_inc(j)
        !! Reverse the increment of column j (0-based) if `y_j + inc` would
        !! violate the constraint on y_j.
        integer(ip),intent(in) :: j  !! column index (0-based)
        real(wp) :: conj  !! constraint flag of y_j
        conj = me%constraints(j+1)
        if (abs(conj) == ONE) then
            if ((yj + inc) * conj < ZERO) inc = -inc
        else if (abs(conj) == TWO) then
            if ((yj + inc) * conj <= ZERO) inc = -inc
        end if
        end subroutine adjust_inc

    end function cvBandPrecDQJac
!*****************************************************************************************

    end submodule fvode_core_bandpre
!*****************************************************************************************

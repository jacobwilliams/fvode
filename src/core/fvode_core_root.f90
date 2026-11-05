!*****************************************************************************************
!>
!  CVODE root finding (SUNDIALS v7.9.0 `cvode.c`: `CVodeRootInit`,
!  `cvRcheck1`, `cvRcheck2`, `cvRcheck3`, `cvRootfind`; `cvode_io.c`:
!  `CVodeSetRootDirection`, `CVodeSetNoInactiveRootWarn`, `CVodeGetRootInfo`).
!
!  During the integration CVODE looks for roots of a set of user functions
!  `g_i(t, y(t))`. After each step, the root functions are evaluated at the
!  end of the step; a sign change locates a root, which is then refined by
!  the Illinois algorithm (a modified secant method) on the interpolated
!  solution, to within `ttol = 100*uround*(|tn| + |h|)`.

    submodule (fvode_core) fvode_core_root

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeRootInit`: initialize the root finding: `nrtfn` root functions,
!  evaluated by `g`. `nrtfn <= 0` turns root finding off.

    module subroutine CVodeRootInit(me, nrtfn, g, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(in) :: nrtfn  !! number of root functions (<= 0: turn root finding off)
    procedure(cv_root_fn),optional :: g  !! the root function
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`

    integer :: nrt  !! number of root functions (>= 0)

    nrt = max(nrtfn, 0)

    ! if rerunning CVodeRootInit() with a different number of root
    ! functions (changing number of gfun components), then free
    ! currently held memory resources
    if ((nrt /= me%nrtfn) .and. (me%nrtfn > 0)) call cvFreeRoots(me)

    ! set variable values in CVode memory block
    ! (if nrtfn = 0, then return)
    if (nrt == 0) then
        me%nrtfn = nrt
        me%gfun => null()
        flag = CV_SUCCESS
        return
    end if

    ! if rerunning CVodeRootInit() with the same number of root functions
    ! (not changing number of gfun components), then check if the root
    ! function argument has changed
    ! (if gfun != NULL, then return)
    if (nrt == me%nrtfn) then
        if (present(g)) then
            me%gfun => g
            flag = CV_SUCCESS
        else
            ! (C frees the arrays but keeps nrtfn, which would leave root
            ! finding on with freed memory; FVODE turns it off)
            call cvFreeRoots(me)
            me%nrtfn = 0
            me%gfun => null()
            call cvProcessError(me, CV_ILL_INPUT, 'CVodeRootInit', 'g = NULL illegal.')
            flag = CV_ILL_INPUT
        end if
        return
    end if

    ! set variable values in CVode memory block
    me%nrtfn = nrt
    if (.not. present(g)) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeRootInit', 'g = NULL illegal.')
        flag = CV_ILL_INPUT
        return
    end if
    me%gfun => g

    ! allocate necessary memory and return
    allocate(me%glo(nrt), me%ghi(nrt), me%grout(nrt), me%iroots(nrt), &
             me%rootdir(nrt), me%gactive(nrt))

    ! set default values for rootdir (both directions)
    me%rootdir = 0

    ! set default values for gactive (all active)
    me%gactive = .true.

    flag = CV_SUCCESS

    end subroutine CVodeRootInit
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetRootDirection`: the direction of zero-crossings to be located
!  for each root function: +1 only increasing, -1 only decreasing, 0 both.

    module subroutine CVodeSetRootDirection(me, rootdir, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,dimension(:),intent(in) :: rootdir  !! direction for each root function
    integer,intent(out) :: flag  !! `CV_SUCCESS` or `CV_ILL_INPUT`

    if (me%nrtfn == 0) then
        call cvProcessError(me, CV_ILL_INPUT, 'CVodeSetRootDirection', &
                            'Rootfinding was not initialized.')
        flag = CV_ILL_INPUT
        return
    end if

    me%rootdir = rootdir(1:me%nrtfn)

    flag = CV_SUCCESS

    end subroutine CVodeSetRootDirection
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeSetNoInactiveRootWarn`: do not warn about root functions that are
!  identically zero at the start of the integration.

    module subroutine CVodeSetNoInactiveRootWarn(me, flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer,intent(out) :: flag  !! `CV_SUCCESS`

    me%mxgnull = 0
    flag = CV_SUCCESS

    end subroutine CVodeSetNoInactiveRootWarn
!*****************************************************************************************

!*****************************************************************************************
!>
!  `CVodeGetRootInfo`: after a `CV_ROOT_RETURN`, which root functions have
!  a root: `rootsfound(i)` = +1 (g_i increasing), -1 (decreasing), 0 (no root).

    module subroutine CVodeGetRootInfo(me, rootsfound)

    class(cvode_solver),intent(in) :: me  !! the CVODE solver
    integer,dimension(:),intent(out) :: rootsfound  !! root information for each root function

    rootsfound(1:me%nrtfn) = me%iroots

    end subroutine CVodeGetRootInfo
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvRcheck1`: called at the first call to `CVode`, to check for zeros of
!  the root functions at t0 (and just after it). Root functions that are
!  zero at t0 are made inactive; if a function is still zero at a slightly
!  later time `t0 + smallh`, it stays inactive, otherwise it is reactivated
!  with its value there.

    module function cvRcheck1(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! `CV_SUCCESS` or `CV_RTFUNC_FAIL`

    integer :: i  !! root function index
    integer :: retval  !! return value of the root function
    real(wp) :: smallh  !! small step from t0 used to evaluate g again
    real(wp) :: hratio  !! smallh / h
    real(wp) :: tplus  !! t0 + smallh
    logical :: zroot  !! is some g_i zero at t0?

    me%iroots = 0
    me%tlo  = me%tn
    me%ttol = (abs(me%tn) + abs(me%h)) * me%uround * HUNDRED

    ! evaluate g at initial t and check for zero values
    retval = me%gfun(me, me%tlo, me%zn(:,0), me%glo)
    me%nge = 1
    if (retval /= 0) then
        flag = CV_RTFUNC_FAIL
        return
    end if

    zroot = .false.
    do i = 1, me%nrtfn
        if (abs(me%glo(i)) == ZERO) then
            zroot = .true.
            me%gactive(i) = .false.
        end if
    end do
    if (.not. zroot) then
        flag = CV_SUCCESS
        return
    end if

    ! some g_i is zero at t0; look at g at t0+(small increment)
    hratio = max(me%ttol / abs(me%h), PT1)
    smallh = hratio * me%h
    tplus  = me%tlo + smallh
    me%y = me%zn(:,0) + hratio * me%zn(:,1)
    retval = me%gfun(me, tplus, me%y, me%ghi)
    me%nge = me%nge + 1
    if (retval /= 0) then
        flag = CV_RTFUNC_FAIL
        return
    end if

    ! we check now only the components of g which were exactly 0.0 at t0
    ! to see if we can 'activate' them
    do i = 1, me%nrtfn
        if (.not. me%gactive(i) .and. abs(me%ghi(i)) /= ZERO) then
            me%gactive(i) = .true.
            me%glo(i) = me%ghi(i)
        end if
    end do

    flag = CV_SUCCESS

    end function cvRcheck1
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvRcheck2`: called at the start of a `CVode` call after a root was
!  returned. Re-evaluates g at the last root `tlo`: if some g_i is zero
!  there, g is evaluated at `tlo + smallh`; a g_i that is zero at both points
!  is a close pair of roots (`CLOSERT`), and a new zero at `tlo + smallh` is
!  returned as a root (`RTFOUND`).

    module function cvRcheck2(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! `CV_SUCCESS`, `RTFOUND`, `CLOSERT` or `CV_RTFUNC_FAIL`

    integer :: i  !! root function index
    integer :: retval  !! return value of the root function
    integer :: ier  !! return flag of CVodeGetDky (not used, as in C)
    real(wp) :: smallh  !! small step from tlo
    real(wp) :: hratio  !! smallh / h
    real(wp) :: tplus  !! tlo + smallh
    logical :: zroot  !! is some active g_i zero?

    if (me%irfnd == 0) then
        flag = CV_SUCCESS
        return
    end if

    call CVodeGetDky(me, me%tlo, 0, me%y, ier)
    retval = me%gfun(me, me%tlo, me%y, me%glo)
    me%nge = me%nge + 1
    if (retval /= 0) then
        flag = CV_RTFUNC_FAIL
        return
    end if

    zroot = .false.
    me%iroots = 0
    do i = 1, me%nrtfn
        if (.not. me%gactive(i)) cycle
        if (abs(me%glo(i)) == ZERO) then
            zroot = .true.
            me%iroots(i) = 1
        end if
    end do
    if (.not. zroot) then
        flag = CV_SUCCESS
        return
    end if

    ! one or more g_i has a zero at tlo. Check g at tlo+smallh
    me%ttol = (abs(me%tn) + abs(me%h)) * me%uround * HUNDRED
    if (me%h > ZERO) then
        smallh = me%ttol
    else
        smallh = -me%ttol
    end if
    tplus = me%tlo + smallh
    if ((tplus - me%tn) * me%h >= ZERO) then
        hratio = smallh / me%h
        me%y = me%y + hratio * me%zn(:,1)
    else
        call CVodeGetDky(me, tplus, 0, me%y, ier)
    end if
    retval = me%gfun(me, tplus, me%y, me%ghi)
    me%nge = me%nge + 1
    if (retval /= 0) then
        flag = CV_RTFUNC_FAIL
        return
    end if

    ! check for close roots (error return), for a new zero at tlo+smallh,
    ! and for a g_i that changed from zero to nonzero
    zroot = .false.
    do i = 1, me%nrtfn
        if (.not. me%gactive(i)) cycle
        if (abs(me%ghi(i)) == ZERO) then
            if (me%iroots(i) == 1) then
                flag = CLOSERT
                return
            end if
            zroot = .true.
            me%iroots(i) = 1
        else
            if (me%iroots(i) == 1) me%glo(i) = me%ghi(i)
        end if
    end do
    if (zroot) then
        flag = RTFOUND
    else
        flag = CV_SUCCESS
    end if

    end function cvRcheck2
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvRcheck3`: called after a step, to look for a root in the interval
!  (tlo, thi], where thi = tn (or tout, in `CV_NORMAL` mode, if tout was
!  passed). On a root, tlo is set to the root location and y to the
!  interpolated solution there.

    module function cvRcheck3(me, tout, itask) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    real(wp),intent(in) :: tout  !! the user's output time
    integer,intent(in) :: itask  !! `CV_NORMAL` or `CV_ONE_STEP`
    integer :: flag  !! `CV_SUCCESS`, `RTFOUND` or `CV_RTFUNC_FAIL`

    integer :: i  !! root function index
    integer :: ier  !! return flag of cvRootfind
    integer :: retval  !! return value of the root function
    integer :: idky  !! return flag of CVodeGetDky (not used, as in C)

    ! set thi = tn or tout, whichever comes first; set y = y(thi)
    if (itask == CV_ONE_STEP) then
        me%thi = me%tn
        me%y = me%zn(:,0)
    end if
    if (itask == CV_NORMAL) then
        if ((tout - me%tn) * me%h >= ZERO) then
            me%thi = me%tn
            me%y = me%zn(:,0)
        else
            me%thi = tout
            call CVodeGetDky(me, me%thi, 0, me%y, idky)
        end if
    end if

    ! set ghi = g(thi) and call cvRootfind to search (tlo,thi) for roots
    retval = me%gfun(me, me%thi, me%y, me%ghi)
    me%nge = me%nge + 1
    if (retval /= 0) then
        flag = CV_RTFUNC_FAIL
        return
    end if

    me%ttol = (abs(me%tn) + abs(me%h)) * me%uround * HUNDRED
    ier = cvRootfind(me)
    if (ier == CV_RTFUNC_FAIL) then
        flag = CV_RTFUNC_FAIL
        return
    end if
    do i = 1, me%nrtfn
        if (.not. me%gactive(i) .and. me%grout(i) /= ZERO) me%gactive(i) = .true.
    end do
    me%tlo = me%trout
    me%glo = me%grout

    ! if no root found, return CV_SUCCESS
    if (ier == CV_SUCCESS) then
        flag = CV_SUCCESS
        return
    end if

    ! if a root was found, interpolate to get y(trout) and return
    call CVodeGetDky(me, me%trout, 0, me%y, idky)
    flag = RTFOUND

    end function cvRcheck3
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvRootfind`: search the interval (tlo, thi] for roots of the active
!  root functions (in the requested directions), given `glo = g(tlo)` and
!  `ghi = g(thi)`.
!
!  Uses the Illinois algorithm, a modified secant method: each iterate is
!  the secant point of the function with the largest relative sign change;
!  if the same end of the bracket is kept twice, the weight `alph` of the
!  other end is halved (or doubled), to guarantee convergence. The iteration
!  stops when the bracket is shorter than `ttol`. On return, `trout` and
!  `grout` hold the root (or thi) and the g values there, and
!  `iroots(i)` = +-1 marks the functions with a root (with the direction).
!
!  Returns `CV_SUCCESS` (no root), `RTFOUND` or `CV_RTFUNC_FAIL`.

    function cvRootfind(me) result(flag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: flag  !! `CV_SUCCESS`, `RTFOUND` or `CV_RTFUNC_FAIL`

    real(wp) :: alph  !! weight of glo in the secant formula (Illinois modification)
    real(wp) :: tmid  !! next trial point
    real(wp) :: gfrac  !! relative size of a sign change |g(thi)/(g(thi) - g(tlo))|
    real(wp) :: maxfrac  !! largest gfrac
    real(wp) :: fracint  !! bracket length in units of ttol
    real(wp) :: fracsub  !! fraction of the bracket used when tmid is too close to an end
    integer :: i  !! root function index
    integer :: retval  !! return value of the root function
    integer :: imax  !! index of the function with the largest gfrac
    integer :: side  !! which end of the bracket was replaced last: 1 = thi, 2 = tlo
    integer :: sideprev  !! previous value of side
    integer :: idky  !! return flag of CVodeGetDky (not used, as in C)
    logical :: zroot  !! is some active g_i exactly zero (with the right direction)?
    logical :: sgnchg  !! does some active g_i change sign (with the right direction)?

    imax = 1

    ! first check for change in sign in ghi or for a zero in ghi
    maxfrac = ZERO
    zroot   = .false.
    sgnchg  = .false.
    do i = 1, me%nrtfn
        if (.not. me%gactive(i)) cycle
        if (abs(me%ghi(i)) == ZERO) then
            if (me%rootdir(i) * me%glo(i) <= ZERO) zroot = .true.
        else
            if (differentsign(me%glo(i), me%ghi(i)) .and. &
                (me%rootdir(i) * me%glo(i) <= ZERO)) then
                gfrac = abs(me%ghi(i) / (me%ghi(i) - me%glo(i)))
                if (gfrac > maxfrac) then
                    sgnchg  = .true.
                    maxfrac = gfrac
                    imax    = i
                end if
            end if
        end if
    end do

    ! if no sign change was found, reset trout and grout. Then return
    ! CV_SUCCESS if no zero was found, or set iroots and return RTFOUND
    if (.not. sgnchg) then
        me%trout = me%thi
        me%grout = me%ghi
        if (.not. zroot) then
            flag = CV_SUCCESS
            return
        end if
        do i = 1, me%nrtfn
            me%iroots(i) = 0
            if (.not. me%gactive(i)) cycle
            if ((abs(me%ghi(i)) == ZERO) .and. (me%rootdir(i) * me%glo(i) <= ZERO)) then
                me%iroots(i) = merge(-1, 1, me%glo(i) > 0)
            end if
        end do
        flag = RTFOUND
        return
    end if

    ! initialize alph to avoid compiler warning
    alph = ONE

    ! a sign change was found. Loop to locate nearest root
    side     = 0
    sideprev = -1
    do  ! looping point

        ! if interval size is already less than tolerance ttol, break
        if (abs(me%thi - me%tlo) <= me%ttol) exit

        ! set weight alph.
        ! On the first two passes, set alph = 1. Thereafter, reset alph
        ! according to the side (low vs high) of the subinterval in which
        ! the sign change was found in the previous two passes.
        ! If the sides were opposite, set alph = 1.
        ! If the sides were the same, then double alph (if high side),
        ! or halve alph (if low side).
        ! The next guess tmid is the secant method value if alph = 1, but
        ! is closer to tlo if alph < 1, and closer to thi if alph > 1.
        if (sideprev == side) then
            if (side == 2) then
                alph = alph * TWO
            else
                alph = alph * HALF
            end if
        else
            alph = ONE
        end if

        ! set next root approximation tmid and get g(tmid).
        ! If tmid is too close to tlo or thi, adjust it inward,
        ! by a fractional distance that is between 0.1 and 0.5.
        tmid = me%thi - (me%thi - me%tlo) * me%ghi(imax) / &
                        (me%ghi(imax) - alph * me%glo(imax))
        if (abs(tmid - me%tlo) < HALF * me%ttol) then
            fracint = abs(me%thi - me%tlo) / me%ttol
            if (fracint > FIVE) then
                fracsub = PT1
            else
                fracsub = HALF / fracint
            end if
            tmid = me%tlo + fracsub * (me%thi - me%tlo)
        end if
        if (abs(me%thi - tmid) < HALF * me%ttol) then
            fracint = abs(me%thi - me%tlo) / me%ttol
            if (fracint > FIVE) then
                fracsub = PT1
            else
                fracsub = HALF / fracint
            end if
            tmid = me%thi - fracsub * (me%thi - me%tlo)
        end if

        call CVodeGetDky(me, tmid, 0, me%y, idky)
        retval = me%gfun(me, tmid, me%y, me%grout)
        me%nge = me%nge + 1
        if (retval /= 0) then
            flag = CV_RTFUNC_FAIL
            return
        end if

        ! check to see in which subinterval g changes sign, and reset imax.
        ! Set side = 1 if sign change is on low side, or 2 if on high side.
        maxfrac  = ZERO
        zroot    = .false.
        sgnchg   = .false.
        sideprev = side
        do i = 1, me%nrtfn
            if (.not. me%gactive(i)) cycle
            if (abs(me%grout(i)) == ZERO) then
                if (me%rootdir(i) * me%glo(i) <= ZERO) zroot = .true.
            else
                if (differentsign(me%glo(i), me%grout(i)) .and. &
                    (me%rootdir(i) * me%glo(i) <= ZERO)) then
                    gfrac = abs(me%grout(i) / (me%grout(i) - me%glo(i)))
                    if (gfrac > maxfrac) then
                        sgnchg  = .true.
                        maxfrac = gfrac
                        imax    = i
                    end if
                end if
            end if
        end do
        if (sgnchg) then
            ! sign change found in (tlo,tmid); replace thi with tmid
            me%thi = tmid
            me%ghi = me%grout
            side = 1
            ! stop at root thi if converged; otherwise loop
            if (abs(me%thi - me%tlo) <= me%ttol) exit
            cycle  ! return to looping point
        end if

        if (zroot) then
            ! no sign change in (tlo,tmid), but g = 0 at tmid; return root tmid
            me%thi = tmid
            me%ghi = me%grout
            exit
        end if

        ! no sign change in (tlo,tmid), and no zero at tmid.
        ! Sign change must be in (tmid,thi). Replace tlo with tmid.
        me%tlo = tmid
        me%glo = me%grout
        side = 2
        ! stop at root thi if converged; otherwise loop back
        if (abs(me%thi - me%tlo) <= me%ttol) exit

    end do  ! end of root-search loop

    ! reset trout and grout, set iroots, and return RTFOUND
    me%trout = me%thi
    do i = 1, me%nrtfn
        me%grout(i)  = me%ghi(i)
        me%iroots(i) = 0
        if (.not. me%gactive(i)) cycle
        if ((abs(me%ghi(i)) == ZERO) .and. (me%rootdir(i) * me%glo(i) <= ZERO)) then
            me%iroots(i) = merge(-1, 1, me%glo(i) > 0)
        end if
        if (differentsign(me%glo(i), me%ghi(i)) .and. (me%rootdir(i) * me%glo(i) <= ZERO)) then
            me%iroots(i) = merge(-1, 1, me%glo(i) > 0)
        end if
    end do
    flag = RTFOUND

    end function cvRootfind
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRdifferentsign(x, y)`: do x and y have different signs? Compares the
!  sign bits (so +0 and -0 have different signs, as in C).

    elemental function differentsign(x, y) result(different)

    real(wp),intent(in) :: x  !! first value
    real(wp),intent(in) :: y  !! second value
    logical :: different  !! true if the signs of x and y differ

    different = sign(ONE, x) /= sign(ONE, y)

    end function differentsign
!*****************************************************************************************

    end submodule fvode_core_root
!*****************************************************************************************

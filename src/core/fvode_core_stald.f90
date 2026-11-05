!*****************************************************************************************
!>
!  CVODE BDF stability limit detection (STALD) (SUNDIALS v7.9.0 `cvode.c`:
!  `cvBDFStab`, `cvSLdet`).
!
!  For BDF orders q >= 3, the step size can be limited by stability rather
!  than accuracy, which shows up as growing oscillations in the higher
!  derivatives. STALD keeps the last five values of three scaled derivative
!  norms (`ssdat`), estimates from them the dominant characteristic root
!  `rr` of the method, and reduces the order when `rr` is close to 1
!  (> 0.98). It is turned on with [[CVodeSetStabLimDet]].

    submodule (fvode_core) fvode_core_stald

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvBDFStab`: update the scaled derivative data `ssdat` and, if the order
!  is not being reduced and enough steps have been taken at the current
!  order (`nscon >= q + 5`), call [[cvSLdet]]; if it detects a stability
!  limit, set `qprime = q - 1` and the step size for order q - 1.

    module subroutine cvBDFStab(me)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver

    integer :: i, k  !! indices into ssdat (i: time level, k: derivative)
    integer :: ldflag  !! return flag of cvSLdet
    integer :: factorial  !! (q-1)!
    real(wp) :: sq, sqm1, sqm2  !! scaled norms of the derivatives of order q+1, q and q-1

    ! if order is 3 or greater, then save scaled derivative data,
    ! push old data down in i, then add current values to top.
    if (me%q >= 3) then
        do k = 1, 3
            do i = 5, 2, -1
                me%ssdat(i,k) = me%ssdat(i-1,k)
            end do
        end do
        factorial = 1
        do i = 1, me%q - 1
            factorial = factorial * i
        end do
        sq   = factorial * me%q * (me%q + 1) * me%acnrm / max(me%tq(5), CV_TINY)
        sqm1 = factorial * me%q * n_vwrmsnorm(me%zn(:,me%q), me%ewt)
        sqm2 = factorial * n_vwrmsnorm(me%zn(:,me%q-1), me%ewt)
        me%ssdat(1,1) = sqm2 * sqm2
        me%ssdat(1,2) = sqm1 * sqm1
        me%ssdat(1,3) = sq * sq
    end if

    if (me%qprime >= me%q) then

        ! if order is 3 or greater, and enough ssdat has been saved,
        ! nscon >= q+5, then call stability limit detection routine.
        if ((me%q >= 3) .and. (me%nscon >= me%q + 5)) then
            ldflag = cvSLdet(me)
            if (ldflag > 3) then
                ! a stability limit violation is indicated by
                ! a return flag of 4, 5, or 6.
                ! Reduce new order.
                me%qprime = me%q - 1
                me%eta    = me%etaqm1
                me%eta    = min(me%eta, me%etamax)
                me%eta    = me%eta / max(ONE, abs(me%h) * me%hmax_inv * me%eta)
                me%hprime = me%h * me%eta
                me%nor    = me%nor + 1
            end if
        end if

    else
        ! otherwise, let order increase happen, and
        ! reset stability limit counter, nscon.
        me%nscon = 0
    end if

    end subroutine cvBDFStab
!*****************************************************************************************

!*****************************************************************************************
!>
!  `cvSLdet`: detect a stability limit from the data in `ssdat`.
!
!  For each of the three derivative sequences (k = 1, 2, 3), the ratios of
!  successive values estimate the characteristic root `rr`: either directly
!  (if the ratios are nearly constant), or by solving the quartic equations
!  for rr given by the data, with up to three Newton corrections. The
!  estimate is checked for consistency with the order q.
!
!  Returns:
!
!  * 1, 2, 3: `rr` found (by ratios, quartic, or corrected quartic), no
!    stability limit;
!  * 4, 5, 6: the same, with `rr > 0.98`: stability limit detected;
!  * < 0: no reliable estimate (-1: data too small, -2/-3: ratios not
!    constant enough, -4: singular system, -5: rr out of range, -6: no
!    convergence, -7/-8/-9: inconsistent data).

    function cvSLdet(me) result(kflag)

    class(cvode_solver),intent(inout) :: me  !! the CVODE solver
    integer :: kflag  !! detection result (see above)

    integer :: i, k, j  !! indices: time level i, derivative k and j
    integer :: it  !! Newton correction iteration
    integer :: kmin  !! derivative with the smallest residual
    real(wp),dimension(5,4) :: rat  !! ratios of successive ssdat values, rat(i,k) (C: [5][4], used as [1..4][1..3])
    real(wp),dimension(4) :: rav  !! average ratio for each k
    real(wp),dimension(4) :: qkr  !! residuals of the quartics at rr
    real(wp),dimension(4) :: sigsq  !! corrected values of ssdat(3,k)
    real(wp),dimension(4) :: smax  !! max of ssdat(:,k)
    real(wp),dimension(4) :: ssmax  !! smax(k)**2
    real(wp),dimension(4) :: drr  !! Newton corrections of rr
    real(wp),dimension(4) :: rrc  !! corrected rr for each k
    real(wp),dimension(4) :: sqmx  !! max scaled residual at rrc(k)
    real(wp),dimension(4,4) :: qjk  !! residual of quartic j at rrc(k)
    real(wp),dimension(5) :: vrat  !! variance of the ratios for each k
    real(wp),dimension(6,4) :: qc  !! quartic coefficients qc(i,k) (C: [6][4], used as [1..5][1..3])
    real(wp),dimension(6,4) :: qco  !! working copy of qc (for the elimination)
    real(wp) :: rr  !! estimated characteristic root
    real(wp) :: rrcut  !! rr above which a stability limit is reported
    real(wp) :: vrrtol, vrrt2  !! tolerances on the variance of the ratios
    real(wp) :: sqtol  !! tolerance on the scaled residual of the quartics
    real(wp) :: rrtol  !! tolerance on the consistency check of rr
    real(wp) :: smink, smaxk  !! min and max of ssdat(:,k)
    real(wp) :: sumrat, sumrsq  !! sum and sum of squares of the ratios
    real(wp) :: vmin, vmax  !! min and max of vrat
    real(wp) :: drrmax, adrr  !! max deviation of rav(k) from rr, and |rav(k) - rr|
    real(wp) :: tem  !! temporary (elimination factor; consistency test value)
    real(wp) :: sqmax  !! max scaled residual of the quartics
    real(wp) :: saqk, saqj  !! scaled residuals |qkr(k)|/ssmax(k), |qjk(j,k)|/ssmax(j)
    real(wp) :: qp  !! derivative of the quartic for the Newton correction
    real(wp) :: s  !! trial rr in the Newton correction
    real(wp) :: sqmaxk  !! max scaled residual at rrc(k)
    real(wp) :: sqmin  !! smallest sqmx(k)
    real(wp) :: rsa, rsb, rsc, rsd  !! ssdat(1..4,k) times powers of rr
    real(wp) :: rd1a, rd1b, rd1c  !! first differences of rsa..rsd
    real(wp) :: rd2a, rd2b  !! second differences
    real(wp) :: rd3a  !! third difference
    real(wp) :: cest1  !! estimate of the error constant
    real(wp) :: corr1  !! correction of ssdat(3,k)
    real(wp) :: ratp, ratm  !! sigsq(3)/sigsq(2) and sigsq(1)/sigsq(2)
    real(wp) :: qfac1, qfac2  !! order-dependent factors of the consistency check
    real(wp) :: bb  !! intermediate value of the consistency check
    real(wp) :: rrb  !! rr predicted by the consistency check

    kmin  = 0
    kflag = 0

    ! the following are cutoffs and tolerances used by this routine
    rrcut  = 0.98_wp
    vrrtol = 1.0e-4_wp
    vrrt2  = 5.0e-4_wp
    sqtol  = 1.0e-3_wp
    rrtol  = 1.0e-2_wp

    rr = ZERO

    ! index k corresponds to the degree of the interpolating polynomial:
    !   k = 1 -> q-1
    !   k = 2 -> q
    !   k = 3 -> q+1
    ! index i is a backward-in-time index, i = 1 -> current time,
    ! i = 2 -> previous step, etc.

    ! get maxima, minima, and variances, and form quartic coefficients
    do k = 1, 3
        smink = me%ssdat(1,k)
        smaxk = ZERO

        do i = 1, 5
            smink = min(smink, me%ssdat(i,k))
            smaxk = max(smaxk, me%ssdat(i,k))
        end do

        if (smink < CV_TINY * smaxk) then
            kflag = -1
            return
        end if
        smax(k)  = smaxk
        ssmax(k) = smaxk * smaxk

        sumrat = ZERO
        sumrsq = ZERO
        do i = 1, 4
            rat(i,k) = me%ssdat(i,k) / me%ssdat(i+1,k)
            sumrat   = sumrat + rat(i,k)
            sumrsq   = sumrsq + rat(i,k) * rat(i,k)
        end do
        rav(k)  = FOURTH * sumrat
        vrat(k) = abs(FOURTH * sumrsq - rav(k) * rav(k))

        qc(5,k) = me%ssdat(1,k) * me%ssdat(3,k) - me%ssdat(2,k) * me%ssdat(2,k)
        qc(4,k) = me%ssdat(2,k) * me%ssdat(3,k) - me%ssdat(1,k) * me%ssdat(4,k)
        qc(3,k) = ZERO
        qc(2,k) = me%ssdat(2,k) * me%ssdat(5,k) - me%ssdat(3,k) * me%ssdat(4,k)
        qc(1,k) = me%ssdat(4,k) * me%ssdat(4,k) - me%ssdat(3,k) * me%ssdat(5,k)

        do i = 1, 5
            qco(i,k) = qc(i,k)
        end do
    end do ! end of k loop

    ! isolate normal or nearly-normal matrix case. The three quartics will
    ! have a common or nearly-common root in this case.
    ! Return a kflag = 1 if this procedure works. If the three roots
    ! differ more than vrrt2, return error kflag = -3.
    vmin = min(vrat(1), min(vrat(2), vrat(3)))
    vmax = max(vrat(1), max(vrat(2), vrat(3)))

    if (vmin < vrrtol * vrrtol) then

        if (vmax > vrrt2 * vrrt2) then
            kflag = -2
            return
        else
            rr = (rav(1) + rav(2) + rav(3)) / THREE
            drrmax = ZERO
            do k = 1, 3
                adrr   = abs(rav(k) - rr)
                drrmax = max(drrmax, adrr)
            end do
            if (drrmax > vrrt2) then
                kflag = -3
                return
            end if
            kflag = 1
            ! can compute characteristic root, drop to next section
        end if

    else

        ! use the quartics to get rr
        if (abs(qco(1,1)) < CV_TINY * ssmax(1)) then
            kflag = -4
            return
        end if

        tem = qco(1,2) / qco(1,1)
        do i = 2, 5
            qco(i,2) = qco(i,2) - tem * qco(i,1)
        end do

        qco(1,2) = ZERO
        tem = qco(1,3) / qco(1,1)
        do i = 2, 5
            qco(i,3) = qco(i,3) - tem * qco(i,1)
        end do
        qco(1,3) = ZERO

        if (abs(qco(2,2)) < CV_TINY * ssmax(2)) then
            kflag = -4
            return
        end if

        tem = qco(2,3) / qco(2,2)
        do i = 3, 5
            qco(i,3) = qco(i,3) - tem * qco(i,2)
        end do

        if (abs(qco(4,3)) < CV_TINY * ssmax(3)) then
            kflag = -4
            return
        end if

        rr = -qco(5,3) / qco(4,3)

        if (rr < CV_TINY .or. rr > HUNDRED) then
            kflag = -5
            return
        end if

        do k = 1, 3
            qkr(k) = qc(5,k) + rr * (qc(4,k) + rr * rr * (qc(2,k) + rr * qc(1,k)))
        end do

        sqmax = ZERO
        do k = 1, 3
            saqk = abs(qkr(k)) / ssmax(k)
            if (saqk > sqmax) sqmax = saqk
        end do

        if (sqmax < sqtol) then
            kflag = 2
            ! can compute characteristic root, drop to "given rr, etc"
        else
            ! do Newton corrections to improve rr
            do it = 1, 3
                do k = 1, 3
                    qp = qc(4,k) + rr * rr * (THREE * qc(2,k) + rr * FOUR * qc(1,k))
                    drr(k) = ZERO
                    if (abs(qp) > CV_TINY * ssmax(k)) drr(k) = -qkr(k) / qp
                    rrc(k) = rr + drr(k)
                end do

                do k = 1, 3
                    s = rrc(k)
                    sqmaxk = ZERO
                    do j = 1, 3
                        qjk(j,k) = qc(5,j) + s * (qc(4,j) + s * s * (qc(2,j) + s * qc(1,j)))
                        saqj = abs(qjk(j,k)) / ssmax(j)
                        if (saqj > sqmaxk) sqmaxk = saqj
                    end do
                    sqmx(k) = sqmaxk
                end do

                sqmin = sqmx(1) + ONE
                do k = 1, 3
                    if (sqmx(k) < sqmin) then
                        kmin  = k
                        sqmin = sqmx(k)
                    end if
                end do
                rr = rrc(kmin)

                if (sqmin < sqtol) then
                    kflag = 3
                    ! can compute characteristic root
                    ! break out of Newton correction loop and drop to "given rr,etc"
                    exit
                else
                    do j = 1, 3
                        qkr(j) = qjk(j,kmin)
                    end do
                end if
            end do ! end of Newton correction loop

            if (sqmin > sqtol) then
                kflag = -6
                return
            end if
        end if ! end of if (sqmax < sqtol) else

    end if ! end of if (vmin < vrrtol*vrrtol) else, quartics to get rr

    ! given rr, find sigsq[k] and verify rr
    ! All positive kflag drop to this section
    do k = 1, 3
        rsa  = me%ssdat(1,k)
        rsb  = me%ssdat(2,k) * rr
        rsc  = me%ssdat(3,k) * rr * rr
        rsd  = me%ssdat(4,k) * rr * rr * rr
        rd1a = rsa - rsb
        rd1b = rsb - rsc
        rd1c = rsc - rsd
        rd2a = rd1a - rd1b
        rd2b = rd1b - rd1c
        rd3a = rd2a - rd2b

        if (abs(rd1b) < CV_TINY * smax(k)) then
            kflag = -7
            return
        end if

        cest1 = -rd3a / rd1b
        if (cest1 < CV_TINY .or. cest1 > FOUR) then
            kflag = -7
            return
        end if
        corr1    = (rd2b / cest1) / (rr * rr)
        sigsq(k) = me%ssdat(3,k) + corr1
    end do

    if (sigsq(2) < CV_TINY) then
        kflag = -8
        return
    end if

    ratp  = sigsq(3) / sigsq(2)
    ratm  = sigsq(1) / sigsq(2)
    qfac1 = FOURTH * (me%q * me%q - ONE)
    qfac2 = TWO / (me%q - ONE)
    bb    = ratp * ratm - ONE - qfac1 * ratp
    tem   = ONE - qfac2 * bb

    if (abs(tem) < CV_TINY) then
        kflag = -8
        return
    end if

    rrb = ONE / tem

    if (abs(rrb - rr) > rrtol) then
        kflag = -9
        return
    end if

    ! check to see if rr is above cutoff rrcut
    if (rr > rrcut) then
        if (kflag == 1) kflag = 4
        if (kflag == 2) kflag = 5
        if (kflag == 3) kflag = 6
    end if

    ! all positive kflag returned at this point

    end function cvSLdet
!*****************************************************************************************

    end submodule fvode_core_stald
!*****************************************************************************************

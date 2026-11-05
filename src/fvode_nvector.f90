!*****************************************************************************************
!>
!  Serial vector operations (SUNDIALS v7.9.0 `N_Vector` replacement).
!
!  CVODE only uses serial vectors, so the generic N_Vector type becomes a plain
!  `real(wp)` array. Each routine is a direct translation of the serial
!  implementation in `src/nvector/serial/nvector_serial.c`, with the same
!  special cases and loop order, so results are bitwise identical to SUNDIALS.
!
!  **Fused and vector-array operations** (`N_VLinearCombination`,
!  `N_VScaleAddMulti`, `N_VDotProdMulti`, `N_VScaleVectorArray`) are
!  *disabled by default* for serial vectors. CVODE then uses the generic
!  fallbacks in `src/sundials/sundials_nvector.c`, which are built from
!  `N_VScale`, `N_VLinearSum` and `N_VDotProd`. Those fallbacks are what is
!  translated here, because they are what a default CVODE run executes.
!
!  **Aliasing.** The C code often passes the same vector as input and output
!  (e.g. `N_VLinearSum(a, x, b, y, y)`). All vector dummy arguments here have the
!  `target` attribute, so this is legal Fortran provided the actual arguments are
!  contiguous targets (e.g. allocatable components of a `target` object). Every
!  loop reads element `i` of the inputs before writing element `i` of the output.
!
!  Arrays of vectors (`N_Vector*`) are arrays of [[nvec_ptr]].

    module fvode_nvector

    use fvode_kinds, only: wp, ip
    use fvode_parameters, only: SUN_BIG_REAL
    use fvode_math, only: sun_min, sun_rsqrt, sun_sqr

    implicit none

    private

    real(wp),parameter :: ZERO   = 0.0_wp
    real(wp),parameter :: HALF   = 0.5_wp
    real(wp),parameter :: ONE    = 1.0_wp
    real(wp),parameter :: ONEPT5 = 1.5_wp

    type,public :: nvec_ptr
        !! Pointer to a vector: the Fortran equivalent of an element of an
        !! `N_Vector*` array.
        real(wp),dimension(:),pointer,contiguous :: v => null()
    end type nvec_ptr

    ! standard vector operations
    public :: n_vlinearsum, n_vconst, n_vprod, n_vdiv, n_vscale, n_vabs, n_vinv
    public :: n_vaddconst, n_vdotprod, n_vmaxnorm, n_vwrmsnorm, n_vwsqrsumlocal
    public :: n_vmin, n_vcompare, n_vinvtest, n_vconstrmask, n_vminquotient
    ! fused / vector array operations (generic fallbacks)
    public :: n_vlinearcombination, n_vscaleaddmulti, n_vdotprodmulti
    public :: n_vscalevectorarray

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VLinearSum`: `z = a*x + b*y`.
!
!  The special cases are kept exactly as in `N_VLinearSum_Serial`, because some
!  of them round differently from `a*x + b*y` (e.g. `a == b` computes
!  `a*(x+y)`). The C code also checks for `z == y` / `z == x` (in-place axpy);
!  those paths give the same bits as the cases below, so no aliasing test is
!  needed.

    subroutine n_vlinearsum(a, x, b, y, z)

    real(wp),intent(in) :: a
    real(wp),dimension(:),intent(in),target :: x
    real(wp),intent(in) :: b
    real(wp),dimension(:),intent(in),target :: y
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i, n
    real(wp) :: c

    n = size(x, kind=ip)

    ! Case: a == b == 1.0
    if (a == ONE .and. b == ONE) then
        do i = 1, n
            z(i) = x(i) + y(i)                          ! VSum
        end do
        return
    end if

    ! Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0
    if (a == ONE .and. b == -ONE) then
        do i = 1, n
            z(i) = x(i) - y(i)                          ! VDiff(x, y)
        end do
        return
    end if
    if (a == -ONE .and. b == ONE) then
        do i = 1, n
            z(i) = y(i) - x(i)                          ! VDiff(y, x)
        end do
        return
    end if

    ! Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0
    if (a == ONE) then
        c = b
        do i = 1, n
            z(i) = (c * y(i)) + x(i)                    ! VLin1(b, y, x)
        end do
        return
    end if
    if (b == ONE) then
        c = a
        do i = 1, n
            z(i) = (c * x(i)) + y(i)                    ! VLin1(a, x, y)
        end do
        return
    end if

    ! Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0
    if (a == -ONE) then
        c = b
        do i = 1, n
            z(i) = (c * y(i)) - x(i)                    ! VLin2(b, y, x)
        end do
        return
    end if
    if (b == -ONE) then
        c = a
        do i = 1, n
            z(i) = (c * x(i)) - y(i)                    ! VLin2(a, x, y)
        end do
        return
    end if

    ! Case: a == b
    if (a == b) then
        do i = 1, n
            z(i) = a * (x(i) + y(i))                    ! VScaleSum
        end do
        return
    end if

    ! Case: a == -b
    if (a == -b) then
        do i = 1, n
            z(i) = a * (x(i) - y(i))                    ! VScaleDiff
        end do
        return
    end if

    ! all other cases
    do i = 1, n
        z(i) = (a * x(i)) + (b * y(i))
    end do

    end subroutine n_vlinearsum
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VConst`: `z = c`

    subroutine n_vconst(c, z)

    real(wp),intent(in) :: c
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(z, kind=ip)
        z(i) = c
    end do

    end subroutine n_vconst
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VProd`: `z = x .* y`

    subroutine n_vprod(x, y, z)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(in),target :: y
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(x, kind=ip)
        z(i) = x(i) * y(i)
    end do

    end subroutine n_vprod
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VDiv`: `z = x ./ y`

    subroutine n_vdiv(x, y, z)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(in),target :: y
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(x, kind=ip)
        z(i) = x(i) / y(i)
    end do

    end subroutine n_vdiv
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VScale`: `z = c*x`.
!
!  `c == 1` copies and `c == -1` negates, as in `N_VScale_Serial` (the in-place
!  `z == x` path, `x *= c`, gives the same bits).

    subroutine n_vscale(c, x, z)

    real(wp),intent(in) :: c
    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i, n

    n = size(x, kind=ip)

    if (c == ONE) then
        do i = 1, n
            z(i) = x(i)                                 ! VCopy
        end do
    else if (c == -ONE) then
        do i = 1, n
            z(i) = -x(i)                                ! VNeg
        end do
    else
        do i = 1, n
            z(i) = c * x(i)
        end do
    end if

    end subroutine n_vscale
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VAbs`: `z = |x|`

    subroutine n_vabs(x, z)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(x, kind=ip)
        z(i) = abs(x(i))
    end do

    end subroutine n_vabs
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VInv`: `z = 1 ./ x`

    subroutine n_vinv(x, z)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(x, kind=ip)
        z(i) = ONE / x(i)
    end do

    end subroutine n_vinv
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VAddConst`: `z = x + b`

    subroutine n_vaddconst(x, b, z)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),intent(in) :: b
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(x, kind=ip)
        z(i) = x(i) + b
    end do

    end subroutine n_vaddconst
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VDotProd`: `sum(x .* y)`, accumulated in index order.

    function n_vdotprod(x, y) result(sum)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(in),target :: y
    real(wp) :: sum

    integer(ip) :: i

    sum = ZERO
    do i = 1, size(x, kind=ip)
        sum = sum + x(i) * y(i)
    end do

    end function n_vdotprod
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VMaxNorm`: `max |x(i)|` (0 for an empty vector).

    function n_vmaxnorm(x) result(max)

    real(wp),dimension(:),intent(in),target :: x
    real(wp) :: max

    integer(ip) :: i

    max = ZERO
    do i = 1, size(x, kind=ip)
        if (abs(x(i)) > max) max = abs(x(i))
    end do

    end function n_vmaxnorm
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VWSqrSumLocal`: `sum((x .* w)**2)`, accumulated in index order.

    function n_vwsqrsumlocal(x, w) result(sum)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(in),target :: w
    real(wp) :: sum

    integer(ip) :: i
    real(wp) :: prodi

    sum = ZERO
    do i = 1, size(x, kind=ip)
        prodi = x(i) * w(i)
        sum = sum + sun_sqr(prodi)
    end do

    end function n_vwsqrsumlocal
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VWrmsNorm`: `sqrt(sum((x .* w)**2) / n)`

    function n_vwrmsnorm(x, w) result(norm)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(in),target :: w
    real(wp) :: norm

    norm = n_vwsqrsumlocal(x, w)
    norm = sun_rsqrt(norm / real(size(x, kind=ip), wp))

    end function n_vwrmsnorm
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VMin`: `min x(i)` (first element wins on ties).

    function n_vmin(x) result(min)

    real(wp),dimension(:),intent(in),target :: x
    real(wp) :: min

    integer(ip) :: i

    min = x(1)
    do i = 2, size(x, kind=ip)
        if (x(i) < min) min = x(i)
    end do

    end function n_vmin
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VCompare`: `z(i) = 1` if `|x(i)| >= c`, else `0`.

    subroutine n_vcompare(c, x, z)

    real(wp),intent(in) :: c
    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(inout),target :: z

    integer(ip) :: i

    do i = 1, size(x, kind=ip)
        if (abs(x(i)) >= c) then
            z(i) = ONE
        else
            z(i) = ZERO
        end if
    end do

    end subroutine n_vcompare
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VInvTest`: `z(i) = 1/x(i)` where `x(i) /= 0` (other elements of `z` are
!  left unchanged). Returns `.false.` if any `x(i) == 0`.

    function n_vinvtest(x, z) result(no_zero_found)

    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(inout),target :: z
    logical :: no_zero_found

    integer(ip) :: i

    no_zero_found = .true.
    do i = 1, size(x, kind=ip)
        if (x(i) == ZERO) then
            no_zero_found = .false.
        else
            z(i) = ONE / x(i)
        end if
    end do

    end function n_vinvtest
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VConstrMask`: checks `x` against the constraints `c`
!  (`c(i)` = 0: none, 1: `x >= 0`, -1: `x <= 0`, 2: `x > 0`, -2: `x < 0`).
!  Sets `m(i) = 1` where a constraint is violated, else `0`. Returns `.true.`
!  if all constraints are satisfied.

    function n_vconstrmask(c, x, m) result(ok)

    real(wp),dimension(:),intent(in),target :: c
    real(wp),dimension(:),intent(in),target :: x
    real(wp),dimension(:),intent(inout),target :: m
    logical :: ok

    integer(ip) :: i
    real(wp) :: temp
    logical :: test

    temp = ZERO
    do i = 1, size(x, kind=ip)
        m(i) = ZERO
        ! continue if no constraints were set for the variable
        if (c(i) == ZERO) cycle
        ! check if a set constraint has been violated
        test = (abs(c(i)) > ONEPT5 .and. x(i) * c(i) <= ZERO) .or. &
               (abs(c(i)) > HALF   .and. x(i) * c(i) <  ZERO)
        if (test) then
            m(i) = ONE
            temp = ONE
        end if
    end do

    ! return false if any constraint was violated
    ok = .not. (temp == ONE)

    end function n_vconstrmask
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VMinQuotient`: `min num(i)/denom(i)` over `denom(i) /= 0`.
!  Returns `SUN_BIG_REAL` if every `denom(i)` is zero.

    function n_vminquotient(num, denom) result(min)

    real(wp),dimension(:),intent(in),target :: num
    real(wp),dimension(:),intent(in),target :: denom
    real(wp) :: min

    integer(ip) :: i
    logical :: not_even_once

    not_even_once = .true.
    min = SUN_BIG_REAL

    do i = 1, size(num, kind=ip)
        if (denom(i) == ZERO) cycle
        if (.not. not_even_once) then
            min = sun_min(min, num(i) / denom(i))
        else
            min = num(i) / denom(i)
            not_even_once = .false.
        end if
    end do

    end function n_vminquotient
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VLinearCombination`: `z = sum_{i=1}^{nvec} c(i) * X(i)`.
!
!  Generic fallback from `sundials_nvector.c` (fused ops are disabled by default):
!  `N_VScale(c[0], X[0], z)`, then `N_VLinearSum(c[i], X[i], 1, z, z)` for `i >= 1`.
!  `z` may be one of the `X(i)`.

    subroutine n_vlinearcombination(nvec, c, X, z)

    integer,intent(in) :: nvec
    real(wp),dimension(:),intent(in) :: c
    type(nvec_ptr),dimension(:),intent(in) :: X
    real(wp),dimension(:),intent(inout),target :: z

    integer :: i

    call n_vscale(c(1), X(1)%v, z)
    do i = 2, nvec
        call n_vlinearsum(c(i), X(i)%v, ONE, z, z)
    end do

    end subroutine n_vlinearcombination
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VScaleAddMulti`: `Z(i) = a(i)*x + Y(i)`, `i = 1..nvec`.
!
!  Generic fallback: `N_VLinearSum(a[i], x, 1, Y[i], Z[i])`.

    subroutine n_vscaleaddmulti(nvec, a, x, Y, Z)

    integer,intent(in) :: nvec
    real(wp),dimension(:),intent(in) :: a
    real(wp),dimension(:),intent(in),target :: x
    type(nvec_ptr),dimension(:),intent(in) :: Y
    type(nvec_ptr),dimension(:),intent(in) :: Z

    integer :: i

    do i = 1, nvec
        call n_vlinearsum(a(i), x, ONE, Y(i)%v, Z(i)%v)
    end do

    end subroutine n_vscaleaddmulti
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VDotProdMulti`: `dotprods(i) = x . Y(i)`, `i = 1..nvec`.
!
!  Generic fallback: one `N_VDotProd` per vector.

    subroutine n_vdotprodmulti(nvec, x, Y, dotprods)

    integer,intent(in) :: nvec
    real(wp),dimension(:),intent(in),target :: x
    type(nvec_ptr),dimension(:),intent(in) :: Y
    real(wp),dimension(:),intent(inout) :: dotprods

    integer :: i

    do i = 1, nvec
        dotprods(i) = n_vdotprod(x, Y(i)%v)
    end do

    end subroutine n_vdotprodmulti
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VScaleVectorArray`: `Z(i) = c(i)*X(i)`, `i = 1..nvec`.
!
!  Generic fallback: one `N_VScale` per vector. `Z(i)` may be `X(i)`.

    subroutine n_vscalevectorarray(nvec, c, X, Z)

    integer,intent(in) :: nvec
    real(wp),dimension(:),intent(in) :: c
    type(nvec_ptr),dimension(:),intent(in) :: X
    type(nvec_ptr),dimension(:),intent(in) :: Z

    integer :: i

    do i = 1, nvec
        call n_vscale(c(i), X(i)%v, Z(i)%v)
    end do

    end subroutine n_vscalevectorarray
!*****************************************************************************************

    end module fvode_nvector
!*****************************************************************************************

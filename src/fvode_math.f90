!*****************************************************************************************
!>
!  Math helpers from `sundials_math.h` / `sundials_math.c` (SUNDIALS v7.9.0).
!
!  These must be used instead of the "obvious" Fortran intrinsics wherever the C
!  code uses the SUNDIALS macro, because the results can differ:
!
!  * [[sun_min]] / [[sun_max]] are the C ternaries `A < B ? A : B` and
!    `A > B ? A : B`. Fortran `min`/`max` treat NaN in a processor-dependent way,
!    and the ternaries return the *second* argument on ties (relevant for
!    signed zeros).
!  * [[sun_rpower_i]] multiplies in a loop like `SUNRpowerI`. Fortran `x**i`
!    uses binary exponentiation, which rounds differently.
!  * [[sun_rsqrt]] returns 0 for `x <= 0` like `SUNRsqrt`.

    module fvode_math

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUN_UNIT_ROUNDOFF, SUN_BIG_REAL

    implicit none

    private

    public :: sun_min, sun_max, sun_sqr, sun_rsqrt, sun_rabs
    public :: sun_ipower_i, sun_rpower_i, sun_rpower_r
    public :: sun_rsamesign, sun_rdifferentsign
    public :: sun_rcompare, sun_rcompare_tol

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNMIN(A, B)`: `((A) < (B) ? (A) : (B))`

    elemental function sun_min(a, b) result(c)

    real(wp),intent(in) :: a
    real(wp),intent(in) :: b
    real(wp) :: c

    if (a < b) then
        c = a
    else
        c = b
    end if

    end function sun_min
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNMAX(A, B)`: `((A) > (B) ? (A) : (B))`

    elemental function sun_max(a, b) result(c)

    real(wp),intent(in) :: a
    real(wp),intent(in) :: b
    real(wp) :: c

    if (a > b) then
        c = a
    else
        c = b
    end if

    end function sun_max
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNSQR(A)`: `((A) * (A))`

    elemental function sun_sqr(a) result(c)

    real(wp),intent(in) :: a
    real(wp) :: c

    c = a * a

    end function sun_sqr
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRsqrt(x)`: `((x) <= 0 ? 0 : sqrt(x))`. NaN input gives NaN.

    elemental function sun_rsqrt(x) result(c)

    real(wp),intent(in) :: x
    real(wp) :: c

    if (x <= 0.0_wp) then
        c = 0.0_wp
    else
        c = sqrt(x)
    end if

    end function sun_rsqrt
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRabs(x)`: `fabs(x)`

    elemental function sun_rabs(x) result(c)

    real(wp),intent(in) :: x
    real(wp) :: c

    c = abs(x)

    end function sun_rabs
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNIpowerI(base, exponent)`: integer power by repeated multiplication
!  (returns 1 for `exponent <= 0`).

    elemental function sun_ipower_i(base, exponent) result(prod)

    integer,intent(in) :: base
    integer,intent(in) :: exponent
    integer :: prod

    integer :: i

    prod = 1
    do i = 1, exponent
        prod = prod * base
    end do

    end function sun_ipower_i
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRpowerI(base, exponent)`: real to an integer power by repeated
!  multiplication. A negative exponent gives the reciprocal of the positive power.

    elemental function sun_rpower_i(base, exponent) result(prod)

    real(wp),intent(in) :: base
    integer,intent(in) :: exponent
    real(wp) :: prod

    integer :: i, expt

    prod = 1.0_wp
    expt = abs(exponent)
    do i = 1, expt
        prod = prod * base
    end do
    if (exponent < 0) prod = 1.0_wp / prod

    end function sun_rpower_i
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRpowerR(base, exponent)`: `pow(base, exponent)`.
!
!  gfortran evaluates a real `**` real power with the C library `pow`, the same
!  function SUNDIALS calls, so this matches bitwise (checked by the
!  `math_powers` test).

    elemental function sun_rpower_r(base, exponent) result(c)

    real(wp),intent(in) :: base
    real(wp),intent(in) :: exponent
    real(wp) :: c

    c = base ** exponent

    end function sun_rpower_r
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRsamesign(x, y)`: `signbit(x) == signbit(y)`.
!
!  `sign(1, x)` copies the sign bit of `x` (including for -0.0 and NaN).

    elemental function sun_rsamesign(x, y) result(same)

    real(wp),intent(in) :: x
    real(wp),intent(in) :: y
    logical :: same

    same = (sign(1.0_wp, x) == sign(1.0_wp, y))

    end function sun_rsamesign
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRdifferentsign(x, y)`: `!SUNRsamesign(x, y)`

    elemental function sun_rdifferentsign(x, y) result(different)

    real(wp),intent(in) :: x
    real(wp),intent(in) :: y
    logical :: different

    different = .not. sun_rsamesign(x, y)

    end function sun_rdifferentsign
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRCompare(a, b)`: `SUNRCompareTol(a, b, 10*SUN_UNIT_ROUNDOFF)`.
!
!  Note the SUNDIALS convention: returns `.false.` if `a` and `b` are
!  (approximately) **equal**, `.true.` if they differ.

    elemental function sun_rcompare(a, b) result(differ)

    real(wp),intent(in) :: a
    real(wp),intent(in) :: b
    logical :: differ

    differ = sun_rcompare_tol(a, b, 10 * SUN_UNIT_ROUNDOFF)

    end function sun_rcompare
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNRCompareTol(a, b, tol)`: returns `.false.` if `a` and `b` are equal to
!  within `tol` (relative to `|a+b|`, or absolute `10*SUN_UNIT_ROUNDOFF` when
!  `|a+b|` is small), `.true.` otherwise. NaN inputs give `.true.`
!  (the C code uses `!isless`; a Fortran `<` with NaN is also false).

    elemental function sun_rcompare_tol(a, b, tol) result(differ)

    real(wp),intent(in) :: a
    real(wp),intent(in) :: b
    real(wp),intent(in) :: tol
    logical :: differ

    real(wp) :: diff, norm

    if (a == b) then
        differ = .false.
        return
    end if

    diff = abs(a - b)
    norm = sun_min(abs(a + b), SUN_BIG_REAL)

    differ = .not. (diff < sun_max(10 * SUN_UNIT_ROUNDOFF, tol * norm))

    end function sun_rcompare_tol
!*****************************************************************************************

    end module fvode_math
!*****************************************************************************************

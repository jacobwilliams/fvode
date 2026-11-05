!*****************************************************************************************
!>
!  Vector operations with nontrivial logic (SUNDIALS v7.9.0 `N_Vector`).
!
!  CVODE only uses serial vectors, so the generic N_Vector type becomes a plain
!  `real(wp)` array, and simple operations are written directly as Fortran
!  array expressions in the translated code:
!
!  | SUNDIALS                     | Fortran                 |
!  |------------------------------|-------------------------|
!  | `N_VLinearSum(a,x,b,y,z)`    | `z = a*x + b*y`         |
!  | `N_VScale(c,x,z)`            | `z = c*x`               |
!  | `N_VConst(c,z)`              | `z = c`                 |
!  | `N_VProd` / `N_VDiv`         | `z = x*y` / `z = x/y`   |
!  | `N_VAbs` / `N_VInv`          | `z = abs(x)` / `z = 1/x`|
!  | `N_VAddConst(x,b,z)`         | `z = x + b`             |
!  | `N_VDotProd(x,y)`            | `dot_product(x,y)`      |
!  | `N_VMaxNorm(x)` / `N_VMin(x)`| `maxval(abs(x))` / `minval(x)` |
!  | `N_VCompare(c,x,z)`          | `z = merge(1, 0, abs(x) >= c)` |
!  | `N_VLinearCombination` etc.  | sums over columns of `zn` |
!
!  This module holds the operations that are more than one expression.
!  Translated from `src/nvector/serial/nvector_serial.c`.

    module fvode_nvector

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUN_BIG_REAL

    implicit none

    private

    public :: n_vwrmsnorm, n_vinvtest, n_vconstrmask, n_vminquotient

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VWrmsNorm`: weighted root-mean-square norm `sqrt(sum((x*w)**2) / n)`.

    pure function n_vwrmsnorm(x, w) result(norm)

    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(in) :: w
    real(wp) :: norm

    norm = sqrt(sum((x*w)**2) / size(x))

    end function n_vwrmsnorm
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VInvTest`: `z(i) = 1/x(i)` where `x(i) /= 0` (other elements of `z` are
!  left unchanged). Returns `.false.` if any `x(i) == 0`.

    function n_vinvtest(x, z) result(no_zero_found)

    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(inout) :: z
    logical :: no_zero_found

    where (x /= 0.0_wp) z = 1.0_wp / x
    no_zero_found = all(x /= 0.0_wp)

    end function n_vinvtest
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VConstrMask`: checks `x` against the constraints `c`
!  (`c(i)` = 0: none, 1: `x >= 0`, -1: `x <= 0`, 2: `x > 0`, -2: `x < 0`).
!  Sets `m(i) = 1` where a constraint is violated, else `0`. Returns `.true.`
!  if all constraints are satisfied.

    function n_vconstrmask(c, x, m) result(ok)

    real(wp),dimension(:),intent(in) :: c
    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(out) :: m
    logical :: ok

    logical,dimension(size(x)) :: violated

    violated = (c /= 0.0_wp) .and. &
               ((abs(c) > 1.5_wp .and. x*c <= 0.0_wp) .or. &
                (abs(c) > 0.5_wp .and. x*c <  0.0_wp))
    m  = merge(1.0_wp, 0.0_wp, violated)
    ok = .not. any(violated)

    end function n_vconstrmask
!*****************************************************************************************

!*****************************************************************************************
!>
!  `N_VMinQuotient`: `min num(i)/denom(i)` over `denom(i) /= 0`.
!  Returns `SUN_BIG_REAL` if every `denom(i)` is zero.

    pure function n_vminquotient(num, denom) result(qmin)

    real(wp),dimension(:),intent(in) :: num
    real(wp),dimension(:),intent(in) :: denom
    real(wp) :: qmin

    if (any(denom /= 0.0_wp)) then
        qmin = minval(num / denom, mask = denom /= 0.0_wp)
    else
        qmin = SUN_BIG_REAL
    end if

    end function n_vminquotient
!*****************************************************************************************

    end module fvode_nvector
!*****************************************************************************************

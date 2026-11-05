!*****************************************************************************************
!>
!  Unit tests for [[fvode_nvector]].

    module test_nvector

    use fvode_kinds, only: wp
    use fvode_parameters, only: SUN_BIG_REAL
    use fvode_nvector
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

    private

    public :: collect_nvector

    real(wp),parameter :: tol = 1.0e-14_wp

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Collect all tests in this suite.

    subroutine collect_nvector(testsuite)

    type(unittest_type),allocatable,intent(out) :: testsuite(:)

    testsuite = [ new_unittest("wrmsnorm",    test_wrmsnorm),    &
                  new_unittest("invtest",     test_invtest),     &
                  new_unittest("constrmask",  test_constrmask),  &
                  new_unittest("minquotient", test_minquotient)  ]

    end subroutine collect_nvector
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_wrmsnorm(error)

    type(error_type),allocatable,intent(out) :: error

    ! sqrt(((1*2)**2 + (-3*1)**2 + (4*0.5)**2) / 3) = sqrt(17/3)
    call check(error, n_vwrmsnorm([1.0_wp, -3.0_wp, 4.0_wp], [2.0_wp, 1.0_wp, 0.5_wp]), &
               sqrt(17.0_wp/3.0_wp), thr=tol, rel=.true.)

    end subroutine test_wrmsnorm
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_invtest(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(3) :: z
    logical :: ok

    z = -7.0_wp
    ok = n_vinvtest([2.0_wp, 0.0_wp, -4.0_wp], z)
    call check(error, .not. ok, 'zero not detected');           if (allocated(error)) return
    call check(error, z(1), 0.5_wp);                            if (allocated(error)) return
    call check(error, z(2), -7.0_wp, 'zero element modified');   if (allocated(error)) return
    call check(error, z(3), -0.25_wp);                          if (allocated(error)) return

    ok = n_vinvtest([2.0_wp, 1.0_wp, -4.0_wp], z)
    call check(error, ok, 'false zero detected')

    end subroutine test_invtest
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_constrmask(error)

    type(error_type),allocatable,intent(out) :: error

    !                        none  >=0   >=0   <=0   >0    >0    <0    <0
    real(wp),dimension(8),parameter :: c = [0.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, &
                                            2.0_wp, 2.0_wp, -2.0_wp, -2.0_wp]
    real(wp),dimension(8),parameter :: x = [-5.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
                                            1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp]
    real(wp),dimension(8),parameter :: expected = [0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
                                                   0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp]
    real(wp),dimension(8) :: m
    logical :: ok

    ok = n_vconstrmask(c, x, m)
    call check(error, .not. ok, 'violation not detected');      if (allocated(error)) return
    call check(error, all(m == expected), 'wrong mask');         if (allocated(error)) return

    ok = n_vconstrmask(c, abs(x) + 1.0_wp, m)
    call check(error, .not. ok, 'violation (c<0) not detected'); if (allocated(error)) return

    ok = n_vconstrmask(c, sign(abs(x) + 1.0_wp, c), m)
    call check(error, ok, 'false violation');                    if (allocated(error)) return
    call check(error, all(m == 0.0_wp), 'mask not cleared')

    end subroutine test_constrmask
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_minquotient(error)

    type(error_type),allocatable,intent(out) :: error

    call check(error, n_vminquotient([1.0_wp, 5.0_wp, -3.0_wp], [2.0_wp, 0.0_wp, 1.0_wp]), &
               -3.0_wp)
    if (allocated(error)) return
    call check(error, n_vminquotient([1.0_wp, 5.0_wp], [0.0_wp, 0.0_wp]), SUN_BIG_REAL)

    end subroutine test_minquotient
!*****************************************************************************************

    end module test_nvector
!*****************************************************************************************

!*****************************************************************************************
!>
!  Bitwise checks of [[fvode_nvector]] and [[fvode_math]] against SUNDIALS.
!
!  `test/data/nvector_serial.txt` is written by `reference/nvector_ops.c`
!  (via `pixi run -e reference ref-run`). It holds the inputs and the results of
!  every serial N_Vector operation CVODE uses, computed by the real SUNDIALS
!  v7.9.0 library. Each test here repeats the operation in Fortran and requires
!  **bitwise** equality.

    module test_nvector

    use fvode_kinds, only: wp
    use fvode_nvector
    use fvode_math, only: sun_rpower_i, sun_rpower_r
    use testdrive, only: new_unittest, unittest_type, error_type, test_failed

    implicit none

    private

    public :: collect_nvector

    character(len=*),parameter :: data_file = 'test/data/nvector_serial.txt'

    type :: block_type
        character(len=:),allocatable :: name
        real(wp),dimension(:),allocatable :: v
    end type block_type

    type(block_type),dimension(:),allocatable,save :: blocks  !! data file contents
    logical,save :: loaded = .false.

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Collect all tests in this suite.

    subroutine collect_nvector(testsuite)

    type(unittest_type),allocatable,intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest("linearsum",        test_linearsum),        &
        new_unittest("scale",            test_scale),            &
        new_unittest("elementwise",      test_elementwise),      &
        new_unittest("invtest",          test_invtest),          &
        new_unittest("constrmask",       test_constrmask),       &
        new_unittest("reductions",       test_reductions),       &
        new_unittest("linearcombination",test_linearcombination),&
        new_unittest("scaleaddmulti",    test_scaleaddmulti),    &
        new_unittest("dotprodmulti",     test_dotprodmulti),     &
        new_unittest("scalevectorarray", test_scalevectorarray), &
        new_unittest("math_powers",      test_math_powers)       ]

    end subroutine collect_nvector
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the data file once (blocks of `name n` followed by `n` values).

    subroutine load(error)

    type(error_type),allocatable,intent(out) :: error

    integer :: iunit, istat, n, i
    character(len=256) :: line, name
    type(block_type) :: b

    if (loaded) return

    open(newunit=iunit, file=data_file, status='old', action='read', iostat=istat)
    if (istat /= 0) then
        call test_failed(error, 'cannot open '//data_file)
        return
    end if

    allocate(blocks(0))
    do
        read(iunit, '(a)', iostat=istat) line
        if (istat /= 0) exit
        read(line, *) name, n
        b%name = trim(name)
        allocate(b%v(n))
        do i = 1, n
            read(iunit, *) b%v(i)
        end do
        blocks = [blocks, b]
        deallocate(b%v)
    end do
    close(iunit)
    loaded = .true.

    end subroutine load
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get a block by name (fails the test if missing).

    subroutine get(error, name, v)

    type(error_type),allocatable,intent(inout) :: error
    character(len=*),intent(in) :: name
    real(wp),dimension(:),allocatable,intent(out) :: v

    integer :: i

    if (allocated(error)) return
    do i = 1, size(blocks)
        if (blocks(i)%name == name) then
            v = blocks(i)%v
            return
        end if
    end do
    call test_failed(error, 'block missing from data file: '//name)

    end subroutine get
!*****************************************************************************************

!*****************************************************************************************
!>
!  Require bitwise equality of `got` with the reference block `name`.

    subroutine expect(error, name, got)

    type(error_type),allocatable,intent(inout) :: error
    character(len=*),intent(in) :: name
    real(wp),dimension(:),intent(in) :: got

    real(wp),dimension(:),allocatable :: ref
    integer :: i
    character(len=200) :: msg

    if (allocated(error)) return
    call get(error, name, ref); if (allocated(error)) return

    if (size(ref) /= size(got)) then
        write(msg, '(a,": size ",i0," vs reference ",i0)') name, size(got), size(ref)
        call test_failed(error, trim(msg))
        return
    end if
    do i = 1, size(ref)
        ! bitwise: compare bit patterns so -0.0 /= 0.0 and NaN == NaN
        if (transfer(got(i), 0_8) /= transfer(ref(i), 0_8)) then
            write(msg, '(a,"(",i0,"): got ",es25.17," reference ",es25.17)') &
                name, i, got(i), ref(i)
            call test_failed(error, trim(msg))
            return
        end if
    end do

    end subroutine expect
!*****************************************************************************************

!*****************************************************************************************
!>
!  Scalar convenience wrapper for [[expect]].

    subroutine expect_s(error, name, got)

    type(error_type),allocatable,intent(inout) :: error
    character(len=*),intent(in) :: name
    real(wp),intent(in) :: got

    call expect(error, name, [got])

    end subroutine expect_s
!*****************************************************************************************

!*****************************************************************************************
!>
!  Logical result as stored in the data file (1/0).

    pure function flag(l) result(r)

    logical,intent(in) :: l
    real(wp) :: r

    r = merge(1.0_wp, 0.0_wp, l)

    end function flag
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_linearsum(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, y, z, xx, yy, coef
    integer :: k
    character(len=32) :: tag

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'in_y', y)
    call get(error, 'linsum_coef', coef)
    if (allocated(error)) return

    allocate(z(size(x)))
    do k = 0, size(coef)/2 - 1
        write(tag, '(i0)') k
        associate (a => coef(2*k+1), b => coef(2*k+2))
            call n_vlinearsum(a, x, b, y, z)
            call expect(error, 'linsum_'//trim(tag), z)
            yy = y
            call n_vlinearsum(a, x, b, yy, yy)          ! z == y
            call expect(error, 'linsum_zy_'//trim(tag), yy)
            xx = x
            call n_vlinearsum(a, xx, b, y, xx)          ! z == x
            call expect(error, 'linsum_zx_'//trim(tag), xx)
        end associate
        if (allocated(error)) return
    end do

    end subroutine test_linearsum
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_scale(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, z, xx, coef
    integer :: k
    character(len=32) :: tag

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'scale_coef', coef)
    if (allocated(error)) return

    allocate(z(size(x)))
    do k = 0, size(coef) - 1
        write(tag, '(i0)') k
        call n_vscale(coef(k+1), x, z)
        call expect(error, 'scale_'//trim(tag), z)
        xx = x
        call n_vscale(coef(k+1), xx, xx)
        call expect(error, 'scale_inplace_'//trim(tag), xx)
        if (allocated(error)) return
    end do

    end subroutine test_scale
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_elementwise(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, y, w, z

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'in_y', y)
    call get(error, 'in_w', w)
    if (allocated(error)) return

    allocate(z(size(x)))
    call n_vconst(0.123_wp, z);      call expect(error, 'const', z)
    call n_vprod(x, y, z);           call expect(error, 'prod', z)
    call n_vdiv(x, w, z);            call expect(error, 'div', z)
    call n_vabs(x, z);               call expect(error, 'abs', z)
    call n_vinv(w, z);               call expect(error, 'inv', z)
    call n_vaddconst(x, 0.25_wp, z); call expect(error, 'addconst', z)
    call n_vcompare(0.5_wp, x, z);   call expect(error, 'compare', z)

    end subroutine test_elementwise
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_invtest(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: xz, w, z
    logical :: ok

    call load(error); if (allocated(error)) return
    call get(error, 'in_xz', xz)
    call get(error, 'in_w', w)
    if (allocated(error)) return

    allocate(z(size(xz)))
    z = -7.0_wp
    ok = n_vinvtest(xz, z)
    call expect_s(error, 'invtest_flag', flag(ok))
    call expect(error, 'invtest', z)

    z = -7.0_wp
    ok = n_vinvtest(w, z)
    call expect_s(error, 'invtest_ok_flag', flag(ok))
    call expect(error, 'invtest_ok', z)

    end subroutine test_invtest
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_constrmask(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: c, x, xs, m
    logical :: ok

    call load(error); if (allocated(error)) return
    call get(error, 'in_c', c)
    call get(error, 'in_x', x)
    call get(error, 'in_xs', xs)
    if (allocated(error)) return

    allocate(m(size(x)))
    ok = n_vconstrmask(c, x, m)
    call expect_s(error, 'constrmask_flag', flag(ok))
    call expect(error, 'constrmask', m)

    ok = n_vconstrmask(c, xs, m)
    call expect_s(error, 'constrmask_ok_flag', flag(ok))
    call expect(error, 'constrmask_ok', m)

    end subroutine test_constrmask
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_reductions(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, y, w, d, zeros

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'in_y', y)
    call get(error, 'in_w', w)
    call get(error, 'in_d', d)
    if (allocated(error)) return

    call expect_s(error, 'dotprod',      n_vdotprod(x, y))
    call expect_s(error, 'maxnorm',      n_vmaxnorm(x))
    call expect_s(error, 'wrmsnorm',     n_vwrmsnorm(x, w))
    call expect_s(error, 'wsqrsumlocal', n_vwsqrsumlocal(x, w))
    call expect_s(error, 'min',          n_vmin(x))
    call expect_s(error, 'minquotient',  n_vminquotient(x, d))
    allocate(zeros(size(x)), source=0.0_wp)
    call expect_s(error, 'minquotient_allzero', n_vminquotient(x, zeros))

    end subroutine test_reductions
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_linearcombination(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, y, w, d, z, xx, coef
    type(nvec_ptr),dimension(4) :: X4

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'in_y', y)
    call get(error, 'in_w', w)
    call get(error, 'in_d', d)
    call get(error, 'lincomb_coef', coef)
    if (allocated(error)) return

    allocate(z(size(x)))
    X4(1)%v => x; X4(2)%v => y; X4(3)%v => w; X4(4)%v => d
    call n_vlinearcombination(4, coef, X4, z)
    call expect(error, 'lincomb', z)

    xx = x
    X4(1)%v => xx                                   ! z == X(1)
    call n_vlinearcombination(4, coef, X4, xx)
    call expect(error, 'lincomb_alias', xx)

    end subroutine test_linearcombination
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_scaleaddmulti(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, y, w, d, r1, r2, r3, coef
    type(nvec_ptr),dimension(3) :: Yv, Zv

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'in_y', y)
    call get(error, 'in_w', w)
    call get(error, 'in_d', d)
    call get(error, 'scaleaddmulti_coef', coef)
    if (allocated(error)) return

    allocate(r1(size(x)), r2(size(x)), r3(size(x)))
    Yv(1)%v => y;  Yv(2)%v => w;  Yv(3)%v => d
    Zv(1)%v => r1; Zv(2)%v => r2; Zv(3)%v => r3
    call n_vscaleaddmulti(3, coef, x, Yv, Zv)
    call expect(error, 'scaleaddmulti_1', r1)
    call expect(error, 'scaleaddmulti_2', r2)
    call expect(error, 'scaleaddmulti_3', r3)

    end subroutine test_scaleaddmulti
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_dotprodmulti(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: x, y, w, d
    type(nvec_ptr),dimension(3) :: Yv
    real(wp),dimension(3) :: dots

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', x)
    call get(error, 'in_y', y)
    call get(error, 'in_w', w)
    call get(error, 'in_d', d)
    if (allocated(error)) return

    Yv(1)%v => y; Yv(2)%v => w; Yv(3)%v => d
    call n_vdotprodmulti(3, x, Yv, dots)
    call expect(error, 'dotprodmulti', dots)

    end subroutine test_dotprodmulti
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_scalevectorarray(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable,target :: a1, a2, a3, coef
    type(nvec_ptr),dimension(3) :: A

    call load(error); if (allocated(error)) return
    call get(error, 'in_x', a1)
    call get(error, 'in_y', a2)
    call get(error, 'in_w', a3)
    call get(error, 'scalevecarray_coef', coef)
    if (allocated(error)) return

    A(1)%v => a1; A(2)%v => a2; A(3)%v => a3
    call n_vscalevectorarray(3, coef, A, A)         ! in place, as in cvRescale
    call expect(error, 'scalevecarray_1', a1)
    call expect(error, 'scalevecarray_2', a2)
    call expect(error, 'scalevecarray_3', a3)

    end subroutine test_scalevectorarray
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_math_powers(error)

    type(error_type),allocatable,intent(out) :: error

    real(wp),dimension(:),allocatable :: pb, pe
    real(wp),dimension(30) :: powi
    real(wp),dimension(20) :: powr
    real(wp),dimension(4),parameter :: rb = [0.3_wp, 1.7_wp, 1.0e-6_wp, 123.456_wp]
    real(wp),dimension(5),parameter :: re = [1.0_wp/3.0_wp, 0.2_wp, 1.0_wp/6.0_wp, &
                                             0.5_wp, 2.0_wp]
    integer :: i, j, k

    call load(error); if (allocated(error)) return
    call get(error, 'powi_base', pb)
    call get(error, 'powi_exp', pe)
    if (allocated(error)) return

    k = 0
    do i = 1, size(pb)
        do j = 1, size(pe)
            k = k + 1
            powi(k) = sun_rpower_i(pb(i), nint(pe(j)))
        end do
    end do
    call expect(error, 'powi', powi)

    k = 0
    do i = 1, size(rb)
        do j = 1, size(re)
            k = k + 1
            powr(k) = sun_rpower_r(rb(i), re(j))
        end do
    end do
    call expect(error, 'powr', powr)

    end subroutine test_math_powers
!*****************************************************************************************

    end module test_nvector
!*****************************************************************************************

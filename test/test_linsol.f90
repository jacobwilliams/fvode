!*****************************************************************************************
!>
!  Unit tests for the matrices and the direct linear solvers.

    module test_linsol

    use fvode_kinds, only: wp, ip
    use fvode_parameters, only: SUN_SUCCESS, SUNLS_LUFACT_FAIL
    use fvode_matrix, only: dense_matrix, band_matrix
    use fvode_linsol_direct, only: sun_linsol_dense, sun_linsol_band
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

    private

    public :: collect_linsol

    integer(ip),parameter :: n = 9, mu = 2, ml = 3  !! matrix size and upper/lower bandwidths of the test matrix

    contains
!*****************************************************************************************

    subroutine collect_linsol(testsuite)
    !! collect all tests in this suite
    type(unittest_type),allocatable,intent(out) :: testsuite(:)  !! the unit tests of this suite
    testsuite = [ new_unittest("dense_lu_solve", test_dense),   &
                  new_unittest("dense_singular", test_singular),&
                  new_unittest("band_lu_solve",  test_band),    &
                  new_unittest("band_layout",    test_band_layout) ]
    end subroutine collect_linsol

!*****************************************************************************************
!>
!  A test matrix with entries in the band -mu..ml: diagonally weak, so
!  partial pivoting is needed.

    pure real(wp) function aij(i, j)
    integer(ip),intent(in) :: i, j  !! row and column index (1-based)
    if (j - i > mu .or. i - j > ml) then
        aij = 0.0_wp
    else
        aij = real(mod(3*i + 7*j, 11) - 5, wp) / 3.0_wp
        if (i == j) aij = 0.1_wp * aij
    end if
    end function aij

!*****************************************************************************************
    subroutine test_dense(error)
    !! test the dense LU solve against a known solution
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(dense_matrix) :: A  !! matrix (overwritten by its LU factors)
    type(sun_linsol_dense) :: S  !! dense linear solver
    real(wp),dimension(n) :: x, b, xtrue  !! computed solution, right-hand side, and exact solution
    real(wp),dimension(n,n) :: Aorig  !! copy of the matrix before factorization
    integer(ip) :: i, j  !! row and column indices
    A = dense_matrix(n, n)
    S = sun_linsol_dense(n)
    do j = 1, n
        do i = 1, n
            A%data(i,j) = aij(i,j) + merge(0.5_wp, 0.0_wp, i == 1 .or. j == n)
        end do
    end do
    Aorig = A%data
    xtrue = [(real(i, wp) - 4.5_wp, i = 1, n)]
    b = matmul(Aorig, xtrue)
    call check(error, S%setup(A), SUN_SUCCESS);      if (allocated(error)) return
    call check(error, S%solve(A, x, b, 0.0_wp), SUN_SUCCESS); if (allocated(error)) return
    call check(error, maxval(abs(x - xtrue)) < 1.0e-12_wp, 'dense solution inaccurate')
    end subroutine test_dense

!*****************************************************************************************
    subroutine test_singular(error)
    !! test that a zero pivot is reported
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(dense_matrix) :: A  !! singular matrix
    type(sun_linsol_dense) :: S  !! dense linear solver
    A = dense_matrix(3_ip, 3_ip)
    S = sun_linsol_dense(3_ip)
    A%data = reshape([1.0_wp, 2.0_wp, 3.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 0.0_wp, 1.0_wp, 1.0_wp], [3,3])
    call check(error, S%setup(A), SUNLS_LUFACT_FAIL);  if (allocated(error)) return
    call check(error, S%last_flag, 2, 'zero pivot should be reported in column 2')
    end subroutine test_singular

!*****************************************************************************************
    subroutine test_band(error)
    !! test the band matvec and LU solve against a known solution
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(band_matrix) :: A  !! band matrix (overwritten by its LU factors)
    type(sun_linsol_band) :: S  !! band linear solver
    real(wp),dimension(n) :: x, b, xtrue  !! computed solution or product, right-hand side, and exact solution
    real(wp),dimension(n,n) :: Afull  !! the same matrix in dense form
    integer(ip) :: i, j  !! row and column indices
    A = band_matrix(n, mu, ml)
    S = sun_linsol_band(n)
    Afull = 0.0_wp
    do j = 1, n
        do i = max(1_ip, j - mu), min(n, j + ml)
            call A%set(i, j, aij(i,j))
            Afull(i,j) = aij(i,j)
        end do
    end do
    xtrue = [(real(i, wp)**2 / 10.0_wp - 1.0_wp, i = 1, n)]
    b = matmul(Afull, xtrue)
    ! band matvec agrees with the dense product
    call A%matvec(xtrue, x)
    call check(error, maxval(abs(x - b)) < 1.0e-13_wp, 'band matvec wrong'); if (allocated(error)) return
    call check(error, S%setup(A), SUN_SUCCESS);      if (allocated(error)) return
    call check(error, S%solve(A, x, b, 0.0_wp), SUN_SUCCESS); if (allocated(error)) return
    call check(error, maxval(abs(x - xtrue)) < 1.0e-12_wp, 'band solution inaccurate')
    end subroutine test_band

!*****************************************************************************************
    subroutine test_band_layout(error)
    !! scale_add_i and copy on the band agree with the dense equivalents
    type(error_type),allocatable,intent(out) :: error  !! error information (allocated if the test fails)
    type(band_matrix) :: A, B  !! band matrix and its copy
    type(dense_matrix) :: D  !! the same matrix in dense form
    integer(ip) :: i, j  !! row and column indices
    integer :: flag  !! return flag of the copy
    A = band_matrix(n, mu, ml)
    D = dense_matrix(n, n)
    do j = 1, n
        do i = max(1_ip, j - mu), min(n, j + ml)
            call A%set(i, j, aij(i,j))
            D%data(i,j) = aij(i,j)
        end do
    end do
    call A%scale_add_i(-0.25_wp)
    call D%scale_add_i(-0.25_wp)
    B = band_matrix(n, mu, ml)
    call A%copy_to(B, flag)
    call check(error, flag, 0);  if (allocated(error)) return
    do j = 1, n
        do i = max(1_ip, j - mu), min(n, j + ml)
            if (B%get(i,j) /= D%data(i,j)) then
                call check(error, .false., 'band scale_add_i/copy mismatch')
                return
            end if
        end do
    end do
    end subroutine test_band_layout

    end module test_linsol
!*****************************************************************************************

!*****************************************************************************************
!>
!  Matrices used by the direct linear solvers (SUNDIALS v7.9.0 `SUNMatrix`):
!  dense and band.
!
!  Translated from `src/sunmatrix/dense/sunmatrix_dense.c` and
!  `src/sunmatrix/band/sunmatrix_band.c`. Only the operations CVODE uses are
!  provided: zero, copy, `A = c*A + I`, and (for Krylov solvers, later) matvec.
!
!### Storage
!
!  * [[dense_matrix]]: `data(1:m, 1:n)`, column-major (as in SUNDIALS).
!  * [[band_matrix]]: SUNDIALS band storage, `data(0:ldim-1, 0:n-1)` with
!    `ldim = smu + ml + 1`. Element (i,j) (0-based, `-mu <= i-j <= ml`) is
!    `data(i - j + smu, j)`. The extra `smu - mu` rows above the band hold the
!    fill-in of the LU factorization. Use [[band_get]]/[[band_set]] with
!    1-based (i,j) to access elements.

    module fvode_matrix

    use fvode_kinds, only: wp, ip

    implicit none

    private

    integer,parameter,public :: SUNMATRIX_DENSE = 0  !! (as in `SUNMatrix_ID`)
    integer,parameter,public :: SUNMATRIX_BAND  = 3  !! (as in `SUNMatrix_ID`)

    type,abstract,public :: sun_matrix
        !! Base type for matrices (`SUNMatrix`).
    contains
        procedure(mat_id_fn),deferred      :: get_id      !! `SUNMatGetID`
        procedure(mat_zero_fn),deferred    :: zero        !! `SUNMatZero`: A = 0
        procedure(mat_copy_fn),deferred    :: copy_to     !! `SUNMatCopy(A, B)`: B = A
        procedure(mat_scaleaddi_fn),deferred :: scale_add_i !! `SUNMatScaleAddI`: A = c*A + I
        procedure(mat_matvec_fn),deferred  :: matvec      !! `SUNMatMatvec`: y = A*x
    end type sun_matrix

    type,extends(sun_matrix),public :: dense_matrix
        !! Dense matrix (`SUNDenseMatrix`).
        integer(ip) :: m = 0  !! number of rows
        integer(ip) :: n = 0  !! number of columns
        real(wp),dimension(:,:),allocatable :: data  !! data(1:m,1:n)
    contains
        procedure :: get_id      => dense_get_id  !! `SUNMatGetID`: returns `SUNMATRIX_DENSE`
        procedure :: zero        => dense_zero  !! `SUNMatZero`: A = 0
        procedure :: copy_to     => dense_copy  !! `SUNMatCopy`: B = A
        procedure :: scale_add_i => dense_scale_add_i  !! `SUNMatScaleAddI`: A = c*A + I
        procedure :: matvec      => dense_matvec  !! `SUNMatMatvec`: y = A*x
    end type dense_matrix

    type,extends(sun_matrix),public :: band_matrix
        !! Band matrix (`SUNBandMatrix`).
        integer(ip) :: n = 0    !! number of rows and columns
        integer(ip) :: mu = 0   !! upper bandwidth
        integer(ip) :: ml = 0   !! lower bandwidth
        integer(ip) :: smu = 0  !! storage upper bandwidth (mu + ml, for LU fill-in)
        integer(ip) :: ldim = 0 !! column storage length (smu + ml + 1)
        real(wp),dimension(:,:),allocatable :: data  !! data(0:ldim-1,0:n-1)
    contains
        procedure :: get_id      => band_get_id  !! `SUNMatGetID`: returns `SUNMATRIX_BAND`
        procedure :: zero        => band_zero  !! `SUNMatZero`: A = 0 (all storage)
        procedure :: copy_to     => band_copy  !! `SUNMatCopy`: B = A (B is enlarged if needed)
        procedure :: scale_add_i => band_scale_add_i  !! `SUNMatScaleAddI`: A = c*A + I on the band
        procedure :: matvec      => band_matvec  !! `SUNMatMatvec`: y = A*x
        procedure :: get         => band_get  !! element (i,j), 1-based (`SM_ELEMENT_B`)
        procedure :: set         => band_set  !! set element (i,j), 1-based (`SM_ELEMENT_B`)
    end type band_matrix

    interface dense_matrix
        module procedure dense_create
    end interface dense_matrix

    interface band_matrix
        module procedure band_create
    end interface band_matrix

    abstract interface
        pure function mat_id_fn(me) result(id)
            !! Matrix type ID (`SUNMatGetID`).
            import :: sun_matrix
            implicit none
            class(sun_matrix),intent(in) :: me  !! the matrix
            integer :: id  !! matrix type ID (`SUNMATRIX_DENSE`, `SUNMATRIX_BAND`)
        end function mat_id_fn
        subroutine mat_zero_fn(me)
            !! Set all entries to zero (`SUNMatZero`).
            import :: sun_matrix
            implicit none
            class(sun_matrix),intent(inout) :: me  !! the matrix, set to zero on return
        end subroutine mat_zero_fn
        subroutine mat_copy_fn(me, b, flag)
            !! Copy A into B (`SUNMatCopy`).
            import :: sun_matrix
            implicit none
            class(sun_matrix),intent(in) :: me  !! the matrix A to copy
            class(sun_matrix),intent(inout) :: b  !! the destination B (same type and size)
            integer,intent(out) :: flag  !! 0 on success, -1 if B is incompatible
        end subroutine mat_copy_fn
        subroutine mat_scaleaddi_fn(me, c)
            !! A = c*A + I (`SUNMatScaleAddI`).
            import :: sun_matrix, wp
            implicit none
            class(sun_matrix),intent(inout) :: me  !! the matrix A, overwritten by c*A + I
            real(wp),intent(in) :: c  !! the scale factor c
        end subroutine mat_scaleaddi_fn
        subroutine mat_matvec_fn(me, x, y)
            !! Matrix-vector product y = A*x (`SUNMatMatvec`).
            import :: sun_matrix, wp
            implicit none
            class(sun_matrix),intent(in) :: me  !! the matrix A
            real(wp),dimension(:),intent(in) :: x  !! the vector x
            real(wp),dimension(:),intent(out) :: y  !! the product y = A*x
        end subroutine mat_matvec_fn
    end interface

    contains
!*****************************************************************************************

!=========================================================================================
!  Dense
!=========================================================================================

    function dense_create(m, n) result(A)
    !! `SUNDenseMatrix(M, N)`: an m x n dense matrix, initialized to zero.
    integer(ip),intent(in) :: m  !! number of rows
    integer(ip),intent(in) :: n  !! number of columns
    type(dense_matrix) :: A  !! the new matrix
    A%m = m
    A%n = n
    allocate(A%data(m, n), source=0.0_wp)
    end function dense_create

    pure function dense_get_id(me) result(id)
    !! `SUNMatGetID_Dense`
    class(dense_matrix),intent(in) :: me  !! the matrix
    integer :: id  !! `SUNMATRIX_DENSE`
    id = SUNMATRIX_DENSE
    end function dense_get_id

    subroutine dense_zero(me)
    !! `SUNMatZero_Dense`: A = 0
    class(dense_matrix),intent(inout) :: me  !! the matrix, set to zero on return
    me%data = 0.0_wp
    end subroutine dense_zero

    subroutine dense_copy(me, b, flag)
    !! `SUNMatCopy_Dense`: b = me (same dimensions required).
    class(dense_matrix),intent(in) :: me  !! the matrix A to copy
    class(sun_matrix),intent(inout) :: b  !! the destination B (a dense matrix of the same size)
    integer,intent(out) :: flag  !! 0 on success, -1 if B is not a dense matrix of the same size
    flag = -1
    select type (b)
    type is (dense_matrix)
        if (b%m /= me%m .or. b%n /= me%n) return
        b%data = me%data
        flag = 0
    end select
    end subroutine dense_copy

    subroutine dense_scale_add_i(me, c)
    !! `SUNMatScaleAddI_Dense`: A = c*A + I
    class(dense_matrix),intent(inout) :: me  !! the matrix A, overwritten by c*A + I
    real(wp),intent(in) :: c  !! the scale factor c
    integer(ip) :: j  !! diagonal index
    me%data = c * me%data
    do j = 1, min(me%m, me%n)
        me%data(j,j) = me%data(j,j) + 1.0_wp
    end do
    end subroutine dense_scale_add_i

    subroutine dense_matvec(me, x, y)
    !! `SUNMatMatvec_Dense`: y = A*x
    class(dense_matrix),intent(in) :: me  !! the matrix A
    real(wp),dimension(:),intent(in) :: x  !! the vector x (length n)
    real(wp),dimension(:),intent(out) :: y  !! the product y = A*x (length m)
    y = matmul(me%data, x)
    end subroutine dense_matvec

!=========================================================================================
!  Band
!=========================================================================================

    function band_create(n, mu, ml) result(A)
    !! `SUNBandMatrix(N, mu, ml)`: an n x n band matrix with upper and lower
    !! bandwidths mu and ml, storage upper bandwidth smu = mu + ml (room for
    !! the LU fill-in), initialized to zero.
    integer(ip),intent(in) :: n  !! number of rows and columns
    integer(ip),intent(in) :: mu  !! upper bandwidth
    integer(ip),intent(in) :: ml  !! lower bandwidth
    type(band_matrix) :: A  !! the new matrix
    A%n    = n
    A%mu   = mu
    A%ml   = ml
    A%smu  = mu + ml
    A%ldim = A%smu + ml + 1
    allocate(A%data(0:A%ldim-1, 0:n-1), source=0.0_wp)
    end function band_create

    pure function band_get_id(me) result(id)
    !! `SUNMatGetID_Band`
    class(band_matrix),intent(in) :: me  !! the matrix
    integer :: id  !! `SUNMATRIX_BAND`
    id = SUNMATRIX_BAND
    end function band_get_id

    pure function band_get(me, i, j) result(aij)
    !! Element (i,j), 1-based (`SM_ELEMENT_B`). Must be inside the band.
    class(band_matrix),intent(in) :: me  !! the matrix
    integer(ip),intent(in) :: i  !! row index (1-based)
    integer(ip),intent(in) :: j  !! column index (1-based)
    real(wp) :: aij  !! the element A(i,j)
    aij = me%data(i - j + me%smu, j - 1)
    end function band_get

    subroutine band_set(me, i, j, aij)
    !! Set element (i,j), 1-based (`SM_ELEMENT_B`). Must be inside the band.
    class(band_matrix),intent(inout) :: me  !! the matrix
    integer(ip),intent(in) :: i  !! row index (1-based)
    integer(ip),intent(in) :: j  !! column index (1-based)
    real(wp),intent(in) :: aij  !! the value to store in A(i,j)
    me%data(i - j + me%smu, j - 1) = aij
    end subroutine band_set

    subroutine band_zero(me)
    !! `SUNMatZero_Band`: A = 0
    class(band_matrix),intent(inout) :: me  !! the matrix, set to zero on return
    me%data = 0.0_wp
    end subroutine band_zero

    subroutine band_copy(me, b, flag)
    !! `SUNMatCopy_Band`: b = me. b is enlarged if me has larger bandwidths.
    class(band_matrix),intent(in) :: me  !! the matrix A to copy
    class(sun_matrix),intent(inout) :: b  !! the destination B (band, same size)
    integer,intent(out) :: flag  !! 0 on success, -1 if B is incompatible
    integer(ip) :: ml, mu, smu  !! bandwidths (lower, upper, storage upper) of the enlarged B
    real(wp),dimension(:,:),allocatable :: newdata  !! storage for the enlarged B
    flag = -1
    select type (b)
    type is (band_matrix)
        if (b%n /= me%n) return
        ! grow b if a's bandwidth is larger
        if (me%mu > b%mu .or. me%ml > b%ml) then
            ml  = max(b%ml, me%ml)
            mu  = max(b%mu, me%mu)
            smu = max(b%smu, me%smu)
            allocate(newdata(0:smu+ml, 0:b%n-1), source=0.0_wp)
            newdata(smu-b%smu:smu-b%smu+b%ldim-1, :) = b%data
            call move_alloc(newdata, b%data)
            b%mu = mu; b%ml = ml; b%smu = smu; b%ldim = smu + ml + 1
        end if
        ! copy the band of a (rows -mu..ml relative to the diagonal)
        b%data(b%smu-me%mu:b%smu+me%ml, :) = me%data(me%smu-me%mu:me%smu+me%ml, :)
        flag = 0
    end select
    end subroutine band_copy

    subroutine band_scale_add_i(me, c)
    !! `SUNMatScaleAddI_Band`: A = c*A + I (on the band)
    class(band_matrix),intent(inout) :: me  !! the matrix A, overwritten by c*A + I on the band
    real(wp),intent(in) :: c  !! the scale factor c
    me%data(me%smu-me%mu:me%smu+me%ml, :) = c * me%data(me%smu-me%mu:me%smu+me%ml, :)
    me%data(me%smu, :) = me%data(me%smu, :) + 1.0_wp
    end subroutine band_scale_add_i

    subroutine band_matvec(me, x, y)
    !! `SUNMatMatvec_Band`: y = A*x
    class(band_matrix),intent(in) :: me  !! the matrix A
    real(wp),dimension(:),intent(in) :: x  !! the vector x (length n)
    real(wp),dimension(:),intent(out) :: y  !! the product y = A*x (length n)
    integer(ip) :: i, j  !! row and column indices (0-based)
    y = 0.0_wp
    do j = 0, me%n - 1
        do i = max(0_ip, j - me%mu), min(me%n - 1, j + me%ml)
            y(i+1) = y(i+1) + me%data(i - j + me%smu, j) * x(j+1)
        end do
    end do
    end subroutine band_matvec

    end module fvode_matrix
!*****************************************************************************************

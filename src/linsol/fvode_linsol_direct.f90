!*****************************************************************************************
!>
!  Direct linear solvers: dense and band LU with partial pivoting
!  (SUNDIALS v7.9.0 `SUNLinSol_Dense`, `SUNLinSol_Band`).
!
!  Translated from `src/sunlinsol/dense/sunlinsol_dense.c`,
!  `src/sunlinsol/band/sunlinsol_band.c`, and the factorization routines
!  `SUNDlsMat_denseGETRF/GETRS` (`src/sundials/sundials_dense.c`) and
!  `SUNDlsMat_bandGBTRF/GBTRS` (`src/sundials/sundials_band.c`). These are
!  SUNDIALS's own LU routines (not LAPACK), so pivoting matches CVODE.

    module fvode_linsol_direct

    use fvode_kinds, only: wp, ip
    use fvode_parameters, only: SUNLINEARSOLVER_DIRECT, SUN_SUCCESS, SUNLS_LUFACT_FAIL, &
                                SUN_ERR_ARG_INCOMPATIBLE
    use fvode_matrix, only: sun_matrix, dense_matrix, band_matrix
    use fvode_linsol, only: sun_linear_solver, ls_problem

    implicit none

    private

    type,extends(sun_linear_solver),public :: sun_linsol_dense
        !! Dense LU solver (`SUNLinSol_Dense`).
        integer(ip),dimension(:),allocatable :: pivots  !! row pivots, 1-based
    contains
        procedure :: get_type => dense_get_type  !! `SUNLinSolGetType_Dense`: direct solver
        procedure :: setup    => dense_setup  !! `SUNLinSolSetup_Dense`: LU factorization
        procedure :: solve    => dense_solve  !! `SUNLinSolSolve_Dense`: solve with the LU factors
    end type sun_linsol_dense

    type,extends(sun_linear_solver),public :: sun_linsol_band
        !! Band LU solver (`SUNLinSol_Band`).
        integer(ip),dimension(:),allocatable :: pivots  !! row pivots, 0-based (as in C)
    contains
        procedure :: get_type => band_get_type  !! `SUNLinSolGetType_Band`: direct solver
        procedure :: setup    => band_setup  !! `SUNLinSolSetup_Band`: band LU factorization
        procedure :: solve    => band_solve  !! `SUNLinSolSolve_Band`: solve with the LU factors
    end type sun_linsol_band

    interface sun_linsol_dense
        module procedure dense_create
    end interface sun_linsol_dense

    interface sun_linsol_band
        module procedure band_create
    end interface sun_linsol_band

    public :: dense_getrf, dense_getrs, band_gbtrf, band_gbtrs

    contains
!*****************************************************************************************

    function dense_create(n) result(S)
    !! `SUNLinSol_Dense(y, A)` for an n x n system.
    integer(ip),intent(in) :: n  !! number of equations
    type(sun_linsol_dense) :: S  !! the new solver
    allocate(S%pivots(n), source=0_ip)
    end function dense_create

    function band_create(n) result(S)
    !! `SUNLinSol_Band(y, A)` for an n x n system.
    integer(ip),intent(in) :: n  !! number of equations
    type(sun_linsol_band) :: S  !! the new solver
    allocate(S%pivots(0:n-1), source=0_ip)
    end function band_create

    pure function dense_get_type(me) result(itype)
    !! `SUNLinSolGetType_Dense`
    class(sun_linsol_dense),intent(in) :: me  !! the linear solver
    integer :: itype  !! `SUNLINEARSOLVER_DIRECT`
    itype = SUNLINEARSOLVER_DIRECT
    end function dense_get_type

    pure function band_get_type(me) result(itype)
    !! `SUNLinSolGetType_Band`
    class(sun_linsol_band),intent(in) :: me  !! the linear solver
    integer :: itype  !! `SUNLINEARSOLVER_DIRECT`
    itype = SUNLINEARSOLVER_DIRECT
    end function band_get_type

!*****************************************************************************************
!>
!  `SUNLinSolSetup_Dense`: LU factorization of A (in place).

    function dense_setup(me, A, prob) result(flag)

    class(sun_linsol_dense),intent(inout) :: me  !! the linear solver (stores the pivots)
    class(sun_matrix),intent(inout),optional :: A  !! a `dense_matrix`, overwritten by its LU factors
    class(ls_problem),intent(inout),optional :: prob  !! the problem (not used by direct solvers)
    integer :: flag  !! `SUN_SUCCESS`, `SUNLS_LUFACT_FAIL` (zero pivot), or < 0

    flag = SUN_ERR_ARG_INCOMPATIBLE
    if (.not. present(A)) return
    select type (A)
    type is (dense_matrix)
        me%last_flag = int(dense_getrf(A%data, me%pivots))
        ! nonzero: this column encountered a zero pivot
        if (me%last_flag > 0) then
            flag = SUNLS_LUFACT_FAIL
        else
            flag = SUN_SUCCESS
        end if
    end select

    end function dense_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSolve_Dense`: solve using the LU factors.

    function dense_solve(me, A, x, b, tol, prob) result(flag)

    class(sun_linsol_dense),intent(inout) :: me  !! the solver (pivots from `setup`)
    class(sun_matrix),intent(inout),optional :: A  !! the LU factors from `setup` (a `dense_matrix`)
    class(ls_problem),intent(inout),optional :: prob  !! the problem (not used by direct solvers)
    real(wp),dimension(:),intent(inout) :: x  !! the solution x
    real(wp),dimension(:),intent(in) :: b  !! the right-hand side b
    real(wp),intent(in) :: tol  !! tolerance (unused by direct solvers)
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_INCOMPATIBLE`

    flag = SUN_ERR_ARG_INCOMPATIBLE
    if (.not. present(A)) return
    select type (A)
    type is (dense_matrix)
        x = b
        call dense_getrs(A%data, me%pivots, x)
        me%last_flag = SUN_SUCCESS
        flag = SUN_SUCCESS
    end select

    end function dense_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNDlsMat_denseGETRF`: LU factorization with partial pivoting of the
!  m x n matrix a (m >= n), overwritten by the factors. On return the
!  multipliers are stored below the diagonal, and `p(k)` is the row swapped
!  with row k. Returns 0 on success, or k > 0 if `U(k,k) = 0`.

    function dense_getrf(a, p) result(info)

    real(wp),dimension(:,:),intent(inout) :: a  !! in: the m x n matrix; out: its LU factors
    integer(ip),dimension(:),intent(out) :: p  !! pivot rows: row k was swapped with row `p(k)`
    integer(ip) :: info  !! 0 on success, k if the k-th pivot is zero

    integer(ip) :: i, j, k, l, m, n  !! row, column and step indices, pivot row, matrix dimensions
    real(wp) :: mult, a_kj  !! row multiplier 1/a(k,k), and the element a(k,j)
    real(wp) :: temp  !! temporary for swapping (C swaps element by element; an automatic
                      !! row here would be heap-allocated by gfortran on every call)

    m = size(a, 1)
    n = size(a, 2)

    ! k-th elimination step number
    do k = 1, n

        ! find l = pivot row number
        l = k
        do i = k + 1, m
            if (abs(a(i,k)) > abs(a(l,k))) l = i
        end do
        p(k) = l

        ! check for zero pivot element
        if (a(l,k) == 0.0_wp) then
            info = k
            return
        end if

        ! swap a(k,:) and a(l,:) if necessary
        if (l /= k) then
            do j = 1, n
                temp = a(l,j)
                a(l,j) = a(k,j)
                a(k,j) = temp
            end do
        end if

        ! scale the elements below the diagonal in column k by 1/a(k,k):
        ! this stores the pivot row multipliers a(i,k)/a(k,k)
        mult = 1.0_wp / a(k,k)
        a(k+1:m,k) = mult * a(k+1:m,k)

        ! row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1
        ! row k is the pivot row after swapping with row l.
        ! The computation is done one column at a time, column j=k+1, ..., n-1.
        do j = k + 1, n
            a_kj = a(k,j)
            ! a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j); a_kj = a(k,j), col_k[i] = a(i,k)/a(k,k)
            if (a_kj /= 0.0_wp) a(k+1:m,j) = a(k+1:m,j) - a_kj * a(k+1:m,k)
        end do

    end do

    ! return 0 to indicate success
    info = 0

    end function dense_getrf
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNDlsMat_denseGETRS`: solve A x = b with the LU factors from
!  [[dense_getrf]]; b is overwritten by x.

    subroutine dense_getrs(a, p, b)

    real(wp),dimension(:,:),intent(in) :: a  !! the LU factors from `dense_getrf`
    integer(ip),dimension(:),intent(in) :: p  !! the pivots from `dense_getrf`
    real(wp),dimension(:),intent(inout) :: b  !! on input the right-hand side b, on output the solution x

    integer(ip) :: k, n, pk  !! step index, matrix dimension, pivot row
    real(wp) :: tmp  !! temporary for swapping

    n = size(a, 2)

    ! permute b, based on pivot information in p
    do k = 1, n
        pk = p(k)
        if (pk /= k) then
            tmp   = b(k)
            b(k)  = b(pk)
            b(pk) = tmp
        end if
    end do

    ! solve Ly = b, store solution y in b
    do k = 1, n - 1
        b(k+1:n) = b(k+1:n) - a(k+1:n,k) * b(k)
    end do

    ! solve Ux = y, store solution x in b
    do k = n, 2, -1
        b(k) = b(k) / a(k,k)
        b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
    end do
    b(1) = b(1) / a(1,1)

    end subroutine dense_getrs
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSetup_Band`: band LU factorization of A (in place).

    function band_setup(me, A, prob) result(flag)

    class(sun_linsol_band),intent(inout) :: me  !! the linear solver (stores the pivots)
    class(sun_matrix),intent(inout),optional :: A  !! a `band_matrix`, overwritten by its LU factors
    class(ls_problem),intent(inout),optional :: prob  !! the problem (not used by direct solvers)
    integer :: flag  !! `SUN_SUCCESS`, `SUNLS_LUFACT_FAIL` (zero pivot), or < 0

    flag = SUN_ERR_ARG_INCOMPATIBLE
    if (.not. present(A)) return
    select type (A)
    type is (band_matrix)
        ! ensure that the storage upper bandwidth is large enough for the fill-in
        if (A%smu < min(A%n - 1, A%mu + A%ml)) return
        me%last_flag = int(band_gbtrf(A%data, A%n, A%mu, A%ml, A%smu, me%pivots))
        ! nonzero: that row encountered a zero pivot
        if (me%last_flag > 0) then
            flag = SUNLS_LUFACT_FAIL
        else
            flag = SUN_SUCCESS
        end if
    end select

    end function band_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNLinSolSolve_Band`: solve using the band LU factors.

    function band_solve(me, A, x, b, tol, prob) result(flag)

    class(sun_linsol_band),intent(inout) :: me  !! the solver (pivots from `setup`)
    class(sun_matrix),intent(inout),optional :: A  !! the band LU factors from `setup` (a `band_matrix`)
    class(ls_problem),intent(inout),optional :: prob  !! the problem (not used by direct solvers)
    real(wp),dimension(:),intent(inout) :: x  !! the solution x
    real(wp),dimension(:),intent(in) :: b  !! the right-hand side b
    real(wp),intent(in) :: tol  !! tolerance (unused by direct solvers)
    integer :: flag  !! `SUN_SUCCESS` or `SUN_ERR_ARG_INCOMPATIBLE`

    flag = SUN_ERR_ARG_INCOMPATIBLE
    if (.not. present(A)) return
    select type (A)
    type is (band_matrix)
        x = b
        call band_gbtrs(A%data, A%n, A%smu, A%ml, me%pivots, x)
        me%last_flag = SUN_SUCCESS
        flag = SUN_SUCCESS
    end select

    end function band_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNDlsMat_bandGBTRF`: LU factorization with partial pivoting of the n x n
!  band matrix stored in `a(0:smu+ml, 0:n-1)` (element (i,j) at
!  `a(i-j+smu, j)`, 0-based). Indices are 0-based, as in C. Returns 0 on
!  success, or k > 0 if `U(k,k) = 0` (1-based k).

    function band_gbtrf(a, n, mu, ml, smu, p) result(info)

    real(wp),dimension(0:,0:),intent(inout) :: a  !! band storage, (i,j) at a(i-j+smu,j); out: LU
    integer(ip),intent(in) :: n    !! number of equations
    integer(ip),intent(in) :: mu   !! upper bandwidth
    integer(ip),intent(in) :: ml   !! lower bandwidth
    integer(ip),intent(in) :: smu  !! storage upper bandwidth (>= min(n-1, mu+ml))
    integer(ip),dimension(0:),intent(out) :: p  !! pivots (0-based): row k swapped with row `p(k)`
    integer(ip) :: info  !! 0 on success, k if the k-th pivot is zero (1-based)

    integer(ip) :: c                       !! column index (zeroing the fill-in rows)
    integer(ip) :: num_rows                !! number of fill-in rows, smu - mu
    integer(ip) :: i, j, k, l              !! row, column, step and pivot row indices
    integer(ip) :: storage_l, storage_k    !! storage rows of elements (l,j) and (k,j)
    integer(ip) :: last_col_k, last_row_k  !! last column and row changed by step k
    real(wp) :: max, temp, mult, a_kj  !! pivot size, swap temp, multiplier, a(k,j)
    logical :: swap  !! rows k and l are swapped in this step

    ! zero out the first smu - mu rows of the rectangular array a
    num_rows = smu - mu
    if (num_rows > 0) then
        do c = 0, n - 1
            a(0:num_rows-1, c) = 0.0_wp
        end do
    end if

    ! k = elimination step number
    do k = 0, n - 2

        ! (diagonal of column k is at a(smu,k); subdiagonal entries below it)
        last_row_k = min(n - 1, k + ml)

        ! find l = pivot row number
        l = k
        max = abs(a(smu, k))
        do i = k + 1, last_row_k
            if (abs(a(i - k + smu, k)) > max) then
                l = i
                max = abs(a(i - k + smu, k))
            end if
        end do
        storage_l = l - k + smu
        p(k) = l

        ! check for zero pivot element
        if (a(storage_l, k) == 0.0_wp) then
            info = k + 1
            return
        end if

        ! swap a(l,k) and a(k,k) if necessary
        swap = (l /= k)
        if (swap) then
            temp = a(storage_l, k)
            a(storage_l, k) = a(smu, k)
            a(smu, k) = temp
        end if

        ! scale the elements below the diagonal in column k by -1.0 / a(k,k).
        ! After the above swap, a(k,k) holds the pivot element. This scaling
        ! stores the pivot row multipliers -a(i,k)/a(k,k) in a(i,k), i=k+1, ...
        mult = -1.0_wp / a(smu, k)
        a(smu+1:smu+last_row_k-k, k) = mult * a(smu+1:smu+last_row_k-k, k)

        ! row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., last_row_k
        ! row k is the pivot row after swapping with row l.
        ! The computation is done one column at a time, column j=k+1, ..., last_col_k
        last_col_k = min(k + smu, n - 1)
        do j = k + 1, last_col_k

            storage_l = l - j + smu
            storage_k = k - j + smu
            a_kj = a(storage_l, j)

            ! swap the elements a(k,j) and a(l,j) if necessary
            if (swap) then
                a(storage_l, j) = a(storage_k, j)
                a(storage_k, j) = a_kj
            end if

            ! a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)
            ! (the multipliers stored are -a(i,k)/a(k,k), so add)
            if (a_kj /= 0.0_wp) then
                a(k+1-j+smu:last_row_k-j+smu, j) = a(k+1-j+smu:last_row_k-j+smu, j) + &
                                                   a_kj * a(smu+1:smu+last_row_k-k, k)
            end if
        end do

    end do

    ! set the last pivot row to be n-1 and check for a zero pivot
    p(n-1) = n - 1
    if (a(smu, n-1) == 0.0_wp) then
        info = n
        return
    end if

    ! return 0 to indicate success
    info = 0

    end function band_gbtrf
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNDlsMat_bandGBTRS`: solve A x = b with the band LU factors from
!  [[band_gbtrf]]; b (1-based array) is overwritten by x.

    subroutine band_gbtrs(a, n, smu, ml, p, bvec)

    real(wp),dimension(0:,0:),intent(in) :: a  !! the band LU factors from `band_gbtrf`
    integer(ip),intent(in) :: n    !! number of equations
    integer(ip),intent(in) :: smu  !! storage upper bandwidth
    integer(ip),intent(in) :: ml   !! lower bandwidth
    integer(ip),dimension(0:),intent(in) :: p  !! the pivots from `band_gbtrf`
    real(wp),dimension(0:),intent(inout) :: bvec  !! on input the right-hand side b, on output the solution x (0-based)

    integer(ip) :: k, l, i                  !! step, pivot row and row indices
    integer(ip) :: first_row_k, last_row_k  !! first and last rows in column k
    real(wp) :: mult  !! multiplier (pivoted b(k), then -x(k))

    ! solve Ly = Pb, store solution y in b
    do k = 0, n - 2
        l = p(k)
        mult = bvec(l)
        if (l /= k) then
            bvec(l) = bvec(k)
            bvec(k) = mult
        end if
        last_row_k = min(n - 1, k + ml)
        do i = k + 1, last_row_k
            bvec(i) = bvec(i) + mult * a(i - k + smu, k)
        end do
    end do

    ! solve Ux = y, store solution x in b
    do k = n - 1, 0, -1
        first_row_k = max(0_ip, k - smu)
        bvec(k) = bvec(k) / a(smu, k)
        mult = -bvec(k)
        do i = first_row_k, k - 1
            bvec(i) = bvec(i) + mult * a(i - k + smu, k)
        end do
    end do

    end subroutine band_gbtrs
!*****************************************************************************************

    end module fvode_linsol_direct
!*****************************************************************************************

!*****************************************************************************************
!>
!  Orthogonalization and QR routines shared by the Krylov solvers
!  (SUNDIALS v7.9.0 `src/sundials/sundials_iterative.c`: `SUNModifiedGS`,
!  `SUNClassicalGS`, `SUNQRfact`, `SUNQRsol`).
!
!  Krylov vectors are the columns `v(:,0:)` of a 2-D array and the
!  Hessenberg matrix is `h(0:,0:)`, with the C (0-based) indices.

    module fvode_iterative

    use fvode_kinds, only: wp

    implicit none

    private

    real(wp),parameter :: ZERO   = 0.0_wp    !! 0
    real(wp),parameter :: ONE    = 1.0_wp    !! 1
    real(wp),parameter :: FACTOR = 1000.0_wp !! reorthogonalization threshold factor

    public :: modified_gs, classical_gs, qr_fact, qr_sol

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNModifiedGS`: orthogonalize `v(:,k)` against the previous `p` vectors
!  `v(:,max(k-p,0):k-1)` by modified Gram-Schmidt, storing the coefficients
!  in `h(i,k-1)`. If the new vector is much shorter than the original
!  (cancellation), it is reorthogonalized once.

    subroutine modified_gs(v, h, k, p, new_vk_norm)

    real(wp),dimension(:,0:),intent(inout) :: v  !! the vectors; v(:,k) is orthogonalized
    real(wp),dimension(0:,0:),intent(inout) :: h  !! Hessenberg matrix (column k-1 is set)
    integer,intent(in) :: k  !! index of the vector to orthogonalize
    integer,intent(in) :: p  !! number of previous vectors to orthogonalize against
    real(wp),intent(out) :: new_vk_norm  !! 2-norm of the orthogonalized v(:,k)

    integer :: i  !! vector index
    integer :: k_minus_1  !! k - 1
    integer :: i0  !! first vector index
    real(wp) :: new_norm_2  !! sum of squares of the reorthogonalization coefficients
    real(wp) :: new_product  !! a reorthogonalization coefficient
    real(wp) :: vk_norm  !! norm of the original v(:,k)
    real(wp) :: temp  !! temporary

    ! get norm of the input vector
    vk_norm = sqrt(dot_product(v(:,k), v(:,k)))

    ! perform modified Gram-Schmidt
    k_minus_1 = k - 1
    i0 = max(k - p, 0)
    do i = i0, k - 1
        h(i,k_minus_1) = dot_product(v(:,i), v(:,k))
        v(:,k) = v(:,k) - h(i,k_minus_1) * v(:,i)
    end do

    ! compute the norm of the new vector at v[k]
    new_vk_norm = sqrt(dot_product(v(:,k), v(:,k)))

    ! if the norm of the new vector at v[k] is less than FACTOR (== 1000)
    ! times unit roundoff times the norm of the input vector v[k], then the
    ! vector will be reorthogonalized in order to ensure that
    ! nonorthogonality is not being masked by a very small vector length.
    temp = FACTOR * vk_norm
    if ((temp + new_vk_norm) /= temp) return

    new_norm_2 = ZERO

    do i = i0, k - 1
        new_product = dot_product(v(:,i), v(:,k))
        temp = FACTOR * h(i,k_minus_1)
        if ((temp + new_product) == temp) cycle
        h(i,k_minus_1) = h(i,k_minus_1) + new_product
        v(:,k) = v(:,k) - new_product * v(:,i)
        new_norm_2 = new_norm_2 + new_product**2
    end do

    if (new_norm_2 /= ZERO) then
        new_product = new_vk_norm**2 - new_norm_2
        if (new_product > ZERO) then
            new_vk_norm = sqrt(new_product)
        else
            new_vk_norm = ZERO
        end if
    end if

    end subroutine modified_gs
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNClassicalGS`: orthogonalize `v(:,k)` against the previous `p` vectors
!  by classical Gram-Schmidt (all coefficients computed at once), with one
!  reorthogonalization if the norm drops by more than a factor 1000.
!
!  @note As in C, the coefficient and vector indices in the update use `i`
!  where `i + i0` would be expected. The Krylov solvers call this with
!  `p = maxl >= k`, so `i0 = 0` and both are the same.

    subroutine classical_gs(v, h, k, p, new_vk_norm, stemp)

    real(wp),dimension(:,0:),intent(inout) :: v  !! the vectors; v(:,k) is orthogonalized
    real(wp),dimension(0:,0:),intent(inout) :: h  !! Hessenberg matrix (column k-1 is set)
    integer,intent(in) :: k  !! index of the vector to orthogonalize
    integer,intent(in) :: p  !! number of previous vectors to orthogonalize against
    real(wp),intent(out) :: new_vk_norm  !! 2-norm of the orthogonalized v(:,k)
    real(wp),dimension(0:),intent(inout) :: stemp  !! work array, size >= k+1

    integer :: i  !! vector index
    integer :: i0  !! first vector index
    integer :: k_minus_1  !! k - 1
    real(wp) :: vk_norm  !! norm of the original v(:,k)

    k_minus_1 = k - 1
    i0 = max(k - p, 0)

    ! compute the dot products with v(:,i0..k): stemp(0:k-i0)
    do i = 0, k - i0
        stemp(i) = dot_product(v(:,k), v(:,i0+i))
    end do

    ! compute the norm of the vector v[k]
    vk_norm = sqrt(stemp(k - i0))

    ! set the coefficients and update v[k] = v[k] - sum_i stemp[i] v[i]
    do i = k - i0 - 1, 0, -1
        h(i,k_minus_1) = stemp(i)
    end do
    do i = 0, k - i0 - 1
        v(:,k) = v(:,k) - stemp(i) * v(:,i)
    end do

    ! compute the norm of the new vector at v[k]
    new_vk_norm = sqrt(dot_product(v(:,k), v(:,k)))

    ! reorthogonalize if necessary
    if ((FACTOR * new_vk_norm) < vk_norm) then
        do i = 1, k - i0
            stemp(i) = dot_product(v(:,k), v(:,i0+i-1))
        end do
        do i = i0, k - 1
            h(i,k_minus_1) = h(i,k_minus_1) + stemp(i - i0 + 1)
        end do
        do i = i0, k - 1
            v(:,k) = v(:,k) - stemp(i - i0 + 1) * v(:,i - i0)
        end do
        new_vk_norm = sqrt(dot_product(v(:,k), v(:,k)))
    end if

    end subroutine classical_gs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Givens rotation (c, s) that zeroes `temp2` in `[c -s; s c] * [temp1; temp2]`
!  (the rotation used by `SUNQRfact`).

    pure subroutine givens_rotation(temp1, temp2, c, s)

    real(wp),intent(in) :: temp1  !! first entry
    real(wp),intent(in) :: temp2  !! second entry (to be eliminated)
    real(wp),intent(out) :: c  !! cosine of the rotation
    real(wp),intent(out) :: s  !! sine of the rotation

    real(wp) :: temp3  !! ratio of the entries

    if (temp2 == ZERO) then
        c = ONE
        s = ZERO
    else if (abs(temp2) >= abs(temp1)) then
        temp3 = temp1 / temp2
        s = -ONE / sqrt(ONE + temp3**2)
        c = -s * temp3
    else
        temp3 = temp2 / temp1
        c = ONE / sqrt(ONE + temp3**2)
        s = -c * temp3
    end if

    end subroutine givens_rotation
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNQRfact`: QR factorization of the (n+1) x n upper Hessenberg matrix h
!  by Givens rotations, stored in `q` (`q(2k)`, `q(2k+1)` = c, s of rotation k).
!
!  * `job = 0`: factor the whole matrix;
!  * otherwise: update the factorization for a new last column
!    (apply the previous rotations to column n-1 and compute a new one).
!
!  Returns 0 on success, or k > 0 if `R(k-1,k-1)` (0-based) is zero.

    function qr_fact(n, h, q, job) result(code)

    integer,intent(in) :: n  !! number of columns
    real(wp),dimension(0:,0:),intent(inout) :: h  !! the Hessenberg matrix, overwritten by R
    real(wp),dimension(0:),intent(inout) :: q  !! the Givens rotations, size 2n
    integer,intent(in) :: job  !! 0: full factorization; otherwise: add the last column
    integer :: code  !! 0, or the (1-based) index of a zero diagonal element of R

    real(wp) :: c, s  !! cosine and sine of a rotation
    real(wp) :: temp1, temp2  !! entries rotated
    integer :: i, j, k  !! indices
    integer :: q_ptr  !! index into q
    integer :: n_minus_1  !! n - 1

    code = 0

    if (job == 0) then

        ! compute a new factorization of H
        do k = 0, n - 1

            ! multiply column k by the previous k-1 Givens rotations
            do j = 0, k - 2
                i = 2 * j
                temp1 = h(j,k)
                temp2 = h(j+1,k)
                c = q(i)
                s = q(i+1)
                h(j,k)   = c * temp1 - s * temp2
                h(j+1,k) = s * temp1 + c * temp2
            end do

            ! compute the Givens rotation components c and s
            q_ptr = 2 * k
            temp1 = h(k,k)
            temp2 = h(k+1,k)
            call givens_rotation(temp1, temp2, c, s)
            q(q_ptr)   = c
            q(q_ptr+1) = s
            h(k,k) = c * temp1 - s * temp2
            if (h(k,k) == ZERO) code = k + 1
        end do

    else

        ! update the factored H to which a new column has been added
        n_minus_1 = n - 1

        ! multiply the new column by the previous n-1 Givens rotations
        do k = 0, n_minus_1 - 1
            i = 2 * k
            temp1 = h(k,n_minus_1)
            temp2 = h(k+1,n_minus_1)
            c = q(i)
            s = q(i+1)
            h(k,n_minus_1)   = c * temp1 - s * temp2
            h(k+1,n_minus_1) = s * temp1 + c * temp2
        end do

        ! compute new Givens rotation and multiply it times the last two
        ! entries in the new column of H. Note that the second entry of
        ! this product will be 0, so it is not necessary to compute it.
        temp1 = h(n_minus_1,n_minus_1)
        temp2 = h(n,n_minus_1)
        call givens_rotation(temp1, temp2, c, s)
        q_ptr = 2 * n_minus_1
        q(q_ptr)   = c
        q(q_ptr+1) = s
        h(n_minus_1,n_minus_1) = c * temp1 - s * temp2
        if (h(n_minus_1,n_minus_1) == ZERO) code = n

    end if

    end function qr_fact
!*****************************************************************************************

!*****************************************************************************************
!>
!  `SUNQRsol`: solve the least-squares problem `min ||b - H x||` with the QR
!  factorization from [[qr_fact]]: `b` (size n+1) is overwritten, the
!  solution is in `b(0:n-1)`. Returns 0, or k > 0 if `R(k-1,k-1)` is zero.

    function qr_sol(n, h, q, b) result(code)

    integer,intent(in) :: n  !! number of columns
    real(wp),dimension(0:,0:),intent(in) :: h  !! the factor R (from qr_fact)
    real(wp),dimension(0:),intent(in) :: q  !! the Givens rotations
    real(wp),dimension(0:),intent(inout) :: b  !! right-hand side (size n+1); solution in b(0:n-1)
    integer :: code  !! 0, or the (1-based) index of a zero diagonal element of R

    real(wp) :: c, s  !! cosine and sine of a rotation
    real(wp) :: temp1, temp2  !! entries rotated
    integer :: i, k  !! indices
    integer :: q_ptr  !! index into q

    code = 0

    ! compute Q*b
    do k = 0, n - 1
        q_ptr = 2 * k
        c = q(q_ptr)
        s = q(q_ptr+1)
        temp1 = b(k)
        temp2 = b(k+1)
        b(k)   = c * temp1 - s * temp2
        b(k+1) = s * temp1 + c * temp2
    end do

    ! solve R*x = Q*b
    do k = n - 1, 0, -1
        if (h(k,k) == ZERO) then
            code = k + 1
            exit
        end if
        b(k) = b(k) / h(k,k)
        do i = 0, k - 1
            b(i) = b(i) - b(k) * h(i,k)
        end do
    end do

    end function qr_sol
!*****************************************************************************************

    end module fvode_iterative
!*****************************************************************************************

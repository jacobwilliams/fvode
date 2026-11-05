!*****************************************************************************************
!>
!  2-D advection-diffusion equation, BDF + Newton + band linear solver with a
!  user-supplied Jacobian.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvAdvDiff_bnd.c`
!  (same problem, settings and output; compare with `cvAdvDiff_bnd.out`).
!
!  The equation is
!
!      du/dt = d2u/dx2 + 0.5 du/dx + d2u/dy2
!
!  on 0 <= x <= 2, 0 <= y <= 1, with zero boundary values and
!  `u(x,y,0) = x(2-x) y(1-y) exp(5xy)`, discretized by central differences
!  on a 10 x 5 interior mesh (50 equations, bandwidths mu = ml = 5).
!
!  It also shows how to pass problem data to the callbacks: the solver type
!  is extended with the grid coefficients (in C, `CVodeSetUserData`).

    module cv_advdiff_bnd_module

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: sun_matrix, band_matrix

    implicit none

    private

    real(wp),parameter,public :: XMAX = 2.0_wp  !! domain boundaries
    real(wp),parameter,public :: YMAX = 1.0_wp  !! upper y boundary of the domain
    integer,parameter,public  :: MX = 10        !! mesh dimensions
    integer,parameter,public  :: MY = 5  !! number of interior mesh points in y
    integer,parameter,public  :: NEQ = MX * MY  !! number of equations

    type,extends(cvode_solver),public :: advdiff_solver
        !! the CVODE solver plus the problem data (C: UserData)
        real(wp) :: dx = 0.0_wp, dy = 0.0_wp  !! mesh spacing in x and y
        real(wp) :: hdcoef = 0.0_wp, hacoef = 0.0_wp, vdcoef = 0.0_wp  !! coefficients of the horizontal diffusion, horizontal advection and vertical diffusion terms
    end type advdiff_solver

    public :: f, jac, ijth

    contains

    pure integer function ijth(i, j)
    !! index of u(i,j) in the solution vector (C: IJth, 1-based i, j)
    integer,intent(in) :: i, j  !! mesh indices in x and y (1-based)
    ijth = j + (i - 1) * MY
    end function ijth

    function f(me, t, y, ydot) result(retval)
    !! right-hand side: central differences, zero boundary values
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y  !! current state vector
    real(wp),dimension(:),intent(out) :: ydot  !! time derivative dy/dt = f(t,y)
    integer :: retval  !! return value: 0 on success
    real(wp) :: uij, udn, uup, ult, urt, hdiff, hadv, vdiff  !! u at (i,j) and its four neighbors, and the diffusion and advection terms
    integer :: i, j  !! mesh indices in x and y
    select type (me)
    class is (advdiff_solver)
        do j = 1, MY
            do i = 1, MX
                ! extract u at x_i, y_j and four neighboring points
                uij = y(ijth(i, j))
                udn = 0.0_wp; if (j /= 1)  udn = y(ijth(i, j - 1))
                uup = 0.0_wp; if (j /= MY) uup = y(ijth(i, j + 1))
                ult = 0.0_wp; if (i /= 1)  ult = y(ijth(i - 1, j))
                urt = 0.0_wp; if (i /= MX) urt = y(ijth(i + 1, j))
                ! set diffusion and advection terms and load into udot
                hdiff = me%hdcoef * (ult - 2.0_wp * uij + urt)
                hadv  = me%hacoef * (urt - ult)
                vdiff = me%vdcoef * (uup - 2.0_wp * uij + udn)
                ydot(ijth(i, j)) = hdiff + hadv + vdiff
            end do
        end do
    end select
    retval = 0
    end function f

    function jac(me, t, y, fy, J, tmp1, tmp2, tmp3) result(retval)
    !! band Jacobian
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current value of the independent variable
    real(wp),dimension(:),intent(in) :: y, fy  !! current state vector and f(t,y)
    class(sun_matrix),intent(inout) :: J  !! Jacobian matrix df/dy (filled in here)
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3  !! work vectors (not used)
    integer :: retval  !! return value: 0 on success
    integer :: i, jj  !! mesh indices in x and y
    integer(ip) :: k  !! index of u(i,j) in the solution vector (the Jacobian column)
    select type (me)
    class is (advdiff_solver)
        select type (J)
        type is (band_matrix)
            do jj = 1, MY
                do i = 1, MX
                    k = ijth(i, jj)
                    call J%set(k, k, -2.0_wp * (me%vdcoef + me%hdcoef))
                    if (i /= 1)   call J%set(k - MY, k, me%hdcoef + me%hacoef)
                    if (i /= MX)  call J%set(k + MY, k, me%hdcoef - me%hacoef)
                    if (jj /= 1)  call J%set(k - 1, k, me%vdcoef)
                    if (jj /= MY) call J%set(k + 1, k, me%vdcoef)
                end do
            end do
        end select
    end select
    retval = 0
    end function jac

    end module cv_advdiff_bnd_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_advdiff_bnd

    use fvode
    use fvode_matrix, only: band_matrix
    use fvode_linsol_direct, only: sun_linsol_band
    use cv_advdiff_bnd_module

    implicit none

    real(wp),parameter :: ATOL  = 1.0e-5_wp  !! scalar absolute tolerance
    real(wp),parameter :: T0    = 0.0_wp     !! initial time
    real(wp),parameter :: T1    = 0.1_wp     !! first output time
    real(wp),parameter :: DTOUT = 0.1_wp     !! output time increment
    integer,parameter  :: NOUT  = 10         !! number of output times

    type(advdiff_solver) :: cvode  !! the CVODE solver with the problem data
    real(wp),dimension(NEQ) :: u  !! solution vector
    real(wp) :: reltol, abstol, t, tout, x, y  !! tolerances, current and output times, and mesh coordinates
    integer :: iout, retval, i, j  !! output counter, return flag and mesh indices

    reltol = 0.0_wp
    abstol = ATOL

    ! set grid coefficients
    cvode%dx = XMAX / (MX + 1)
    cvode%dy = YMAX / (MY + 1)
    cvode%hdcoef = 1.0_wp / (cvode%dx * cvode%dx)
    cvode%hacoef = 0.5_wp / (2.0_wp * cvode%dx)
    cvode%vdcoef = 1.0_wp / (cvode%dy * cvode%dy)

    ! initial profile
    do j = 1, MY
        y = j * cvode%dy
        do i = 1, MX
            x = i * cvode%dx
            u(ijth(i, j)) = x * (XMAX - x) * y * (YMAX - y) * exp(5.0_wp * x * y)
        end do
    end do

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, T0, u, retval)
    call cvode%ss_tolerances(reltol, abstol, retval)

    ! band matrix and linear solver (the storage bandwidth for the LU
    ! fill-in is set automatically to mu + ml)
    call cvode%set_linear_solver(sun_linsol_band(int(NEQ, ip)), &
                                 band_matrix(int(NEQ, ip), int(MY, ip), int(MY, ip)), retval)
    call cvode%set_jac_fn(jac, retval)

    ! in loop over output points: call CVode, print results, test for errors
    write(*,'(/a)') '2-D Advection-Diffusion Equation'
    write(*,'(a,i0,a,i0)') 'Mesh dimensions = ', MX, ' X ', MY
    write(*,'(a,i0)') 'Total system size = ', NEQ
    write(*,'(a/)') 'Tolerance parameters: reltol = 0   abstol = 1e-05'
    write(*,'(a,es14.6,a)') 'At t = 0      max.norm(u) =', maxval(abs(u)), ' '

    tout = T1
    do iout = 1, NOUT
        call cvode%solve(tout, u, t, CV_NORMAL, retval)
        write(*,'(a,f4.2,a,es14.6,a,i4)') 'At t = ', t, '   max.norm(u) =', &
            maxval(abs(u)), '   nst = ', cvode%get_num_steps()
        if (retval < 0) exit
        tout = tout + DTOUT
    end do

    write(*,'(/a)') 'Final Statistics:'
    write(*,'(5a)') 'nst = ', lj(cvode%get_num_steps()), ' nfe  = ', lj(cvode%get_num_rhs_evals()), &
        ' nsetups = '//lj(cvode%get_num_lin_solv_setups())//' nfeLS = '// &
        lj(cvode%get_num_lin_rhs_evals())//' nje = '//itoa(cvode%get_num_jac_evals())
    write(*,'(5a)') 'nni = ', lj(cvode%get_num_nonlin_solv_iters()), ' ncfn = ', &
        lj(cvode%get_num_nonlin_solv_conv_fails()), ' netf = '//itoa(cvode%get_num_err_test_fails())

    call cvode%destroy()

    contains

        function itoa(i) result(s)
        !! integer to string (no blanks)
        integer(li),intent(in) :: i  !! the integer to convert
        character(len=:),allocatable :: s  !! the integer as a string
        character(len=24) :: buf  !! write buffer
        write(buf,'(i0)') i
        s = trim(buf)
        end function itoa

        function lj(i) result(s)
        !! left-justified in 6 characters (C: %-6ld)
        integer(li),intent(in) :: i  !! the integer to format
        character(len=6) :: s  !! the integer, left-justified
        s = itoa(i)
        end function lj

    end program cv_advdiff_bnd
!*****************************************************************************************

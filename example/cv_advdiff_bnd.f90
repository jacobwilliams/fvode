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
    real(wp),parameter,public :: YMAX = 1.0_wp
    integer,parameter,public  :: MX = 10        !! mesh dimensions
    integer,parameter,public  :: MY = 5
    integer,parameter,public  :: NEQ = MX * MY  !! number of equations

    type,extends(cvode_solver),public :: advdiff_solver
        !! the CVODE solver plus the problem data (C: UserData)
        real(wp) :: dx = 0.0_wp, dy = 0.0_wp
        real(wp) :: hdcoef = 0.0_wp, hacoef = 0.0_wp, vdcoef = 0.0_wp
    end type advdiff_solver

    public :: f, jac, ijth

    contains

    pure integer function ijth(i, j)
    !! index of u(i,j) in the solution vector (C: IJth, 1-based i, j)
    integer,intent(in) :: i, j
    ijth = j + (i - 1) * MY
    end function ijth

    function f(me, t, y, ydot) result(retval)
    !! right-hand side: central differences, zero boundary values
    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    real(wp),dimension(:),intent(in) :: y
    real(wp),dimension(:),intent(out) :: ydot
    integer :: retval
    real(wp) :: uij, udn, uup, ult, urt, hdiff, hadv, vdiff
    integer :: i, j
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
    class(cvode_solver),intent(inout) :: me
    real(wp),intent(in) :: t
    real(wp),dimension(:),intent(in) :: y, fy
    class(sun_matrix),intent(inout) :: J
    real(wp),dimension(:),intent(inout) :: tmp1, tmp2, tmp3
    integer :: retval
    integer :: i, jj
    integer(ip) :: k
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

    type(advdiff_solver) :: cvode
    real(wp),dimension(NEQ) :: u
    real(wp) :: reltol, abstol, t, tout, x, y
    integer :: iout, retval, i, j

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
        integer(li),intent(in) :: i
        character(len=:),allocatable :: s
        character(len=24) :: buf
        write(buf,'(i0)') i
        s = trim(buf)
        end function itoa

        function lj(i) result(s)
        !! left-justified in 6 characters (C: %-6ld)
        integer(li),intent(in) :: i
        character(len=6) :: s
        s = itoa(i)
        end function lj

    end program cv_advdiff_bnd
!*****************************************************************************************

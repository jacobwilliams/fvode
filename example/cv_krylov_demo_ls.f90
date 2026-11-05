!*****************************************************************************************
!>
!  2-species diurnal kinetics advection-diffusion PDE in 2-D, solved four
!  times with BDF + Newton and the Krylov solvers SPGMR, SPFGMR, SPBCGS and
!  SPTFQMR, each with the block-diagonal part of the Newton matrix as a
!  left preconditioner (user-supplied preconditioner setup and solve). A
!  copy of the block-diagonal part of the Jacobian is saved and
!  conditionally reused in the setup.
!
!  FVODE port of the SUNDIALS example `examples/cvode/serial/cvKrylovDemo_ls.c`
!  (same problem, settings and output, except the workspace sizes; compare
!  with `cvKrylovDemo_ls.out`). The optional first command-line argument
!  selects the norm factor as in C (0: default sqrt(N), 1: sqrt(N) given
!  explicitly, 2: computed from a dot product); the residual monitoring of
!  the C example (second argument) is not ported.
!
!  The system is
!
!      dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
!                 + Ri(c1,c2,t)      for i = 1,2,   where
!      R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
!      R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
!      Kv(y) = Kv0*exp(y/5) ,
!
!  on 0 <= x <= 20, 30 <= y <= 50 (km) with homogeneous Neumann boundary
!  conditions, for 0 <= t <= 86400 s (1 day), on a uniform 10 x 10 mesh.

    module cv_krylov_demo_ls_module

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_linsol_direct, only: dense_getrf, dense_getrs

    implicit none

    private

    integer,parameter,public :: NUM_SPECIES = 2  !! number of species
    real(wp),parameter :: KH  = 4.0e-6_wp   !! horizontal diffusivity Kh
    real(wp),parameter :: VEL = 0.001_wp    !! advection velocity V
    real(wp),parameter :: KV0 = 1.0e-8_wp   !! coefficient in Kv(y)
    real(wp),parameter :: Q1  = 1.63e-16_wp !! coefficient q1
    real(wp),parameter :: Q2  = 4.66e-16_wp !! coefficient q2
    real(wp),parameter :: C3  = 3.7e16_wp   !! coefficient c3
    real(wp),parameter :: A3  = 22.62_wp    !! coefficient in the expression for q3(t)
    real(wp),parameter :: A4  = 7.601_wp    !! coefficient in the expression for q4(t)
    real(wp),parameter :: C1_SCALE = 1.0e6_wp   !! scale of the initial profile of c1
    real(wp),parameter :: C2_SCALE = 1.0e12_wp  !! scale of the initial profile of c2
    real(wp),parameter :: HALFDAY = 4.32e4_wp   !! number of seconds in a half day
    real(wp),parameter :: PI   = 3.1415926535898_wp  !! pi (as in the C example)
    real(wp),parameter :: XMIN = 0.0_wp   !! lower grid boundary in x
    real(wp),parameter :: XMAX = 20.0_wp  !! upper grid boundary in x
    real(wp),parameter :: YMIN = 30.0_wp  !! lower grid boundary in y
    real(wp),parameter :: YMAX = 50.0_wp  !! upper grid boundary in y
    real(wp),parameter :: XMID = 10.0_wp  !! grid midpoint in x
    real(wp),parameter :: YMID = 40.0_wp  !! grid midpoint in y
    integer,parameter,public :: MX = 10  !! number of x mesh points
    integer,parameter,public :: MY = 10  !! number of y mesh points
    integer,parameter,public :: NSMX = NUM_SPECIES * MX  !! number of unknowns in a mesh row
    integer,parameter,public :: MM = MX * MY  !! number of mesh points
    integer,parameter,public :: NEQ = NUM_SPECIES * MM  !! number of equations

    type,extends(cvode_solver),public :: krylov_demo_solver
        !! the CVODE solver plus the problem data (C: UserData)
        real(wp),dimension(NUM_SPECIES,NUM_SPECIES,0:MX-1,0:MY-1) :: P = 0.0_wp  !! preconditioner blocks (LU factors of I - gamma*Jbd)
        real(wp),dimension(NUM_SPECIES,NUM_SPECIES,0:MX-1,0:MY-1) :: Jbd = 0.0_wp  !! saved diagonal Jacobian blocks
        integer(ip),dimension(NUM_SPECIES,0:MX-1,0:MY-1) :: pivot = 0  !! pivots of the LU factors of the blocks
        real(wp) :: q4 = 0.0_wp  !! current value of q4(t)
        real(wp) :: om = 0.0_wp  !! diurnal frequency, pi/halfday
        real(wp) :: dx = 0.0_wp  !! mesh spacing in x
        real(wp) :: dy = 0.0_wp  !! mesh spacing in y
        real(wp) :: hdco = 0.0_wp  !! horizontal diffusion coefficient Kh/dx**2
        real(wp) :: haco = 0.0_wp  !! horizontal advection coefficient V/(2 dx)
        real(wp) :: vdco = 0.0_wp  !! vertical diffusion coefficient Kv0/dy**2
        integer :: linsolver = 0  !! the linear solver in use (0: SPGMR, 1: SPFGMR, 2: SPBCGS, 3: SPTFQMR)
    end type krylov_demo_solver

    public :: f, precond, psolve, ijkth, init_user_data, set_initial_profiles

    contains

    pure integer function ijkth(i, jx, jy)
    !! index of species i at mesh point (jx,jy) in the solution vector
    !! (C: `IJKth`, with 1 <= i <= 2, 0 <= jx <= MX-1, 0 <= jy <= MY-1)
    integer,intent(in) :: i  !! species (1-based)
    integer,intent(in) :: jx  !! mesh index in x (0-based)
    integer,intent(in) :: jy  !! mesh index in y (0-based)
    ijkth = i + jx * NUM_SPECIES + jy * NSMX
    end function ijkth

    subroutine init_user_data(me)
    !! load the problem constants (C: `InitUserData`)
    type(krylov_demo_solver),intent(inout) :: me  !! the solver with the problem data
    me%om   = PI / HALFDAY
    me%dx   = (XMAX - XMIN) / (MX - 1)
    me%dy   = (YMAX - YMIN) / (MY - 1)
    me%hdco = KH / me%dx**2
    me%haco = VEL / (2.0_wp * me%dx)
    me%vdco = (1.0_wp / me%dy**2) * KV0
    end subroutine init_user_data

    subroutine set_initial_profiles(u, dx, dy)
    !! polynomial initial profiles of c1 and c2 (C: `SetInitialProfiles`)
    real(wp),dimension(:),intent(out) :: u  !! the solution vector
    real(wp),intent(in) :: dx  !! mesh spacing in x
    real(wp),intent(in) :: dy  !! mesh spacing in y
    integer :: jx, jy  !! mesh indices
    real(wp) :: x, y, cx, cy  !! mesh coordinates and profile factors in x and y
    do jy = 0, MY - 1
        y  = YMIN + jy * dy
        cy = (0.1_wp * (y - YMID))**2
        cy = 1.0_wp - cy + 0.5_wp * cy**2
        do jx = 0, MX - 1
            x  = XMIN + jx * dx
            cx = (0.1_wp * (x - XMID))**2
            cx = 1.0_wp - cx + 0.5_wp * cx**2
            u(ijkth(1, jx, jy)) = C1_SCALE * cx * cy
            u(ijkth(2, jx, jy)) = C2_SCALE * cx * cy
        end do
    end do
    end subroutine set_initial_profiles

    function f(me, t, u, udot) result(retval)
    !! right-hand side f(t,u)
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: u  !! current solution
    real(wp),dimension(:),intent(out) :: udot  !! du/dt = f(t,u)
    integer :: retval  !! 0 on success
    real(wp) :: q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt  !! q3(t), the concentrations at the point and its lower and upper neighbors in y, left neighbor in x
    real(wp) :: c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2  !! right neighbor in x, vertical diffusion coefficients, horizontal diffusion and advection terms
    real(wp) :: qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup  !! kinetic terms, sin(om*t), vertical diffusion terms and the y values at jy -+ 1/2
    real(wp) :: q4coef, dely, verdco, hordco, horaco  !! local copies of the problem data
    integer :: jx, jy, idn, iup, ileft, iright  !! mesh indices and the neighbor offsets
    select type (me)
    class is (krylov_demo_solver)
        ! set diurnal rate coefficients
        s = sin(me%om * t)
        if (s > 0.0_wp) then
            q3    = exp(-A3 / s)
            me%q4 = exp(-A4 / s)
        else
            q3    = 0.0_wp
            me%q4 = 0.0_wp
        end if

        ! make local copies of problem variables, for efficiency
        q4coef = me%q4
        dely   = me%dy
        verdco = me%vdco
        hordco = me%hdco
        horaco = me%haco

        ! loop over all grid points
        do jy = 0, MY - 1
            ! set vertical diffusion coefficients at jy +- 1/2
            ydn  = YMIN + (jy - 0.5_wp) * dely
            yup  = ydn + dely
            cydn = verdco * exp(0.2_wp * ydn)
            cyup = verdco * exp(0.2_wp * yup)
            idn = -1; if (jy == 0) idn = 1
            iup = 1;  if (jy == MY - 1) iup = -1
            do jx = 0, MX - 1
                ! extract c1 and c2, and set kinetic rate terms
                c1    = u(ijkth(1, jx, jy))
                c2    = u(ijkth(2, jx, jy))
                qq1   = Q1 * c1 * C3
                qq2   = Q2 * c1 * c2
                qq3   = q3 * C3
                qq4   = q4coef * c2
                rkin1 = -qq1 - qq2 + 2.0_wp * qq3 + qq4
                rkin2 = qq1 - qq2 - qq4

                ! set vertical diffusion terms
                c1dn   = u(ijkth(1, jx, jy + idn))
                c2dn   = u(ijkth(2, jx, jy + idn))
                c1up   = u(ijkth(1, jx, jy + iup))
                c2up   = u(ijkth(2, jx, jy + iup))
                vertd1 = cyup * (c1up - c1) - cydn * (c1 - c1dn)
                vertd2 = cyup * (c2up - c2) - cydn * (c2 - c2dn)

                ! set horizontal diffusion and advection terms
                ileft  = -1; if (jx == 0) ileft = 1
                iright = 1;  if (jx == MX - 1) iright = -1
                c1lt   = u(ijkth(1, jx + ileft, jy))
                c2lt   = u(ijkth(2, jx + ileft, jy))
                c1rt   = u(ijkth(1, jx + iright, jy))
                c2rt   = u(ijkth(2, jx + iright, jy))
                hord1  = hordco * (c1rt - 2.0_wp * c1 + c1lt)
                hord2  = hordco * (c2rt - 2.0_wp * c2 + c2lt)
                horad1 = horaco * (c1rt - c1lt)
                horad2 = horaco * (c2rt - c2lt)

                ! load all terms into udot
                udot(ijkth(1, jx, jy)) = vertd1 + hord1 + horad1 + rkin1
                udot(ijkth(2, jx, jy)) = vertd2 + hord2 + horad2 + rkin2
            end do
        end do
    end select
    retval = 0
    end function f

    function precond(me, t, y, fy, jok, jcur, gamma) result(retval)
    !! preconditioner setup: the 2x2 diagonal blocks of `I - gamma*J`
    !! (Jacobian blocks recomputed unless `jok`), LU-factored in place
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    logical,intent(in) :: jok  !! may the saved Jacobian blocks be reused?
    logical,intent(out) :: jcur  !! were the Jacobian blocks recomputed?
    real(wp),intent(in) :: gamma  !! the scalar gamma in P = I - gamma*J
    integer :: retval  !! 0 on success, 1 if a block is singular
    real(wp) :: c1, c2, cydn, cyup, diag, ydn, yup, q4coef, dely, verdco, hordco  !! concentrations, vertical diffusion coefficients, the diagonal diffusion term, the y values at jy -+ 1/2, and local copies of the problem data
    integer :: jx, jy, i  !! mesh indices and the diagonal index
    retval = 0
    select type (me)
    class is (krylov_demo_solver)
        if (jok) then
            ! jok = true: copy Jbd to P
            me%P = me%Jbd
            jcur = .false.
        else
            ! jok = false: generate Jbd from scratch and copy to P

            ! make local copies of problem variables, for efficiency
            q4coef = me%q4
            dely   = me%dy
            verdco = me%vdco
            hordco = me%hdco

            ! compute 2x2 diagonal Jacobian blocks (using q4 values
            ! computed on the last f call); load into P
            do jy = 0, MY - 1
                ydn  = YMIN + (jy - 0.5_wp) * dely
                yup  = ydn + dely
                cydn = verdco * exp(0.2_wp * ydn)
                cyup = verdco * exp(0.2_wp * yup)
                diag = -(cydn + cyup + 2.0_wp * hordco)
                do jx = 0, MX - 1
                    c1 = y(ijkth(1, jx, jy))
                    c2 = y(ijkth(2, jx, jy))
                    me%Jbd(1,1,jx,jy) = (-Q1 * C3 - Q2 * c2) + diag
                    me%Jbd(1,2,jx,jy) = -Q2 * c1 + q4coef
                    me%Jbd(2,1,jx,jy) = Q1 * C3 - Q2 * c2
                    me%Jbd(2,2,jx,jy) = (-Q2 * c1 - q4coef) + diag
                    me%P(:,:,jx,jy) = me%Jbd(:,:,jx,jy)
                end do
            end do
            jcur = .true.
        end if

        ! scale by -gamma
        me%P = -gamma * me%P

        ! add identity matrix and do LU decompositions on blocks in place
        do jx = 0, MX - 1
            do jy = 0, MY - 1
                do i = 1, NUM_SPECIES
                    me%P(i,i,jx,jy) = me%P(i,i,jx,jy) + 1.0_wp
                end do
                if (dense_getrf(me%P(:,:,jx,jy), me%pivot(:,jx,jy)) /= 0) then
                    retval = 1
                    return
                end if
            end do
        end do
    end select
    end function precond

    function psolve(me, t, y, fy, r, z, gamma, delta, lr) result(retval)
    !! preconditioner solve: the block-diagonal system P z = r, with the LU
    !! factors from [[precond]]
    class(cvode_solver),intent(inout) :: me  !! the CVODE solver (gives access to problem data)
    real(wp),intent(in) :: t  !! current time
    real(wp),dimension(:),intent(in) :: y  !! current solution
    real(wp),dimension(:),intent(in) :: fy  !! f(t, y)
    real(wp),dimension(:),intent(in) :: r  !! the right-hand side
    real(wp),dimension(:),intent(out) :: z  !! the solution
    real(wp),intent(in) :: gamma  !! the scalar gamma in P = I - gamma*J
    real(wp),intent(in) :: delta  !! tolerance (not used: direct solve)
    integer,intent(in) :: lr  !! left or right preconditioner (only left is used)
    integer :: retval  !! 0 on success
    integer :: jx, jy  !! mesh indices
    select type (me)
    class is (krylov_demo_solver)
        z = r
        ! solve the block-diagonal system Px = r using LU factors stored
        ! in P and pivot data in pivot, and return the solution in z
        do jx = 0, MX - 1
            do jy = 0, MY - 1
                call dense_getrs(me%P(:,:,jx,jy), me%pivot(:,jx,jy), &
                                 z(ijkth(1, jx, jy):ijkth(NUM_SPECIES, jx, jy)))
            end do
        end do
    end select
    retval = 0
    end function psolve

    end module cv_krylov_demo_ls_module
!*****************************************************************************************

!*****************************************************************************************
    program cv_krylov_demo_ls

    use fvode
    use fvode_linsol, only: sun_linear_solver
    use fvode_linsol_spgmr, only: sun_linsol_spgmr
    use fvode_linsol_spfgmr, only: sun_linsol_spfgmr
    use fvode_linsol_spbcgs, only: sun_linsol_spbcgs
    use fvode_linsol_sptfqmr, only: sun_linsol_sptfqmr
    use cv_krylov_demo_ls_module

    implicit none

    real(wp),parameter :: RTOL  = 1.0e-5_wp      !! scalar relative tolerance
    real(wp),parameter :: FLOOR = 100.0_wp       !! value of c1 or c2 at which tolerances change from relative to absolute
    real(wp),parameter :: ATOL  = RTOL * FLOOR   !! scalar absolute tolerance
    real(wp),parameter :: T0    = 0.0_wp         !! initial time
    real(wp),parameter :: TWOHR = 7200.0_wp      !! number of seconds in two hours
    integer,parameter  :: NOUT  = 12             !! number of output times
    integer,parameter  :: USE_SPGMR   = 0  !! linear solver: SPGMR
    integer,parameter  :: USE_SPFGMR  = 1  !! linear solver: SPFGMR
    integer,parameter  :: USE_SPBCG   = 2  !! linear solver: SPBCGS
    integer,parameter  :: USE_SPTFQMR = 3  !! linear solver: SPTFQMR

    type(krylov_demo_solver) :: cvode  !! the CVODE solver with the problem data
    real(wp),dimension(NEQ) :: u  !! solution vector
    real(wp) :: t, tout  !! current and output times
    real(wp) :: nrmfac  !! LS norm conversion factor
    integer :: retval, iout, linsolver  !! return flag, output counter and linear solver
    integer :: nrmfactor  !! LS norm conversion factor flag (command-line argument)
    character(len=32) :: arg  !! command-line argument

    ! retrieve the command-line options
    nrmfactor = 0
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read(arg, *) nrmfactor
    end if

    call init_user_data(cvode)
    call set_initial_profiles(u, cvode%dx, cvode%dy)

    call cvode%create(CV_BDF, retval)
    call cvode%init(f, T0, u, retval)
    call cvode%ss_tolerances(RTOL, ATOL, retval)
    ! (CVodeInit attaches a Newton nonlinear solver, as SUNNonlinSol_Newton
    ! + CVodeSetNonlinearSolver do in the C example)

    do linsolver = USE_SPGMR, USE_SPTFQMR

        if (linsolver /= 0) then
            call init_user_data(cvode)
            call set_initial_profiles(u, cvode%dx, cvode%dy)
            ! re-initialize CVode for the solution of the same problem, but
            ! using a different linear solver module
            call cvode%reinit(T0, u, retval)
            if (check_retval(retval, 'CVodeReInit')) stop 1
        end if

        cvode%linsolver = linsolver

        select case (linsolver)
        case (USE_SPGMR)
            write(*,'(a)') ' ------- ', '| SPGMR |', ' -------'
            call attach(sun_linsol_spgmr(NEQ, SUN_PREC_LEFT, 0))
        case (USE_SPFGMR)
            write(*,'(a)') ' --------- ', '| SPFGMR |', ' ---------'
            call attach(sun_linsol_spfgmr(NEQ, SUN_PREC_LEFT, 0))
        case (USE_SPBCG)
            write(*,'(a)') ' ------- ', '| SPBCGS |', ' -------'
            call attach(sun_linsol_spbcgs(NEQ, SUN_PREC_LEFT, 0))
        case (USE_SPTFQMR)
            write(*,'(a)') ' --------- ', '| SPTFQMR |', ' ---------'
            call attach(sun_linsol_sptfqmr(NEQ, SUN_PREC_LEFT, 0))
        end select

        ! set preconditioner setup and solve routines precond and psolve
        call cvode%set_preconditioner(precond, psolve, retval)
        if (check_retval(retval, 'CVodeSetPreconditioner')) stop 1

        select case (nrmfactor)
        case (1)
            nrmfac = sqrt(real(NEQ, wp))
        case (2)
            nrmfac = -1.0_wp
        case default
            nrmfac = 0.0_wp
        end select
        call cvode%set_ls_norm_factor(nrmfac, retval)
        if (check_retval(retval, 'CVodeSetLSNormFactor')) stop 1

        write(*,'(a/)') ' ', '2-species diurnal advection-diffusion problem'

        ! in loop over output points, call CVode, print results, test for error
        tout = TWOHR
        do iout = 1, NOUT
            call cvode%solve(tout, u, t, CV_NORMAL, retval)
            call print_output(t)
            if (check_retval(retval, 'CVode')) exit
            tout = tout + TWOHR
        end do

        call print_stats(linsolver)

    end do

    call cvode%destroy()

    contains

        subroutine attach(LS)
        !! attach the linear solver (C: `CVodeSetLinearSolver(cvode_mem, LS, NULL)`)
        class(sun_linear_solver),intent(in) :: LS  !! the linear solver
        call cvode%set_linear_solver(LS, flag=retval)
        if (check_retval(retval, 'CVodeSetLinearSolver')) stop 1
        end subroutine attach

        subroutine print_output(t)
        !! print current t, step count, order, step size, and sampled c1,c2 values
        real(wp),intent(in) :: t  !! current time
        integer,parameter :: mxh = MX / 2 - 1  !! x index of the middle point
        integer,parameter :: myh = MY / 2 - 1  !! y index of the middle point
        integer,parameter :: mx1 = MX - 1  !! x index of the top right point
        integer,parameter :: my1 = MY - 1  !! y index of the top right point
        write(*,'(a,es8.2,a,i0,a,i0,a,es8.2)') 't = ', t, '   no. steps = ', cvode%get_num_steps(), &
            '   order = ', cvode%get_last_order(), '   stepsize = ', cvode%get_last_step()
        write(*,'(a,es12.3,2x,es12.3,2x,es12.3)') 'c1 (bot.left/middle/top rt.) = ', &
            u(ijkth(1, 0, 0)), u(ijkth(1, mxh, myh)), u(ijkth(1, mx1, my1))
        write(*,'(a,es12.3,2x,es12.3,2x,es12.3/)') 'c2 (bot.left/middle/top rt.) = ', &
            u(ijkth(2, 0, 0)), u(ijkth(2, mxh, myh)), u(ijkth(2, mx1, my1))
        end subroutine print_output

        subroutine print_stats(linsolver)
        !! print the final statistics (except the workspace sizes)
        integer,intent(in) :: linsolver  !! the linear solver in use
        write(*,'(/a/)') 'Final Statistics.. '
        write(*,'(a,i5)') 'nst     = ', cvode%get_num_steps()
        write(*,'(a,i5,a,i5)') 'nfe     = ', cvode%get_num_rhs_evals(), &
            '     nfeLS   = ', cvode%get_num_lin_rhs_evals()
        write(*,'(a,i5,a,i5)') 'nni     = ', cvode%get_num_nonlin_solv_iters(), &
            '     nli     = ', cvode%get_num_lin_iters()
        write(*,'(a,i5,a,i5)') 'nsetups = ', cvode%get_num_lin_solv_setups(), &
            '     netf    = ', cvode%get_num_err_test_fails()
        write(*,'(a,i5,a,i5)') 'npe     = ', cvode%get_num_prec_evals(), &
            '     nps     = ', cvode%get_num_prec_solves()
        write(*,'(a,i5,a,i5/)') 'ncfn    = ', cvode%get_num_nonlin_solv_conv_fails(), &
            '     ncfl    = ', cvode%get_num_lin_conv_fails()
        if (linsolver < 2) write(*,'(a/)') repeat('=', 70)
        end subroutine print_stats

        logical function check_retval(retval, funcname)
        !! report a negative return flag (C: `check_retval(&retval, funcname, 1)`)
        integer,intent(in) :: retval  !! the return flag
        character(len=*),intent(in) :: funcname  !! name of the called function
        check_retval = retval < 0
        if (check_retval) write(*,'(/3a,i0/)') 'SUNDIALS_ERROR: ', funcname, '() failed with retval = ', retval
        end function check_retval

    end program cv_krylov_demo_ls
!*****************************************************************************************

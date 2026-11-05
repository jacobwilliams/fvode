!*****************************************************************************************
!>
!  Test problems with known solutions for the solver comparison
!  ([[test_solver_comparison]]): the right-hand sides, initial
!  conditions, final times and exact solutions, and a counter of the
!  right-hand side evaluations shared by all solvers.

    module comparison_problems

    use fvode_kinds, only: wp

    implicit none

    private

    integer,parameter,public :: PROB_KEPLER     = 1  !! elliptical Kepler orbit (non-stiff)
    integer,parameter,public :: PROB_OSCILLATOR = 2  !! harmonic oscillator over many periods (non-stiff)
    integer,parameter,public :: PROB_ANALYTIC   = 3  !! y' = lambda*(y - atan t) + 1/(1+t^2), lambda = -100 (stiff)
    integer,parameter,public :: PROB_STIFF_LINEAR = 4  !! rotated linear system with eigenvalues -1 and -1000 (stiff)
    integer,parameter,public :: NUM_PROBLEMS    = 4  !! number of problems

    real(wp),parameter :: PI = acos(-1.0_wp)  !! pi
    real(wp),parameter :: ECC = 0.7_wp  !! eccentricity of the Kepler orbit (mu = 1, a = 1)
    real(wp),parameter :: LAMBDA = -100.0_wp  !! stiffness parameter of the analytic problem
    real(wp),dimension(2),parameter :: EIGS = [-1.0_wp, -1000.0_wp]  !! eigenvalues of the stiff linear system

    integer,public :: iprob = PROB_KEPLER  !! the current problem
    integer,public :: nfe = 0  !! right-hand side evaluations (reset before each solve)

    public :: problem_name, problem_neq, problem_stiff, problem_tf, problem_y0, rhs, exact

    contains
!*****************************************************************************************

    function problem_name(ip) result(name)
    !! short name of a problem (used in file names)
    integer,intent(in) :: ip  !! the problem
    character(len=:),allocatable :: name  !! the name
    select case (ip)
    case (PROB_KEPLER);       name = 'kepler'
    case (PROB_OSCILLATOR);   name = 'oscillator'
    case (PROB_ANALYTIC);     name = 'analytic_stiff'
    case (PROB_STIFF_LINEAR); name = 'stiff_linear'
    end select
    end function problem_name

    pure integer function problem_neq(ip)
    !! number of equations
    integer,intent(in) :: ip  !! the problem
    select case (ip)
    case (PROB_KEPLER);   problem_neq = 4
    case (PROB_ANALYTIC); problem_neq = 1
    case default;         problem_neq = 2
    end select
    end function problem_neq

    pure logical function problem_stiff(ip)
    !! is the problem stiff (solved with BDF + Newton by fvode and dvode)?
    integer,intent(in) :: ip  !! the problem
    problem_stiff = (ip == PROB_ANALYTIC .or. ip == PROB_STIFF_LINEAR)
    end function problem_stiff

    pure real(wp) function problem_tf(ip)
    !! final time (the initial time is 0)
    integer,intent(in) :: ip  !! the problem
    select case (ip)
    case (PROB_KEPLER);       problem_tf = 5.3_wp * 2.0_wp * PI  ! 5.3 orbits (period 2 pi)
    case (PROB_OSCILLATOR);   problem_tf = 100.0_wp
    case (PROB_ANALYTIC);     problem_tf = 10.0_wp
    case (PROB_STIFF_LINEAR); problem_tf = 1.0_wp
    end select
    end function problem_tf

    function problem_y0(ip) result(y0)
    !! initial condition (the exact solution at t = 0)
    integer,intent(in) :: ip  !! the problem
    real(wp),dimension(:),allocatable :: y0  !! the initial condition
    allocate(y0(problem_neq(ip)))
    call exact(ip, 0.0_wp, y0)
    end function problem_y0

    subroutine rhs(t, y, ydot)
    !! the right-hand side of the current problem (counts the evaluations)
    real(wp),intent(in) :: t  !! time
    real(wp),dimension(:),intent(in) :: y  !! state
    real(wp),dimension(:),intent(out) :: ydot  !! dy/dt
    real(wp) :: r3  !! |r|^3
    nfe = nfe + 1
    select case (iprob)
    case (PROB_KEPLER)
        r3 = norm2(y(1:2))**3
        ydot(1:2) = y(3:4)
        ydot(3:4) = -y(1:2) / r3
    case (PROB_OSCILLATOR)
        ydot(1) = y(2)
        ydot(2) = -y(1)
    case (PROB_ANALYTIC)
        ydot(1) = LAMBDA * (y(1) - atan(t)) + 1.0_wp / (1.0_wp + t * t)
    case (PROB_STIFF_LINEAR)
        ! A = Q diag(EIGS) Q^T, Q = rotation by 45 degrees
        ydot(1) = 0.5_wp * (EIGS(1) + EIGS(2)) * y(1) + 0.5_wp * (EIGS(1) - EIGS(2)) * y(2)
        ydot(2) = 0.5_wp * (EIGS(1) - EIGS(2)) * y(1) + 0.5_wp * (EIGS(1) + EIGS(2)) * y(2)
    end select
    end subroutine rhs

    subroutine exact(ip, t, y)
    !! the exact solution of a problem at time t
    integer,intent(in) :: ip  !! the problem
    real(wp),intent(in) :: t  !! time
    real(wp),dimension(:),intent(out) :: y  !! the exact state
    real(wp) :: E  !! eccentric anomaly
    real(wp) :: w1, w2  !! components in the eigenvector basis
    select case (ip)
    case (PROB_KEPLER)
        ! start at periapsis; mean motion n = 1, so M = t
        E = kepler_E(t)
        y(1) = cos(E) - ECC
        y(2) = sqrt(1.0_wp - ECC**2) * sin(E)
        y(3) = -sin(E) / (1.0_wp - ECC * cos(E))
        y(4) = sqrt(1.0_wp - ECC**2) * cos(E) / (1.0_wp - ECC * cos(E))
    case (PROB_OSCILLATOR)
        y(1) = cos(t)
        y(2) = -sin(t)
    case (PROB_ANALYTIC)
        y(1) = atan(t)
    case (PROB_STIFF_LINEAR)
        ! y(0) = [1, 0]: components along the eigenvectors (1,1)/sqrt2 and (1,-1)/sqrt2
        w1 = sqrt(0.5_wp) * exp(EIGS(1) * t)
        w2 = sqrt(0.5_wp) * exp(EIGS(2) * t)
        y(1) = sqrt(0.5_wp) * (w1 + w2)
        y(2) = sqrt(0.5_wp) * (w1 - w2)
    end select
    end subroutine exact

    pure real(wp) function kepler_E(M)
    !! solve Kepler's equation `E - e sin E = M` (Newton's method)
    real(wp),intent(in) :: M  !! mean anomaly
    real(wp) :: dE  !! Newton correction
    integer :: iter  !! iteration counter
    kepler_E = M + ECC * sin(M)
    do iter = 1, 100
        dE = (kepler_E - ECC * sin(kepler_E) - M) / (1.0_wp - ECC * cos(kepler_E))
        kepler_E = kepler_E - dE
        if (abs(dE) <= 4.0_wp * epsilon(1.0_wp) * max(1.0_wp, abs(kepler_E))) exit
    end do
    end function kepler_E

    end module comparison_problems
!*****************************************************************************************

!*****************************************************************************************
!>
!  Run one problem with FVODE, DVODE or DDEABM to its final time at a
!  given tolerance (rtol = atol = tol), returning the solution, the number
!  of right-hand side evaluations and the wall-clock time per solve.

    module comparison_runners

    use fvode
    use fvode_core, only: cvode_solver
    use fvode_matrix, only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense
    use fvode_nonlinsol_fixedpoint, only: sun_nonlinsol_fixedpoint
    use dvode_module, only: dvode_t
    use dvode_kinds_module, only: dvode_wp
    use ddeabm_module, only: ddeabm_class, ddeabm_rk
    use comparison_problems

    implicit none

    private

    integer,parameter,public :: SOLVER_FVODE  = 1  !! FVODE (this library)
    integer,parameter,public :: SOLVER_DVODE  = 2  !! DVODE (Modern Fortran edition)
    integer,parameter,public :: SOLVER_DDEABM = 3  !! DDEABM (Adams-Bashforth-Moulton)
    integer,parameter,public :: NUM_SOLVERS   = 3  !! number of solvers

    public :: solver_name, run_solver

    contains
!*****************************************************************************************

    function solver_name(is) result(name)
    !! name of a solver
    integer,intent(in) :: is  !! the solver
    character(len=:),allocatable :: name  !! the name
    select case (is)
    case (SOLVER_FVODE);  name = 'FVODE'
    case (SOLVER_DVODE);  name = 'DVODE'
    case (SOLVER_DDEABM); name = 'DDEABM'
    end select
    end function solver_name

    subroutine run_solver(is, iprb, tol, nrep, y, nevals, seconds, ok)
    !! solve problem iprb with solver is (nrep times, for the timing)
    integer,intent(in) :: is  !! the solver
    integer,intent(in) :: iprb  !! the problem
    real(wp),intent(in) :: tol  !! relative and absolute tolerance
    integer,intent(in) :: nrep  !! number of repetitions (timing)
    real(wp),dimension(:),intent(out) :: y  !! solution at the final time
    integer,intent(out) :: nevals  !! right-hand side evaluations of one solve
    real(wp),intent(out) :: seconds  !! wall-clock time of one solve
    logical,intent(out) :: ok  !! did the solver reach the final time?
    integer(8) :: c0, c1, rate  !! clock counts and rate
    integer :: irep  !! repetition counter
    iprob = iprb
    call system_clock(c0, rate)
    do irep = 1, nrep
        nfe = 0
        select case (is)
        case (SOLVER_FVODE);  call run_fvode(iprb, tol, y, ok)
        case (SOLVER_DVODE);  call run_dvode(iprb, tol, y, ok)
        case (SOLVER_DDEABM); call run_ddeabm(iprb, tol, y, ok)
        end select
        if (.not. ok) exit
    end do
    call system_clock(c1)
    nevals = nfe
    seconds = real(c1 - c0, wp) / real(rate, wp) / nrep
    end subroutine run_solver

    subroutine run_fvode(iprb, tol, y, ok)
    !! FVODE: Adams + fixed point (non-stiff) or BDF + Newton with the
    !! dense difference-quotient Jacobian (stiff), CV_NORMAL to tf
    integer,intent(in) :: iprb  !! the problem
    real(wp),intent(in) :: tol  !! tolerance
    real(wp),dimension(:),intent(out) :: y  !! solution at tf
    logical,intent(out) :: ok  !! success?
    type(cvode_solver) :: cv  !! the solver
    real(wp) :: t  !! return time
    integer :: n, flag  !! problem size and return flag
    n = problem_neq(iprb)
    y = problem_y0(iprb)
    if (problem_stiff(iprb)) then
        call cv%create(CV_BDF, flag)
    else
        call cv%create(CV_ADAMS, flag)
    end if
    call cv%init(fv_rhs, 0.0_wp, y, flag)
    call cv%ss_tolerances(tol, tol, flag)
    call cv%set_max_num_steps(10000000_li, flag)
    if (problem_stiff(iprb)) then
        call cv%set_linear_solver(sun_linsol_dense(int(n, ip)), dense_matrix(int(n, ip), int(n, ip)), flag)
    else
        call cv%set_nonlinear_solver(sun_nonlinsol_fixedpoint(n, 0), flag)
    end if
    call cv%solve(problem_tf(iprb), y, t, CV_NORMAL, flag)
    ok = (flag >= 0)
    call cv%destroy()
    end subroutine run_fvode

    function fv_rhs(me, t, y, ydot) result(retval)
    !! FVODE right-hand side wrapper
    class(cvode_solver),intent(inout) :: me  !! the solver
    real(wp),intent(in) :: t  !! time
    real(wp),dimension(:),intent(in) :: y  !! state
    real(wp),dimension(:),intent(out) :: ydot  !! dy/dt
    integer :: retval  !! 0
    call rhs(t, y, ydot)
    retval = 0
    end function fv_rhs

    subroutine run_dvode(iprb, tol, y, ok)
    !! DVODE: mf = 10 (Adams, functional iteration) for non-stiff problems,
    !! mf = 22 (BDF, internal difference-quotient full Jacobian) for stiff
    !! ones, itask = 1 to tf
    integer,intent(in) :: iprb  !! the problem
    real(wp),intent(in) :: tol  !! tolerance
    real(wp),dimension(:),intent(out) :: y  !! solution at tf
    logical,intent(out) :: ok  !! success?
    type(dvode_t) :: s  !! the solver
    real(dvode_wp) :: t, tout  !! time and output time
    real(dvode_wp),dimension(1) :: rtol, atol  !! tolerances
    real(dvode_wp),dimension(:),allocatable :: rwork  !! real work array
    integer,dimension(:),allocatable :: iwork  !! integer work array
    integer :: n, mf, istate, lrw, liw  !! problem size, method flag, state, work array sizes
    n = problem_neq(iprb)
    y = problem_y0(iprb)
    mf = merge(22, 10, problem_stiff(iprb))
    lrw = 22 + 16 * n + 2 * n * n
    liw = 30 + n
    allocate(rwork(lrw), source=0.0_dvode_wp)
    allocate(iwork(liw), source=0)
    iwork(6) = 10000000  ! maximum number of steps
    rtol = tol
    atol = tol
    t = 0.0_dvode_wp
    tout = problem_tf(iprb)
    istate = 1
    call s%initialize(f=dv_rhs)
    call s%solve(n, y, t, tout, 1, rtol, atol, 1, istate, 1, rwork, lrw, iwork, liw, mf)
    ok = (istate == 2)
    end subroutine run_dvode

    subroutine dv_rhs(me, neq, t, y, ydot)
    !! DVODE right-hand side wrapper
    class(dvode_t),intent(inout) :: me  !! the solver
    integer :: neq  !! number of equations
    real(dvode_wp) :: t  !! time
    real(dvode_wp) :: y(neq)  !! state
    real(dvode_wp) :: ydot(neq)  !! dy/dt
    call rhs(t, y, ydot)
    end subroutine dv_rhs

    subroutine run_ddeabm(iprb, tol, y, ok)
    !! DDEABM: variable-order Adams-Bashforth-Moulton PECE, to tf
    integer,intent(in) :: iprb  !! the problem
    real(wp),intent(in) :: tol  !! tolerance
    real(wp),dimension(:),intent(out) :: y  !! solution at tf
    logical,intent(out) :: ok  !! success?
    type(ddeabm_class) :: s  !! the solver
    real(ddeabm_rk) :: t  !! time
    integer :: idid  !! return flag
    y = problem_y0(iprb)
    call s%initialize(problem_neq(iprb), maxnum=100000000, df=dd_rhs, rtol=[tol], atol=[tol])
    t = 0.0_ddeabm_rk
    call s%integrate(t, y, problem_tf(iprb), idid=idid)
    ok = (idid > 0)
    call s%destroy()
    end subroutine run_ddeabm

    subroutine dd_rhs(me, t, x, xdot)
    !! DDEABM right-hand side wrapper
    class(ddeabm_class),intent(inout) :: me  !! the solver
    real(ddeabm_rk),intent(in) :: t  !! time
    real(ddeabm_rk),dimension(:),intent(in) :: x  !! state
    real(ddeabm_rk),dimension(:),intent(out) :: xdot  !! dx/dt
    call rhs(t, x, xdot)
    end subroutine dd_rhs

    end module comparison_runners
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solver comparison: FVODE vs DVODE and DDEABM on problems with known
!  solutions (an elliptical Kepler orbit, a harmonic oscillator, and two
!  stiff problems). For each problem and solver, rtol = atol is swept from
!  1e-3 to 1e-12; the error at the final time (max norm, against the
!  exact solution), the number of right-hand side evaluations and the
!  time per solve are tabulated, and work-precision plots (function
!  evaluations vs error) are written to `test_output/` with
!  pyplot-fortran (if python with matplotlib is available).
!
!  Checks: every solver reaches the final time at every tolerance; FVODE's
!  error goes down with the tolerance and is small at tight tolerances;
!  and FVODE's cost is comparable to DVODE's (CVODE's predecessor).

    program test_solver_comparison

    use iso_fortran_env, only: error_unit
    use fvode_kinds, only: wp
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, new_unittest, &
                         unittest_type, error_type, check
    use pyplot_module, only: pyplot
    use comparison_problems
    use comparison_runners

    implicit none

    integer,parameter :: NTOL = 10  !! number of tolerances (1e-3 ... 1e-12)
    integer :: stat  !! number of failed tests
    type(testsuite_type),allocatable :: testsuites(:)  !! the test suite

    stat = 0
    testsuites = [ new_testsuite("solver_comparison", collect) ]
    write(error_unit, '("#", *(1x, a))') "Testing:", testsuites(1)%name
    call run_testsuite(testsuites(1)%collect, error_unit, stat)
    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

    contains

    subroutine collect(testsuite)
    !! collect the tests (one per problem)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)  !! the unit tests
    testsuite = [ new_unittest("kepler_orbit",   test_kepler),     &
                  new_unittest("oscillator",     test_oscillator), &
                  new_unittest("analytic_stiff", test_analytic),   &
                  new_unittest("stiff_linear",   test_stiff_linear) ]
    end subroutine collect

    subroutine test_kepler(error)
    !! elliptical orbit (e = 0.7) against the Kepler solution
    type(error_type),allocatable,intent(out) :: error  !! error information
    call compare_problem(PROB_KEPLER, 'Elliptical orbit (e = 0.7), 5.3 periods', error)
    end subroutine test_kepler

    subroutine test_oscillator(error)
    !! harmonic oscillator over about 16 periods
    type(error_type),allocatable,intent(out) :: error  !! error information
    call compare_problem(PROB_OSCILLATOR, 'Harmonic oscillator, t = 0 ... 100', error)
    end subroutine test_oscillator

    subroutine test_analytic(error)
    !! stiff scalar problem with solution atan(t)
    type(error_type),allocatable,intent(out) :: error  !! error information
    call compare_problem(PROB_ANALYTIC, 'Stiff analytic problem (lambda = -100)', error)
    end subroutine test_analytic

    subroutine test_stiff_linear(error)
    !! stiff linear system with eigenvalues -1 and -1000
    type(error_type),allocatable,intent(out) :: error  !! error information
    call compare_problem(PROB_STIFF_LINEAR, 'Stiff linear system (eigenvalues -1, -1000)', error)
    end subroutine test_stiff_linear

    subroutine compare_problem(ip, title, error)
    !! sweep the tolerances for all solvers, tabulate, plot and check
    integer,intent(in) :: ip  !! the problem
    character(len=*),intent(in) :: title  !! plot title
    type(error_type),allocatable,intent(out) :: error  !! error information
    real(wp),dimension(NTOL) :: tols  !! the tolerances
    real(wp),dimension(NTOL,NUM_SOLVERS) :: err  !! max-norm error at tf
    real(wp),dimension(NTOL,NUM_SOLVERS) :: secs  !! time per solve
    integer,dimension(NTOL,NUM_SOLVERS) :: nevals  !! right-hand side evaluations
    logical,dimension(NTOL,NUM_SOLVERS) :: ok  !! success
    real(wp),dimension(:),allocatable :: y, ytrue  !! computed and exact solution
    integer :: it, is, nrep  !! tolerance and solver indices, repetitions
    integer,parameter :: nrep_min = 3  !! minimum repetitions for the timing
    real(wp) :: ratio  !! FVODE / DVODE evaluations

    allocate(y(problem_neq(ip)), ytrue(problem_neq(ip)))
    call exact(ip, problem_tf(ip), ytrue)
    do it = 1, NTOL
        tols(it) = 10.0_wp**(-(it + 2))
    end do

    do is = 1, NUM_SOLVERS
        do it = 1, NTOL
            ! more repetitions for cheap runs, so the timing is meaningful
            nrep = nrep_min
            call run_solver(is, ip, tols(it), 1, y, nevals(it,is), secs(it,is), ok(it,is))
            if (ok(it,is) .and. nevals(it,is) < 20000) then
                nrep = max(nrep_min, 20000 / max(1, nevals(it,is)))
                call run_solver(is, ip, tols(it), nrep, y, nevals(it,is), secs(it,is), ok(it,is))
            end if
            err(it,is) = maxval(abs(y - ytrue))
        end do
    end do

    call print_table(ip, tols, err, nevals, secs, ok)
    call make_plot(ip, title, err, nevals, ok)

    ! every solver reaches tf at every tolerance
    do is = 1, NUM_SOLVERS
        call check(error, all(ok(:,is)), solver_name(is)//' failed on '//problem_name(ip))
        if (allocated(error)) return
    end do

    ! FVODE converges: the error decreases by orders of magnitude over the
    ! tolerance range and is small at the tight tolerances
    call check(error, err(NTOL,SOLVER_FVODE) < 1.0e-4_wp * err(1,SOLVER_FVODE), &
               'FVODE error does not decrease with the tolerance')
    if (allocated(error)) return
    call check(error, maxval(err(NTOL-2:NTOL,SOLVER_FVODE)) < 1.0e-6_wp, &
               'FVODE error too large at tight tolerances')
    if (allocated(error)) return

    ! FVODE's cost is comparable to DVODE's (same family of methods)
    do it = 1, NTOL
        ratio = real(nevals(it,SOLVER_FVODE), wp) / real(max(1, nevals(it,SOLVER_DVODE)), wp)
        call check(error, ratio < 3.0_wp, 'FVODE needs more than 3x the f evaluations of DVODE')
        if (allocated(error)) return
    end do

    end subroutine compare_problem

    subroutine print_table(ip, tols, err, nevals, secs, ok)
    !! print the results of one problem
    integer,intent(in) :: ip  !! the problem
    real(wp),dimension(NTOL),intent(in) :: tols  !! the tolerances
    real(wp),dimension(NTOL,NUM_SOLVERS),intent(in) :: err  !! errors
    integer,dimension(NTOL,NUM_SOLVERS),intent(in) :: nevals  !! evaluations
    real(wp),dimension(NTOL,NUM_SOLVERS),intent(in) :: secs  !! time per solve
    logical,dimension(NTOL,NUM_SOLVERS),intent(in) :: ok  !! success
    integer :: it, is  !! indices
    write(*,'(/a)') 'Problem: '//problem_name(ip)
    write(*,'(a7,3(3x,a8,a10,a11))') 'tol', ('solver', 'nfe', 'error', is = 1, NUM_SOLVERS)
    do it = 1, NTOL
        write(*,'(es7.0)',advance='no') tols(it)
        do is = 1, NUM_SOLVERS
            if (ok(it,is)) then
                write(*,'(3x,a8,i10,es11.2)',advance='no') solver_name(is), nevals(it,is), err(it,is)
            else
                write(*,'(3x,a8,a21)',advance='no') solver_name(is), 'FAILED'
            end if
        end do
        write(*,'(a)') ''
    end do
    write(*,'(a)') 'time per solve [microseconds]:'
    do it = 1, NTOL
        write(*,'(es7.0,3(3x,a8,f10.1))') tols(it), (solver_name(is), 1.0e6_wp * secs(it,is), is = 1, NUM_SOLVERS)
    end do
    end subroutine print_table

    subroutine make_plot(ip, title, err, nevals, ok)
    !! work-precision plot: function evaluations vs error (log-log)
    integer,intent(in) :: ip  !! the problem
    character(len=*),intent(in) :: title  !! plot title
    real(wp),dimension(NTOL,NUM_SOLVERS),intent(in) :: err  !! errors
    integer,dimension(NTOL,NUM_SOLVERS),intent(in) :: nevals  !! evaluations
    logical,dimension(NTOL,NUM_SOLVERS),intent(in) :: ok  !! success
    type(pyplot) :: plt  !! the plot
    integer :: is, istat  !! solver index and status
    character(len=*),parameter :: dir = 'test_output'  !! output directory
    character(len=*),dimension(NUM_SOLVERS),parameter :: styles = ['o--', 's- ', '^- ']  !! line styles
    integer,dimension(NUM_SOLVERS),parameter :: order = [SOLVER_DVODE, SOLVER_DDEABM, SOLVER_FVODE]  !! drawing order (FVODE on top: it often coincides with DVODE)
    integer :: k  !! drawing index
    logical,dimension(NTOL) :: mask  !! points to plot (successful, nonzero error)
    call execute_command_line('mkdir -p '//dir, exitstat=istat)
    call plt%initialize(grid=.true., xlabel='max error at final time', &
                        ylabel='number of function evaluations', title=title, &
                        legend=.true., figsize=[8,6], use_oo_api=.true., tight_layout=.true.)
    do k = 1, NUM_SOLVERS
        is = order(k)
        mask = ok(:,is) .and. err(:,is) > 0.0_wp
        call plt%add_plot(pack(err(:,is), mask), real(pack(nevals(:,is), mask), wp), &
                          label=solver_name(is), linestyle=trim(styles(is)), &
                          markersize=merge(9, 6, is == SOLVER_FVODE), &
                          linewidth=merge(1, 3, is == SOLVER_FVODE), &
                          xscale='log', yscale='log', istat=istat)
    end do
    call plt%savefig(dir//'/'//problem_name(ip)//'.png', &
                     pyfile=dir//'/'//problem_name(ip)//'.py', istat=istat)
    if (istat /= 0) then
        write(*,'(a)') 'note: the plot was not generated (python/matplotlib not available?)'
    else
        write(*,'(a)') 'wrote '//dir//'/'//problem_name(ip)//'.png'
    end if
    end subroutine make_plot

    end program test_solver_comparison
!*****************************************************************************************

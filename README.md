FVODE: Modern Fortran port of CVODE

**Status:** early development. See [PLAN.md](PLAN.md) for the conversion plan and progress.

FVODE is a native Modern Fortran translation of the CVODE ODE integrator from
[SUNDIALS](https://github.com/LLNL/sundials) v7.9.0. The goal is a faithful
translation: given the same problem and options, FVODE should take exactly the same
steps as CVODE. This is checked by comparing step-by-step traces against the C library.

### Building

Everything runs through [pixi](https://pixi.sh):

```
pixi run build      # build the library (fpm)
pixi run test       # run the unit tests
pixi run docs       # FORD documentation in ./doc
```

### Example

The Robertson chemical kinetics problem, a classic stiff test case, solved with
BDF, a Newton iteration and a dense direct linear solver. The Jacobian is
approximated by difference quotients.

```fortran
program robertson

    use fvode                                     ! kinds and constants
    use fvode_core,          only: cvode_solver
    use fvode_matrix,        only: dense_matrix
    use fvode_linsol_direct, only: sun_linsol_dense

    implicit none

    type(cvode_solver) :: cvode
    real(wp) :: y(3), t, tout
    integer  :: i, flag

    y = [1.0_wp, 0.0_wp, 0.0_wp]                  ! initial condition at t = 0

    call cvode%create(CV_BDF, flag)               ! BDF method for a stiff problem
    call cvode%init(rhs, 0.0_wp, y, flag)
    call cvode%sv_tolerances(1.0e-4_wp, [1.0e-8_wp, 1.0e-14_wp, 1.0e-6_wp], flag)
    call cvode%set_linear_solver(sun_linsol_dense(3_ip), dense_matrix(3_ip, 3_ip), flag)

    tout = 0.4_wp
    do i = 1, 6
        call cvode%solve(tout, y, t, CV_NORMAL, flag)
        if (flag < 0) error stop 'integration failed'
        write(*,'(a,es10.3,a,3es14.6)') 't = ', t, '   y =', y
        tout = tout * 10.0_wp
    end do

    write(*,'(a,i0,a,i0)') 'steps: ', cvode%get_num_steps(), &
                           '   rhs evals: ', cvode%get_num_rhs_evals()
    call cvode%destroy()

contains

    function rhs(me, t, y, ydot) result(retval)
        !! Robertson chemical kinetics: ydot = f(t, y)
        class(cvode_solver), intent(inout) :: me
        real(wp), intent(in)  :: t, y(:)
        real(wp), intent(out) :: ydot(:)
        integer :: retval
        ydot(1) = -0.04_wp*y(1) + 1.0e4_wp*y(2)*y(3)
        ydot(3) =  3.0e7_wp*y(2)**2
        ydot(2) = -ydot(1) - ydot(3)
        retval = 0                                ! 0 = success
    end function rhs

end program robertson
```

Output:

```
t =  4.000E-01   y =  9.851641E-01  3.386242E-05  1.480205E-02
t =  4.000E+00   y =  9.055097E-01  2.240338E-05  9.446793E-02
t =  4.000E+01   y =  7.157956E-01  9.183493E-06  2.841952E-01
t =  4.000E+02   y =  4.505294E-01  3.223044E-06  5.494674E-01
t =  4.000E+03   y =  1.831814E-01  8.940625E-07  8.168177E-01
t =  4.000E+04   y =  3.897848E-02  1.621552E-07  9.610214E-01
steps: 389   rhs evals: 586
```

To use FVODE in your own fpm project, add it as a dependency in `fpm.toml`:

```toml
[dependencies]
fvode = { git = "https://github.com/jacobwilliams/fvode" }
```

For non-stiff problems, use `CV_ADAMS` and attach the fixed-point nonlinear
solver instead of a linear solver:
`call cvode%set_nonlinear_solver(sun_nonlinsol_fixedpoint(n, 0), flag)` from
`fvode_nonlinsol_fixedpoint` (see
[example/cv_analytic_fp.f90](example/cv_analytic_fp.f90)). To pass problem data to
the right-hand side, extend `cvode_solver` with your own components and use
`select type` inside the function. The programs in [example/](example/) are
ports of the SUNDIALS C examples and show this and the other solver options.

### Verifying against the C reference

The `reference` environment builds the pinned C SUNDIALS and runs the trace drivers
in [reference/](reference/):

```
pixi run -e reference ref-build   # fetch + build SUNDIALS v7.9.0 (once)
pixi run -e reference ref-run     # run the C trace drivers, regenerate test data
pixi run -e reference compare     # compare FVODE traces against the C traces
```

### License

FVODE is released under the [MIT License](LICENSE).

FVODE is derived from SUNDIALS, which is Copyright (c) 2002-2026, Lawrence Livermore
National Security, Southern Methodist University, University of Maryland Baltimore
County, and the SUNDIALS contributors, and is distributed under the BSD 3-Clause
License. The SUNDIALS license and notice are reproduced in
[LICENSE-SUNDIALS](LICENSE-SUNDIALS).

### References

* Original CVODE: https://github.com/LLNL/sundials
* A. C. Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and Differential/Algebraic
  Equation Solvers," ACM Trans. Math. Softw. 31(3), 2005.
* D. J. Gardner et al., "Enabling new flexibility in the SUNDIALS suite of nonlinear
  and differential/algebraic equation solvers," ACM Trans. Math. Softw. 48(3), 2022.

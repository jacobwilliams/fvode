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

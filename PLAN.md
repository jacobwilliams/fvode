# FVODE Conversion Plan: CVODE (C) → Modern Fortran

## 1. Goal and guiding principles

Port the CVODE integrator from SUNDIALS to idiomatic Modern Fortran (F2008/F2018).
Build it with **fpm** and run every tool from the **pixi** environment.

* **Faithful algorithm.** Every numerical decision should match CVODE: step size and order
  selection, error test, Nordsieck history updates, nonlinear solver convergence test,
  Jacobian/LS setup heuristics, root finding, constraints, stability limit detection.
  The translation must be **mathematically equivalent**, not bit-for-bit: the compiler may
  optimize (fused multiply-add, vectorization), so real values drift from the C values by
  rounding. Given the same problem, options and tolerances, FVODE should take the same steps
  with the same integrator statistics (`nst`, `nfe`, `nsetups`, `netf`, `nni`, `ncfn`,
  `qlast`, …), with step sizes and solutions that agree to within that rounding drift.
* **Pinned reference.** Translate from one fixed upstream version:
  **SUNDIALS v7.9.0** (commit `312fc0f`, `src/cvode/`). Every translated routine gets a
  provenance comment (`!! Translated from cvode.c: cvStep (SUNDIALS v7.9.0)`), so a later
  upstream diff can be mapped back onto our code.
* **Modern Fortran, not transliterated C.** Use derived types, allocatable arrays, abstract
  interfaces, type-bound procedures, `pure`/`elemental` where possible, and no global state.
  The *algorithm* stays the same. The *plumbing* (SUNContext, N_Vector op tables, void
  pointers, error-handler stacks, loggers, profilers) is replaced by native Fortran
  constructs.
* **Licensing.** SUNDIALS is BSD-3-Clause. We must keep its copyright notice and license
  text (add `LICENSE-SUNDIALS` and a NOTICE section in the README), even though FVODE
  itself is MIT.

## 2. Scope

### In scope (serial, single-node CVODE)

| Upstream file | Lines | Content | FVODE target |
|---|---|---|---|
| `cvode.c` | 5003 | Main integrator: `CVodeCreate/Init/ReInit`, `CVode`, `cvStep`, `cvHin`, `cvPredict`, `cvSet*`, `cvAdjust*`, `cvDoErrorTest`, `cvCompleteStep`, `cvPrepareNextStep`, `cvSetEta`, `cvChooseEta`, `cvBDFStab`/`cvSLdet`, root finding (`cvRcheck1-3`, `cvRootfind`), `CVodeGetDky`, `cvEwtSet*`, constraints | `fvode_core.f90`, `fvode_root.f90`, `fvode_stab.f90` |
| `cvode_impl.h` | 800 | `CVodeMemRec` (157 `cv_` fields), constants (`L_MAX`, `ETA_*_DEFAULT`, `CORTES`, …) | `fvode_types.f90`, `fvode_internal_parameters.f90` |
| `cvode_io.c` | 1841 | `CVodeSet*` / `CVodeGet*` / `CVodePrintAllStats` | type-bound setters/getters in `fvode_io.f90` |
| `cvode_nls.c` | 490 | Interface between CVODE and the nonlinear solver (`cvNlsResid`, `cvNlsLSetup`, `cvNlsConvTest`, …) | `fvode_nls.f90` |
| `cvode_ls.c` / `cvode_ls_impl.h` | 1976 | Generic linear solver interface: DQ Jacobian (dense/band), Jacobian-times-vector DQ, `cvLsSetup/Solve`, `msbj`/`dgmax` heuristics | `fvode_ls.f90` |
| `cvode_diag.c` | 460 | Diagonal approximate Jacobian solver | `fvode_diag.f90` |
| `cvode_bandpre.c` | 596 | Banded preconditioner for Krylov solvers | `fvode_bandpre.f90` |
| `cvode_proj.c` | 456 | Projection onto constraint manifold | `fvode_proj.f90` (phase 6) |
| `cvode_resize.c` | 511 | `CVodeResizeHistory` | optional, phase 6 |

The **SUNDIALS shared components** that CVODE depends on algorithmically are also in scope:

| Upstream | Content | FVODE target |
|---|---|---|
| `nvector/serial` | vector ops (`N_VLinearSum`, `N_VWrmsNorm`, `N_VLinearCombination`, `N_VConstrMask`, `N_VMinQuotient`, …) | Fortran array expressions; `fvode_nvector.f90` for the few non-trivial ones |
| `sundials_math.c` | `SUNRsqrt`, `SUNRpowerR`, `SUNMIN`, … | Fortran intrinsics |
| `sundials_dense.c`, `sundials_band.c` | `denseGETRF/GETRS`, `bandGBTRF/GBTRS` (SUNDIALS's own LU, **not** LAPACK) | `fvode_dense.f90`, `fvode_band.f90` |
| `sunmatrix/dense`, `sunmatrix/band` | matrix storage, `ScaleAddI`, `MatVec` | `fvode_matrix*.f90` |
| `sunlinsol/dense`, `band` | direct solvers | `fvode_linsol_dense.f90`, `fvode_linsol_band.f90` |
| `sunlinsol/spgmr`, `spfgmr`, `spbcgs`, `sptfqmr`, `pcg` + `sundials_iterative.c` | Krylov solvers, `SUNModifiedGS`, `SUNClassicalGS`, `SUNQRfact/QRsol` | `fvode_linsol_sp*.f90`, `fvode_iterative.f90` |
| `sunnonlinsol/newton`, `fixedpoint` | Newton iteration, fixed point with Anderson acceleration | `fvode_nonlinsol_newton.f90`, `fvode_nonlinsol_fixedpoint.f90` |

### Out of scope (initially)

* MPI/parallel/GPU vectors, `cvode_bbdpre.c` (parallel preconditioner), `cvode_fused_*`.
* `SUNContext`, logger, profiler, adiak metadata, `cvode_cli.c` (command-line options),
  SUNAdaptController and adjoint/checkpointing infrastructure.
* Sparse matrix and third-party solvers (KLU, SuperLU, MAGMA, oneMKL, PETSc).
  Sparse + an optional LAPACK backend may come later (see phase 7).
* The existing SUNDIALS Fortran 2003 bindings (`fmod_int*`). FVODE is a native
  implementation, not a binding.
* CVODES (sensitivities). The architecture should not block a later port.

## 3. Key design decisions

### 3.1 Kinds and precision

* `fvode_kinds` module: `integer, parameter :: wp = real64` (default), `ip = int64` for
  `sunindextype` equivalents (or `int32` to match a 32-bit build).
* Use fpm's `[preprocess]` (cpp) with `-DREAL32` / `-DREAL128` to select precision, as
  SUNDIALS does with `SUNDIALS_SINGLE/DOUBLE/EXTENDED_PRECISION`.
* `UNIT_ROUNDOFF` = `epsilon(1.0_wp)`, which is the same as SUNDIALS `SUN_UNIT_ROUNDOFF` (DBL_EPSILON).

### 3.2 N_Vector → plain arrays and array syntax

CVODE only ever uses serial vectors through the N_Vector op table, so the abstraction adds
nothing algorithmically. Vectors are `real(wp)` arrays, and the vector operations are written
as ordinary Fortran array expressions in the translated code:

| SUNDIALS | Fortran |
|---|---|
| `N_VLinearSum(a,x,b,y,z)` | `z = a*x + b*y` |
| `N_VScale(c,x,z)`, `N_VConst(c,z)` | `z = c*x`, `z = c` |
| `N_VProd`, `N_VDiv`, `N_VAbs`, `N_VInv`, `N_VAddConst` | `x*y`, `x/y`, `abs(x)`, `1/x`, `x + b` |
| `N_VDotProd`, `N_VMaxNorm`, `N_VMin` | `dot_product`, `maxval(abs(x))`, `minval` |
| `N_VLinearCombination`, `N_VScaleAddMulti`, `N_VScaleVectorArray` | loops over columns of `zn` |

Only the operations with real logic are functions in `fvode_nvector`: `n_vwrmsnorm`,
`n_vinvtest`, `n_vconstrmask`, `n_vminquotient`.

The Nordsieck array `cv_zn[L_MAX]` (an array of N_Vectors) becomes
`real(wp), allocatable :: zn(:,:)` with bounds `zn(1:n, 0:L_MAX-1)`. **Keep zero-based
second dimension** so indices like `zn[j]`, `l[j]`, `tq[1..5]`, `tau[1..L_MAX]` map 1:1
onto the C source. This removes a large class of off-by-one translation errors.

### 3.3 Solver object and user callbacks

* `type :: cvode_solver` carries all `CVodeMemRec` fields as components. Pointer fields
  become `allocatable` arrays, and `sunbooleantype` becomes `logical`.
* User callbacks (`f`, `g` root function, `jac`, `jtimes`, `psetup`, `psolve`, `ewtfn`,
  `projfn`) are **procedure pointer components** with abstract interfaces whose first
  argument is `class(cvode_solver), intent(inout) :: me`. Users pass data by **extending**
  `cvode_solver` (no `void* user_data`). An optional `class(*), pointer :: user_data` can
  stay for C-style use.
* The public API is type-bound and mirrors the C API names so the CVODE docs still apply:

  ```fortran
  call solver%create(lmm=CV_BDF)                    ! CVodeCreate
  call solver%init(f=rhs, t0=t0, y0=y0, flag=ierr)  ! CVodeInit
  call solver%set_tolerances(reltol, abstol)        ! CVodeSStolerances / SVtolerances
  call solver%set_linear_solver(dense_linsol())     ! CVodeSetLinearSolver
  call solver%root_init(nrtfn, g)                   ! CVodeRootInit
  call solver%solve(tout, y, t, CV_NORMAL, ierr)    ! CVode
  call solver%get_dky(t, k, dky, ierr)              ! CVodeGetDky
  call solver%print_all_stats(unit)                 ! CVodePrintAllStats
  ```

### 3.4 Linear and nonlinear solver abstraction

These *are* algorithmic interfaces (`lsetup`, `lsolve`, `ninit`, convergence tests), so
keep them polymorphic:

* `type, abstract :: sun_linear_solver` with deferred `setup`, `solve`, `get_type`
  (direct / iterative / matrix-iterative), `set_atimes`, `set_preconditioner`,
  `set_scaling_vectors`, `num_iters`, `res_norm`, …
* `type, abstract :: sun_matrix` with deferred `zero`, `copy`, `scale_add_i`, `matvec`.
  Concrete types are `dense_matrix` (column-major `a(n,n)`) and `band_matrix` (keep
  SUNDIALS's `smu = mu + ml` storage layout so `bandGBTRF` ports line for line).
* `type, abstract :: sun_nonlinear_solver` with deferred `setup`, `solve`, and
  `set_sys_fn`/`set_lsetup_fn`/`set_lsolve_fn`/`set_conv_test_fn`. The CVODE→NLS function
  hooks are procedure pointers bound to the owning `cvode_solver`.

### 3.5 Return codes and errors

* Keep **every C return code as a named constant with the same integer value**
  (`CV_SUCCESS = 0`, `CV_TSTOP_RETURN = 1`, `CV_ROOT_RETURN = 2`, `CV_WARNING = 99`,
  `CV_TOO_MUCH_WORK = -1`, …, `SUN_NLS_CONV_RECVR`, `SUNLS_*`). The control flow in
  `cvHandleNFlag`/`cvHandleFailure` depends on these values.
* Error messages come from a single `process_error` (replacing `cvProcessError`) that
  writes to a user-selectable unit and keeps the same message text, plus an optional
  user error-handler procedure pointer.
* No `stop`/`error stop` inside the library.

### 3.6 Coding conventions

| C | Fortran |
|---|---|
| `sunrealtype` | `real(wp)` |
| `long int` counters | `integer(int64)` |
| `SUNRabs`, `SUNMAX`, `SUNMIN`, `SUNRsqrt` | `abs`, `max`, `min`, `sqrt` |
| `SUNRpowerR(x,y)`, `SUNRpowerI(x,i)` | `x**y`, `x**i` |
| `cv_mem->cv_X` | `me%X`; renames for case clashes: `cv_L` → `LL`, `M[]` → `mm` |
| `#define` constants | `parameter`s with the same names: public codes in `fvode_parameters`, internal ones (`ONE`, `CORTES`, `ETA_MAX_FS_DEFAULT`, …) in `fvode_internal_parameters` (not re-exported). Only rename: `TINY` → `CV_TINY` (clashes with the intrinsic) |
| `static` functions | `private` module procedures |
| `goto`-free loops with `break` | `do … exit`; named blocks/`cycle` where C uses `continue` |
| `0`-based arrays | explicit lower bounds (`0:`) where the algorithm indexes from 0 |

**Documentation** (FORD): every procedure has a `!>` header (what it does and which C function
it translates), or a `!!` docstring for small contained procedures. Every derived type,
argument, local variable, parameter, type component, procedure-pointer component and
type-bound binding has a `!!` comment. `pixi run doc-audit`
([tools/doc_audit.py](tools/doc_audit.py)) checks this and must report 0 findings. Don't start
a `!!` line with `word:`, because FORD reads that as a metadata key.

FVODE is compiled with normal optimizations. For debugging, an FVODE build with
`-ffp-contract=off` reproduces the C reference (which is built that way) exactly; see §6.

## 4. Repository layout

```
fvode/
├── fpm.toml
├── pixi.toml
├── LICENSE, LICENSE-SUNDIALS
├── src/
│   ├── fvode.f90                 ! umbrella module (public API)
│   ├── fvode_kinds.F90           ! wp / ip / li kinds (cpp-selectable precision)
│   ├── fvode_parameters.f90      ! public constants + return codes
│   ├── fvode_internal_parameters.f90 ! internal CVODE constants
│   ├── fvode_nvector.f90         ! the few non-trivial vector ops
│   ├── core/  fvode_types.f90 fvode_core.f90 fvode_root.f90 fvode_stab.f90
│   │          fvode_io.f90 fvode_nls.f90 fvode_ls.f90 fvode_diag.f90
│   │          fvode_bandpre.f90 fvode_proj.f90
│   ├── matrix/   fvode_matrix.f90 fvode_matrix_dense.f90 fvode_matrix_band.f90
│   │             fvode_dense.f90 fvode_band.f90
│   ├── linsol/   fvode_linsol.f90 fvode_linsol_dense.f90 fvode_linsol_band.f90
│   │             fvode_iterative.f90 fvode_linsol_spgmr.f90 … _spfgmr _spbcgs _sptfqmr _pcg
│   └── nonlinsol/ fvode_nonlinsol.f90 fvode_nonlinsol_newton.f90 fvode_nonlinsol_fixedpoint.f90
├── test/          ! unit tests (test-drive) + regression vs C reference
├── example/       ! Fortran ports of examples/cvode/serial/*.c
├── reference/     ! C harness + scripts to build/run pinned SUNDIALS
└── doc/           ! FORD output
```

## 5. Build and environment (pixi + fpm)

The real configuration is in [pixi.toml](pixi.toml) and [fpm.toml](fpm.toml). This section
summarises it.

**Default environment** (gfortran 16, fpm, python, FORD):

| Task | What it does |
|---|---|
| `pixi run build` | `fpm build --profile release` |
| `pixi run test` | `fpm test --profile debug` (test-drive suites in `test/`) |
| `pixi run examples` | `fpm run --example --all` |
| `pixi run lint` | debug build with `-Wall -Wextra -std=f2018 -fcheck=all` |
| `pixi run docs` | FORD docs into `doc/` (registers graphviz plugins first: pixi skips conda post-link scripts) |
| `pixi run doc-audit` | check that everything is documented ([tools/doc_audit.py](tools/doc_audit.py)) |

**`reference` environment** (cmake, ninja, git, plus a C compiler: `gcc 16` on Linux,
because the `c-compiler` metapackage pins gcc 15 and clashes with gfortran 16; `c-compiler`
(clang) on macOS; Windows is not supported for this env):

| Task | What it does |
|---|---|
| `ref-build` | [reference/build_reference.sh](reference/build_reference.sh): clone tag `v7.9.0`, build static CVODE + serial components with `-ffp-contract=off` into `reference/install` (runs once) |
| `ref-run` | [reference/run_reference.sh](reference/run_reference.sh): regenerate `test/data/sundials_constants.txt`, then build and run every `reference/trace_*.c` into `reference/out/trace_*.csv` |
| `fvode-traces` | [reference/run_fvode_traces.sh](reference/run_fvode_traces.sh): run every FVODE trace driver `reference/fortran/trace_*.f90` (fpm examples) into `reference/out/fvode_*.csv` |
| `compare-exact` | diagnostic: FVODE trace drivers built with `-ffp-contract=off`; reports how long each trace stays bit-for-bit identical to C, then applies the tolerance checks ([reference/compare.py](reference/compare.py) `--exact`) |
| `compare-examples` | [reference/compare_examples.py](reference/compare_examples.py): run each C example and its FVODE port (no FMA) and require identical output |
| `compare` | sanity check of the normal optimized build (tolerances; see §6) |

We build the reference from source at the pinned tag, not from the conda-forge `sundials`
package, so the version matches exactly.

**fpm.toml:** `test-drive` v0.5.0 is a dev-dependency. Kinds come from `fvode_kinds.F90`
(preprocessed, `-DREAL32` / `-DREAL128` / `-DINT32`). The FORD settings moved to
`[extra.ford]`: FORD ignored the old `[extra.fvode]` section, and `front_page_file`
is not a FORD 7 option.

## 6. Verification strategy

1. **Unit tests per ported kernel** (test-drive): vector ops, dense/band LU vs known
   factors, Krylov solvers on small SPD/non-symmetric systems, Newton/fixed point on
   scalar problems.
2. **Step-trace comparison (the main fidelity check).** Each C driver
   `reference/trace_<name>.c` runs a problem in `CV_ONE_STEP` mode and writes one CSV row
   per internal step through the public `CVodeGet*` API (see
   [reference/trace_common.h](reference/trace_common.h)):
   `nst,t,hu,qu,hcur,qcur,nfe,nni,nncf,netf,nsetups,njev,nfels,flag,y1..yN`.
   This needs no changes to SUNDIALS. The matching FVODE program writes the same columns to
   `reference/out/fvode_<name>.csv`. `reference/compare.py` compares the two traces in two
   modes:
   * **`compare-exact` (diagnostic).** FVODE is built with `-ffp-contract=off`, like the C
     reference, and the report says how many rows are bit-for-bit identical. Where the
     arithmetic matches C's, a faithful translation stays identical, so the first differing
     row points to the code that ran differently. One intended difference: integer powers
     (`x**i` squares, `SUNRpowerI` multiplies in a loop). It first matters at an order
     increase to q >= 4 with a changed step size. After the exact report, the tolerance
     checks decide pass/fail.
   * **`compare-examples`.** Each FVODE example port, built without FMA, must print exactly
     what the C example built here prints (workspace lines excluded).
   * **`compare` (sanity check of the optimized build).** If the integer columns (steps,
     orders, counters) match, the real columns must agree to a mixed tolerance (default
     `1e-3` of the value or of the column maximum). On sensitive problems a rounding
     difference can flip a near-threshold decision, after which the step sequences differ
     (`analytic_fp` diverges at its second output). That is reported as DIVERGED, not as
     a failure, if the final solution agrees to `1e-4` of the column scale and the step
     count differs by at most 10%.
3. **Example regression.** Port the upstream serial examples and compare against the
   shipped `.out` files. Solution values should agree to the printed digits. Statistics
   only for problems that aren't rounding-sensitive: for `cv_analytic_fp` the `.out` file
   reports 2065 steps, the pinned C build here takes 2160, and optimized FVODE 2024,
   while the exact-mode trace matches the C build bit for bit. Planned examples:
   `cvRoberts_dns` (BDF + Newton + dense, with roots), `cvRoberts_dns_uw` (user `ewt`),
   `cvRoberts_dns_constraints`, `cvRoberts_dns_negsol`, `cvAdvDiff_bnd` (band),
   `cvDirectDemo_ls` (Adams/BDF × functional/dense/band/diag matrix), `cvDiurnal_kry`
   (SPGMR + user preconditioner), `cvDiurnal_kry_bp` (bandpre), `cvKrylovDemo_ls`
   (all Krylov solvers), `cvKrylovDemo_prec`, `cvDisc_dns` (discontinuities/ReInit),
   `cvAnalytic_mels`, `cvVdp_auto_nls`, `cvPendulum_dns` (projection), `cvRocket_dns`
   (root finding + stop), `cvParticle_dns`.
4. **CI** with GitHub Actions running `pixi run test`, and in the reference env
   `compare`, `compare-exact` and `compare-examples`, on linux-64 and osx-arm64.

## 7. Phased work plan

Each phase ends with green tests and a commit. The order follows dependencies, so every
phase can be verified on its own.

### Phase 0: Scaffolding ✅ done
Delivered: `LICENSE` (MIT + derivation note), `LICENSE-SUNDIALS`, `fvode_kinds`,
`fvode_parameters`, `fvode_internal_parameters`, and the `fvode` umbrella module. There is a
unit test that checks all 207 upstream constants against
`test/data/sundials_constants.txt`, which is generated from the C headers by
`reference/constants.c` plus the `.c`-local `#define`s parsed by `reference/parse_defines.py`.
The reference build and trace harness exist, with two traces so far:
`vdp_adams_fp` (Adams + fixed point, 217 steps; Phase 2 target) and
`roberts_bdf_dense` (BDF + Newton + dense DQ Jacobian, 648 steps; Phase 3 target).

Original task list:
* Update `pixi.toml` (tasks, reference feature) and `fpm.toml` (test-drive, keywords).
* Add `LICENSE-SUNDIALS`, `fvode_kinds`, `fvode_parameters` (all constants and return codes
  copied from `cvode.h`, `cvode_impl.h`, `cvode_ls.h`, `sundials_*` headers).
* Build the C reference (`pixi run -e reference ref-build`) and write the C trace harness.
* **Exit criterion:** `pixi run build`, `pixi run test` (empty), and `pixi run -e reference ref-build` all succeed.

### Phase 1: Vector operations ✅ done
* First done as bitwise-exact ports of the serial C kernels and `sundials_math` helpers.
  Then simplified when we decided that mathematical equivalence is the goal: simple
  operations are now array expressions, `fvode_math` was removed, and `fvode_nvector` keeps
  `n_vwrmsnorm`, `n_vinvtest`, `n_vconstrmask` and `n_vminquotient`, with unit tests.
* Finding kept from the exact port: serial vectors have fused operations **disabled** by
  default, so CVODE uses the generic `N_VLinearCombination`/`N_VScaleAddMulti` fallbacks.
  Mathematically these are plain sums, which is what the array code does.

### Phase 2: Core integrator, explicit path (functional/fixed-point iteration) ✅ done
(Since Phase 3, `compare-exact` is a diagnostic. Its exact match holds until the first
order increase past 3 with a changed step size, because of `x**i`; see §6.)
`fvode_core` holds the solver type, the step logic and the options/statistics
(`cvode.c`, `cvode_nls.c`, most of `cvode_io.c`). `fvode_nonlinsol` holds the nonlinear
solver interface, and `fvode_nonlinsol_fixedpoint` the fixed-point solver with Anderson
acceleration and damping. In v7.9.0 the fixed-point solver does its own QR (modified
Gram-Schmidt plus Givens downdating); `SUNQRAdd_*` is not used.

Verified with four traces. All are bit-for-bit identical to C in `compare-exact`:

| Trace | Problem | Steps |
|---|---|---|
| `vdp_adams_fp` | van der Pol, Adams + fixed point | 217 |
| `vdp_adams_aa` | same, Anderson m = 2, up to 5 iterations (subspace recycling) | 270 |
| `vdp_adams_aa_damped` | same, Anderson m = 3, beta = 0.75 | 241 |
| `analytic_fp` | `cv_analytic_fp` problem, 10 CV_NORMAL outputs | 2160 |

There is also a user example, [example/cv_analytic_fp.f90](example/cv_analytic_fp.f90),
a port of the SUNDIALS F2003 example, and unit tests `test/test_core.f90` (accuracy,
`get_dky`, stop time, input checks).

Not yet ported and marked `[not yet ported]` in the code: linear-solver hooks (Phase 3),
roots/constraints/STALD (Phase 4), projection and history resize (Phase 6). Also left for
later: `CVodeSetUserData` (replaced by extending the type), monitor functions, and
`CVodeGetWorkSpace`.

Port in this order. Each item is a direct translation of the named C function(s):
1. `cvode_solver` type ← `CVodeMemRec`; `create`, `init`, `reinit`, `free` (finalizer).
2. Tolerances and weights: `cvEwtSet`, `cvEwtSetSS`, `cvEwtSetSV`, user `ewt` function.
3. `cvInitialSetup`, `cvHin`, `cvUpperBoundH0`, `cvYddNorm`.
4. Method coefficients: `cvSet`, `cvSetAdams`, `cvAdamsStart`, `cvAdamsFinish`, `cvAltSum`,
   `cvSetBDF`, `cvSetTqBDF`.
5. History: `cvPredict`, `cvRescale`, `cvRestore`, `cvAdjustParams`, `cvAdjustOrder`,
   `cvAdjustAdams`, `cvAdjustBDF`, `cvIncreaseBDF`, `cvDecreaseBDF`.
6. Step control: `cvStep`, `cvDoErrorTest`, `cvCompleteStep`, `cvPrepareNextStep`,
   `cvSetEta`, `cvChooseEta`, `cvComputeEtaqm1/qp1`, `cvHandleFailure`, `cvHandleNFlag`.
7. Driver: `CVode` (`CV_NORMAL`/`CV_ONE_STEP`, `tstop`, `mxstep`, `hmin/hmax`, warnings),
   `CVodeGetDky`.
8. `fvode_nls` (`cvNls`, `cvNlsInit`, `cvNlsResid*`, `cvNlsFPFunction`, `cvNlsConvTest`) and
   the fixed-point nonlinear solver (including Anderson acceleration).
9. `fvode_io`: all `Set*`/`Get*` and `print_all_stats` (table and CSV formats).
* **Exit criterion:** nonstiff Adams + fixed-point examples (`cvAnalytic` style, a
  `cvDirectDemo_ls` problem 1 with functional iteration) match the C step trace (same steps
  and counters; exact with `compare-exact`).

### Phase 3: Newton and direct linear solvers ✅ done
Delivered:
* `fvode_matrix`: `dense_matrix`, `band_matrix` (SUNDIALS band storage, `smu = mu + ml`).
* `fvode_linsol`: the abstract `sun_linear_solver`; `fvode_linsol_direct`: dense and band
  LU (`denseGETRF/GETRS`, `bandGBTRF/GBTRS`).
* `fvode_nonlinsol_newton`. `init` attaches a Newton solver by default, as in C.
* Submodule `fvode_core_ls` (`cvode_ls.c`: `jbad` heuristic, saved Jacobian, `A = I - gamma J`,
  `gamrat` solution scaling, DQ dense/band Jacobians, user Jacobian). Submodule
  `fvode_core_diag` (`cvode_diag.c`). An `ls_kind` field selects the attached interface,
  replacing C's `cv_linit`/`cv_lsetup`/... function pointers.
* Examples `cv_direct_demo` (port of `cvDirectDemo_ls`: 16 runs covering dense user/DQ, band
  user/DQ and diagonal, with Adams and BDF) and `cv_advdiff_bnd` (port of `cvAdvDiff_bnd`, user
  data by extending `cvode_solver`). Unit tests `test/test_linsol.f90`.

Verification:
* `compare-examples` (new): both example ports, built without FMA, print **exactly** the same
  output as the C examples built here. That covers every direct-solver path above.
* `roberts_bdf_dense` trace: bit-for-bit for 278 of 648 steps, i.e. up to the first order
  increase to 4. Integer powers there differ by design (`x**i` vs the `SUNRpowerI` loop; the
  user chose idiomatic `x**i`). The run then stays within tolerance (final solution within
  1e-8).
* The optimized `cv_advdiff_bnd` reproduces the upstream `.out` statistics exactly.

Upstream quirk: `CVodePrintAllStats` reads the `CVLsMem` statistics whenever any linear solver
is attached, including `CVDiag` (wrong memory). FVODE prints them only for the generic interface.

Original task list:
* Dense and band matrices, `denseGETRF/GETRS`, `bandGBTRF/GBTRS`, `ScaleAddI`.
* Dense and band linear solvers.
* Newton nonlinear solver.
* `fvode_ls` ← `cvode_ls.c`: `cvLsInitialize`, `cvLsSetup` (the `gamma`/`gamrat`/`DGMAX_LSETUP`/
  `MSBP` logic and `jbad`/`jcur` handling), `cvLsSolve` (including `gamrat` scaling of
  the solution), `cvLsDenseDQJac`, `cvLsBandDQJac`, `cvLsLinSys`, user `jac`.
* `fvode_diag` ← `cvode_diag.c`.
* **Exit criterion:** `cvRoberts_dns` (without roots), `cvAdvDiff_bnd`, all of `cvDirectDemo_ls`
  match the traces and `.out` files.

### Phase 4: Root finding, constraints, stability limit detection
* `cvRcheck1`, `cvRcheck2`, `cvRcheck3`, `cvRootfind` (Illinois algorithm), `CVodeRootInit`,
  `SetRootDirection`, `SetNoInactiveRootWarn`, `GetRootInfo`.
* `cvCheckConstraints` and `CVodeSetConstraints`.
* `cvBDFStab`, `cvSLdet` (STALD).
* **Exit criterion:** `cvRoberts_dns`, `cvRoberts_dns_constraints`, `cvRoberts_dns_negsol`,
  `cvRocket_dns`, `cvDisc_dns` match.

### Phase 5: Iterative linear solvers and preconditioning
* `fvode_iterative` (`SUNModifiedGS`, `SUNClassicalGS`, `SUNQRfact`, `SUNQRsol`).
* SPGMR, SPFGMR, SPBCGS, SPTFQMR, PCG (left/right/both preconditioning, scaling vectors).
* `cvLsATimes`, `cvLsDQJtimes`, `cvLsPSetup`, `cvLsPSolve`, `eplifac`/`nrmfac` handling,
  user `jtimes`/`jtsetup`/`psetup`/`psolve`.
* `fvode_bandpre` ← `cvode_bandpre.c`.
* **Exit criterion:** `cvDiurnal_kry`, `cvDiurnal_kry_bp`, `cvKrylovDemo_ls` (all solver variants),
  `cvKrylovDemo_prec` match.

### Phase 6: Remaining features
* `cvode_proj.c` (projection: `CVodeSetProjFn`, error projection, `cvDoProjection`), verified with `cvPendulum_dns`.
* `CVodeSetEtaFixedStepBounds`, `SetEtaMaxFirstStep/EarlyStep/General`, `SetEtaMin*`,
  `SetDeltaGammaMaxLSetup`, `SetNonlinConvCoef`, and other newer tuning knobs from `cvode_io.c`.
* `CVodeSetMaxHnilWarns`, user-supplied nonlinear system function (`cvAnalytic_mels`, `cvVdp_auto_nls`).
* Optional: `cvode_resize.c` (`CVodeResizeHistory`), checked against `cvParticle_dns` and the resize tests.

### Phase 7: Polish and optional extensions
* FORD documentation for every public procedure, with links to the corresponding CVODE
  User Guide sections. `README` usage example.
* Performance pass: `contiguous`, avoid temporaries in hot loops, `-O3 -march=native`
  release profile, and a benchmark against C CVODE. Optimisations must keep `compare` passing.
* Optional: LAPACK dense/band backends (`dgetrf`/`dgbtrf` via conda-forge `liblapack`),
  sparse CSC matrix + a native sparse LU, or a KLU binding. These come after fidelity is
  established and are clearly marked as non-reference paths.
* Optional: precision variants (`real32`, `real128`) tested in CI.

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Rounding drift (FMA, vectorization) hides or mimics translation errors | `compare-exact` (`-ffp-contract=off` on both sides, bit-for-bit) is the fidelity check; tolerance-based `compare` only sanity-checks the optimized build |
| Off-by-one errors from 0-based C indexing | Zero-based lower bounds on `zn`, `l`, `tq`, `tau`, `ssdat` and similar arrays; step-trace comparison detects divergence on the first step |
| Large, stateful `cvStep` with many interacting flags | Port `cvode.c` routines nearly line for line with provenance comments; refactor only after the trace matches |
| Hidden dependencies on SUNDIALS internals (e.g., `N_VGetLength`, `SUNNonlinSol` callbacks reaching into `cv_mem`) | Grep each in-scope file for every external symbol before porting it, and record the mapping table in `doc/mapping.md` |
| Upstream evolves after v7.9.0 | Provenance comments + pinned tag; later, a scripted `git diff v7.9.0..vX src/cvode` review |
| Callback design limits users (no `void*`) | Extension of `cvode_solver` *plus* optional `class(*)` user data |

## 9. Deliverables checklist

- [ ] `pixi run build|test|examples|docs` all green on linux-64 and osx-arm64
- [ ] `pixi run -e reference compare` passes (identical steps and counters) for all ported examples
- [ ] Every upstream serial CVODE example in §6.3 has a Fortran port under `example/`
- [ ] `doc/mapping.md`: C function → Fortran procedure table
- [ ] `LICENSE-SUNDIALS` + NOTICE in the README
- [ ] FORD docs published

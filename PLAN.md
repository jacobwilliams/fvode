# FVODE Conversion Plan: CVODE (C) → Modern Fortran

## 1. Goal and guiding principles

Port the CVODE integrator from SUNDIALS to idiomatic Modern Fortran (F2008/F2018).
Build it with **fpm** and run every tool from the **pixi** environment.

* **Faithful algorithm.** Every numerical decision should match CVODE: step size and order
  selection, error test, Nordsieck history updates, nonlinear solver convergence test,
  Jacobian/LS setup heuristics, root finding, constraints, stability limit detection.
  Given the same problem, options and tolerances, FVODE should produce **identical
  integrator statistics** (`nst`, `nfe`, `nsetups`, `netf`, `nni`, `ncfn`, `qlast`, `hlast`, …)
  and solutions that agree to near round-off.
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
| `cvode_impl.h` | 800 | `CVodeMemRec` (157 `cv_` fields), constants (`L_MAX`, `ETAMX*`, `CORTES`, …) | `fvode_types.f90`, `fvode_parameters.f90` |
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
| `nvector/serial` | vector ops (`N_VLinearSum`, `N_VWrmsNorm`, `N_VWL2Norm`, `N_VLinearCombination`, `N_VScaleAddMulti`, `N_VConstrMask`, `N_VMinQuotient`, …) | `fvode_nvector.f90` (procedures on `real(wp)` arrays) |
| `sundials_math.c` | `SUNRsqrt`, `SUNRpowerR`, `SUNRCompare` | `fvode_math.f90` |
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

### 3.2 N_Vector → plain contiguous arrays

CVODE only ever uses serial vectors through the N_Vector op table, so the abstraction adds
nothing algorithmically. Replace it with `real(wp), dimension(:), contiguous` arrays and a
module of vector kernels (`fvode_nvector`) that **replicate the serial implementation
exactly**, including floating-point operation order. Examples:

* `N_VWrmsNorm`: `sqrt(sum((x*w)**2)/n)`, accumulated in the same loop order. Do **not**
  use `norm2`: it rescales internally and would change rounding.
* `N_VLinearSum`: keep the special cases (`a==1, b==1`, `a==1, b==-1`, `a==b`, …), because
  they change which operations run (for example `y - x` instead of `1*y + (-1)*x`).
* Fused ops (`N_VLinearCombination`, `N_VScaleAddMulti`) are used in `cvPredict`,
  `cvRescale` and `CVodeGetDky`. Port the serial loop order exactly.

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
| `SUNRabs`, `SUNMAX`, `SUNMIN` | `abs`, `max`, `min` |
| `SUNRpowerR(x,y)`, `SUNRpowerI` | `x**y` (real), `x**i` (integer). Check against `pow` for edge cases |
| `SUNRsqrt(x)` | `merge(sqrt(x), 0, x>0)`, the same guard as SUNDIALS |
| `#define` constants | `real(wp), parameter` in `fvode_parameters` (`ONE`, `ETAMX1 = 10000`, `CORTES = 0.1`, …), same names |
| `static` functions | `private` module procedures |
| `goto`-free loops with `break` | `do … exit`; named blocks/`cycle` where C uses `continue` |
| `0`-based arrays | explicit lower bounds (`0:`) where the algorithm indexes from 0 |

Compile with `-ffp-contract=off` (no FMA contraction) in both FVODE and the C reference
during verification, so that operation-order fidelity can be checked.

## 4. Repository layout

```
fvode/
├── fpm.toml
├── pixi.toml
├── LICENSE, LICENSE-SUNDIALS
├── src/
│   ├── fvode.f90                 ! umbrella module (public API)
│   ├── fvode_kinds.f90
│   ├── fvode_parameters.f90      ! constants + return codes
│   ├── fvode_math.f90
│   ├── fvode_nvector.f90
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

### 5.1 `pixi.toml` additions

```toml
[dependencies]
gfortran = "16.*"
fpm = ">=0.13.0,<0.14"
python = "3.13.*"
ford = ">=7.0"

[feature.reference.dependencies]      # pinned C SUNDIALS for verification
c-compiler = "*"
cmake = ">=3.24"
ninja = "*"
numpy = "*"

[environments]
default   = { solver-group = "main" }
reference = { features = ["reference"], solver-group = "main" }

[tasks]
build    = "fpm build --profile release"
test     = "fpm test --profile debug"
examples = "fpm run --example --all --profile release"
docs     = "ford fpm.toml"
lint     = "fpm build --profile debug --flag '-Wall -Wextra -std=f2018 -fcheck=all'"

[feature.reference.tasks]
ref-fetch = "git clone --depth 1 --branch v7.9.0 https://github.com/LLNL/sundials.git reference/sundials"
ref-build = { cmd = "cmake -S reference/sundials -B reference/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-ffp-contract=off -DEXAMPLES_ENABLE_C=ON -DBUILD_ARKODE=OFF -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF && cmake --build reference/build", depends-on = ["ref-fetch"] }
compare  = { cmd = "python reference/compare.py", depends-on = ["ref-build", "examples"] }
```

We build the reference from source at the pinned tag, not from the conda-forge `sundials`
package, so that the version matches exactly and we can add per-step trace output.

### 5.2 `fpm.toml` changes

* `[dev-dependencies] test-drive = { git = "https://github.com/fortran-lang/test-drive", tag = "v0.5.0" }`
* `[preprocess] [preprocess.cpp] macros = []` for precision selection.
* Fix `keywords` (`["ode", "cvode", "bdf", "adams", "stiff", "integrator"]`).
* Keep `auto-tests`/`auto-examples`, and add `[[example]]` entries if names need control.

## 6. Verification strategy

1. **Unit tests per ported kernel** (test-drive): vector ops, dense/band LU vs known
   factors, Krylov solvers on small SPD/non-symmetric systems, Newton/fixed point on
   scalar problems.
2. **Step-trace comparison (the main fidelity check).** Add a small C harness in
   `reference/` that runs a problem through C CVODE and logs, for every internal step,
   `nst, t, h, q, hu, qu, nfe, nni, netf, ncfn, nsetups, acnrm, ‖zn[0]‖`. FVODE gets an
   equivalent optional trace hook (`solver%set_step_trace(unit)`). `reference/compare.py`
   diffs the two traces:
   * integer columns (`nst`, `q`, counters) must match **exactly**,
   * real columns must match to `≤ 1e-12` relative (target: bitwise with `-ffp-contract=off`).
   The first divergent step points straight to the routine that was mistranslated.
3. **Example regression.** Port the upstream serial examples and compare against the
   shipped `.out` files (final stats and solution values):
   `cvRoberts_dns` (BDF + Newton + dense, with roots), `cvRoberts_dns_uw` (user `ewt`),
   `cvRoberts_dns_constraints`, `cvRoberts_dns_negsol`, `cvAdvDiff_bnd` (band),
   `cvDirectDemo_ls` (Adams/BDF × functional/dense/band/diag matrix), `cvDiurnal_kry`
   (SPGMR + user preconditioner), `cvDiurnal_kry_bp` (bandpre), `cvKrylovDemo_ls`
   (all Krylov solvers), `cvKrylovDemo_prec`, `cvDisc_dns` (discontinuities/ReInit),
   `cvAnalytic_mels`, `cvVdp_auto_nls`, `cvPendulum_dns` (projection), `cvRocket_dns`
   (root finding + stop), `cvParticle_dns`.
4. **CI** with GitHub Actions running `pixi run test` and `pixi run compare` on linux-64 and
   osx-arm64.

## 7. Phased work plan

Each phase ends with green tests and a commit. The order follows dependencies, so every
phase can be verified on its own.

### Phase 0: Scaffolding (½ day)
* Update `pixi.toml` (tasks, reference feature) and `fpm.toml` (test-drive, keywords).
* Add `LICENSE-SUNDIALS`, `fvode_kinds`, `fvode_parameters` (all constants and return codes
  copied from `cvode.h`, `cvode_impl.h`, `cvode_ls.h`, `sundials_*` headers).
* Build the C reference (`pixi run -e reference ref-build`) and write the C trace harness.
* **Exit criterion:** `pixi run build`, `pixi run test` (empty), and `pixi run -e reference ref-build` all succeed.

### Phase 1: Math and vector kernels
* `fvode_math`, then `fvode_nvector` with every serial N_Vector op that CVODE and the solvers call
  (grep `N_V` across the in-scope files to get the exact list).
* **Exit criterion:** unit tests compare each op bitwise with C outputs, produced by a tiny C driver.

### Phase 2: Core integrator, explicit path (functional/fixed-point iteration)
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
   the fixed-point nonlinear solver (including Anderson acceleration and `SUNQRAdd_*` from
   `sundials_iterative.c`).
9. `fvode_io`: all `Set*`/`Get*` and `print_all_stats` (table and CSV formats).
* **Exit criterion:** nonstiff Adams + fixed-point examples (`cvAnalytic` style, a
  `cvDirectDemo_ls` problem 1 with functional iteration) match the C step trace exactly.

### Phase 3: Newton and direct linear solvers
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
  release profile, and a benchmark against C CVODE. Optimisations must not change the step trace.
* Optional: LAPACK dense/band backends (`dgetrf`/`dgbtrf` via conda-forge `liblapack`),
  sparse CSC matrix + a native sparse LU, or a KLU binding. These come after fidelity is
  established and are clearly marked as non-reference paths.
* Optional: precision variants (`real32`, `real128`) tested in CI.

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Floating-point differences from operation reordering, FMA, compiler vectorization | `-ffp-contract=off` for both builds during verification; port the exact loop order of serial N_Vector kernels; avoid whole-array intrinsics where they could reorder sums (`sum` is fine for WRMS only if it matches serial order; otherwise use explicit loops) |
| Off-by-one errors from 0-based C indexing | Zero-based lower bounds on `zn`, `l`, `tq`, `tau`, `ssdat` and similar arrays; step-trace comparison detects divergence on the first step |
| Large, stateful `cvStep` with many interacting flags | Port `cvode.c` routines nearly line for line with provenance comments; refactor only after the trace matches |
| Hidden dependencies on SUNDIALS internals (e.g., `N_VGetLength`, `SUNNonlinSol` callbacks reaching into `cv_mem`) | Grep each in-scope file for every external symbol before porting it, and record the mapping table in `doc/mapping.md` |
| Upstream evolves after v7.9.0 | Provenance comments + pinned tag; later, a scripted `git diff v7.9.0..vX src/cvode` review |
| Callback design limits users (no `void*`) | Extension of `cvode_solver` *plus* optional `class(*)` user data |

## 9. Deliverables checklist

- [ ] `pixi run build|test|examples|docs` all green on linux-64 and osx-arm64
- [ ] `pixi run -e reference compare` reports identical integer stats for all ported examples
- [ ] Every upstream serial CVODE example in §6.3 has a Fortran port under `example/`
- [ ] `doc/mapping.md`: C function → Fortran procedure table
- [ ] `LICENSE-SUNDIALS` + NOTICE in the README
- [ ] FORD docs published

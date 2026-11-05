/* Step trace: stability limit detection (STALD). A damped fast oscillation
 * (eigenvalues -1 +- 100i) plus a slow decay (-0.01), BDF + Newton + dense DQ
 * Jacobian, rtol = 1e-6, atol = 1e-10, with CVodeSetStabLimDet on. Once the
 * oscillation has decayed, the steps grow until hλ approaches the imaginary
 * axis, where BDF orders 3-5 are unstable: STALD reduces the order. */
#include <sundials/sundials_context.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "trace_common.h"

#define T0     SUN_RCONST(0.0)
#define TFINAL SUN_RCONST(500.0)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  sunrealtype* dd = N_VGetArrayPointer(ydot);
  dd[0] = -yd[0] + SUN_RCONST(100.0) * yd[1];
  dd[1] = -SUN_RCONST(100.0) * yd[0] - yd[1];
  dd[2] = -SUN_RCONST(0.01) * yd[2];
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  TRACE_CHECK(SUNContext_Create(SUN_COMM_NULL, &sunctx));
  N_Vector y = N_VNew_Serial(3, sunctx);
  N_VGetArrayPointer(y)[0] = SUN_RCONST(1.0);
  N_VGetArrayPointer(y)[1] = SUN_RCONST(0.0);
  N_VGetArrayPointer(y)[2] = SUN_RCONST(1.0);

  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  TRACE_CHECK(CVodeInit(cvode_mem, f, T0, y));
  TRACE_CHECK(CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-10)));
  SUNMatrix A        = SUNDenseMatrix(3, 3, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  TRACE_CHECK(CVodeSetLinearSolver(cvode_mem, LS, A));
  TRACE_CHECK(CVodeSetStabLimDet(cvode_mem, SUNTRUE));
  TRACE_CHECK(CVodeSetMaxNumSteps(cvode_mem, 100000));

  trace_header(3);
  sunrealtype t = T0;
  while (t < TFINAL)
  {
    int flag = CVode(cvode_mem, TFINAL, y, &t, CV_ONE_STEP);
    trace_row(cvode_mem, flag, t, y, 1);
    if (flag < 0) { break; }
  }
  long int nor;
  TRACE_CHECK(CVodeGetNumStabLimOrderReds(cvode_mem, &nor));
  fprintf(stderr, "steps: %s, order reductions by STALD: %ld\n", "see trace", nor);

  N_VDestroy(y);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&sunctx);
  return 0;
}

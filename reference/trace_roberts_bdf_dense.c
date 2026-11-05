/* Step trace: Robertson kinetics (cvRoberts_dns without root finding), BDF +
 * Newton + dense linear solver with the internal difference-quotient Jacobian,
 * vector absolute tolerances. Exercises the Phase 3 path. */
#include <sundials/sundials_context.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "trace_common.h"

#define T0     SUN_RCONST(0.0)
#define TFINAL SUN_RCONST(4.0e10)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  sunrealtype* dd = N_VGetArrayPointer(ydot);
  sunrealtype yd1, yd3;
  yd1 = dd[0] = SUN_RCONST(-0.04) * yd[0] + SUN_RCONST(1.0e4) * yd[1] * yd[2];
  yd3 = dd[2] = SUN_RCONST(3.0e7) * yd[1] * yd[1];
  dd[1]       = -yd1 - yd3;
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  TRACE_CHECK(SUNContext_Create(SUN_COMM_NULL, &sunctx));

  N_Vector y = N_VNew_Serial(3, sunctx);
  N_VGetArrayPointer(y)[0] = SUN_RCONST(1.0);
  N_VGetArrayPointer(y)[1] = SUN_RCONST(0.0);
  N_VGetArrayPointer(y)[2] = SUN_RCONST(0.0);

  N_Vector abstol = N_VNew_Serial(3, sunctx);
  N_VGetArrayPointer(abstol)[0] = SUN_RCONST(1.0e-8);
  N_VGetArrayPointer(abstol)[1] = SUN_RCONST(1.0e-14);
  N_VGetArrayPointer(abstol)[2] = SUN_RCONST(1.0e-6);

  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  TRACE_CHECK(CVodeInit(cvode_mem, f, T0, y));
  TRACE_CHECK(CVodeSVtolerances(cvode_mem, SUN_RCONST(1.0e-4), abstol));
  SUNMatrix A        = SUNDenseMatrix(3, 3, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  TRACE_CHECK(CVodeSetLinearSolver(cvode_mem, LS, A));
  TRACE_CHECK(CVodeSetMaxNumSteps(cvode_mem, 100000));

  trace_header(3);
  sunrealtype t = T0;
  while (t < TFINAL)
  {
    int flag = CVode(cvode_mem, TFINAL, y, &t, CV_ONE_STEP);
    trace_row(cvode_mem, flag, t, y, 1);
    if (flag < 0) { break; }
  }

  N_VDestroy(y);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&sunctx);
  return 0;
}

/* Step trace: van der Pol (cvDirectDemo_ls problem 1), Adams + fixed-point
 * iteration, scalar tolerances (reltol = 0, abstol = 1e-6).
 * Exercises the Phase 2 path: no linear solver. */
#include <sundials/sundials_context.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include "trace_common.h"

#define ETA    SUN_RCONST(3.0)
#define T0     SUN_RCONST(0.0)
#define TFINAL (SUN_RCONST(1.39283880203) + SUN_RCONST(3.0) * SUN_RCONST(2.214773875))

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  sunrealtype* dd = N_VGetArrayPointer(ydot);
  dd[0] = yd[1];
  dd[1] = (SUN_RCONST(1.0) - yd[0] * yd[0]) * ETA * yd[1] - yd[0];
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  TRACE_CHECK(SUNContext_Create(SUN_COMM_NULL, &sunctx));

  N_Vector y = N_VNew_Serial(2, sunctx);
  N_VGetArrayPointer(y)[0] = SUN_RCONST(2.0);
  N_VGetArrayPointer(y)[1] = SUN_RCONST(0.0);

  void* cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  TRACE_CHECK(CVodeInit(cvode_mem, f, T0, y));
  TRACE_CHECK(CVodeSStolerances(cvode_mem, SUN_RCONST(0.0), SUN_RCONST(1.0e-6)));
  SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
  TRACE_CHECK(CVodeSetNonlinearSolver(cvode_mem, NLS));
  TRACE_CHECK(CVodeSetMaxNumSteps(cvode_mem, 100000));

  trace_header(2);
  sunrealtype t = T0;
  while (t < TFINAL)
  {
    int flag = CVode(cvode_mem, TFINAL, y, &t, CV_ONE_STEP);
    trace_row(cvode_mem, flag, t, y, 0);
    if (flag < 0) { break; }
  }

  N_VDestroy(y);
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  SUNContext_Free(&sunctx);
  return 0;
}

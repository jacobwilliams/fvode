/* Step trace: analytic problem dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t),
 * lambda = -100 (examples/cvode/F2003_serial/cv_analytic_fp_f2003.f90),
 * Adams + fixed point, reltol = 1e-6, abstol = 1e-10, CV_NORMAL outputs at
 * t = 1, 2, ..., 10 (one row per output; the row is written after each
 * CVode call, like the example). */
#include <math.h>
#include <sundials/sundials_context.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include "trace_common.h"

#define LAMBDA SUN_RCONST(-100.0)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  sunrealtype* dd = N_VGetArrayPointer(ydot);
  dd[0] = LAMBDA * yd[0] + SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t * t) - LAMBDA * atan(t);
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  TRACE_CHECK(SUNContext_Create(SUN_COMM_NULL, &sunctx));
  N_Vector y = N_VNew_Serial(1, sunctx);
  N_VGetArrayPointer(y)[0] = SUN_RCONST(0.0);

  void* cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  TRACE_CHECK(CVodeInit(cvode_mem, f, SUN_RCONST(0.0), y));
  TRACE_CHECK(CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-10)));
  SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
  TRACE_CHECK(CVodeSetNonlinearSolver(cvode_mem, NLS));

  trace_header(1);
  sunrealtype t = SUN_RCONST(0.0), tout = SUN_RCONST(0.0);
  for (int i = 1; i <= 10; i++)
  {
    tout += SUN_RCONST(1.0);
    int flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    trace_row(cvode_mem, flag, t, y, 0);
    if (flag < 0) { break; }
  }

  N_VDestroy(y);
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  SUNContext_Free(&sunctx);
  return 0;
}

/* Shared driver for the CVodeResizeHistory step traces, after the SUNDIALS
 * unit test test/unit_tests/cvode/CXX_serial/cv_test_resize_history.cpp:
 * y' = lambda*(y - atan(t)) + 1/(1+t^2), lambda = -1, y(0) = 0 (solution
 * atan(t)), rtol = 1e-5, atol = 1e-8, 60 steps in CV_ONE_STEP mode towards
 * tf = 10. After every step the (t, y, f) history is saved and passed to
 * CVodeResizeHistory. Include after defining:
 *   RESIZE_MODE    1: resize without changing the size; 2: grow by one
 *                  element each step
 *   RESIZE_METHOD  CV_ADAMS or CV_BDF
 *   RESIZE_NEWTON  0: fixed point (m = 2, 10 iterations); 1: Newton with a
 *                  band (diagonal) linear solver and a user Jacobian
 * One CSV row per step; only y1 is written (all components are equal). */
#include <math.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include "trace_common.h"

static sunrealtype lambda = SUN_RCONST(-1.0);

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  sunrealtype phi_t       = atan(t);
  sunrealtype phi_prime_t = SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t * t);
  sunrealtype* yd         = N_VGetArrayPointer(y);
  sunrealtype* fd         = N_VGetArrayPointer(f);
  for (sunindextype i = 0; i < N_VGetLength(y); i++)
  {
    fd[i] = lambda * (yd[i] - phi_t) + phi_prime_t;
  }
  return 0;
}

static int ode_rhs_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                       void* user_data, N_Vector tmp1, N_Vector tmp2,
                       N_Vector tmp3)
{
  sunrealtype* Jd = SUNBandMatrix_Data(J);
  for (sunindextype i = 0; i < N_VGetLength(y); i++) { Jd[i] = lambda; }
  return 0;
}

static void save_history(sunrealtype t_n, N_Vector y_n, sunrealtype* t_hist,
                         N_Vector* y_hist, N_Vector* f_hist, int hist_size,
                         int step)
{
  int i_start = hist_size - 2;
  if (step < hist_size) { i_start = step - 1; }
  for (int i = i_start; i >= 0; i--) { t_hist[i + 1] = t_hist[i]; }
  t_hist[0] = t_n;
  for (int i = i_start; i >= 0; i--)
  {
    N_VScale(SUN_RCONST(1.0), y_hist[i], y_hist[i + 1]);
    N_VScale(SUN_RCONST(1.0), f_hist[i], f_hist[i + 1]);
  }
  N_VScale(SUN_RCONST(1.0), y_n, y_hist[0]);
  ode_rhs(t_n, y_n, f_hist[0], NULL);
}

static void resize_history(sunrealtype* t_hist, N_Vector* y_hist,
                           N_Vector* f_hist, int hist_size, SUNContext sunctx)
{
  int new_size = N_VGetLength(y_hist[0]) + 1;
  for (int i = 0; i < hist_size; i++)
  {
    N_Vector new_vec      = N_VNew_Serial(new_size, sunctx);
    sunrealtype* old_data = N_VGetArrayPointer(y_hist[i]);
    sunrealtype* new_data = N_VGetArrayPointer(new_vec);
    for (int j = 0; j < new_size; j++) { new_data[j] = old_data[0]; }
    N_VDestroy(y_hist[i]);
    y_hist[i] = new_vec;
  }
  for (int i = 0; i < hist_size; i++)
  {
    N_VDestroy(f_hist[i]);
    f_hist[i] = N_VClone(y_hist[i]);
    ode_rhs(t_hist[i], y_hist[i], f_hist[i], NULL);
  }
}

static void attach_solvers(void* cvode_mem, N_Vector y, SUNContext sunctx,
                           SUNNonlinearSolver* NLS, SUNLinearSolver* LS,
                           SUNMatrix* A)
{
  if (RESIZE_NEWTON == 0)
  {
    if (*NLS) { SUNNonlinSolFree(*NLS); }
    *NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
    TRACE_CHECK(CVodeSetNonlinearSolver(cvode_mem, *NLS));
    TRACE_CHECK(CVodeSetMaxNonlinIters(cvode_mem, 10));
  }
  else
  {
    if (*A) { SUNMatDestroy(*A); }
    if (*LS) { SUNLinSolFree(*LS); }
    *A  = SUNBandMatrix(N_VGetLength(y), 0, 0, sunctx);
    *LS = SUNLinSol_Band(y, *A, sunctx);
    TRACE_CHECK(CVodeSetLinearSolver(cvode_mem, *LS, *A));
    TRACE_CHECK(CVodeSetJacFn(cvode_mem, ode_rhs_jac));
  }
}

int main(void)
{
  SUNContext sunctx;
  TRACE_CHECK(SUNContext_Create(SUN_COMM_NULL, &sunctx));

  N_Vector y  = N_VNew_Serial(1, sunctx);
  N_Vector y1 = N_VNew_Serial(1, sunctx); /* y(1) for the trace rows */
  N_VConst(SUN_RCONST(0.0), y);           /* atan(0) */

  void* cvode_mem = CVodeCreate(RESIZE_METHOD, sunctx);
  TRACE_CHECK(CVodeInit(cvode_mem, ode_rhs, SUN_RCONST(0.0), y));
  TRACE_CHECK(CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-8)));

  SUNNonlinearSolver NLS = NULL;
  SUNLinearSolver LS     = NULL;
  SUNMatrix A            = NULL;
  attach_solvers(cvode_mem, y, sunctx, &NLS, &LS, &A);

  int hist_size        = (RESIZE_METHOD == CV_BDF) ? 5 : 12;
  sunrealtype t_hist[12] = {0};
  N_Vector* y_hist     = N_VCloneVectorArray(hist_size, y);
  N_Vector* f_hist     = N_VCloneVectorArray(hist_size, y);
  N_VScale(SUN_RCONST(1.0), y, y_hist[0]);
  ode_rhs(t_hist[0], y_hist[0], f_hist[0], NULL);

  trace_header(1);
  sunrealtype t_ret;
  for (int i = 1; i <= 60; i++)
  {
    int flag = CVode(cvode_mem, SUN_RCONST(10.0), y, &t_ret, CV_ONE_STEP);
    N_VGetArrayPointer(y1)[0] = N_VGetArrayPointer(y)[0];
    trace_row(cvode_mem, flag, t_ret, y1, RESIZE_NEWTON);
    if (flag < 0) { break; }

    save_history(t_ret, y, t_hist, y_hist, f_hist, hist_size, i);
    if (RESIZE_MODE == 2) { resize_history(t_hist, y_hist, f_hist, hist_size, sunctx); }
    int n_hist = (i < hist_size) ? i + 1 : hist_size;
    TRACE_CHECK(CVodeResizeHistory(cvode_mem, t_hist, y_hist, f_hist, n_hist, n_hist));
    if (RESIZE_MODE == 2)
    {
      N_VDestroy(y);
      y = N_VClone(y_hist[0]);
      attach_solvers(cvode_mem, y, sunctx, &NLS, &LS, &A);
    }
  }

  N_VDestroy(y);
  N_VDestroy(y1);
  N_VDestroyVectorArray(y_hist, hist_size);
  N_VDestroyVectorArray(f_hist, hist_size);
  CVodeFree(&cvode_mem);
  if (NLS) { SUNNonlinSolFree(NLS); }
  if (LS) { SUNLinSolFree(LS); }
  if (A) { SUNMatDestroy(A); }
  SUNContext_Free(&sunctx);
  return 0;
}

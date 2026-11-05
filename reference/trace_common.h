/* Shared step-trace output for the reference drivers (trace_*.c).
 *
 * After every internal step (CV_ONE_STEP mode) one CSV row is written with the
 * integrator state that FVODE must reproduce. The FVODE trace drivers write the
 * same columns, and reference/compare.py diffs the two files.
 *
 * Integer columns must match exactly; real columns are printed with 17
 * significant digits so a bitwise comparison is possible. */
#ifndef FVODE_TRACE_COMMON_H
#define FVODE_TRACE_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#define TRACE_CHECK(call)                                              \
  do {                                                                 \
    int trace_rv_ = (call);                                            \
    if (trace_rv_ < 0) {                                               \
      fprintf(stderr, "%s failed with %d\n", #call, trace_rv_);        \
      exit(1);                                                         \
    }                                                                  \
  } while (0)

static void trace_header(int neq)
{
  printf("nst,t,hu,qu,hcur,qcur,nfe,nni,nncf,netf,nsetups,njev,nfels,flag");
  for (int i = 0; i < neq; i++) { printf(",y%d", i + 1); }
  printf("\n");
}

static void trace_row(void* cvode_mem, int flag, sunrealtype t, N_Vector y,
                      int has_ls)
{
  long int nst, nfe, nni, nncf, netf, nsetups, njev = 0, nfels = 0;
  int qu, qcur;
  sunrealtype hu, hcur;

  TRACE_CHECK(CVodeGetNumSteps(cvode_mem, &nst));
  TRACE_CHECK(CVodeGetLastStep(cvode_mem, &hu));
  TRACE_CHECK(CVodeGetLastOrder(cvode_mem, &qu));
  TRACE_CHECK(CVodeGetCurrentStep(cvode_mem, &hcur));
  TRACE_CHECK(CVodeGetCurrentOrder(cvode_mem, &qcur));
  TRACE_CHECK(CVodeGetNumRhsEvals(cvode_mem, &nfe));
  TRACE_CHECK(CVodeGetNumNonlinSolvIters(cvode_mem, &nni));
  TRACE_CHECK(CVodeGetNumNonlinSolvConvFails(cvode_mem, &nncf));
  TRACE_CHECK(CVodeGetNumErrTestFails(cvode_mem, &netf));
  TRACE_CHECK(CVodeGetNumLinSolvSetups(cvode_mem, &nsetups));
  if (has_ls)
  {
    TRACE_CHECK(CVodeGetNumJacEvals(cvode_mem, &njev));
    TRACE_CHECK(CVodeGetNumLinRhsEvals(cvode_mem, &nfels));
  }

  printf("%ld,%.17e,%.17e,%d,%.17e,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%d", nst, t,
         hu, qu, hcur, qcur, nfe, nni, nncf, netf, nsetups, njev, nfels, flag);
  sunrealtype* yd = N_VGetArrayPointer(y);
  for (sunindextype i = 0; i < N_VGetLength(y); i++) { printf(",%.17e", yd[i]); }
  printf("\n");
}

#endif

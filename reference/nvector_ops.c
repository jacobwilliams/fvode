/* Bitwise reference for FVODE's vector operations (fvode_nvector).
 *
 * Builds deterministic, awkward inputs (mixed magnitudes, zeros, constraint
 * values in {0, +-1, +-2}), applies every serial N_Vector operation that CVODE
 * uses -- through the public N_V* API, so the default (non-fused) code paths
 * are exercised -- and prints inputs and results to stdout.
 *
 * Output format, one block per array:
 *     <name> <n>
 *     <value>          (n lines, %.17e; booleans are 1/0)
 *
 * test/test_nvector.f90 reads the inputs, repeats each operation and requires
 * bitwise equality. Built and run by reference/run_reference.sh. */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>

#define N 257 /* long enough for summation order to matter */

static SUNContext ctx;

static unsigned long long lcg_state = 12345u;
static double lcg(void) /* uniform in [-1, 1) */
{
  lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
  return (double)(lcg_state >> 11) / 4503599627370496.0 - 1.0;
}

static void put(const char* name, N_Vector v)
{
  double* d = N_VGetArrayPointer(v);
  printf("%s %d\n", name, N);
  for (int i = 0; i < N; i++) { printf("%.17e\n", d[i]); }
}

static void put_s(const char* name, double s)
{
  printf("%s 1\n%.17e\n", name, s);
}

static N_Vector vec(void) { return N_VNew_Serial(N, ctx); }

static N_Vector copy(N_Vector x)
{
  N_Vector v = vec();
  N_VScale(1.0, x, v);
  return v;
}

int main(void)
{
  SUNContext_Create(SUN_COMM_NULL, &ctx);

  /* ---- inputs ---- */
  N_Vector x = vec(), y = vec(), w = vec(), c = vec(), d = vec(), xz = vec();
  N_Vector xs = vec(); /* satisfies the constraints c */
  double *xd = N_VGetArrayPointer(x), *yd = N_VGetArrayPointer(y);
  double *wd = N_VGetArrayPointer(w), *cd = N_VGetArrayPointer(c);
  double *dd = N_VGetArrayPointer(d), *xzd = N_VGetArrayPointer(xz);
  double* xsd = N_VGetArrayPointer(xs);
  for (int i = 0; i < N; i++)
  {
    xd[i]  = lcg() * pow(10.0, (double)((i % 9) - 4));
    yd[i]  = lcg() * pow(10.0, (double)((i % 7) - 3));
    wd[i]  = 1.0 / (1.0e-4 * fabs(lcg()) + 1.0e-8);
    cd[i]  = (double)((i % 5) - 2); /* -2,-1,0,1,2 */
    dd[i]  = (i % 4 == 0) ? 0.0 : lcg();
    xzd[i] = (i % 6 == 0) ? 0.0 : xd[i];
    /* satisfy c: c=2 -> >0, c=1 -> >=0, c=-1 -> <=0, c=-2 -> <0 */
    xsd[i] = (cd[i] == 0.0) ? xd[i] : (cd[i] > 0 ? 1.0 : -1.0) * (fabs(xd[i]) + 1.0e-3);
    if (i % 10 == 3 && fabs(cd[i]) == 1.0) { xsd[i] = 0.0; } /* boundary, allowed */
  }
  put("in_x", x); put("in_y", y); put("in_w", w); put("in_c", c);
  put("in_d", d); put("in_xz", xz); put("in_xs", xs);

  /* ---- N_VLinearSum: every special case, plain and in place ---- */
  const double ab[][2] = {{1, 1},     {1, -1},   {-1, 1},     {1, 0.3},
                          {0.3, 1},   {-1, 0.7}, {0.7, -1},   {0.6, 0.6},
                          {0.6, -0.6}, {0.4, 1.7}, {0, 0.5},  {0.5, 0},
                          {0, 0},     {-1, -1},  {-0.3, -0.3}};
  const int nab = sizeof(ab) / sizeof(ab[0]);
  printf("linsum_coef %d\n", 2 * nab);
  for (int k = 0; k < nab; k++) { printf("%.17e\n%.17e\n", ab[k][0], ab[k][1]); }
  for (int k = 0; k < nab; k++)
  {
    char name[64];
    N_Vector z = vec(), yy = copy(y), xx = copy(x);
    N_VLinearSum(ab[k][0], x, ab[k][1], y, z);
    snprintf(name, sizeof name, "linsum_%d", k); put(name, z);
    N_VLinearSum(ab[k][0], x, ab[k][1], yy, yy); /* z == y */
    snprintf(name, sizeof name, "linsum_zy_%d", k); put(name, yy);
    N_VLinearSum(ab[k][0], xx, ab[k][1], y, xx); /* z == x */
    snprintf(name, sizeof name, "linsum_zx_%d", k); put(name, xx);
    N_VDestroy(z); N_VDestroy(yy); N_VDestroy(xx);
  }

  /* ---- N_VScale ---- */
  const double sc[] = {1.0, -1.0, 0.37, 0.0};
  const int nsc = sizeof(sc) / sizeof(sc[0]);
  printf("scale_coef %d\n", nsc);
  for (int k = 0; k < nsc; k++) { printf("%.17e\n", sc[k]); }
  for (int k = 0; k < nsc; k++)
  {
    char name[64];
    N_Vector z = vec(), xx = copy(x);
    N_VScale(sc[k], x, z);
    snprintf(name, sizeof name, "scale_%d", k); put(name, z);
    N_VScale(sc[k], xx, xx);
    snprintf(name, sizeof name, "scale_inplace_%d", k); put(name, xx);
    N_VDestroy(z); N_VDestroy(xx);
  }

  /* ---- elementwise ops ---- */
  N_Vector z = vec();
  N_VConst(0.123, z);        put("const", z);
  N_VProd(x, y, z);          put("prod", z);
  N_VDiv(x, w, z);           put("div", z);
  N_VAbs(x, z);              put("abs", z);
  N_VInv(w, z);              put("inv", z);
  N_VAddConst(x, 0.25, z);   put("addconst", z);
  N_VCompare(0.5, x, z);     put("compare", z);

  N_VConst(-7.0, z); /* elements where xz == 0 must keep this value */
  put_s("invtest_flag", N_VInvTest(xz, z) ? 1.0 : 0.0); put("invtest", z);
  N_VConst(-7.0, z);
  put_s("invtest_ok_flag", N_VInvTest(w, z) ? 1.0 : 0.0); put("invtest_ok", z);

  put_s("constrmask_flag", N_VConstrMask(c, x, z) ? 1.0 : 0.0); put("constrmask", z);
  put_s("constrmask_ok_flag", N_VConstrMask(c, xs, z) ? 1.0 : 0.0);
  put("constrmask_ok", z);

  /* ---- reductions ---- */
  put_s("dotprod", N_VDotProd(x, y));
  put_s("maxnorm", N_VMaxNorm(x));
  put_s("wrmsnorm", N_VWrmsNorm(x, w));
  put_s("wsqrsumlocal", N_VWSqrSumLocal(x, w));
  put_s("min", N_VMin(x));
  put_s("minquotient", N_VMinQuotient(x, d));
  N_Vector zeros = vec(); N_VConst(0.0, zeros);
  put_s("minquotient_allzero", N_VMinQuotient(x, zeros));

  /* ---- fused / vector array ops (generic fallbacks) ---- */
  double cc[4] = {1.0, -0.731, -1.0, 2.5e-3};
  N_Vector X[4] = {x, y, w, d};
  printf("lincomb_coef 4\n");
  for (int k = 0; k < 4; k++) { printf("%.17e\n", cc[k]); }
  N_VLinearCombination(4, cc, X, z); put("lincomb", z);
  {
    N_Vector xx = copy(x);
    N_Vector Xa[4] = {xx, y, w, d}; /* z == X[0] (as in cvPredict-style calls) */
    N_VLinearCombination(4, cc, Xa, xx); put("lincomb_alias", xx);
    N_VDestroy(xx);
  }

  double aa[3] = {0.25, 1.0, -1.0};
  N_Vector Y[3] = {y, w, d};
  N_Vector Z[3] = {vec(), vec(), vec()};
  printf("scaleaddmulti_coef 3\n");
  for (int k = 0; k < 3; k++) { printf("%.17e\n", aa[k]); }
  N_VScaleAddMulti(3, aa, x, Y, Z);
  put("scaleaddmulti_1", Z[0]); put("scaleaddmulti_2", Z[1]); put("scaleaddmulti_3", Z[2]);

  double dots[3];
  N_VDotProdMulti(3, x, Y, dots);
  printf("dotprodmulti 3\n%.17e\n%.17e\n%.17e\n", dots[0], dots[1], dots[2]);

  double sv[3] = {0.5, -1.0, 3.0};
  printf("scalevecarray_coef 3\n");
  for (int k = 0; k < 3; k++) { printf("%.17e\n", sv[k]); }
  {
    N_Vector A[3] = {copy(x), copy(y), copy(w)}; /* in place, as in cvRescale */
    N_VScaleVectorArray(3, sv, A, A);
    put("scalevecarray_1", A[0]); put("scalevecarray_2", A[1]); put("scalevecarray_3", A[2]);
  }

  /* ---- sundials_math ---- */
  const double pb[] = {1.1, 0.97, -1.3, 1.0e-3, 7.0};
  const int pe[]    = {0, 1, 5, 13, -3, -1};
  printf("powi_base 5\n");
  for (int i = 0; i < 5; i++) { printf("%.17e\n", pb[i]); }
  printf("powi_exp 6\n");
  for (int j = 0; j < 6; j++) { printf("%d\n", pe[j]); }
  printf("powi 30\n");
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 6; j++) { printf("%.17e\n", SUNRpowerI(pb[i], pe[j])); }
  const double rb[] = {0.3, 1.7, 1.0e-6, 123.456};
  const double re[] = {1.0 / 3.0, 0.2, 1.0 / 6.0, 0.5, 2.0};
  printf("powr 20\n");
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 5; j++) { printf("%.17e\n", SUNRpowerR(rb[i], re[j])); }

  return 0;
}

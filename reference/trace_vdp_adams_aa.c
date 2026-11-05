/* Step trace: van der Pol, Adams + Anderson-accelerated fixed point (m = 2,
 * up to 5 iterations, so the acceleration subspace fills and is recycled). */
#define VDP_AA_M     2
#define VDP_MAXITERS 5
#define VDP_BETA     SUN_RCONST(1.0)
#include "vdp_adams_driver.h"

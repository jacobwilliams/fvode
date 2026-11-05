/* Step trace: van der Pol, Adams + damped Anderson acceleration (m = 3,
 * up to 5 iterations, beta = 0.75). */
#define VDP_AA_M     3
#define VDP_MAXITERS 5
#define VDP_BETA     SUN_RCONST(0.75)
#include "vdp_adams_driver.h"

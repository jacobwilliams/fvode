/* Step trace: CVodeResizeHistory after every step (resize_driver.h):
 * grow by one element each step, CV_ADAMS, fixed point. */
#define RESIZE_MODE   2
#define RESIZE_METHOD CV_ADAMS
#define RESIZE_NEWTON 0
#include "resize_driver.h"

/* Step trace: CVodeResizeHistory after every step (resize_driver.h):
 * grow by one element each step, CV_BDF, Newton + band. */
#define RESIZE_MODE   2
#define RESIZE_METHOD CV_BDF
#define RESIZE_NEWTON 1
#include "resize_driver.h"

/* Step trace: CVodeResizeHistory after every step (resize_driver.h):
 * same size, CV_BDF, Newton + band. */
#define RESIZE_MODE   1
#define RESIZE_METHOD CV_BDF
#define RESIZE_NEWTON 1
#include "resize_driver.h"

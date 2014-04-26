#ifndef __SEAM_LCI
#define __SEAM_LCI

#include "utils.h"

/* local cubic interpolation */

void lci(int nx,       /* output length                     */
	 int ny,       /* input length                      */
	 float dy,     /* input step,                       */
	 float oy,     /* input origin                      */
	 float * yx,   /* abscissa map (length = nx)        */
	 float * fy,   /* input data, fcn y (length = ny)   */
	 float * fx    /* output data, fcn x (length = nx)  */
	 );

#endif

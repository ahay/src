#ifndef _vp_stdplot_h
#define _vp_stdplot_h

#include <rsf.h>

void vp_stdplot_init (float min1, float max1, float min2, float max2,
		      bool transp, bool xreverse, bool yreverse, bool pad);
void vp_frame_init (sf_file in, const char *where);
void vp_frame(void);
void vp_minmax(float *x1, float* y1, float* x2, float* y2);

#endif

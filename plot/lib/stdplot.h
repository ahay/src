#ifndef _vp_stdplot_h
#define _vp_stdplot_h

#include <rsf.h>

void vp_stdplot_init (float min1, float max1, float min2, float max2,
		      bool transp, bool xreverse, bool yreverse, bool pad);
void vp_frame_init (sf_file in, const char *where);
void vp_barframe_init (float min, float max);
void vp_simpleframe(void);
void vp_frame(void);
void vp_barraster (int nbuf, unsigned char** buf);
void vp_simplebarframe (void);
void vp_barframe(void);
void vp_barline (int nc, float *c, float cmin, float cmax);
void vp_framenum(float num);

#endif

/* 	$Id: stdplot.h,v 1.11 2003/10/18 18:21:45 fomels Exp $	 */

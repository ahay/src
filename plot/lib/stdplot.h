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
void vp_cubeplot_init (int n1pix, int n2pix, int n1front, int n2front,
		       bool flat);
void vp_cuberaster(int n1, int n2, unsigned char** buf,
		   int f1, int f2, int f3);

#endif

/* 	$Id: stdplot.h,v 1.12 2004/03/29 08:00:12 fomels Exp $	 */

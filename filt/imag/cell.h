#ifndef _cell__h
#define _cell__h

#include <rsf.h>

bool cell_snap (float *z, int *iz, float eps);
float cell_p2a (float* p);

/* second-order symplectic */
void cell_intersect (float a, float x, float dy, float p, 
		     float *sx, int *jx);
float cell_update1 (int dim, float s, float v, float *p, const float *g);
float cell_update2 (int dim, float s, float v, float *p, const float *g);

/* first-order symplectic */
void cell1_intersect (float a, float x, float dy, float p, float *sx, int *jx);
float cell1_update1 (int dim, float s, float v, float *p, const float *g);
float cell1_update2 (int dim, float s, float v, float *p, const float *g);

/* first-order nonsymplectic */

void cell11_intersect2 (float a, float da, 
			const float *p, const float *g, float *sp, int *jp);
float cell11_update1 (int dim, float s, float v, float *p, const float *g);
float cell11_update2 (int dim, float s, float v, float *p, const float *g);


#endif

/* 	$Id: cell.h,v 1.3 2003/09/30 14:30:52 fomels Exp $	 */

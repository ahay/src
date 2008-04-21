/* Cell ray tracing */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>

#include "_bool.h"
/*^*/

#include "cell.h"
#include "quadratic.h"
#include "c99.h"
#include "_defs.h"

static float pg;

void sf_cell1_intersect (float a, float x, float dy, float p, 
			 float *sx, int *jx)
/*< intersecting a straight ray with cell boundaries >*/
{
    float si; 
    int i;
    
    *sx = SF_HUGE;
    if (fabsf(p) > SF_EPS) {
	for (i = 0; i < 3; i++) {
	    si = dy*(1-i-x)/p;
	    if (si > 0. && si < *sx) {
		*sx = si;
		*jx = 1-i;
	    }
	}
    }
}

float sf_cell1_update1 (int dim, float s, float v, float *p, const float *g) 
/*< symplectic first-order: step 1 >*/
{
    int i;
    
    for (i=0; i < dim; i++) {
	p[i] = v*p[i];
    }
    
    return (0.5*v*v*s*(1. + s*pg));
}


float sf_cell1_update2 (int dim, float s, float v, float *p, const float *g) 
/*< symplectic first-order: step 2 >*/
{
    int i;
    float d;
    
    d = 0.;
    for (i=0; i < dim; i++) {
    	p[i] += g[i]*s*v;
	d += p[i]*p[i];
    }
    d = sqrt(d);
    pg = 0.;
    for (i=0; i < dim; i++) {
	p[i] /= d;
	pg += p[i]*g[i];
    }
    pg /= 3.;
    
    return (0.5*v*v*s*(1. - s*pg));
}

void sf_cell11_intersect2 (float a, float da, 
			const float* p, const float* g, 
			float *sp, int *jp)
/*< intersecting a straight ray with cell boundaries >*/
{
    float den, s1, p1[2];
    
    den = g[0]*g[0]+g[1]*g[1];
    
    *sp = SF_HUGE;
    
    if (den < SF_EPS) return;
    
    p1[0] = -cosf(a+da)-p[0];
    p1[1] = sinf(a+da)-p[1];
    s1 = (p1[0]*g[0]+p1[1]*g[1])/den;
    if (s1 > 0. && s1 < *sp) {
	*sp = s1;
	*jp = 1;
    }
    
    p1[0] = -cosf(a-da)-p[0];
    p1[1] = sinf(a-da)-p[1];
    s1 = (p1[0]*g[0]+p1[1]*g[1])/den;
    if (s1 > 0. && s1 < *sp) {
	*sp = s1;
	*jp = -1;
    }
}

float sf_cell11_update1 (int dim, float s, float v, float *p, const float *g) 
/*< nonsymplectic first-order: step 1 >*/
{
    int i;
    
    for (i=0; i < dim; i++) {
	p[i] = v*(p[i] + g[i]*s);
    }
    
    return (0.5*v*v*s*(1. + s*pg));
}

float sf_cell11_update2 (int dim, float s, float v, float *p, const float *g) 
/*< nonsymplectic first-order: step 2 >*/
{
    int i;
    float d;
    
    d = 0.;
    for (i=0; i < dim; i++) {
	d += p[i]*p[i];
    }
    d = sqrt(d);
    pg = 0.;
    for (i=0; i < dim; i++) {
	p[i] /= d;
	pg += p[i]*g[i];
    }
    pg /= 3.;
    
    return (0.5*v*v*s*(1. - s*pg));
}

void sf_cell_intersect (float a, float x, float dy, float p, 
			float *sx, int *jx)
/*< intersecting a parabolic ray with cell boundaries >*/
{
    float si; 
    int i;
    
    *sx = SF_HUGE;
    for (i = 0; i < 3; i++) {
	si = sf_quadratic_solve (0.5*a,0.5*p,dy*(x+i-1));
	if (si < *sx) {
	    *sx = si;
	    *jx = 1-i;
	}
    }
}

bool sf_cell_snap (float *z, int *iz, float eps)
/*< round to the nearest boundary >*/
{
    if (*z > 1.-eps) {
	*z = 0.;
	(*iz)++; 
	return true;
    } else if (fabsf(*z) < eps) {
	*z = 0.;
	return true;
    } else if (*z < eps-1.) {
	*z = 0.;
	(*iz)--;
	return true;
    } else if (*z < 0.) {
	(*z) += 1.;
	(*iz)--;
	return false;
    } else {
	return false;
    }
}

float sf_cell_update1 (int dim, float s, float v, float *p, const float *g) 
/*< symplectic second-order: step 1 >*/
{
    int i;
    
    for (i=0; i < dim; i++) {
	p[i] = v*(p[i] + 0.5*g[i]*s);
    }
    
    return (0.5*v*v*s*(1. + s*pg));
}


float sf_cell_update2 (int dim        /* number of dimensions */, 
		       float s        /* sigma */, 
		       float v        /* slowness */, 
		       float *p       /* in - ?, out - direction */, 
		       const float *g /* slowness gradient */) 
/*< symplectic second-order: step 2 >*/
{
    int i;
    float d;
    
    d = 0.;
    for (i=0; i < dim; i++) {
    	p[i] += 0.5*g[i]*s*v;
	d += p[i]*p[i];
    }
    d = sqrt(d);
    pg = 0.;
    for (i=0; i < dim; i++) {
	p[i] /= d;
	pg += p[i]*g[i];
    }
    pg /= 3.;
    
    return (0.5*v*v*s*(1. - s*pg));
}

float sf_cell_p2a (float* p)
/*< convert ray parameter to angle >*/
{
    float a;
    
    if (p[0] <= 0.) {
	if (p[1] >= 1.-10.*SF_EPS) {
	    a = asinf(1.);
	} else if (p[1] <= -1.+10.*SF_EPS) {
	    a = asinf(-1.);
	} else {
	    a = asinf(p[1]);
	}
    } else {
	if (p[1] >= 0.) {
	    if (p[1] >= 1.-10.*SF_EPS) {
		a = asinf(1.);
	    } else {
		a = SF_PI - asinf(p[1]);
	    }
	} else {
	    if (p[1] <= -1.+10.*SF_EPS) {
		a = asinf(-1.);
	    } else {
		a = -SF_PI - asinf(p[1]);
	    }
	}
    }
    
    return a;
}

/* 	$Id$	 */

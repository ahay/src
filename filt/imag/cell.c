#include <math.h>

#include <rsf.h>

#include "cell.h"
#include "quadratic.h"

static float pg;

void cell1_intersect (float a, float x, float dy, float p, 
		      float *sx, int *jx)
{
    float si; 
    int i;
    
    *sx = HUG;
    if (fabsf(p) > EPS) {
	for (i = 0; i < 3; i++) {
	    si = dy*(1-i-x)/p;
	    if (si > 0. && si < *sx) {
		*sx = si;
		*jx = 1-i;
	    }
	}
    }
}


/* symplectic first-order */

float cell1_update1 (int dim, float s, float v, float *p, const float *g) 
{
    int i;
    
    for (i=0; i < dim; i++) {
	p[i] = v*p[i];
    }
    
    return (0.5*v*v*s*(1. + s*pg));
}


float cell1_update2 (int dim, float s, float v, float *p, const float *g) 
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

void cell11_intersect2 (float a, float da, 
			const float* p, const float* g, 
			float *sp, int *jp)
{
    float den, s1, p1[2];
    
    den = g[0]*g[0]+g[1]*g[1];
    
    *sp = HUG;
    
    if (den < EPS) return;
    
    p1[0] = -cosf(a+da)-p[0];
    p1[1] = sinf(a+da)-p[1];
    s1 = (p1[0]*g[0]+p1[1]*g[1])/den;
    if (s1 > 0. && s1 < *sp) {
	*sp = s1;
	*jp = 1.;
    }
    
    p1[0] = -cosf(a-da)-p[0];
    p1[1] = sinf(a-da)-p[1];
    s1 = (p1[0]*g[0]+p1[1]*g[1])/den;
    if (s1 > 0. && s1 < *sp) {
	*sp = s1;
	*jp = -1.;
    }
}


/* nonsymplectic first-order */

float cell11_update1 (int dim, float s, float v, float *p, const float *g) 
{
    int i;
    
    for (i=0; i < dim; i++) {
	p[i] = v*(p[i] + g[i]*s);
    }
    
    return (0.5*v*v*s*(1. + s*pg));
}


float cell11_update2 (int dim, float s, float v, float *p, const float *g) 
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

void cell_intersect (float a, float x, float dy, float p, 
		     float *sx, int *jx)
{
    float si; 
    int i;
    
    *sx = HUG;
    for (i = 0; i < 3; i++) {
	si = quadratic_solve (0.5*a,0.5*p,dy*(x+i-1));
	if (si < *sx) {
	    *sx = si;
	    *jx = 1-i;
	}
    }
}

bool cell_snap (float *z, int *iz, float eps)
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

float cell_update1 (int dim, float s, float v, float *p, const float *g) 
{
    int i;
    
    for (i=0; i < dim; i++) {
	p[i] = v*(p[i] + 0.5*g[i]*s);
    }
    
    return (0.5*v*v*s*(1. + s*pg));
}


float cell_update2 (int dim, float s, float v, float *p, const float *g) 
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

float cell_p2a (float* p)
{
    float a;
    
    if (p[0] <= 0.) {
	if (p[1] >= 1.) {
	    a = asinf(1.);
	} else if (p[1] <= -1.) {
	    a = asinf(-1.);
	} else {
	    a = asinf(p[1]);
	}
    } else {
	if (p[1] >= 0.) {
	    if (p[1] >= 1.) {
		a = asinf(1.);
	    } else {
		a = SF_PI - asinf(p[1]);
	    }
	} else {
	    if (p[1] <= -1.) {
		a = asinf(-1.);
	    } else {
		a = -SF_PI - asinf(p[1]);
	    }
	}
    }
    
    return a;
}

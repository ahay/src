/* Inner part of fast marching. */
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

#include <float.h>
#include <math.h>

#include <rsf.h>
/*^*/

#include "neighborsvti.h"

struct Upd {
    double stencil, stencil2, value;
    double delta;
};

static float approx(float cos2, float sx, float sz, float q);
static int update (float value, float value2, int i);
static float qsolve(int i, float *res2); 
static void stencil (float t, float t2, struct Upd *x); 
static bool updaten (int m, float* res, float* res2, struct Upd *v[]);
static void grid (int *i, const int *n);

static int *in, *n, s[3], order;
static float *ttime, *itime, *vx, *vz, *q, rdx[3];
static double vx1, vz1, q1;
static const float big_value = FLT_MAX;
static struct Upd x[3];

void neighbors_init (int *in1     /* status flag [n[0]*n[1]*n[2]] */, 
		     float *rdx1  /* grid sampling [3] */, 
		     int *n1      /* grid samples [3] */, 
		     int order1   /* accuracy order */, 
		     float *time1 /* traveltime [n[0]*n[1]*n[2]] */,
		     float *itime1 /* isotropic traveltime [n[0]*n[1]*n[2]] */)
/*< Initialize >*/
{
    in = in1; ttime = time1; itime = itime1;
    n = n1; order = order1;
    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];
    rdx[0] = 1./(rdx1[0]*rdx1[0]);
    rdx[1] = 1./(rdx1[1]*rdx1[1]);
    rdx[2] = 1./(rdx1[2]*rdx1[2]);
}

int  neighbours(int i) 
/*< Update neighbors of gridpoint i, return number of updated points >*/
{
    int j, k, ix, npoints;
    float res, res2;
    
    npoints = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[k] != SF_IN) {
		res = qsolve(k,&res2);
		npoints += update(res,res2,k);
	    }
	}
	if (ix-1 >= 0  ) {
	    k = i-s[j];
	    if (in[k] != SF_IN) {
		res = qsolve(k,&res2);
		npoints += update(res,res2,k);
	    }
	}
    }
    return npoints;
}

static int update (float value, float value2, int i)
/* update gridpoint i with new value */
{
    if (value < itime[i]) {
	itime[i] = value;
	ttime[i] = value2;
	if (in[i] == SF_OUT) { 
	    in[i] = SF_FRONT;      
	    sf_pqueue_insert (itime+i);
	    return 1;
	}
/*	sf_pqueue_update (&(ttime+i)); */
    }
    
    return 0;
}

static float qsolve(int i, float *res2)
/* find new traveltime at gridpoint i */
{
    int j, k1, k2, ix;
    float a, b, a2, b2, t, res;
    struct Upd *v[3], *xj;

    for (j=0; j<3; j++) {
	ix = (i/s[j])%n[j];
	
	if (ix > 0) { 
	    k1 = i-s[j];
	    a = itime[k1];
	    a2 = ttime[k1];
	} else {
	    a = big_value;
	    a2 = big_value;
	}

	if (ix < n[j]-1) {
	    k2 = i+s[j];
	    b = itime[k2];
	    b2 = ttime[k2];
	} else {
	    b = big_value;
	    b2 = big_value;
	}

	xj = x+j;
	xj->delta = rdx[j];

	if (a < b) {
	    xj->stencil = xj->value = a;
	    xj->stencil2 = a2;
	} else {
	    xj->stencil = xj->value = b;
	    xj->stencil2 = b2;
	}

	if (order > 1) {
	    if (a < b  && ix-2 >= 0) { 
		k1 = i-2*s[j];
		if (in[k1] != SF_OUT && a >= (t=itime[k1]))
		    stencil(t,ttime[k1],xj);
	    }
	    if (a > b && ix+2 <= n[j]-1) { 
		k2 = i+2*s[j];
		if (in[k2] != SF_OUT && b >= (t=itime[k2]))
		    stencil(t,ttime[k2],xj);
	    }
	}
    }

    if (x[0].value <= x[1].value) {
	if (x[1].value <= x[2].value) {
	    v[0] = x; v[1] = x+1; v[2] = x+2;
	} else if (x[2].value <= x[0].value) {
	    v[0] = x+2; v[1] = x; v[2] = x+1;
	} else {
	    v[0] = x; v[1] = x+2; v[2] = x+1;
	}
    } else {
	if (x[0].value <= x[2].value) {
	    v[0] = x+1; v[1] = x; v[2] = x+2;
	} else if (x[2].value <= x[1].value) {
	    v[0] = x+2; v[1] = x+1; v[2] = x;
	} else {
	    v[0] = x+1; v[1] = x+2; v[2] = x;
	}
    }
    
    vx1=vx[i];
    vz1=vz[i];
    q1=q[i];

    if(v[2]->value < big_value) {   /* ALL THREE DIRECTIONS CONTRIBUTE */
	if (updaten(3, &res, res2, v) || 
	    updaten(2, &res, res2, v) || 
	    updaten(1, &res, res2, v)) return res;
    } else if(v[1]->value < big_value) { /* TWO DIRECTIONS CONTRIBUTE */
	if (updaten(2, &res, res2, v) || 
	    updaten(1, &res, res2, v)) return res;
    } else if(v[0]->value < big_value) { /* ONE DIRECTION CONTRIBUTES */
	if (updaten(1, &res, res2, v)) return res;
    }
	
    return big_value;
}

static void stencil (float t, float t2, struct Upd *x)
/* second-order stencil */
{
    x->delta *= 2.25;
    x->stencil  = (4.0*x->value - t)/3.0;
    x->stencil2 = (4.0*x->stencil2 - t2)/3.0;
}

static bool updaten (int m, float* res, float* res2, struct Upd *v[]) 
/* updating */
{
    double a, b, c, discr, t;
    float cos2, v1;
    int j;

    a = b = c = 0.;
    for (j=0; j<m; j++) {
	a += v[j]->delta;
	b += v[j]->stencil*v[j]->delta;
	c += v[j]->stencil*v[j]->stencil*v[j]->delta;
    }
    b /= a;

    discr=b*b+(vz1-c)/a;

    if (discr < 0.) return false;
    
    t = b +sqrt(discr);
    if (t <= v[m-1]->value) return false;

    *res = t;
    
    cos2 = 0.;
    for (j=0; j<m; j++) {
    	if (x == v[j]) {
	    cos2 = t-v[j]->stencil;
	    cos2 = cos2*cos2*v[j]->delta/vz1;
	    break;
	}
    }

    v1 = approx(cos2, vx1, vz1, q1);

    b = c = 0.;
    for (j=0; j<m; j++) {
	b += v[j]->stencil2*v[j]->delta;
	c += v[j]->stencil2*v[j]->stencil2*v[j]->delta;
    }
    b /= a;

    discr=b*b+(v1-c)/a;

    if (discr < 0.) discr=0.;
    
    t = b +sqrt(discr);
   
    *res2 = t;

    return true;
}

static void grid (int *i, const int *n)
/* restrict i[3] to the grid n[3] */
{ 
    int j;

    for (j=0; j < 3; j++) {
	if (i[j] < 0) {
	    i[j]=0;
	} else if (i[j] >= n[j]) {
	    i[j]=n[j]-1;
	}
    }
}

static float approx(float cos2, float sx, float sz, float q) 
/* approximation of slowness squared */
{
    float s2;

    /* elliptical part */
    s2 = cos2*sz + (1.-cos2)*sx;

    if (q != 1.) { /* anelliptical part */
	s2 = 0.5/(1.+q)*((1.+2.*q)*s2 + 
			 sqrtf(s2*s2 + 4.*(q*q-1.)*cos2*(1.-cos2)*sx*sz));
    }
    
    return s2;
}

int nearsource(float* xs   /* source location [3] */, 
	       int* b      /* constant-velocity box around it [3] */, 
	       float* d    /* grid sampling [3] */, 
	       float* vx2  /* horizontal slowness squared [n[0]*n[1]*n[2]] */,
	       float* vz2  /* vertical slowness squared */,
	       float* q2   /* non-ellipticity */,
	       bool *plane /* if plane-wave source */)
/*< initialize the source >*/
{
    int npoints, ic, i, j, is, start[3], endx[3], ix, iy, iz;
    float cos2, vxi, vzi, qi;
    double delta[3], delta2;

    /* initialize everywhere */
    for (i=0; i < n[0]*n[1]*n[2]; i++) {
	in[i] = SF_OUT;
	itime[i] = big_value;
    }

    vx = vx2;
    vz = vz2;
    q  = q2;

    /* Find index of the source location and project it to the grid */
    for (j=0; j < 3; j++) {
	is = xs[j]/d[j]+0.5;
	start[j] = is-b[j]; 
	endx[j]  = is+b[j];
    } 
    
    grid(start, n);
    grid(endx, n);
    
    ic = (start[0]+endx[0])/2 + 
	n[0]*((start[1]+endx[1])/2 +
	      n[1]*(start[2]+endx[2])/2);
    
    vxi = vx[ic];
    vzi = vz[ic];
    qi  = q[ic];

    /* loop in a small box around the source */
    npoints = n[0]*n[1]*n[2];
    for (ix=start[2]; ix <= endx[2]; ix++) {
	for (iy=start[1]; iy <= endx[1]; iy++) {
	    for (iz=start[0]; iz <= endx[0]; iz++) {
		npoints--;
		i = iz + n[0]*(iy + n[1]*ix);

		delta[0] = xs[0]-iz*d[0];
		delta[1] = xs[1]-iy*d[1];
		delta[2] = xs[2]-ix*d[2];

		delta2 = 0.;
		for (j=0; j < 3; j++) {
		    if (!plane[j]) delta2 += delta[j]*delta[j];
		}

		/* analytical formula (Euclid) for isotropic time */ 
		itime[i] = sqrtf(vzi*delta2);
		if (delta2 > 0.) {
		    cos2 = delta[0]*delta[0]/delta2;
		} else {
		    cos2 = 1.;
		}
		/* anisotropic time */
		ttime[i] = sqrtf(approx(cos2,vxi,vzi,qi)*delta2);
		in[i] = SF_IN;
		
		if ((n[0] > 1 && (iz == start[0] || iz == endx[0])) ||
		    (n[1] > 1 && (iy == start[1] || iy == endx[1])) ||
		    (n[2] > 1 && (ix == start[2] || ix == endx[2]))) {
		    sf_pqueue_insert (itime+i);
		}
	    }
	}
    }
    
    return npoints;
}

/* 	$Id: neighbors.c 1140 2005-04-25 06:23:30Z fomels $	 */


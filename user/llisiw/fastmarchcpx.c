/* Fast marching interface for marching from the surface. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include <rsf.h>

#include <assert.h>
#include "fastmarchcpx.h"

struct Upd {
    double stencil, value, delta;
};

static float *o, *v, *d;
static int *n, *in, s[3], order;
static float **x, **xn, **x1;
static int *offsets;
static float *t;

int neighbors_mask(float* tref, bool* known);
void pqueue_insert(float* v1);
float* pqueue_extract(void);
void pqueue_update(int index);
int neighbours(float* time, int i);
int update(float value, float* time, int i);
float qsolve(float* time, int i);
bool updaten(int i, int m, float* res, struct Upd *vv[]);

void fastmarchcpx_init(int *n_in    /* length */, 
		       float *o_in  /* origin */,
		       float *d_in  /* sampling */,
		       int order_in /* order */) 
/*< initailize >*/
{
    int maxband;

    n = n_in;
    o = o_in;
    d = d_in;

    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];

    maxband = 0;
    /*
    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];
    */
    /* NOTE: maxband above might be insufficient for cpxeikonal source wave-front */
    maxband = n[0]*n[1]*n[2];
    
    x = (float **) sf_alloc ((10*maxband+1),sizeof (float *));
    in = sf_intalloc(n[0]*n[1]*n[2]);

    offsets = (int *) sf_alloc (n[0]*n[1]*n[2],sizeof (int));

    order = order_in;
}

void fastmarchcpx(float* time /* time */,
		  float* t0   /* fixed traveltime */,
		  bool* m     /* known mask */,
		  float* v_in /* slowness squared */)
/*< Run fast marching eikonal solver >*/
{
    float *p;
    int npoints, i;

    t = time;
    v = v_in;
    
    xn = x;
    x1 = x+1;

    /* initialize from boundary */
    for (npoints =  neighbors_mask(t0,m);
	 npoints > 0;
	 npoints -= neighbours(t,i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p-t;

	in[i] = SF_IN;
    }
}    

void fastmarchcpx_close(void)
/*< free allocated storage >*/
{
    free(x);
    free(in);
}

int neighbors_mask(float* tref /* reference traveltime */,
		   bool* known /* known mask */)
/* initialize the source using a mask */
{
    int npoints, i, nxy;

    /* total number of points */
    nxy = n[0]*n[1]*n[2];
    npoints = nxy;

    for (i=0; i < nxy; i++) {
	if (known[i]) {
	    in[i] = SF_IN;
	    t[i] = tref[i];

	    pqueue_insert(t+i);

	    npoints--;
	} else {
	    in[i] = SF_OUT;
	    t[i] = SF_HUGE;
	    offsets[i] = -1;
	}
    }
    
    return npoints;
}

void pqueue_insert(float* v1)
/* insert an element (smallest first) */
{
    int newOffset, tIndex;
    float **xi, **xq;
    unsigned int q;
    
    xi = ++xn;

    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x + q;
	if (*v1 > **xq) break;
	*xi = *xq;

	/* now that it moved down, update its offset */
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;

	xi = xq;
    }
    *xi = v1; 
    
    /* now that we moved it far enough up, record the offset */
    newOffset = xi-x;
    tIndex = v1-t;
    offsets[tIndex] = newOffset;
}

float* pqueue_extract(void)
/* extract the smallest element */
{
    int newOffset, tIndex;
    unsigned int c;
    int nn;
    float *vv, *formerlyLast;
    float **xi, **xc;
    
    /* check if queue is empty */
    nn = (int) (xn-x);
    if (nn == 0) return NULL;

    vv = *x1;
    /* label vv to be out by set its offset as -1 */
    tIndex = vv-t;
    offsets[tIndex] = -1;

    *(xi = x1) = formerlyLast = *(xn--);
    /* NOTE: formerlyLast's offset is inconsistent until the very last step */
    nn--;
    for (c = 2; c <= (unsigned int) nn; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) nn && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*formerlyLast <= **xc) break;
	*xi = *xc; 

	/* now that it moved up, update its offset */
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;

	xi = xc;
    }
    *xi = formerlyLast;

    /* now that we moved it far enough down, record the offset */
    newOffset = xi-x;
    tIndex = *xi-t;
    offsets[tIndex] = newOffset;

    return vv;
}

void pqueue_update(int index)
/* restore the heap */
{
    int newOffset, tIndex;
    unsigned int c;
    float **xc, **xi;
    
    c = offsets[index];
    xi = x+c;

    for (c >>= 1; c > 0; c >>= 1) {
	xc = x + c;
	if (t[index] > **xc) break;
	*xi = *xc;
	
	/* now that it moved down, update its offset */
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;
	
	xi = xc; 
    }
    *xi = t+index; 

    /* now that we moved it far enough up, record the offset */
    newOffset = xi-x;
    tIndex = *xi-t;
    offsets[tIndex] = newOffset;
}

int neighbours(float* time, int i) 
/* update neighbors of gridpoint i, return number of updated points */
{
    int j, k, ix, np;
    
    np = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];

	/* try both directions */
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[k] != SF_IN) np += update(qsolve(time,k),time,k);
	}
	if (ix-1 >= 0 ) {
	    k = i-s[j];
	    if (in[k] != SF_IN) np += update(qsolve(time,k),time,k);
	}
    }
    return np;
}

int update(float value, float* time, int i)
/* update gridpoint i with new value and modify wave front */
{
    /* only update when smaller than current value */
    if (value < time[i]) {
	time[i] = value;
	if (in[i] == SF_OUT) { 
	    in[i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
	} else {
	    assert(in[i] == SF_FRONT);
	    pqueue_update(i);
	}
    }
    
    return 0;
}

float qsolve(float* time, int i)
/* find new traveltime at gridpoint i */
{
    int j, k, ix;
    float a, b, res;
    struct Upd *vv[3], xx[3], *xj;

    for (j=0; j<3; j++) {
	ix = (i/s[j])%n[j];
	
	if (ix > 0) { 
	    k = i-s[j];
	    a = time[k];
	} else {
	    a = SF_HUGE;
	}
	
	if (ix < n[j]-1) {
	    k = i+s[j];
	    b = time[k];
	} else {
	    b = SF_HUGE;
	}
	
	xj = xx+j;
	xj->delta = 1./(d[j]*d[j]);

	if (a < b) {
	    xj->stencil = xj->value = a;
	} else {
	    xj->stencil = xj->value = b;
	}

	/* second order local upwind stencil */
	if (order > 1) {
	    if (a < b  && ix-2 >= 0) { 
		k = i-2*s[j];
		if (in[k] != SF_OUT && a >= time[k]) {
		    xj->delta *= 2.25;
		    xj->stencil = (4.0*xj->value - time[k])/3.0;
		}
	    }
	    if (a > b && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (in[k] != SF_OUT && b >= time[k]) {
		    xj->delta *= 2.25;
		    xj->stencil = (4.0*xj->value - time[k])/3.0;
		}
	    }
	}
    }

    if (xx[0].value <= xx[1].value) {
	if (xx[1].value <= xx[2].value) {
	    vv[0] = xx; vv[1] = xx+1; vv[2] = xx+2;
	} else if (xx[2].value <= xx[0].value) {
	    vv[0] = xx+2; vv[1] = xx; vv[2] = xx+1;
	} else {
	    vv[0] = xx; vv[1] = xx+2; vv[2] = xx+1;
	}
    } else {
	if (xx[0].value <= xx[2].value) {
	    vv[0] = xx+1; vv[1] = xx; vv[2] = xx+2;
	} else if (xx[2].value <= xx[1].value) {
	    vv[0] = xx+2; vv[1] = xx+1; vv[2] = xx;
	} else {
	    vv[0] = xx+1; vv[1] = xx+2; vv[2] = xx;
	}
    }

    if(vv[2]->value < SF_HUGE) {   /* update from all three directions */
	if (updaten(i,3,&res,vv) || 
	    updaten(i,2,&res,vv) || 
	    updaten(i,1,&res,vv)) return res;
    } else if(vv[1]->value < SF_HUGE) { /* update from two directions */
	if (updaten(i,2,&res,vv) || 
	    updaten(i,1,&res,vv)) return res;
    } else if(vv[0]->value < SF_HUGE) { /* update from only one direction */
	if (updaten(i,1,&res,vv)) return res;
    }
	
    return SF_HUGE;
}

bool updaten(int i, int m, float* res, struct Upd *vv[])
/* calculate new traveltime */
{
    double a, b, c, discr, t;
    int j;

    if (m == 1) {
	/* special coded to ensure that a one-sided update is always successful */
	t = vv[0]->stencil+sqrt((double)v[i]/vv[0]->delta);
    
    } else{
	/* solve quadratic equation */
	a = b = c = 0.;

	for (j=0; j<m; j++) {
	    a += vv[j]->delta;
	    b += vv[j]->stencil*vv[j]->delta;
	    c += vv[j]->stencil*vv[j]->stencil*vv[j]->delta;
	}
	b /= a;

	discr=b*b+(((double)v[i])-c)/a;

	if (discr < 0.) return false;
    
	t = b + sqrt(discr);
	/* NOTE: change from <= to < for complex eikonal */
	if (t < vv[m-1]->value) return false;
    }

    *res = t;
    return true;
}

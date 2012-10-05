/* Fast marching interface */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "fastmarchtd.h"

struct Upd {
    double stencil, value, delta;
};

static float *o, *v, *vd, *vdd, *d, dt0, dx0, ot0, ox0;
static int *n, *tin, *xin, s[3], order, nt0, nx0;
static float **x, **xn, **x1;
static int *offsets;
static float *t, *t00;

int mask_t0(float* tref, bool* known);
int mask_x0(float* xref, bool* known);
void pqueue_insert(float* v1);
float* pqueue_extract(void);
void pqueue_update(int index);
int neighbours_t0(float* time, int i);
int neighbours_x0(float* time, int i);
int updatet(float value, float* time, int i);
int updatex(float value, float* time, int i);
float qsolve(float* time, int i);
float nsolve(float* time, int i);
bool updaten(int i, int m, float* res, struct Upd *vv[]);
bool newton(int i, int m, float* res, struct Upd *vv[], float res0);
bool bisect(int i, int m, float* res, struct Upd *vv[], float res0);
float bilinear(float tin, float xin, const float* tbl);

void fastmarch_init(int *n_in    /* length */, 
		    float *o_in  /* origin */,
		    float *d_in  /* sampling */,
		    int order_in /* accuracy order */,
		    int nt0_in, int nx0_in,
		    float dt0_in, float dx0_in,
		    float ot0_in, float ox0_in) 
/*< initailize >*/
{
    int maxband;

    /* size of depth-domain */
    n = n_in;
    o = o_in;
    d = d_in;

    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];

    /* size of heap */
    maxband = 0;

    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];
    
    x = (float **) sf_alloc ((10*maxband+1),sizeof (float *));

    tin = sf_intalloc(n[0]*n[1]*n[2]);
    xin = sf_intalloc(n[0]*n[1]*n[2]);

    offsets = (int *) sf_alloc (n[0]*n[1]*n[2],sizeof (int));

    order = order_in;

    /* size of time-domain */
    nt0 = nt0_in; nx0 = nx0_in;
    dt0 = dt0_in; dx0 = dx0_in;
    ot0 = ot0_in; ox0 = ox0_in;
}

void fastmarch_t0(float* t0   /* t0 */,
		  float* tt0  /* known value */,
		  bool* maskt /* known mask */,
		  float* s    /* slowness squared */)
/*< fast marching t0 >*/
{
    float *p;
    int npoints, i;

    t = t0;
    v = s;
    
    xn = x;
    x1 = x+1;

    /* initialize from boundary */
    for (npoints =  mask_t0(tt0,maskt);
	 npoints > 0;
	 npoints -= neighbours_t0(t,i)) {
	/* eikonal with known right-hand side */

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p-t;

	/* break if exceeds t0 range */
	if (t[i] > ot0+(nt0-1)*dt0) break;

	tin[i] = SF_IN;
    }
}

void fastmarch_x0(float* x0    /* x0 */,
		  float* xx0   /* known value */,
		  bool* maskx  /* known mask */,
		  float *t0    /* t0 */,
		  float* s     /* slowness squared */,
		  float* vdix  /* Dix velocity */,
		  float* vdixd /* Dix velocity derivative in x0 */)
/*< fast marching x0 >*/
{
    float *p;
    int npoints, i;

    t = x0;
    t00 = t0;
    v = s;
    vd = vdix;
    vdd = vdixd;

    xn = x;
    x1 = x+1;

    /* initialize from boundary */
    for (npoints =  mask_x0(xx0,maskx);
	 npoints > 0;
	 npoints -= neighbours_x0(t,i)) {
	/* eikonal with unknown right-hand side */

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p-t;

	/* break if exceeds x0 range */
	if (t[i] > ox0+(nx0-1)*dx0) break;

	xin[i] = SF_IN;
    }
}

void fastmarch_close(void)
/*< free allocated storage >*/
{
    free(x);
    free(tin);
    free(xin);
}

int mask_t0(float* tref /* reference traveltime */,
	    bool* known /* known mask */)
/* initialize t0 */
{
    int npoints, i, nxy;
    float t0max;

    t0max = ot0+(nt0-1)*dt0;

    /* total number of points */
    nxy = n[0]*n[1]*n[2];
    npoints = nxy;

    for (i=0; i < nxy; i++) {
	if (known[i] && tref[i] < t0max) {
	    tin[i] = SF_IN;
	    t[i] = tref[i];

	    pqueue_insert(t+i);

	    npoints--;
	} else {
	    tin[i] = SF_OUT;
	    t[i] = SF_HUGE;
	    offsets[i] = -1;
	}
    }
    
    return npoints;
}

int mask_x0(float* xref /* reference traveltime */,
	    bool* known /* known mask */)
/* initialize x0 */
{
    int npoints, i, nxy;
    float t0max;

    t0max = ot0+(nt0-1)*dt0;

    /* total number of points */
    nxy = n[0]*n[1]*n[2];
    npoints = nxy;

    for (i=0; i < nxy; i++) {
	if (t00[i] <= t0max) {
	    if (known[i]) {
		xin[i] = SF_IN;
		t[i] = xref[i];

		pqueue_insert(t+i);

		npoints--;
	    } else {
		xin[i] = SF_OUT;
		t[i] = SF_HUGE;
		offsets[i] = -1;
	    }
	} else {
	    xin[i] = SF_IN;
	    t[i] = SF_HUGE;
	    offsets[i] = -1;

	    npoints--;
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

	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;

	xi = xq;
    }
    *xi = v1; 
    
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
    tIndex = vv-t;
    offsets[tIndex] = -1;

    *(xi = x1) = formerlyLast = *(xn--);

    nn--;
    for (c = 2; c <= (unsigned int) nn; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) nn && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*formerlyLast <= **xc) break;
	*xi = *xc; 

	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;

	xi = xc;
    }
    *xi = formerlyLast;

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
	
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;
	
	xi = xc; 
    }
    *xi = t+index; 

    newOffset = xi-x;
    tIndex = *xi-t;
    offsets[tIndex] = newOffset;
}

int neighbours_t0(float* time, int i) 
/* update neighbors of gridpoint i, return number of updated points */
{
    int j, k, ix, np;
    
    np = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];

	/* try both directions */
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (tin[k] != SF_IN) np += updatet(qsolve(time,k),time,k);
	}
	if (ix-1 >= 0 ) {
	    k = i-s[j];
	    if (tin[k] != SF_IN) np += updatet(qsolve(time,k),time,k);
	}
    }
    return np;
}

int neighbours_x0(float* time, int i) 
/* update neighbors of gridpoint i, return number of updated points */
{
    int j, k, ix, np;
    
    np = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];

	/* try both directions */
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (xin[k] != SF_IN) np += updatex(nsolve(time,k),time,k);
	}
	if (ix-1 >= 0 ) {
	    k = i-s[j];
	    if (xin[k] != SF_IN) np += updatex(nsolve(time,k),time,k);
	}
    }
    return np;
}

int updatet(float value, float* time, int i)
/* update gridpoint i with new value and modify wave front */
{
    /* only update when smaller than current value */
    if (value < time[i]) {
	time[i] = value;
	if (tin[i] == SF_OUT) { 
	    tin[i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
	} else {
	    pqueue_update(i);
	}
    }
    
    return 0;
}

int updatex(float value, float* time, int i)
/* update gridpoint i with new value and modify wave front */
{
    /* only update when smaller than current value */
    if (value < time[i]) {
	time[i] = value;
	if (xin[i] == SF_OUT) { 
	    xin[i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
	} else {
	    pqueue_update(i);
	}
    }
    
    return 0;
}

float qsolve(float* time, int i)
/* eikonal with known right-hand side */
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
		if (tin[k] != SF_OUT && a >= time[k]) {
		    xj->delta *= 2.25;
		    xj->stencil = (4.0*xj->value - time[k])/3.0;
		}
	    }
	    if (a > b && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (tin[k] != SF_OUT && b >= time[k]) {
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
/* solve quadratic equation (analytical) */
{
    double a, b, c, discr, t;
    int j;

    if (m == 1) {
	/* one-sided update is always successful */
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

	if (t <= vv[m-1]->value) return false;
    }

    *res = t;
    return true;
}

float nsolve(float* time, int i)
/* eikonal with unknown right-hand side */
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
		if (xin[k] != SF_OUT && a >= time[k]) {
		    xj->delta *= 2.25;
		    xj->stencil = (4.0*xj->value - time[k])/3.0;
		}
	    }
	    if (a > b && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (xin[k] != SF_OUT && b >= time[k]) {
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
	if (bisect(i,3,&res,vv,vv[2]->value) || 
	    bisect(i,2,&res,vv,vv[1]->value) || 
	    bisect(i,1,&res,vv,vv[0]->value)) return res;
    } else if(vv[1]->value < SF_HUGE) { /* update from two directions */
	if (bisect(i,2,&res,vv,vv[1]->value) || 
	    bisect(i,1,&res,vv,vv[0]->value)) return res;
    } else if(vv[0]->value < SF_HUGE) { /* update from only one direction */
	if (bisect(i,1,&res,vv,vv[0]->value)) return res;
    }
	
    return SF_HUGE;
}

bool newton(int i, int m, float* res, struct Upd *vv[], float res0)
/* solve quadratic equation (Newton's method) */
{
    int j;
    float temp, val, grad;

    val = -res0;
    grad = 1.;
    temp = 0.;

    do {
	/* Newton step */
	temp -= val/grad;

	if (temp < ox0 || temp >= ox0+(nx0-1)*dx0) return false;

	/* evaluate function value */
	val = 0.;
	for (j=0; j < m; j++) {
	    val += vv[j]->delta*(temp-vv[j]->stencil)*(temp-vv[j]->stencil);
	}
	val -= v[i]*bilinear(t00[i],temp,vd)*bilinear(t00[i],temp,vd);
	
	/* evaluate derivative */
	grad = 0.;
	for (j=0; j < m; j++) {
	    grad += 2.*vv[j]->delta*(temp-vv[j]->stencil);
	}
	grad -= 2.*v[i]*bilinear(t00[i],temp,vd)*bilinear(t00[i],temp,vdd);
	
    } while (fabs(val) > 1.e-2 && fabs(grad) > 1.e-2);

    if (fabs(val) <= 1.e-2) {
	*res = temp;
	return true;
    } else {
	return false;
    }
}

bool bisect(int i, int m, float* res, struct Upd *vv[], float res0)
/* solve quadratic equation (bisection method) */
{
    int j, k;
    float lef, mid, rht;
    float vlef, vmid, vrht;

    /* initial left */
    lef = res0; 

    vlef = 0.;
    for (j=0; j < m; j++) {
	vlef += vv[j]->delta*(lef-vv[j]->stencil)*(lef-vv[j]->stencil);
    }
    vlef -= v[i]*bilinear(t00[i],lef,vd)*bilinear(t00[i],lef,vd);

    /* find right (brute search) */
    for (k=0; k < nx0; k++) {
	rht = k*dx0+ox0;

	if (rht <= res0) continue;

	vrht = 0.;
	for (j=0; j < m; j++) {
	    vrht += vv[j]->delta*(rht-vv[j]->stencil)*(rht-vv[j]->stencil);
	}
	vrht -= v[i]*bilinear(t00[i],rht,vd)*bilinear(t00[i],rht,vd);

	if (vlef*vrht > 0.) 
	    continue;
	else
	    break;
    }

    if (vlef*vrht > 0.) return false;

    /* bisection */
    do {
	mid = (lef+rht)/2.;
	
	vmid = 0.;
	for (j=0; j < m; j++) {
	    vmid += vv[j]->delta*(mid-vv[j]->stencil)*(mid-vv[j]->stencil);
	}
	vmid -= v[i]*bilinear(t00[i],mid,vd)*bilinear(t00[i],mid,vd);

	if (fabs(vmid) <= 5.e-3) break;

	if (vlef*vmid > 0.) {
	    lef = mid;
	    vlef = vmid;
	} else {
	    rht = mid;
	    vrht = vmid;
	}

    } while (rht-lef >= 0.1*dx0);
    
    *res = mid;
    return true;
}

float bilinear(float tin, float xin, const float* tbl)
/* bilinear interpolation */
{
    int neart, nearx;
    float distt, distx, val;

    /* find position on grid */
    neart = (tin-ot0)/dt0;
    nearx = (xin-ox0)/dx0;

    if (neart < 0 || neart > nt0-1) sf_error("Interpolation beyond t0 scope.");
    if (nearx < 0 || nearx > nx0-1) sf_error("Interpolation beyond x0 scope.");

    if (neart == nt0-1) {
	neart = nt0-2;
	distt = 1.;
    } else {
	distt = (tin-(ot0+((float)neart)*dt0))/dt0;
    }
    
    if (nearx == nx0-1) {
	nearx = nx0-2;
	distx = 1.;
    } else {
	distx = (xin-(ox0+((float)nearx)*dx0))/dx0;
    }    

    val = (1.-distt)*(1.-distx)*tbl[nearx*nt0+neart]
	 +distt     *(1.-distx)*tbl[nearx*nt0+neart+1]
	 +(1.-distt)*distx     *tbl[(nearx+1)*nt0+neart]
	 +distt     *distx     *tbl[(nearx+1)*nt0+neart+1];
    
    return val;
}

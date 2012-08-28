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

#include "irays.h"

struct Upd {
    double tstencil, tvalue;
    double xstencil, xvalue;
    double delta;
    int direct;
};

static float *o, *v, *d;
static int *n, *in, s[3], order;
static float **x, **xn, **x1;
static int *offsets;
static float *t_0, *x_0;
static int *f_0;

int init_surface(void);
void pqueue_insert(float* v1);
float* pqueue_extract(void);
void pqueue_update(int index);
/*
int neighbours(float* time, float* xpos, int i);
int update(float tvalue, float xvalue, float* time, float* xpos, int i);
void qsolve(float* time, float* xpos, int i, float* tval, float* xval);
bool updaten(int i, int m, float* tres, float* xres, struct Upd *vv[]);
*/
int neighbours(float* time, float* xpos, int* nupd, int i);
int update(float tvalue, float xvalue, int fvalue, float* time, float* xpos, int* nupd, int i);
void qsolve(float* time, float* xpos, int i, float* tval, float* xval, int* fval);
bool updaten(int i, int m, float* tres, float* xres, int* fres, struct Upd *vv[]);

void fastmarch_init(int *n_in    /* length */, 
		    float *o_in  /* origin */,
		    float *d_in  /* sampling */,
		    int order_in /* accuracy order */)
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

    in = sf_intalloc(n[0]*n[1]*n[2]);
    offsets = (int *) sf_alloc (n[0]*n[1]*n[2],sizeof (int));

    order = order_in;
}

void fastmarch(float* t_in /* t0 */,
	       float* x_in /* x0 */,
	       int* f_in   /* f0 */,
	       float* v_in /* slowness squared */)
/*< fast marching >*/
{
    float *p;
    int npoints, i;

    t_0 = t_in;
    x_0 = x_in;
    f_0 = f_in;
    v = v_in;
    
    xn = x;
    x1 = x+1;

    /* initialize from boundary */
    for (npoints =  init_surface();
	 npoints > 0;
	 npoints -= neighbours(t_0,x_0,f_0,i)) {
	/* decrease by number of updated points */

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p-t_0;

	in[i] = SF_IN;
    }
}

void fastmarch_close(void)
/*< free allocated storage >*/
{
    free(x);
    free(in);
    free(offsets);
}

int init_surface(void)
/* initialize surface source */
{
    int npoints, i, j, k, nxy;

    /* total number of points */
    nxy = n[0]*n[1]*n[2];
    npoints = nxy;

    for (k=0; k < n[2]; k++) {
	for (j=0; j < n[1]; j++) {
	    in[k*s[2]+j*s[1]] = SF_IN;
	    t_0[k*s[2]+j*s[1]] = 0.;
	    x_0[k*s[2]+j*s[1]] = o[1]+j*d[1];
	    f_0[k*s[2]+j*s[1]] = 0;

	    pqueue_insert(t_0+k*s[2]+j*s[1]);	    
	    npoints--;

	    for (i=1; i < n[0]; i++) {
		in[k*s[2]+j*s[1]+i] = SF_OUT;
		t_0[k*s[2]+j*s[1]+i] = SF_HUGE;
		x_0[k*s[2]+j*s[1]+i] = SF_HUGE;
		f_0[k*s[2]+j*s[1]+i] = -1;
		
		offsets[k*s[2]+j*s[1]+i] = -1;
	    }	    
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
	tIndex = (*xi)-t_0;
	offsets[tIndex] = newOffset;

	xi = xq;
    }
    *xi = v1; 
    
    newOffset = xi-x;
    tIndex = v1-t_0;
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
    tIndex = vv-t_0;
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
	tIndex = (*xi)-t_0;
	offsets[tIndex] = newOffset;

	xi = xc;
    }
    *xi = formerlyLast;

    newOffset = xi-x;
    tIndex = *xi-t_0;
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
	if (t_0[index] > **xc) break;
	*xi = *xc;
	
	newOffset = xi-x;
	tIndex = (*xi)-t_0;
	offsets[tIndex] = newOffset;
	
	xi = xc; 
    }
    *xi = t_0+index; 

    newOffset = xi-x;
    tIndex = *xi-t_0;
    offsets[tIndex] = newOffset;
}

int neighbours(float* time, float* xpos, int* nupd, int i) 
/* update neighbors of gridpoint i, return number of updated points */
{
    int j, k, ix, np;
    float tupdate, xupdate;
    int fupdate;

    np = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];

	/* try both directions */
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[k] != SF_IN) {
		qsolve(time,xpos,k,&tupdate,&xupdate,&fupdate);
		np += update(tupdate,xupdate,fupdate,time,xpos,nupd,k);
	    }
	}
	if (ix-1 >= 0 ) {
	    k = i-s[j];
	    if (in[k] != SF_IN) {
		qsolve(time,xpos,k,&tupdate,&xupdate,&fupdate);
		np += update(tupdate,xupdate,fupdate,time,xpos,nupd,k);
	    }
	}
    }
    return np;
}

int update(float tvalue, float xvalue, int fvalue, float* time, float* xpos, int* nupd, int i)
/* update gridpoint i with new value and modify wave front */
{
    /* only update when smaller than current value */
    if (tvalue < time[i]) {
	time[i] = tvalue;
	xpos[i] = xvalue;
	nupd[i] = fvalue;
	if (in[i] == SF_OUT) { 
	    in[i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
	} else {
	    pqueue_update(i);
	}
    }
    
    return 0;
}

void qsolve(float* time, float* xpos, int i, float* tval, float* xval, int* fval)
/* solve new traveltime */
{
    int j, k, ix;
    float ta, tb, xa, xb, tres, xres;
    int fres;
    struct Upd *vv[3], xx[3], *xj;

    for (j=0; j<3; j++) {
	ix = (i/s[j])%n[j];
	
	if (ix > 0) { 
	    k = i-s[j];
	    ta = time[k];
	    xa = xpos[k];
	} else {
	    ta = SF_HUGE;
	    xa = SF_HUGE;
	}
	
	if (ix < n[j]-1) {
	    k = i+s[j];
	    tb = time[k];
	    xb = xpos[k];
	} else {
	    tb = SF_HUGE;
	    xb = SF_HUGE;
	}
	
	xj = xx+j;
	xj->delta = 1./(d[j]*d[j]);
	xj->direct = j;

	if (ta < tb) {
	    xj->tstencil = xj->tvalue = ta;
	    xj->xstencil = xj->xvalue = xa;
	} else {
	    xj->tstencil = xj->tvalue = tb;
	    xj->xstencil = xj->xvalue = xb;
	}

	/* NOTE: be careful with second-order upwind */
	if (order > 1) {
	    if (ta < tb  && ix-2 >= 0) { 
		k = i-2*s[j];
		if (in[k] != SF_OUT && ta >= time[k]) {
		    xj->delta *= 2.25;
		    xj->tstencil = (4.0*xj->tvalue - time[k])/3.0;
		    xj->xstencil = (4.0*xj->xvalue - xpos[k])/3.0;
		}
	    }
	    if (ta > tb && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (in[k] != SF_OUT && tb >= time[k]) {
		    xj->delta *= 2.25;
		    xj->tstencil = (4.0*xj->tvalue - time[k])/3.0;
		    xj->xstencil = (4.0*xj->xvalue - xpos[k])/3.0;
		}
	    }
	}
    }

    if (xx[0].tvalue <= xx[1].tvalue) {
	if (xx[1].tvalue <= xx[2].tvalue) {
	    vv[0] = xx; vv[1] = xx+1; vv[2] = xx+2;
	} else if (xx[2].tvalue <= xx[0].tvalue) {
	    vv[0] = xx+2; vv[1] = xx; vv[2] = xx+1;
	} else {
	    vv[0] = xx; vv[1] = xx+2; vv[2] = xx+1;
	}
    } else {
	if (xx[0].tvalue <= xx[2].tvalue) {
	    vv[0] = xx+1; vv[1] = xx; vv[2] = xx+2;
	} else if (xx[2].tvalue <= xx[1].tvalue) {
	    vv[0] = xx+2; vv[1] = xx+1; vv[2] = xx;
	} else {
	    vv[0] = xx+1; vv[1] = xx+2; vv[2] = xx;
	}
    }

    if(vv[2]->tvalue < SF_HUGE) {   /* update from all three directions */
	if (updaten(i,3,&tres,&xres,&fres,vv) || 
	    updaten(i,2,&tres,&xres,&fres,vv) || 
	    updaten(i,1,&tres,&xres,&fres,vv)) {
	    *tval = tres;
	    *xval = xres;
	    *fval = fres;
	    return;
	}
    } else if(vv[1]->tvalue < SF_HUGE) { /* update from two directions */
	if (updaten(i,2,&tres,&xres,&fres,vv) || 
	    updaten(i,1,&tres,&xres,&fres,vv)) {
	    *tval = tres;
	    *xval = xres;
	    *fval = fres;
	    return;
	}
    } else if(vv[0]->tvalue < SF_HUGE) { /* update from only one direction */
	if (updaten(i,1,&tres,&xres,&fres,vv)) {
	    *tval = tres;
	    *xval = xres;
	    *fval = fres;
	    return;
	}
    }

    *tval = SF_HUGE;
    *xval = SF_HUGE;
    *fval = -1;
    return;
}

bool updaten(int i, int m, float* tres, float* xres, int* fres, struct Upd *vv[])
/* solve quadratic equation (analytical) */
{
    double a, b, c, discr, ttemp, xtemp;
    int ftemp;
    int j;

    if (m == 1) {
	/* one-sided update is always successful */
	ttemp = vv[0]->tstencil+sqrt((double)v[i]/vv[0]->delta);
	xtemp = vv[0]->xvalue;
	ftemp = vv[0]->direct;
    } else{
	/* solve quadratic equation */
	a = b = c = 0.;

	for (j=0; j<m; j++) {
	    a += vv[j]->delta;
	    b += vv[j]->tstencil*vv[j]->delta;
	    c += vv[j]->tstencil*vv[j]->tstencil*vv[j]->delta;
	}
	b /= a;

	discr=b*b+(((double)v[i])-c)/a;

	if (discr < 0.) return false;
    
	ttemp = b + sqrt(discr);

	if (ttemp <= vv[m-1]->tvalue+1.e-7) return false;

	a = b = 0.;

	for (j=0; j<m; j++) {
	    a += vv[j]->delta*(ttemp-vv[j]->tstencil);
	    b += vv[j]->delta*(ttemp-vv[j]->tstencil)*vv[j]->xstencil;
	}

	xtemp = b/a;

	ftemp = -1;
    }

    *tres = ttemp;
    *xres = xtemp;
    *fres = ftemp;
    return true;
}

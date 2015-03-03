/* Fast marching interface (OMP) */
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
/*^*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <assert.h>
#include "upgradomp.h"
#include "fastmarchomp.h"

struct Upd {
    double stencil, value, delta;
};

static int neighbors_nearsource(float* time, float* xs);
void pqueue_insert(float* v1);
float* pqueue_extract(void);
void pqueue_update(int index);
int neighbours(float* time, int i);
int update(float value, float* time, int i);
float qsolve(float* time, int i);
bool updaten(int i, int m, float* res, struct Upd *vv[]);

static float ***x, ***xn, ***x1;
static float **t;
static int **offsets;
static int **in, *n, s[3], order;
static float *v, *o, *d;

void fastmarch_init(int *n1    /* grid samples [3] */, 
		    float *o1  /* grid origin [3] */,
		    float *d1  /* grid sampling [3] */,
		    int order1 /* accuracy order */)
/*< initialize model dimensions and upwind order >*/
{
    int its, mts;
    int maxband;

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    /* model dimensions */
    n = n1; order = order1; o = o1; d = d1;
    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];

    /* allocate shared memory */
    in = sf_intalloc2(n[0]*n[1]*n[2],mts);

    x  = (float ***) sf_alloc(mts,sizeof (float **));
    xn = (float ***) sf_alloc(mts,sizeof (float **));
    x1 = (float ***) sf_alloc(mts,sizeof (float **));

    offsets = (int **) sf_alloc(mts,sizeof (int *));
    t = (float **) sf_alloc(mts,sizeof (float *));

    maxband = 0;
    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];

    for (its=0; its < mts; its++) {
	x[its] = (float **) sf_alloc ((10*maxband+1),sizeof (float *));
	offsets[its] = (int *) sf_alloc (n[0]*n[1]*n[2],sizeof (int));
    }
}

void fastmarch_set(float *v1  /* slowness squared */)
/*< set velocity model (slowness squared) >*/
{
    v = v1;
}

int fastmarch(float *time   /* time */,
	      float *source /* source */,
	      int *list     /* list */,
	      int *mask     /* recv */,
	      float *data   /* reco */,
	      float *rhs    /* rhs */,
	      upgrad upg    /* stencil */)
/*< run fast marching eikonal solver >*/
{
    int its;
    float xs[3], *p;
    int npoints, i, j, k, count, length;
 
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    /* point t to time for pqueue operations */
    t[its] = time;

    /* source distance from origin */
    xs[0] = source[0]-o[0];
    xs[1] = source[1]-o[1];
    xs[2] = source[2]-o[2];
    
    /* initialize priority queue */
    xn[its] = x[its];
    x1[its] = x[its]+1;
    
    count = 0;
    length = 0;
    
    /* fast marching */
    for (npoints =  neighbors_nearsource(time,xs);
	 npoints > 0;
	 npoints -= neighbours(time,i)) {

	/* smallest value in queue */
	p = pqueue_extract();
	
	if (p == NULL) {
	    sf_warning("%s: shot (%g,%g,%g) heap exausted!",__FILE__,source[0],source[1],source[2]);
	    break;
	}
	
	i = p-time;

	/* update wave front */
	in[its][i] = SF_IN;

	/* update stencil */
	upgrad_set(upg,time,i,in[its],length);
	length++;
	
	/* update rhs */
	k = list[0];
	for (j=0; j < list[1]; j++) {
	    if (i == mask[j]) {
		rhs[k+j] = data[j]-time[i];
		count++;
		break;
	    }
	}
	
	/* break for limited acquisition */
	if (count == list[1]) break;
    }
    
    return(length);
}

void fastmarch_close(void)
/*< free allocated memory >*/
{
    int its, mts;

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    for (its=0; its < mts; its++)
	free(x[its]);
    
    free(x);
    free(xn);
    free(x1);
}

static int neighbors_nearsource(float* time /* time */,
			 float* xs   /* source location [3] */)
/* initialize point source */
{
    int its;
    int ic, i, j, is[3];
    double delta[3], delta2;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    /* initialize everywhere */
    for (i=0; i < n[0]*n[1]*n[2]; i++) {
	in[its][i] = SF_OUT;
	time[i] = SF_HUGE;
	offsets[its][i] = -1;
    }

    /* Find index of the source location and project it to the grid */
    delta2 = 0.;

    for (j=0; j < 3; j++) {
	is[j] = xs[j]/d[j]+0.5;

	if (is[j] < 0) {
	    is[j]=0;
	} else if (is[j] >= n[j]) {
	    is[j]=n[j]-1;
	}

	delta[j] = xs[j]-is[j]*d[j];
	delta2 += delta[j]*delta[j];
    }
    
    /* source index */
    ic = is[0]+n[0]*(is[1]+n[1]*is[2]);
    
    /* initialize source */
    time[ic] = sqrtf(((double)v[ic])*delta2);
/*  TEMP: the next line seems to be redundant; to be checked later...*/
    in[its][ic] = SF_IN;
    pqueue_insert(time+ic);
    
    return (n[0]*n[1]*n[2]-1);
}

void pqueue_insert(float* v1)
/* insert an element (smallest first) */
{
    int its, newOffset, tIndex;
    float **xi, **xq;
    unsigned int q;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif


    xi = ++xn[its];
/*  TEMP: the next line seems to be redundant; to be checked later...*/
    *xi = v1;    
    q = (unsigned int) (xn[its]-x[its]);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x[its] + q;
	if (*v1 > **xq) break;
	*xi = *xq;

	/* now that it moved down, update its offset */
	newOffset = xi-x[its];
	tIndex = (*xi)-t[its];
	offsets[its][tIndex] = newOffset;

	xi = xq;
    }
    *xi = v1; 
    
    /* now that we moved it far enough up, record the offset */
    newOffset = xi-x[its];
    tIndex = v1-t[its];
    offsets[its][tIndex] = newOffset;
}

float* pqueue_extract(void)
/* extract the smallest element */
{
    int its, newOffset, tIndex;
    unsigned int c;
    int nn;
    float *vv, *formerlyLast;
    float **xi, **xc;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    /* check if queue is empty */
    nn = (int) (xn[its]-x[its]);
    if (nn == 0) return NULL;

    vv = *(x1[its]);
    /* label vv to be out by set its offset as -1 */
    tIndex = vv-t[its];
    offsets[its][tIndex] = -1;

    *(xi = x1[its]) = formerlyLast = *(xn[its]--);
    /* note: formerlyLast's offset is inconsistent until the very last step */
    nn--;
    for (c = 2; c <= (unsigned int) nn; c <<= 1) {
	xc = x[its] + c;
	if (c < (unsigned int) nn && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*formerlyLast <= **xc) break;
	*xi = *xc; 

	/* now that it moved up, update its offset */
	newOffset = xi-x[its];
	tIndex = (*xi)-t[its];
	offsets[its][tIndex] = newOffset;

	xi = xc;
    }
    *xi = formerlyLast;

    /* now that we moved it far enough down, record the offset */
    newOffset = xi-x[its];
    tIndex = *xi-t[its];
    offsets[its][tIndex] = newOffset;

    return vv;
}

void pqueue_update(int index)
/* restore the heap */
{
    int its, newOffset, tIndex;
    unsigned int c;
/*
    int n;
*/
    float **xc, **xi;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    c = offsets[its][index];
    xi = x[its]+c;
    /* do not need to visit its children */
/*
    n = (int) (xn[its]-x[its]); c = (unsigned int) (xi-x[its]);
    for (c <<= 1; c <= (unsigned int) n; c <<= 1) {
	xc = x[its] + c;
	if (c < (unsigned int) n && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (**vv <= **xc) break;
	*xi = *xc; xi = xc;
    }
*/
    for (c >>= 1; c > 0; c >>= 1) {
	xc = x[its] + c;
	if (t[its][index] > **xc) break;
	*xi = *xc;
	
	/* now that it moved down, update its offset */
	newOffset = xi-x[its];
	tIndex = (*xi)-t[its];
	offsets[its][tIndex] = newOffset;
	
	xi = xc; 
    }
    *xi = t[its]+index; 

    /* now that we moved it far enough up, record the offset */
    newOffset = xi-x[its];
    tIndex = *xi-t[its];
    offsets[its][tIndex] = newOffset;
}

int neighbours(float* time, int i) 
/* update neighbors of gridpoint i, return number of updated points */
{
    int its;
    int j, k, ix, np;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif    

    np = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];

	/* try both directions */
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[its][k] != SF_IN) np += update(qsolve(time,k),time,k);
	}
	if (ix-1 >= 0  ) {
	    k = i-s[j];
	    if (in[its][k] != SF_IN) np += update(qsolve(time,k),time,k);
	}
    }
    return np;
}

int update(float value, float* time, int i)
/* update gridpoint i with new value and modify wave front */
{
    int its;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    /* only update when smaller than current value */
    if (value < time[i]) {
	time[i] = value;
	if (in[its][i] == SF_OUT) { 
	    in[its][i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
	} else {
	    assert(in[its][i] == SF_FRONT);
	    pqueue_update(i);
	}
    }

    return 0;
}

float qsolve(float* time, int i)
/* find new traveltime at gridpoint i */
{
    int its;
    int j, k, ix;
    float a, b, res;
    struct Upd *vv[3], xx[3], *xj;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif   

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
	xj->delta = 1./(d[j]*d[j]);;

	if (a < b) {
	    xj->stencil = xj->value = a;
	} else {
	    xj->stencil = xj->value = b;
	}

	/* second order local upwind stencil */
	if (order > 1) {
	    if (a < b  && ix-2 >= 0) { 
		k = i-2*s[j];
		if (in[its][k] != SF_OUT && a >= time[k]) {
		    xj->delta *= 2.25;
		    xj->stencil = (4.0*xj->value - time[k])/3.0;
		}
	    }
	    if (a > b && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (in[its][k] != SF_OUT && b >= time[k]) {
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

    *res = t;
    return true;
}

/* 	$Id: fastmarch.c 5686 2010-04-07 16:33:34Z llisiw $	 */

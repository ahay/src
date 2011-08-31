/* Fast marching main interface (OMP) */
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

#include <rsf.h>
/*^*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fastmarchomp.h"

struct Upd {
    double stencil, value, delta;
};

int neighbors_nearsource(float* time, float* xs);
void pqueue_insert(float* v1);
float* pqueue_extract(void);
int neighbours(float* time, int i);
int update(float value, float* time, int i);
float qsolve(float* time, int i);
bool updaten (int i, int m, float* res, struct Upd *vv[]);

static float ***x, ***xn, ***x1;
static int **in, *n, s[3], order;
static float *v, *o, *d;

void fastmarch_init (int *n1    /* grid samples [3] */, 
		     float *o1  /* grid origin [3] */,
		     float *d1  /* grid sampling [3] */,
		     int order1 /* accuracy order */)
/*< Initialize >*/
{
    int its, mts;
    int maxband;

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    n = n1; order = order1; o = o1; d = d1;
    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];

    in = sf_intalloc2(n[0]*n[1]*n[2],mts);
    
    x  = (float ***) sf_alloc (mts,sizeof (float **));
    xn = (float ***) sf_alloc (mts,sizeof (float **));
    x1 = (float ***) sf_alloc (mts,sizeof (float **));

    maxband = 0;
    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];

    for (its=0; its < mts; its++)
	x[its] = (float **) sf_alloc ((10*maxband+1),sizeof (float *));
}

void fastmarch_set (float *v1  /* slowness squared */)
/*< set velocity model >*/
{
    v = v1;
}

void fastmarch (float* time                /* time */,
		float s0,float s1,float s2 /* source */)
/*< Run fast marching eikonal solver >*/
{
    int its;
    float xs[3], *p;
    int npoints, i;
 
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
   
    xs[0] = s0-o[0]; xs[1] = s1-o[1]; xs[2] = s2-o[2];

    xn[its] = x[its];
    x1[its] = x[its]+1;

    for (npoints =  neighbors_nearsource (time,xs);
	 npoints > 0;
	 npoints -= neighbours(time,i)) {

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p - time;

	in[its][i] = SF_IN;
    }
}

void fastmarch_close (void)
/*< Free allocated storage >*/
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

int neighbors_nearsource(float* time /* time */,
			 float* xs   /* source location [3] */)
/*< initialize the source >*/
{
    int its;
    int np, ic, i, j, is, start[3], endx[3], ix, iy, iz;
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
    }

    /* Find index of the source location and project it to the grid */
    for (j=0; j < 3; j++) {
	is = xs[j]/d[j]+0.5;
	start[j] = is-1; 
	endx[j]  = is+1;
    } 
    
    for (j=0; j < 3; j++) {
	if (start[j] < 0) {
	    start[j]=0;
	} else if (start[j] >= n[j]) {
	    start[j]=n[j]-1;
	}
    }
    for (j=0; j < 3; j++) {
	if (endx[j] < 0) {
	    endx[j]=0;
	} else if (endx[j] >= n[j]) {
	    endx[j]=n[j]-1;
	}
    }
    
    ic = (start[0]+endx[0])/2 + 
	n[0]*((start[1]+endx[1])/2 +
	      n[1]*(start[2]+endx[2])/2);
    
    /* loop in a small box around the source */
    np = n[0]*n[1]*n[2];
    for (ix=start[2]; ix <= endx[2]; ix++) {
	for (iy=start[1]; iy <= endx[1]; iy++) {
	    for (iz=start[0]; iz <= endx[0]; iz++) {
		np--;
		i = iz + n[0]*(iy + n[1]*ix);

		delta[0] = xs[0]-iz*d[0];
		delta[1] = xs[1]-iy*d[1];
		delta[2] = xs[2]-ix*d[2];

		delta2 = 0.;
		for (j=0; j < 3; j++) {
		    delta2 += delta[j]*delta[j];
		}

		/* analytical formula (Euclid) */ 
		time[i] = sqrtf(((double)v[ic])*delta2);
		in[its][i] = SF_IN;

		if ((n[0] > 1 && (iz == start[0] || iz == endx[0])) ||
		    (n[1] > 1 && (iy == start[1] || iy == endx[1])) ||
		    (n[2] > 1 && (ix == start[2] || ix == endx[2]))) {
		    pqueue_insert(time+i);
		}
	    }
	}
    }
    
    return np;
}

void pqueue_insert(float* v1)
/*< Insert an element (smallest first) >*/
{
    int its;
    float **xi, **xq;
    unsigned int q;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    xi = ++xn[its];
    *xi = v1;
    q = (unsigned int) (xn[its]-x[its]);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x[its] + q;
	if (*v1 > **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v1; 
}

float* pqueue_extract(void)
/*< Extract the smallest element >*/
{
    int its;
    unsigned int c;
    int nn;
    float *vv, *t;
    float **xi, **xc;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    vv = *(x1[its]);
    *(xi = x1[its]) = t = *(xn[its]--);
    nn = (int) (xn[its]-x[its]);
    if (nn < 0) return NULL;
    for (c = 2; c <= (unsigned int) nn; c <<= 1) {
	xc = x[its] + c;
	if (c < (unsigned int) nn && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*t <= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return vv;
}

int neighbours(float* time, int i) 
/*< Update neighbors of gridpoint i, return number of updated points >*/
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
/* update gridpoint i with new value */
{
    int its;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    if (value < time[i]) {
	time[i] = value;
	if (in[its][i] == SF_OUT) { 
	    in[its][i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
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

    if(vv[2]->value < SF_HUGE) {   /* ALL THREE DIRECTIONS CONTRIBUTE */
	if (updaten(i,3,&res,vv) || 
	    updaten(i,2,&res,vv) || 
	    updaten(i,1,&res,vv)) return res;
    } else if(vv[1]->value < SF_HUGE) { /* TWO DIRECTIONS CONTRIBUTE */
	if (updaten(i,2,&res,vv) || 
	    updaten(i,1,&res,vv)) return res;
    } else if(vv[0]->value < SF_HUGE) { /* ONE DIRECTION CONTRIBUTES */
	if (updaten(i,1,&res,vv)) return res;
    }
	
    return SF_HUGE;
}

bool updaten (int i, int m, float* res, struct Upd *vv[]) 
/* updating */
{
    double a, b, c, discr, t;
    int j;

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

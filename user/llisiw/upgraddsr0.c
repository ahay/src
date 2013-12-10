/* Dot product with upwind gradient */
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

#include "upgraddsr0.h"

#ifndef _upgraddsr0_h

typedef struct Upgrad *upgrad;
/* abstract data type */
/*^*/

#endif

struct Upgrad {
    long *order;
    int **pos;
    unsigned char **update;
    double **ww, **qq;
};

static int ndim, ss[3];
static long nt, ntt;
static const int *nn;
static double dd[3];
static const float *t0, *w0;

static int fermat(const void *a, const void *b)
/* comparison for traveltime sorting from small to large */
{
    float ta, tb;

    ta = t0[*(long *)a];
    tb = t0[*(long *)b];

    if (ta >  tb) return 1;
    if (ta == tb) return 0;
    return -1;
}

void line2cart(int dim       /* number of dimensions */, 
	       const int* nn /* box size [dim] */, 
	       long i        /* line coordinate */, 
	       int* ii       /* cartesian coordinates [dim] */)
/* convert line to Cartesian */
{
    int axis;
 
    for (axis = 0; axis < dim; axis++) {
      ii[axis] = i%((long) nn[axis]);
      i /= (long) nn[axis];
    }
}

long* upgrad_order(upgrad upg)
/*< upwind order >*/
{
    return (upg->order);
}

upgrad upgrad_init(int ndim_in      /* number of dimensions */,
		   const int *nn_in /* [dim] data size */,
		   const float *d   /* [dim] data sampling */)
/*< initialize >*/
{
    upgrad upg;
    int i;
    long it;

    /* for 2D model only */
    if (ndim_in > 3) sf_error("%s: dim=%d > 3",__FILE__,ndim_in);

    ndim = ndim_in;
    nn = nn_in;

    nt = nn[1]*(nn[1]+1)/2;
    for (i=0; i < ndim; i++) {
	ss[i] = nt;
	dd[i] = 1./(d[i]*d[i]);
    }

    upg = (upgrad) sf_alloc(1,sizeof(*upg));

    /* traveltime sort */
    upg->order = (long *) sf_alloc(nt, sizeof(long));
    
    /* source-receiver location */
    upg->pos = sf_intalloc2(2,nt);

    /* upwind indicator */
    upg->update = sf_ucharalloc2(2,nt);

    /* upwind stencil for traveltime */
    upg->ww = (double **) sf_alloc(nt, sizeof(double *));
    for (it=0; it < nt; it++) {
	upg->ww[it] = (double *) sf_alloc(ndim+1, sizeof(double));
    }

    /* upwind stencil for slowness-squared */
    upg->qq = (double **) sf_alloc(nt, sizeof(double *));
    for (it=0; it < nt; it++) {
	upg->qq[it] = (double *) sf_alloc(2, sizeof(double));
    }
    
    return upg;
}

void upgrad_set(upgrad upg      /* upwind stencil */,
		const float *r0 /* reference time */,
		const float *s0 /* reference slowness-squared */,
		const int *f0   /* reference flag */)
/*< supply reference >*/
{
    bool source;
    int i, m, ii[3];
    long it, jt, a, b;
    unsigned char *up;
    double t, t2, dt[3], ws, wr;

    t0 = r0;
    w0 = s0;

    for (it = 0; it < nt; it++) {
	jt = upg->order[it];

	line2cart(ndim,nn,jt,ii);

	up = upg->update[it];
	up[0] = up[1] = 0;
	t = t0[jt];
	upg->ww[it][ndim] = 0.;
	upg->qq[jt][0] = 0.;
	upg->qq[jt][1] = 0.;

	if (t == SF_HUGE) {
	    ntt = it;
	    break;
	}

	/* source-receiver position */
	upg->pos[jt][0] = ii[2]*ss[1]+ii[0];
	upg->pos[jt][1] = ii[1]*ss[1]+ii[0];

	/* slowness-squared */
	ws = w0[upg->pos[jt][0]];
	wr = w0[upg->pos[jt][1]];

	source = true;
	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    a = jt-ss[i];
	    b = jt+ss[i];
	    if ((ii[i] == 0) || 
		(ii[i] != nn[i]-1 && 1==fermat(&a,&b))) {
		up[1] |= m;
		t2 = t0[b];
	    } else {
		t2 = t0[a];
	    }

	    if (t2 < t) {
		up[0] |= m;
		dt[i] = t-t2;
		source = false;
	    } else {
		dt[i] = 0.;
	    }
	}

	/* throw away h < 0 */
	if (upg->pos[jt][0] > upg->pos[jt][1]) source = true;

	if (source) continue;

	/* one-sided */
	if ((f0 == NULL && dt[1] == 0. && dt[2] == 0.) || (f0 != NULL && f0[jt] == 1)) {
	    upg->ww[it][0] = sqrt(dd[0]);
	    upg->ww[it][1] = 0.;
	    upg->ww[it][2] = 0.;
	    
	    upg->ww[it][ndim] = upg->ww[it][0];
	    
	    upg->qq[jt][0] = 0.5/sqrt(ws);
	    upg->qq[jt][1] = 0.5/sqrt(wr);
	    
	    continue;
	}

	if ((f0 == NULL && dt[0] == 0. && dt[2] == 0.) || (f0 != NULL && f0[jt] == 2)) {
	    upg->ww[it][0] = 0.;
	    upg->ww[it][1] = sqrt(dd[1]);
	    upg->ww[it][2] = 0.;
	    
	    upg->ww[it][ndim] = upg->ww[it][1];
	    
	    upg->qq[jt][0] = 0.;
	    upg->qq[jt][1] = 0.5/sqrt(wr);

	    continue;
	}

	if ((f0 == NULL && dt[0] == 0. && dt[1] == 0.) || (f0 != NULL && f0[jt] == 3)) {
	    upg->ww[it][0] = 0.;
	    upg->ww[it][1] = 0.;
	    upg->ww[it][2] = sqrt(dd[2]);
	    
	    upg->ww[it][ndim] = upg->ww[it][2];
	    
	    upg->qq[jt][0] = 0.5/sqrt(ws);
	    upg->qq[jt][1] = 0.;

	    continue;
	}

	/* two-sided */
	if ((f0 == NULL && dt[0] == 0.) || (f0 != NULL && f0[jt] == 4)) {
	    sf_error("flag[%d] = 4",jt);
	    
	    continue;
	}

	if ((f0 == NULL && dt[1] == 0.) || (f0 != NULL && f0[jt] == 5)) {
	    upg->ww[it][0] = sqrt(dd[0]);
	    upg->ww[it][1] = 0.;
	    upg->ww[it][2] = dt[2]*dd[2]/sqrt(ws-dt[2]*dt[2]*dd[2]);
	    
	    upg->ww[it][ndim] = upg->ww[it][0]+upg->ww[it][2];
	    
	    upg->qq[jt][0] = 0.5/sqrt(ws-dt[2]*dt[2]*dd[2]);
	    upg->qq[jt][1] = 0.5/sqrt(wr);

	    continue;
	}

	if ((f0 == NULL && dt[2] == 0.) || (f0 != NULL && f0[jt] == 6)) {
	    upg->ww[it][0] = sqrt(dd[0]);
	    upg->ww[it][1] = dt[1]*dd[1]/sqrt(wr-dt[1]*dt[1]*dd[1]);
	    upg->ww[it][2] = 0.;
	    
	    upg->ww[it][ndim] = upg->ww[it][0]+upg->ww[it][1];
	    
	    upg->qq[jt][0] = 0.5/sqrt(ws);
	    upg->qq[jt][1] = 0.5/sqrt(wr-dt[1]*dt[1]*dd[1]);

	    continue;
	}

	/* three-sided */
	if ((f0 == NULL) || (f0 != NULL && f0[jt] == 7)) {
	    upg->ww[it][0] = sqrt(dd[0]);
	    upg->ww[it][1] = dt[1]*dd[1]/sqrt(wr-dt[1]*dt[1]*dd[1]);
	    upg->ww[it][2] = dt[2]*dd[2]/sqrt(ws-dt[2]*dt[2]*dd[2]);
	    
	    upg->ww[it][ndim] = upg->ww[it][0]+upg->ww[it][1]+upg->ww[it][2];
	    
	    upg->qq[jt][0] = 0.5/sqrt(ws-dt[2]*dt[2]*dd[2]);
	    upg->qq[jt][1] = 0.5/sqrt(wr-dt[1]*dt[1]*dd[1]);
	}
    }

    ntt = nt;
}

void upgrad_close(upgrad upg)
/*< free allocated memory >*/
{
    free(upg->order);
    free(*(upg->pos));
    free(upg->pos);
    free(*(upg->ww));
    free(upg->ww);
    free(*(upg->qq));
    free(upg->qq);
    free(*(upg->update));
    free(upg->update);
    free(upg);
}

void upgrad_solve(upgrad upg,
		  const float *rhs /* right-hand side */, 
		  float *x         /* solution */,
		  const float *x0  /* initial solution */)
/*< inv(alpha) >*/
{
    int i, m, k;
    long it, jt, j;
    unsigned char *up;
    double num, den;
   
    for (it = 0; it < ntt; it++) {
	jt = upg->order[it];

	num = rhs[jt];
	up = upg->update[it];
	den = upg->ww[it][ndim];

	if (den == 0.) { /* at the source, use boundary conditions */
	    x[jt] = (NULL != x0)? x0[jt]: 0.;
	    continue;
	}

	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    if (up[0] & m) {
		j = (up[1] & m)? jt+ss[i]:jt-ss[i];	
		num += upg->ww[it][i]*x[j];
	    }
	}

	x[jt] = num/den;
    }

    /* copy */
    for (k=1; k < nn[2]; k++) {
	for (i=0; i < k; i++) {
	    for (m=0; m < nn[0]; m++) {
		x[(long) k*ss[2]+i*ss[1]+m] = x[(long) i*ss[2]+k*ss[1]+m];
	    }
	}
    }
}

void upgrad_inverse(upgrad upg,
		    float *rhs      /* right-hand side */,
		    const float *x  /* solution */,
		    const float *x0 /* initial solution */)
/*< adjoint inv(alpha) >*/
{
    int i, m, ii[3];
    long it, jt, j;
    unsigned char *up;
    double den, w;

    for (it = 0; it < nt; it++) {
	rhs[it] = 0.;
    }
   
    for (it = ntt-1; it >= 0; it--) {
	jt = upg->order[it];

	line2cart(ndim,nn,jt,ii);

	/* paste */
	rhs[jt] += x[jt]+x[(long) ii[1]*ss[2]+ii[2]*ss[1]+ii[0]];

	up = upg->update[it];
	den = upg->ww[it][ndim];

	if (den == 0.) { /* at the source, use boundary conditions */
	    rhs[jt] = (NULL != x0)? x0[jt]: 0.;
	} else {
	    rhs[jt] = rhs[jt]/den;
	}
	
	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    if (up[0] & m) {
		j = (up[1] & m)? jt+ss[i]:jt-ss[i];
		w = upg->ww[it][i]*rhs[jt];
		rhs[j] += w;
	    }
	}
    }
}

void upgrad_collect(upgrad upg,
		    const float *rhs /* right-hand side */, 
		    float *x         /* solution */)
/*< beta >*/
{
    long it;

    for (it = 0; it < ntt; it++) {
	x[it] = upg->qq[it][0]*rhs[upg->pos[it][0]]
	    +upg->qq[it][1]*rhs[upg->pos[it][1]];
    }
}

void upgrad_spread(upgrad upg,
		   float *rhs     /* right-hand side */,
		   const float *x /* solution */)
/*< adjoint beta >*/
{
    long it;

    for (it = 0; it < ntt; it++) {
	rhs[upg->pos[it][0]] += upg->qq[it][0]*x[it];
	rhs[upg->pos[it][1]] += upg->qq[it][1]*x[it];
    }
}

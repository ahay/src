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

#include "upgraddsr.h"

#ifndef _upgraddsr_h

typedef struct Upgrad *upgrad;
/* abstract data type */
/*^*/

#endif

struct Upgrad {
    int *order, **pos;
    unsigned char **update;
    double **ww, **qq;
    bool *source;
};

static int ndim, nt, ss[3];
static const int *nn;
static double dd[3];
static const float *t0, *w0;

static int fermat(const void *a, const void *b)
/* comparison for traveltime sorting from small to large */
{
    float ta, tb;

    ta = t0[*(int *)a];
    tb = t0[*(int *)b];

    if (ta >  tb) return 1;
    if (ta == tb) return 0;
    return -1;
}

upgrad upgrad_init(int mdim        /* number of dimensions */,
		   const int *mm   /* [dim] data size */,
		   const float *d  /* [dim] data sampling */)
/*< initialize >*/
{
    upgrad upg;
    int i;

    if (mdim > 3) sf_error("%s: dim=%d > 3",__FILE__,mdim);

    ndim = mdim;
    nn = mm;

    nt = 1;
    for (i=0; i < ndim; i++) {
	ss[i] = nt;
	nt *= nn[i];
	dd[i] = 1.0/(d[i]*d[i]);
    }

    upg = (upgrad) sf_alloc(1,sizeof(*upg));

    upg->update = sf_ucharalloc2(2,nt);

    upg->ww = (double **) sf_alloc(nt, sizeof(double *));
    for (i=0; i < nt; i++) {
	upg->ww[i] = (double *) sf_alloc(ndim+1, sizeof(double));
    }

    upg->qq = (double **) sf_alloc(nt, sizeof(double *));
    for (i=0; i < nt; i++) {
	upg->qq[i] = (double *) sf_alloc(2, sizeof(double));
    }

    upg->order = sf_intalloc(nt);
    upg->source = sf_boolalloc(nt);
    upg->pos = sf_intalloc2(2,nt);

    return upg;
}

void upgrad_set(upgrad upg      /* upwind stencil */,
		const float *r0 /* reference time */,
		const float *s0 /* reference slowness-squared */,
		const int *f0   /* reference flag */)
/*< supply reference >*/
{
    int i, m, it, jt, ii[3], a, b;
    unsigned char *up;
    double t, t2, dt[3], ws, wr;

    t0 = r0;
    w0 = s0;

    /* sort from small to large traveltime */
    for (it = 0; it < nt; it++) {
	upg->order[it] = it;
    }
    qsort(upg->order, nt, sizeof(int), fermat);
     
    for (it = 0; it < nt; it++) {
	jt = upg->order[it];

	sf_line2cart(ndim,nn,jt,ii);

	up = upg->update[it];
	up[0] = up[1] = 0;
	t = t0[jt];
	upg->ww[it][ndim] = 0.;
	upg->source[it] = true;

	/* position */
	upg->pos[jt][0] = ii[2]*ss[1]+ii[0];
	upg->pos[jt][1] = ii[1]*ss[1]+ii[0];

	/* slowness-squared */
	ws = w0[upg->pos[jt][0]];
	wr = w0[upg->pos[jt][1]];

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
		upg->source[it] = false;
		dt[i] = t-t2;
	    } else {
		dt[i] = 0.;
	    }
	}

	/* throw away h <= 0  */
	if (upg->pos[jt][0] > upg->pos[jt][1]) upg->source[it] = true;

	if (upg->source[it]) continue;

	/* one-sided */
	if ((f0 == NULL && dt[1] == 0. && dt[2] == 0.) || (f0 != NULL && f0[jt] == 1)) {
	    upg->ww[it][0] = 2.*dt[0]*dd[0]*(dt[0]*dt[0]*dd[0]-ws-wr);
	    upg->ww[it][1] = 0.;
	    upg->ww[it][2] = 0.;
	    
	    upg->ww[it][ndim] = upg->ww[it][0];
	    
	    upg->qq[jt][0] = dt[0]*dt[0]*dd[0]-ws+wr;
	    upg->qq[jt][1] = dt[0]*dt[0]*dd[0]+ws-wr;
	    
	    continue;
	}

	if ((f0 == NULL && dt[0] == 0. && dt[2] == 0.) || (f0 != NULL && f0[jt] == 2)) {
	    upg->ww[it][0] = 0.;
	    upg->ww[it][1] = 2.*dt[1]*dd[1]*(dt[1]*dt[1]*dd[1]+ws-wr);
	    upg->ww[it][2] = 0.;
	    
	    upg->ww[it][ndim] = upg->ww[it][1];
	    
	    upg->qq[jt][0] = -dt[1]*dt[1]*dd[1]-ws+wr;
	    upg->qq[jt][1] = dt[1]*dt[1]*dd[1]+ws-wr;

	    continue;
	}

	if ((f0 == NULL && dt[0] == 0. && dt[1] == 0.) || (f0 != NULL && f0[jt] == 3)) {
	    upg->ww[it][0] = 0.;
	    upg->ww[it][1] = 0.;
	    upg->ww[it][2] = 2.*dt[2]*dd[2]*(dt[2]*dt[2]*dd[2]-ws+wr);
	    
	    upg->ww[it][ndim] = upg->ww[it][2];
	    
	    upg->qq[jt][0] = dt[2]*dt[2]*dd[2]-ws+wr;
	    upg->qq[jt][1] = -dt[2]*dt[2]*dd[2]+ws-wr;

	    continue;
	}

	/* two-sided */
	if ((f0 == NULL && dt[0] == 0.) || (f0 != NULL && f0[jt] == 4)) {
	    upg->ww[it][0] = 0.;
	    upg->ww[it][1] = sqrt(dd[1]);
	    upg->ww[it][2] = sqrt(dd[2]);
	    
	    upg->ww[it][ndim] = upg->ww[it][1]+upg->ww[it][2];
	    
	    upg->qq[jt][0] = 1./(2.*sqrt(ws));
	    upg->qq[jt][1] = 1./(2.*sqrt(wr));
	    
	    continue;
	}

	if ((f0 == NULL && dt[1] == 0.) || (f0 != NULL && f0[jt] == 5)) {
	    upg->ww[it][0] = 2.*dt[0]*dd[0]*(dt[0]*dt[0]*dd[0]+dt[2]*dt[2]*dd[2]-ws-wr);
	    upg->ww[it][1] = 0.;
	    upg->ww[it][2] = 2.*dt[2]*dd[2]*(dt[2]*dt[2]*dd[2]+dt[0]*dt[0]*dd[0]-ws+wr);
	    
	    upg->ww[it][ndim] = upg->ww[it][0]+upg->ww[it][2];
	    
	    upg->qq[jt][0] = dt[2]*dt[2]*dd[2]+dt[0]*dt[0]*dd[0]-ws+wr;
	    upg->qq[jt][1] = -dt[2]*dt[2]*dd[2]+dt[0]*dt[0]*dd[0]+ws-wr;

	    continue;
	}

	if ((f0 == NULL && dt[2] == 0.) || (f0 != NULL && f0[jt] == 6)) {
	    upg->ww[it][0] = 2.*dt[0]*dd[0]*(dt[0]*dt[0]*dd[0]+dt[1]*dt[1]*dd[1]-ws-wr);
	    upg->ww[it][1] = 2.*dt[1]*dd[1]*(dt[1]*dt[1]*dd[1]+dt[0]*dt[0]*dd[0]+ws-wr);
	    upg->ww[it][2] = 0.;
	    
	    upg->ww[it][ndim] = upg->ww[it][0]+upg->ww[it][1];
	    
	    upg->qq[jt][0] = -dt[1]*dt[1]*dd[1]+dt[0]*dt[0]*dd[0]-ws+wr;
	    upg->qq[jt][1] = dt[1]*dt[1]*dd[1]+dt[0]*dt[0]*dd[0]+ws-wr;

	    continue;
	}

	/* three-sided */
	if ((f0 == NULL) || (f0 != NULL && f0[jt] == 7)) {
	    upg->ww[it][0] = 2.*dt[0]*dd[0]*(dt[0]*dt[0]*dd[0]+dt[2]*dt[2]*dd[2]+dt[1]*dt[1]*dd[1]-ws-wr);
	    upg->ww[it][1] = 2.*dt[1]*dd[1]*(dt[1]*dt[1]*dd[1]-dt[2]*dt[2]*dd[2]+dt[0]*dt[0]*dd[0]+ws-wr);
	    upg->ww[it][2] = 2.*dt[2]*dd[2]*(dt[2]*dt[2]*dd[2]-dt[1]*dt[1]*dd[1]+dt[0]*dt[0]*dd[0]-ws+wr);
	    
	    upg->ww[it][ndim] = upg->ww[it][0]+upg->ww[it][1]+upg->ww[it][2];
	    
	    upg->qq[jt][0] = dt[2]*dt[2]*dd[2]-dt[1]*dt[1]*dd[1]+dt[0]*dt[0]*dd[0]-ws+wr;
	    upg->qq[jt][1] = dt[1]*dt[1]*dd[1]-dt[2]*dt[2]*dd[2]+dt[0]*dt[0]*dd[0]+ws-wr;
	}
    }
}

void upgrad_close(upgrad upg)
/*< free allocated storage >*/
{
    free(*(upg->ww));
    free(upg->ww);
    free(*(upg->qq));
    free(upg->qq);
    free(*(upg->update));
    free(upg->update);
    free(upg->order);
    free(upg->source);
    free(*(upg->pos));
    free(upg->pos);
    free(upg);
}

void upgrad_solve(upgrad upg,
		  const float *rhs /* right-hand side */, 
		  float *x         /* solution */,
		  const float *x0  /* initial solution */)
/*< inv(alpha) >*/
{
    int it, jt, i, m, j;
    unsigned char *up;
    double num, den;
   
    for (it = 0; it < nt; it++) {
	jt = upg->order[it];

	num = rhs[jt];
	up = upg->update[it];
	den = upg->ww[it][ndim];

	if (upg->source[it]) { /* at the source, use boundary conditions */
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
}

void upgrad_inverse(upgrad upg,
		    float *rhs      /* right-hand side */,
		    const float *x  /* solution */,
		    const float *x0 /* initial solution */)
/*< adjoint inv(alpha) >*/
{
    int it, jt, i, m, j;
    unsigned char *up;
    double den, w;

    for (it = 0; it < nt; it++) {
	rhs[it] = 0.;
    }
   
    for (it = nt-1; it >= 0; it--) {
	jt = upg->order[it];

	rhs[jt] += x[jt];

	up = upg->update[it];
	den = upg->ww[it][ndim];

	if (upg->source[it]) { /* at the source, use boundary conditions */
	    rhs[jt] = (NULL != x0)? x0[jt]: 0.;
	} else {
	    rhs[jt] = rhs[jt]/den;
	}
	/*
	if (rhs[jt] >= 100.) sf_error("rhs[%d] = %g, it = %d",jt,rhs[jt],it);
	*/
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
    int it;

    for (it = 0; it < nt; it++) {
	x[it] = upg->qq[it][0]*rhs[upg->pos[it][0]]+upg->qq[it][1]*rhs[upg->pos[it][1]];
    }
}

void upgrad_spread(upgrad upg,
		   float *rhs     /* right-hand side */,
		   const float *x /* solution */)
/*< adjoint beta >*/
{
    int it;

    for (it = 0; it < nt; it++) {
	/*
	if (upg->pos[it][0] == 6315)
	    sf_warning("rhs = %g, qq[%d][0] = %g, x[%d] = %g",rhs[upg->pos[it][0]],it,upg->qq[it][0],it,x[it]);
	if (upg->pos[it][1] == 6315)
	    sf_warning("rhs = %g, qq[%d][0] = %g, x[%d] = %g",rhs[upg->pos[it][1]],it,upg->qq[it][1],it,x[it]);
	*/

	rhs[upg->pos[it][0]] += upg->qq[it][0]*x[it];
	rhs[upg->pos[it][1]] += upg->qq[it][1]*x[it];
    }
}

void upgrad_copy(float *time)
/*< gama >*/
{
    int i, j, k;

    for (k=1; k < nn[1]; k++) {
	for (j=0; j < k; j++) {
	    for (i=0; i < nn[0]; i++) {
		time[k*ss[2]+j*ss[1]+i] = time[j*ss[2]+k*ss[1]+i];
	    }
	}
    }
}

void upgrad_paste(float *time)
/*< adjoint gama >*/
{
    int i, j, k;

    for (j=1; j < nn[1]; j++) {
	for (k=0; k < j; k++) {
	    for (i=0; i < nn[0]; i++) {
		time[k*ss[2]+j*ss[1]+i] += time[j*ss[2]+k*ss[1]+i];
	    }
	}
    }
}

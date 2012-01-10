/* Dot product with upwind gradient */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#ifndef _upgrad_h

typedef struct Upgrad *upgrad;
/* abstract data type */
/*^*/

#endif

struct Upgrad {
    int *order;
    unsigned char **update;
    float **ww;
};

static int ndim, nt, ss[3];
static const int *nn;
static float dd[3];

void upgrad_init(int mdim        /* number of dimensions */,
		 const int *mm   /* [dim] data size */,
		 const float *d  /* [dim] data sampling */)
/*< initialize >*/
{
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
}

upgrad upgrad_alloc()
/*< allocate memory for stencil >*/
{
    upgrad upg;

    upg = (upgrad) sf_alloc(1,sizeof(*upg));

    upg->update = sf_ucharalloc2(2,nt);
    upg->ww = sf_floatalloc2(ndim+1,nt);
    upg->order = sf_intalloc(nt);

    return upg;
}

void upgrad_set(upgrad upg /* stencil */,
		float *r0  /* reference */,
		int index  /* index */,
		int *flag  /* flag */,
		int length /* length */)
/*< supply reference >*/
{
    int i, m, ii[3], a, b, c;
    unsigned char *up;
    float t, t2;

    upg->order[length] = index;

    sf_line2cart(ndim,nn,index,ii);
    up = upg->update[length];
    up[0] = up[1] = 0;
    t = r0[index];
    upg->ww[length][ndim] = 0.;

    if (length == 0) return;

    for (i=0, m=1; i < ndim; i++, m <<= 1) {
	a = index-ss[i];
	b = index+ss[i];
	if ((ii[i] == 0) || 
	    (ii[i] != nn[i]-1 && r0[a] > r0[b])) {
	    up[1] |= m;
	    t2 = r0[b];
	    c = b;
	} else {
	    t2 = r0[a];
	    c = a;
	}
	
	/* accept only from upwind direction */
	if ((t2 < t) && (flag[c] == SF_IN)) {
	    up[0] |= m;
	    upg->ww[length][i] = (t-t2)*dd[i];
	    upg->ww[length][ndim] += upg->ww[length][i];
	}	    
    }
}

void upgrad_close(upgrad upg)
/*< free allocated storage >*/
{
    free(*(upg->ww));
    free(upg->ww);
    free(*(upg->update));
    free(upg->update);
    free(upg->order);
    free(upg);
}

void upgrad_solve(int length       /* length */,
		  upgrad upg       /* stencil */,
		  const float *rhs /* right-hand side */, 
		  float *x         /* solution */,
		  const float *x0  /* initial solution */)
/*< inverse operator >*/
{
    int it, jt, i, m, j;
    unsigned char *up;
    float num, den;
   
    for (it = 0; it < length; it++) {
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
}

void upgrad_inverse(int length       /* length */,
		    upgrad upg       /* stencil */,
		    float *rhs       /* right-hand side */,
		    const float *x   /* solution */,
		    const float *x0  /* initial solution */)
/*< adjoint of inverse operator >*/
{
    int it, jt, i, m, j;
    unsigned char *up;
    float den, w;

    for (it = 0; it < nt; it++) {
	rhs[it] = 0.;
    }
   
    for (it = length-1; it >= 0; it--) {
	jt = upg->order[it];

	rhs[jt] += x[jt];

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

void upgrad_forw(int length     /* length */,
		 upgrad upg     /* stencil */,
		 const float *x /* solution */,
		 float *rhs     /* right-hand side */)
/*< forward operator >*/
{
    int it, jt, i, m, j;
    unsigned char *up;
    float num, x2;
   
    for (it = 0; it < length; it++) {
	jt = upg->order[it];

	x2 = x[jt];
	up = upg->update[it];
	num = 0.;

	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    if (up[0] & m) {
		j = (up[1] & m)? jt+ss[i]:jt-ss[i];		
		num += upg->ww[it][i]*(x2-x[j]);		
	    }
	}
	
	rhs[jt] = num;
    }
}

void upgrad_adj(int length       /* length */,
		upgrad upg       /* stencil */,
		float *x         /* solution */,
		const float *rhs /* right-hand side */)
/*< adjoint operator >*/
{
    int it, jt, i, m, j;
    unsigned char *up;
    float w;

    for (it = 0; it < nt; it++) {
	x[it] = 0.;
    }
   
    for (it = length-1; it >= 0; it--) {
	jt = upg->order[it];
	up = upg->update[it];

	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    if (up[0] & m) {
		j = (up[1] & m)? jt+ss[i]:jt-ss[i];
		w = upg->ww[it][i]*rhs[jt];

		x[jt] += w;
		x[j]  -= w;
	    }
	}
    }
}

void upgrad_ray(int length       /* length */,
		upgrad upg       /* stencil */,
		float *ray       /* ray */)
/*< extract ray density >*/
{
    int it, jt, i, m, j;
    unsigned char *up;

    for (it = length-1; it >= 0; it--) {
	jt = upg->order[it];

	if (ray[jt] != 0.) {
	    up = upg->update[it];
	    for (i=0, m=1; i < ndim; i++, m <<= 1) {
		if (up[0] & m) {
		    j = (up[1] & m)? jt+ss[i]:jt-ss[i];
		    ray[j] = 1.;
		}
	    }
	}
    }
}

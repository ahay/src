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
static const float *t0;

static int fermat2(const void *a, const void *b)
/* comparison for traveltime sorting from small to large */
{
    float ta, tb;

    ta = t0[*(int *)a];
    tb = t0[*(int *)b];

    if (ta >  tb) return 1;
    if (ta == tb) return 0;
    return -1;
}

upgrad upgrad_init2(int mdim        /* number of dimensions */,
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

    upg->update = sf_ucharalloc2(4,nt);
    upg->ww = sf_floatalloc2(ndim+1,nt);
    upg->order = sf_intalloc(nt);

    return upg;
}

void upgrad_set2(upgrad upg, const float *r0 /* reference */)
/*< supply reference >*/
{
    int i, m, it, jt, ii[3], a, b, c;
    unsigned char *up;
    float t, t2, t3;

    t0 = r0;

    /* sort from small to large traveltime */
    for (it = 0; it < nt; it++) {
	upg->order[it] = it;
    }
    qsort(upg->order, nt, sizeof(int), fermat2);
     
    for (it = 0; it < nt; it++) {
	jt = upg->order[it];

	sf_line2cart(ndim,nn,jt,ii);
	up = upg->update[it];
	up[0] = up[1] = 0;
	up[2] = up[3] = 0;
	t = t0[jt];
	upg->ww[it][ndim] = 0.;
	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    a = jt-ss[i];
	    b = jt+ss[i];
	    if ((ii[i] == 0) || 
		(ii[i] != nn[i]-1 && 1==fermat2(&a,&b))) {
		up[1] |= m;
		t2 = t0[b];
		
		c = jt+2*ss[i];
		if (ii[i] < nn[i]-2 && 1==fermat2(&b,&c)) {
		    up[2] |= m;
		    up[3] |= m;
		    t3 = t0[c];
		}
	    } else {
		t2 = t0[a];
		
		c = jt-2*ss[i];
		if (ii[i] > 1 && 1==fermat2(&a,&c)) {
		    up[2] |= m;
		    t3 = t0[c];
		}
	    }

	    if (t2 < t) {
		up[0] |= m;
		
		if (up[2] & m) {
		    upg->ww[it][i] = (3.*t-4.*t2+t3)/4.*dd[i];
		    upg->ww[it][ndim] += 3.*upg->ww[it][i];
		} else {
		    upg->ww[it][i] = (t-t2)*dd[i];
		    upg->ww[it][ndim] += upg->ww[it][i];
		}
	    }	    
	}
    }
}

void upgrad_close2(upgrad upg)
/*< free allocated storage >*/
{
    free(*(upg->ww));
    free(upg->ww);
    free(*(upg->update));
    free(upg->update);
    free(upg->order);
    free(upg);
}

void upgrad_solve2(upgrad upg,
		   const float *rhs /* right-hand side */, 
		   float *x         /* solution */,
		   const float *x0  /* initial solution */)
/*< inverse operator >*/
{
    int it, jt, i, m, j, k;
    unsigned char *up;
    float num, den;
   
    for (it = 0; it < nt; it++) {
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

		if (up[2] & m) {
		    k = (up[3] & m)? jt+2*ss[i]:jt-2*ss[i];
		    num += 4.*upg->ww[it][i]*x[j];
		    num -= upg->ww[it][i]*x[k];
		} else {
		    num += upg->ww[it][i]*x[j];
		}
	    }
	}
	
	x[jt] = num/den;
    }
}

void upgrad_inverse2(upgrad upg,
		     float *rhs       /* right-hand side */,
		     const float *x   /* solution */,
		     const float *x0  /* initial solution */)
/*< adjoint of inverse operator >*/
{
    int it, jt, i, m, j, k;
    unsigned char *up;
    float den, w;

    for (it = 0; it < nt; it++) {
	rhs[it] = 0.;
    }
   
    for (it = nt-1; it >= 0; it--) {
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

		if (up[2] & m) {
		    k = (up[3] & m)? jt+2*ss[i]:jt-2*ss[i];
		    w = upg->ww[it][i]*rhs[jt];

		    rhs[j] += 4.*w;
		    rhs[k] -= w;
		} else {
		    w = upg->ww[it][i]*rhs[jt];

		    rhs[j] += w;
		}
	    }
	}
    }
}

void upgrad_forw2(upgrad upg,
		  const float *x /* solution */,
		  float *rhs     /* right-hand side */)
/*< forward operator >*/
{
    int it, jt, i, m, j, k;
    unsigned char *up;
    float num, x2;
   
    for (it = 0; it < nt; it++) {
	jt = upg->order[it];

	x2 = x[jt];
	up = upg->update[it];
	num = 0.;

	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    if (up[0] & m) {
		j = (up[1] & m)? jt+ss[i]:jt-ss[i];

		if (up[2] & m) {
		    k = (up[3] & m)? jt+2*ss[i]:jt-2*ss[i];
		    num += upg->ww[it][i]*(3.*x2-4.*x[j]+x[k]);
		} else {
		    num += upg->ww[it][i]*(x2-x[j]);		
		}
	    }
	}
	
	rhs[jt] = num;
    }
}

void upgrad_adj2(upgrad upg,
		 float *x         /* solution */,
		 const float *rhs /* right-hand side */)
/*< adjoint operator >*/
{
    int it, jt, i, m, j, k;
    unsigned char *up;
    float w;

    for (it = 0; it < nt; it++) {
	x[it] = 0.;
    }
   
    for (it = nt-1; it >= 0; it--) {
	jt = upg->order[it];
	up = upg->update[it];

	for (i=0, m=1; i < ndim; i++, m <<= 1) {
	    if (up[0] & m) {
		j = (up[1] & m)? jt+ss[i]:jt-ss[i];

		if (up[2] & m) {
		    k = (up[3] & m)? jt+2*ss[i]:jt-2*ss[i];
		    w = upg->ww[it][i]*rhs[jt];
		
		    x[jt] += 3.*w;
		    x[j]  -= 4.*w;;
		    x[k]  += w;
		} else {
		    w = upg->ww[it][i]*rhs[jt];

		    x[jt] += w;
		    x[j]  -= w;
		}
	    }
	}
    }
}

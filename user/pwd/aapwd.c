/* Amplitude-adjusted PWD */
/*
  Copyright (C) 2021 University of Texas at Austin
  
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
#include "allp3.h"

static int n, *m;
static float *s, *p, *a, *x, *tmp;
static allpass ap;

void aapwd_init(int n1,
		int n2     /* data size */, 
		int nw     /* filter order */,
		bool drift /* if shift filter */,
		float *s1  /* [n] input */,
		int   *m1  /* [n] amplitude mask */,
		float *p1  /* [n] dip */,
		float *a1  /* [n] amplitude */,
    		float *x1  /* [n] intermediate */)
/*< initialize >*/
{
    n = n1*n2;
    s = s1;
    m = m1;
    p = p1;
    a = a1;
    x = x1;

    tmp = sf_floatalloc(n);

    ap = allpass_init (nw,1,n1,n2,1,drift,p);
}

void aapwd_close(void)
/*< free allocated storage >*/
{
    free(tmp);
    free(ap);
}

void aapwd_apply(const float *x2 /* target */,
		 float *y        /* RHS */)
/*< apply the matrix operator >*/
{
    int i;

    allpass1(false, false, ap, x, tmp);
    for (i=0; i < n; i++) {
	/* t - P x */
	y[i] = x2[i] - tmp[i];
    }
    for (i=0; i < n; i++) {
	/* x - A s */
	y[n+i] = x[i] - a[i]*s[i];
    }
}

void aapwd(float *y)
/*< apply the chain >*/
{
    int i;

    for (i=0; i < n; i++) {
	tmp[i] = a[i]*s[i];
    }
    allpass1(false, false, ap, tmp, y);
}

void aapwd_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< linear operator >*/
{
    int i;
    float *pi, *ai, *xi;

    if (nx != 3*n || ny != 2*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    pi = xx;
    ai = xx+n;
    xi = xx+2*n;

    if (adj) {
	/* delta x */ 
	allpass1t(false, false, ap, tmp, yy);
	for (i=0; i < n; i++) {
	    xi[i] -= tmp[i];
	}
	for (i=0; i < n; i++) {
	    xi[i] += yy[i+n];
	}
	/* delta A */
	allpass1(false, true, ap, x, tmp);
	for (i=0; i < n; i++) {
	    /* - P' dp x */
	    pi[i] -= tmp[i]*yy[i];
	    /* - da s */
	    if (m[i]) ai[i] -= s[i]*yy[i+n];
	}
    } else {
	/* delta x */
	allpass1(false, false, ap, xi, tmp);	    
	for (i=0; i < n; i++) {
	    /* - P dx */
	    yy[i] -= tmp[i];
	} 
	for (i=0; i < n; i++) {
	    /* dx */
	    yy[i+n] += xi[i];
	}
	/* delta A */
	allpass1(false, true, ap, x, tmp);
	for (i=0; i < n; i++) {
	    /* - P' dp x */
	    yy[i] -= tmp[i]*pi[i];
	    /* - da s */
	    if (m[i]) yy[i+n] -= s[i]*ai[i];
	}
    }
} 

void aapwdx_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< linear operator for x only >*/
{
    int i;

    if (nx != n || ny != 2*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
	allpass1t(false, false, ap, tmp, yy);
	for (i=0; i < n; i++) {
	    xx[i] -= tmp[i];
	}
	for (i=0; i < n; i++) {
	    xx[i] += yy[i+n];
	} 
    } else {
	allpass1(false, false, ap, xx, tmp);
	for (i=0; i < n; i++) {
	    /* - P dx */
	    yy[i] -= tmp[i];
	}
	for (i=0; i < n; i++) {
	    /* dx */
	    yy[i+n] += xx[i];
	} 
    }
} 

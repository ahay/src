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

static int n;
static float *s, *p, *a, *x, *tmp;
static allpass ap;

void aapwd_init(int n1,
		int n2     /* data size */, 
		int nw     /* filter order */,
		bool drift /* if shift filter */,
		float *s1  /* [n] input */,
		float *p1  /* [n] dip */,
		float *a1  /* [n] amplitude */,
    		float *x1  /* [n] intermediate */)
/*< initialize >*/
{
    n = n1*n2;
    s = s1;
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

void aapwd_apply(const float *x2, float *y)
/*< apply the matrix operator >*/
{
    int i;

    allpass1(false, false, ap, x, y);
    for (i=0; i < n; i++) {
	y[i] = x2[i] - y[i];
    }
    for (i=0; i < n; i++) {
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

/*
void pwdchain_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
[*< linear operator >*]
{
    int ic, i, j, k;

    if (nx != (2*nc-1)*n || ny != nc*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	if (nc > 1) {
	    allpass1t(false, false, ap[0], tmp, y);
	    for (i=0; i < n; i++) {
		x[i+n*nc] -= tmp[i];
	    }
	    for (ic=1; ic < nc-1; ic++) {
		allpass1t(false, false, ap[ic], tmp, y);
		for (i=0; i < n; i++) {		    
		    j = ic*n+i;
		    k = j+n*nc;
		    x[k-n] += y[j];
		    x[k] -= tmp[i];
		}
	    } 
	    for (i=0; i < n; i++) {
		j = (nc-1)*n+i;
		k = j+n*nc;
		x[k-n] += y[j];
	    } 
	} 
	[* delta A *] 
	for (ic=0; ic < nc-1; ic++) {
	    allpass1(false, true, ap[ic], xn+ic*n, tmp);
	    for (i=0; i < n; i++) {
		j = ic*n+i;
		x[j] += y[j]*tmp[i];
	    }
	} 
	allpass1(false, true, ap[nc-1], s, tmp);
	for (i=0; i < n; i++) {
	    j = (nc-1)*n+i;
	    x[j] += y[j]*tmp[i];
	} 
    } else {
	if (nc > 1) {
	    allpass1(false, false, ap[0], x+nc*n, tmp);
	    for (i=0; i < n; i++) {
		y[i] -= tmp[i];
	    }
	    for (ic=1; ic < nc-1; ic++) {
		allpass1(false, false, ap[ic], x+(nc+ic)*n, tmp);
		for (i=0; i < n; i++) {
		    j = ic*n+i;
		    k = j+n*nc;
		    y[j]  += x[k-n] - tmp[i];
		}
	    } 
	    for (i=0; i < n; i++) {
		j = (nc-1)*n+i;
		k = j+n*nc;
		y[j] += x[k-n];
	    } 
	} 
	[* delta A *]
	allpass1(false, true, ap, x, tmp);	    
	for (i=0; i < n; i++) {
		y[j] += x[j]*tmp[i];
	    }
	} 
	allpass1(false, true, ap[nc-1], s, tmp);
	for (i=0; i < n; i++) {
	    j = (nc-1)*n+i;
	    y[j] += x[j]*tmp[i];
	} 
    }
} */

void aapwdx_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator for x only >*/
{
    int i;

    if (nx != n || ny != 2*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	allpass1t(false, false, ap, tmp, y);
	for (i=0; i < n; i++) {
	    x[i] -= tmp[i];
	}
	for (i=0; i < n; i++) {
	    x[i] += y[i+n];
	} 
    } else {
	allpass1(false, false, ap, x, tmp);
	for (i=0; i < n; i++) {
	    y[i] -= tmp[i];
	}
	for (i=0; i < n; i++) {
	    y[i+n] += x[i];
	} 
    }
} 

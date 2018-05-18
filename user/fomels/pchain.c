/* Chain of 1-D PEFs */
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

static int n, nc;
static sf_complex *x0, *xn, *an;

void pchain_init(int n1  /* data size */, 
		 int nc1 /* number of chain elements */,
		 sf_complex *x1  /* [n] input */,
		 sf_complex *an1 /* [n*nc] weights */,
		 sf_complex *xn1 /* [n*(nc-1)] layers */)
/*< initialize >*/
{
    n = n1;
    nc = nc1;
    x0 = x1;
    xn = xn1;
    an = an1;
}

void pchain_apply(const sf_complex *x2, sf_complex *y)
/*< apply the matrix operator >*/
{
    int ic, i, j;

    if (nc > 1) {
	y[0] = x2[0];
	for (i=1; i < n; i++) {
	    y[i] = x2[i] - xn[i]+an[i]*xn[i-1];
	}
	for (ic=1; ic < nc-1; ic++) {
	    y[ic*n] = xn[(ic-1)*n];
	    for (i=1; i < n; i++) {
		j = ic*n+i;
		y[j] = xn[j-n] - xn[j]+an[j]*xn[j-1];
	    }
	}
	y[(nc-1)*n] = xn[(nc-2)*n];
	for (i=1; i < n; i++) {
	    j = (nc-1)*n+i;
	    y[j] = xn[j-n] - x0[i]+an[j]*x0[i-1];
	}
    } else {
	y[0] = x2[0];
	for (i=1; i < n; i++) {
	    y[i] = x2[i] - x0[i]+an[i]*x0[i-1];
	}
    }
}

void pchain(sf_complex *y)
/*< apply the chain >*/
{
    int i, ic;
    sf_complex *tmp;

    tmp = sf_complexalloc(n);

     for (i=0; i < n; i++) {
	tmp[i] = x0[i];
    }
    for (ic=nc-1; ic >= 0; ic--) {
	y[0] = sf_cmplx(0.0f,0.0f);
	for (i=1; i < n; i++) {
	    y[i] = tmp[i] - an[ic*n+i]*tmp[i-1];
	}
	for (i=0; i < n; i++) {
	    tmp[i] = y[i];
	}
    }
    
    free(tmp);
}

void pchain_lop (bool adj, bool add, int nx, int ny, sf_complex* x, sf_complex* y) 
/*< linear operator >*/
{
    int ic, i, j, k;

    if (nx != (2*nc-1)*n || ny != nc*n) sf_error("%s: Wrong size",__FILE__);

    sf_cadjnull(adj,add,nx,ny,x,y);

    if (adj) {
	if (nc > 1) {
	    for (i=1; i < n; i++) {
		x[i+n*nc] -= y[i];
		x[i+n*nc-1] += conjf(an[i])*y[i];
	    }
	    for (ic=1; ic < nc-1; ic++) {
		x[(ic+nc-1)*n] += y[ic*n];
		for (i=1; i < n; i++) {
		    j = ic*n+i;
		    k = j+n*nc;
		    x[k-n] += y[j];
		    x[k]   -= y[j];
		    x[k-1] += conjf(an[j])*y[j]; 
		}
	    } 
	    for (i=0; i < n; i++) {
		j = (nc-1)*n+i;
		k = j+n*nc;
		x[k-n] += y[j];
	    } 
	} 
	/* delta A */
	for (ic=0; ic < nc-1; ic++) {
	    for (i=1; i < n; i++) {
		j = ic*n+i;
		x[j] += y[j]*conjf(xn[j-1]);
	    }
	}
	for (i=1; i < n; i++) {
	    j = (nc-1)*n+i;
	    x[j] += y[j]*conjf(x0[i-1]);
	}
    } else {
	if (nc > 1) {
	    for (i=1; i < n; i++) {
		y[i] -= x[i+n*nc]-an[i]*x[i+n*nc-1];
	    }
	    for (ic=1; ic < nc-1; ic++) {
		y[ic*n] += x[(ic+nc-1)*n];
		for (i=1; i < n; i++) {
		    j = ic*n+i;
		    k = j+n*nc;
		    y[j] += x[k-n] - x[k]+an[j]*x[k-1]; 
		}
	    } 
	    for (i=0; i < n; i++) {
		j = (nc-1)*n+i;
		k = j+n*nc;
		y[j] += x[k-n];
	    } 
	} 
	/* delta A */
	for (ic=0; ic < nc-1; ic++) {
	    for (i=1; i < n; i++) {
		j = ic*n+i;
		y[j] += x[j]*xn[j-1];
	    }
	}
	for (i=1; i < n; i++) {
	    j = (nc-1)*n+i;
	    y[j] += x[j]*x0[i-1];
	}
    }
}

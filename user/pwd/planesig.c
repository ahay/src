/* Signal separation with plane-wave destruction */
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

#include "allp2.h"

static int n, nc;
static float *eps, *tmp;
static allpas2 *sig;

void planesig_init (int nc1          /* number of dip components */,
		    int nw           /* filter size */, 
		    int nj           /* dealiasing stretch */, 
		    int nx, int ny   /* data size */, 
		    bool drift       /* if shift filter */,
		    float ***ss      /* signal slopes */, 
		    float *eps1      /* regularization parameter */)
/*< initialize >*/
{
    int ic;

    nc = nc1;
    eps = eps1;

    sig = (allpas2 *) sf_alloc(nc,sizeof(allpas2));

    for (ic=0; ic < nc; ic++) {
	sig[ic] = allpass2_init(nw,nj,nx,ny,drift,ss[ic]);
    }

    n = nx*ny;
    tmp = sf_floatalloc(n);
}

void planesig_rhs(float *d, float *dd)
/*< compute right-hand side >*/
{
    int ic, i;
    
    for (ic=0; ic < nc; ic++) {
	allpass22_init(sig[ic]);
	allpass21_lop (false, false, n, n, d, dd+ic*n);
    }
    for (i=0; i < nc*n; i++) {
	dd[nc*n+i] = 0.0f;
    }
}

void planesig_lop (bool adj, bool add, int ns, int nd, float *s, float *d)
/*< linear operator >*/
{
    int i, ic, jc;

    if (ns != nc*n || nd != 2*nc*n) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,ns,nd,s,d);

    for (ic=0; ic < nc; ic++) {		
	allpass22_init(sig[ic]);
	for (jc=0; jc < nc; jc++) {
	    if (jc != ic) {		
		allpass21_lop (adj, true, n, n, s+jc*n, d+ic*n);
	    }
	}
	if (adj) {
	    for (i=0; i < n; i++) {
		tmp[i] = eps[ic]*d[(nc+ic)*n+i];
	    }
	    allpass21_lop (true, true, n, n, s+ic*n, tmp);
	} else {
	    allpass21_lop (false, false, n, n, s+ic*n, tmp);
	    for (i=0; i < n; i++) {
		d[(nc+ic)*n+i] += eps[ic]*tmp[i];
	    }
	}
    }
}

void planesig_close(void)
/*< free allocated storage >*/
{
    int ic;

    for (ic=0; ic < nc; ic++) {
	free(sig[ic]);
    }
    free(sig);
    free(tmp);
}

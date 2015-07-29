/* Testing out-of-core patching */
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

#include <stdio.h>
/*^*/

#include <rsf.h>
/*^*/

#include "ocparcel.h"
#include "ocpatch.h"

static int np, nw;
static float *tmp;

void ocparcel_init(int dim     /* number of dimensions */, 
		   int *npatch /* number of patches [dim] */, 
		   int *nwall  /* data size [dim] */, 
		   int *nwind  /* patch size [dim] */)
/*< initialize >*/
{
    int i;

    np = 1;
    nw = 1;  
    for (i=0; i < dim; i++) {
	np *= npatch[i]; /* compute number of patches */
	nw *= nwind[i];  /* compute window size */
    }
    
    ocpatch_init (dim, nw, np, npatch, nwall, nwind);
    tmp = sf_floatalloc(nw);
}

void ocparcel_close(void)
/*< free allocated storage >*/
{
    free (tmp);
}

void ocparcel_lop(bool adj    /* adjoint flag */, 
		  int n       /* total data size */, 
		  int mw      /* total patch size */, 
		  FILE* wall  /* out-of-core weighting */, 
		  float* wind /* patch [mw] */)
/*< parcel operator >*/
{
    int ip, iw;

    if (mw != np*nw) sf_error("%s: wrong dimensions",__FILE__);
  
    for (ip=0; ip < np; ip++) {
	if (adj) {
	    ocpatch_lop (ip, false, wall, tmp);
	    for (iw=0; iw < nw; iw++) {
		tmp[iw] += wind[iw+ip*nw];
	    }
	    ocpatch_lop (ip, true, wall, tmp);
	} else {
	    ocpatch_lop (ip, false, wall, wind+ip*nw);
	}
    }
}

/* Out-of-core  wall weighting from window weighting */
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

#include "ocmkwallwt.h"
#include "ocpatch.h"
#include "oc.h"

void ocmkwallwt(bool inv      /* if compute 1/weight */, 
		int dim       /* number of dimensions */, 
		int* npatch   /* number of patches [dim] */, 
		int* nwall    /* data size [dim] */, 
		int* nwind    /* patch size [dim] */, 
		float* windwt /* window weighting */, 
		FILE* wallwt  /* out-of-core wall weighting */)
/*< make wall weight >*/
{
    int j, iw, ip, np, nw;
    off_t n;
    float *tmp;
    
    np = 1; 
    n = sizeof(float);
    nw = 1;

    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
    }

    tmp = sf_floatalloc(nw);

    oc_zero (n, wallwt);
    ocpatch_init(dim, nw, np, npatch, nwall, nwind);
  
    for (ip=0; ip < np; ip++) {
	ocpatch_lop (ip, false, wallwt, tmp);
	for (iw=0; iw < nw; iw++) {
	    tmp[iw] += windwt[iw];
	}
	ocpatch_lop(ip, true, wallwt, tmp);
    }

    if (inv) oc_invert(n,wallwt);
}

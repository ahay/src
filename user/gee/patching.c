/* Apply a linear operator in patches */
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

#include "patching.h"
#include "patch.h"
#include "mkwallwt.h"

void patching(sf_operator oper /* operator */, 
	      float* modl      /* input */, 
	      float* data      /* output */, 
	      int dim          /* number of dimensions */, 
	      int* npatch      /* number of patches [dim] */, 
	      int* nwall       /* data size [dim] */, 
	      int* nwind       /* patch size [dim] */, 
	      float* windwt    /* window weight */)
/*< patch an operator >*/
{
    float *winmodl, *windata, *wallwt;
    int i, j, iw, ip, np, n, nw;

    np = n = nw = 1; 
    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
    }
  
    winmodl = sf_floatalloc(nw);
    windata = sf_floatalloc(nw);
    wallwt  = sf_floatalloc(n);

    for (i=0; i < n; i++) data[i] = 0.;

    patch_init(dim, npatch, nwall, nwind);
    for (ip = 0; ip < np; ip++) {
	/* modl -> winmodl */
	patch_lop(false, false, n, nw, modl, winmodl);
	/* winmodl -> windata */
	oper(false, false, nw, nw, winmodl, windata);
	/* apply window weighting */
	for (iw=0; iw < nw; iw++) windata[iw] *= windwt[iw];
	/* data <- windata */
	patch_lop(true, true, n, nw, data, windata);
	patch_close();
    }

    /* windwt -> wallwt */
    mkwallwt(dim, npatch, nwall, nwind, windwt, wallwt);

    /* apply wall weighting */
    for (i=0; i < n; i++) data[i] *= wallwt[i];

    free (winmodl);
    free (windata);
    free (wallwt);
}


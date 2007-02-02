/* Local prediction-error filter estimation in patches */
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

#include "lopef.h"
#include "patch.h"
#include "misinput.h"
#include "pef.h"

void find_lopef(int dim      /* number of dimensions */, 
		float *wall  /* whole data */, 
		sf_filter aa /* PEF */, 
		int *npatch  /* number of patches [dim] */, 
		int *nwall   /* data size [dim] */, 
		int *nwind   /* patch size [dim] */, 
		float *mask  /* mask for known data */)
/*< estimate local PEF >*/ 
{
    float *windata, *winmask;
    int *known;
    int mis, ih, ip, iw, j, nw, np, n, nh;
    sf_filter bb;

    nw = np = n = 1;
    for (j=0; j < dim; j++) {
	n *= nwall[j];
	nw *= nwind[j];
	np *= npatch[j];
    }
    nh = aa->nh;

    windata = sf_floatalloc(nw); 
    if (NULL != mask) {
	winmask = sf_floatalloc(nw);
	known = sf_intalloc(nw);
    } else {
	winmask = NULL;
	known = NULL;
    }

    patch_init(dim, npatch, nwall, nwind);
    for (ip=0; ip < np; ip++) {
	bb = aa+ip;

	patch_lop(false, false, n, nw, wall, windata);
	if (NULL != mask) {
	    patch_lop(false, false, n, nw, mask, winmask);
	    for (iw=0; iw < nw; iw++) {
		known[iw] = (winmask[iw] != 0.);
	    }
	    find_mask(nw, known, bb);
	}
	for (mis=iw=0; iw < nw; iw++) {
	    if (!bb->mis[iw]) mis++;
	}
	if (mis > nh) { /* enough equations */
	    find_pef(nw, windata, bb, nh);
	} else if (ip > 1) { /* use last PEF */
	    for (ih=0; ih < nh; ih++) {
		bb->flt[ih] = (bb-1)->flt[ih];
	    }
	}
	patch_close();
    }

    free (windata);
    if (NULL != mask) {
	free (winmask);
	free (known);
    }
}


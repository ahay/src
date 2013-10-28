/* Local matched-filter estimation in patches */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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
#include <rsfgee.h>

#include "match.h"
#include "matchL1.h"

#include "lomatch.h"

void find_lomatch(int dim      /* number of dimensions */,
		float *wall  /* whole dataset [master]*/,
		float *mwall /* whole match dataset [slave]*/,
		sf_filter aa /* matched filter */,
		int *npatch  /* number of patches [dim] */, 
		int *nwall   /* data size [dim] */, 
		int *nwind   /* patch size [dim] */, 
		float *mask  /* mask for known data */,
		float eps    /* regularization parameter */,
		int liter    /* number of linear iteration [L2]*/,
		int niter    /* number of POCS iterations [L1]*/,
		float perc   /* percentage for sharpening [L1]*/,
		bool verb    /* verbosity flag */)
/*< estimate local matched filter>*/
{
    float *windata, *winmatch,*winmask;
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
    winmatch = sf_floatalloc(nw);

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
	patch_lop(false, false, n, nw, mwall, winmatch);

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
	    if (niter==1) {
	    	find_match(nw, winmatch, windata, bb, liter,eps);
	    }
	    else {
	    	find_matchL1(nw, winmatch, windata, bb, liter,niter,perc,eps);
	    }
	} else if (ip > 1) { /* use last matched filter */
	    for (ih=0; ih < nh; ih++) {
		bb->flt[ih] = (bb-1)->flt[ih];
	    }
	}
	patch_close();
	if (verb)
	sf_warning("\r\t\t\t\t\t\t\t\t\t %3.2f%% ",(float)100*(ip+1)/np);
    }

    free (windata);
    free (winmatch);
    if (NULL != mask) {
	free (winmask);
	free (known);
    }
}


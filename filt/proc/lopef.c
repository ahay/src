#include <rsf.h>

#include "lopef.h"
#include "patch.h"
#include "misinput.h"
#include "pef.h"

void find_lopef(int dim, float *wall, filter aa, 
		int *npatch, int *nwall, int *nwind, float *mask) 
{
    float *windata, *winmask;
    int *known;
    int mis, ih, ip, iw, j, nw, np, n, nh;
    filter bb;

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


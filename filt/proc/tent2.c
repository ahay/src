#include <math.h>

#include <rsf.h>

#include "tent2.h"

/* tent2
   -----
   cosine window weighting function - an alternative to tent
   dim - window dimensionality
   nwind[dim] - window size
   windwt[product(nwind)] - window weight (output) */ 
void tent2 (int dim, const int* nwind, float* windwt)
{
    int i, j, nw, x[SF_MAX_DIM];
    double w;

    /* compute window size */
    nw = 1;
    for (j=0; j < dim; j++) {
	nw *= nwind[j];
    }

    /* loop inside the windoe */
    for (i=0; i < nw; i++) { 
	sf_line2cart(dim, nwind, i, x);
    
	windwt[i] = 1.;
	for (j=0; j < dim; j++) {
	    if (nwind[j] > 1) {
		w = cosf(2.*SF_PI*(x[j]+1.)/(nwind[j] + 1.));
		w = 0.5*(1.-w);
		if (w > 0.) 
		    windwt[i] *= w;
		else
		    windwt[i] = 0.;
	    }
	}
    }
}



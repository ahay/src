#include <math.h>

#include <rsf.h>

#include "tent.h"

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

/* tent
   -----
   tent-like window weighting function 
   dim - window dimensionality
   nwind[dim] - window size
   windwt[product(nwind)] - window weight (output) */
void tent (int dim, const int* nwind, const int* center, const int* a, 
	   float* windwt)
{
    int i, j, nw, start[SF_MAX_DIM], end[SF_MAX_DIM], x[SF_MAX_DIM];
    float w, mid[SF_MAX_DIM], wid[SF_MAX_DIM];

    nw = 1;
    for (j=0; j < dim; j++) {
	start[j] = a[j]-center[j];
	end[j] = nwind[j]-center[j];
	mid[j]= (end[j]+start[j])/2.;
	wid[j]= (end[j]-start[j])/2.;
	nw *= nwind[j]; /* compute window size */
    }

    /* loop in the window */
    for (i=0; i < nw; i++) {
	sf_line2cart(dim, nwind, i, x);
    
	windwt[i] = 1.;
	for (j=0; j < dim; j++) {
	    if (x[j] >= start[j] && x[j] <= end[j]) {
		w = (x[j]-mid[j])/wid[j];
		windwt[i] *= MAX(0.,1.-fabs(w));
	    }	else {
		windwt[i] = 0.;
	    }
	}
    }
}

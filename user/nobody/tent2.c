#include <math.h>

#include "cartesian.h"
#include "tent2.h"

/* tent2
   -----
   cosine window weighting function - an alternative to tent
   dim - window dimensionality
   nwind[dim] - window size
   windwt[product(nwind)] - window weight (output) */ 
void tent2 (int dim, const int* nwind, float* windwt)
{
  int i, j, nw, x[MAX_DIM];
  double w, pi;

  /* compute window size */
  nw = 1;
  for (j=0; j < dim; j++) {
    nw *= nwind[j];
  }

  pi = acos(-1.);
  /* loop inside the windoe */
  for (i=0; i < nw; i++) { 
    line2cart(dim, nwind, i, x);
    
    windwt[i] = 1.;
    for (j=0; j < dim; j++) {
      if (nwind[j] > 1) {
	w = cos(2.*pi*(x[j]+1.)/(nwind[j] + 1.));
	w = 0.5*(1.-w);
	if (w > 0.) 
	  windwt[i] *= w;
	else
	  windwt[i] = 0.;
      }
    }
  }
}



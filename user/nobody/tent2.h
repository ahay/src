#ifndef _tent2_h
#define _tent2_h

/* tent2
   -----
   cosine window weighting function - an alternative to tent
   dim - window dimensionality
   nwind[dim] - window size
   windwt[product(nwind)] - window weight (output) */ 
void tent2 (int dim, const int* nwind, float* windwt);

#endif

#include "interp.h"

void bin_int (float x, int n, float* w) 
{
  int i;

  w[0] = 1.;
  
  for (i = 1; i < n; i++)
      w[i] = 0.;
}
    
void lin_int (float x, int n, float* w) 
{
    int i;
    
    if (n == 1) {
	w[0] = 1.;
    } else {
	w[1] = x;
	w[0] = 1. - x;
	for (i = 2; i < n; i++)
	    w[i] = 0.;
    }
}

void lg_int (float x, int n, float* w) 
{
    int i, j, nc;
    float f, xi;

    nc = (n-1)*0.5;
    for (i=0; i < n; i++) {
	f = 1.;
	xi = x + nc - i;
	for (j=0; j < n; j++) {
	    if (i != j) f *= (1. + xi / (i - j));
	}
	w[i] = f;
    }
}

  
/* 	$Id$	 */


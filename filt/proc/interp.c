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
  
/* 	$Id: interp.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */


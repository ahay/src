#include "cstd.h"

void acd_3d_2(float *** uc, 
	      float *** up, 
	      float *** csq, 
	      int * s, 
	      int * e, 
	      float c0, 
	      float * c1) {

  int i0, i1, i2;
  
  for (i2=s[2]; i2<=e[2]; i2++) {
    for (i1=s[1]; i1<=e[1]; i1++) {
      for (i0=s[0]; i0<=e[0]; i0++) {
	up[i2][i1][i0] = 2.0*uc[i2][i1][i0] - up[i2][i1][i0] +
	  csq[i2][i1][i0] * 
	  ( c0*uc[i2][i1][i0] +
	    c1[0]*(uc[i2][i1][i0+1] + uc[i2][i1][i0-1]) +
	    c1[1]*(uc[i2][i1+1][i0] + uc[i2][i1-1][i0]) +
	    c1[2]*(uc[i2+1][i1][i0] + uc[i2-1][i1][i0]) 
	    );
      }
    }
  }

  /* boundary condns: up[-1][i]=0 etc. these do not need to be enforced as they are 
  // a by-product of initialization and involve only non-updated allocated array 
  // entries */
}


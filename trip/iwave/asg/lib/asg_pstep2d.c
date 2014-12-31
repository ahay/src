#include "utils.h"

void asg_pstep2d(float ** restrict bulk, 
		 float ** restrict p0, float ** restrict p1,
		 float ** restrict v0, float ** restrict v1,
		 float ** restrict ep, float ** restrict epp,
		 float * restrict sdiv,
		 int * gsc, int * gec, 
		 int maxoff, float * restrict c[RARR_MAX_NDIM]) {
  
  int i0, i1;
  int ioff;

  for (i1=gsc[1]; i1 <= gec[1]; i1++) {
    for (i0=gsc[0]; i0 <= gec[0]; i0++) sdiv[i0]=0.0f;
    for (ioff = 0; ioff<maxoff; ioff++) {
      for (i0=gsc[0]; i0 <= gec[0]; i0++) {
	sdiv[i0] += 
	  (c[0][ioff]*(v0[i1][i0+ioff]-v0[i1][i0-ioff-1]) +
	   c[1][ioff]*(v1[i1+ioff][i0]-v1[i1-ioff-1][i0]));
      }
    }
    for (i0=gsc[0]; i0 <= gec[0]; i0++) {
      p0[i1][i0] = (p0[i1][i0]*ep[0][i0] + bulk[i1][i0]*sdiv[i0])*epp[0][i0];
      p1[i1][i0] = (p1[i1][i0]*ep[1][i1] + bulk[i1][i0]*sdiv[i0])*epp[1][i1];
    }
  }
}

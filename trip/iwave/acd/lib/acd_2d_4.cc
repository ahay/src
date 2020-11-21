#include "cstd.h"

void acd_2d_4(float ** uc, 
	      float ** up, 
	      float ** csq, 
	      int * s, 
	      int * e, 
	      float c0, 
	      float * c1,
	      float * c2,
	      int * lbc,
	      int * rbc) {

  int i0, i1;
  int s0=s[0];
  int e0=e[0];
  for (i1=s[1]; i1<=e[1]; i1++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
    for (i0=s0; i0<=e0; i0++) {
      float lap = c0*uc[i1][i0] +
        c1[0]*(uc[i1][i0+1] + uc[i1][i0-1]) +
        c1[1]*(uc[i1+1][i0] + uc[i1-1][i0]) +
        c2[0]*(uc[i1][i0+2] + uc[i1][i0-2]) +
        c2[1]*(uc[i1+2][i0] + uc[i1-2][i0]);
      up[i1][i0] = 2.0*uc[i1][i0] - up[i1][i0] +
        csq[i1][i0] * lap;
    }
  }
  /* boundary conditions - note that uc[-1][i]=0 etc. */
  if (lbc[1]) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
    for (i0=s0;i0<=e0;i0++) {
      up[s[1]-2][i0]=-up[s[1]][i0];
    }
  }
  if (rbc[1]) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
    for (i0=s0;i0<=e0;i0++) {
      up[e[1]+2][i0]=-up[e[1]][i0];
    }
  }
  if (lbc[0]) {
    for (i1=s[1];i1<=e[1];i1++) {
      up[i1][s[0]-2]=-up[i1][s[0]];
    }
  }
  if (rbc[0]) {
    for (i1=s[1];i1<=e[1];i1++) {
      up[i1][e[0]+2]=-up[i1][e[0]];
    }
  }
}



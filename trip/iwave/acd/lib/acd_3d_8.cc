#include "cstd.h"

void acd_3d_8(float *** uc, 
	      float *** up, 
	      float *** csq, 
	      int * s, 
	      int * e, 
	      float c0, 
	      float * c1,
	      float * c2,
	      float * c3,
	      float * c4,
	      int * lbc,
	      int * rbc) {

  int i0, i1, i2;

  for (i2=s[2]; i2<=e[2]; i2++) {
    for (i1=s[1]; i1<=e[1]; i1++) {
      for (i0=s[0]; i0<=e[0]; i0++) {
	up[i2][i1][i0] = 2.0*uc[i2][i1][i0] - up[i2][i1][i0] +
	  csq[i2][i1][i0] * 
	  ( c0*uc[i2][i1][i0]	+
	    c1[0]*(uc[i2][i1][i0+1] + uc[i2][i1][i0-1]) +
	    c1[1]*(uc[i2][i1+1][i0] + uc[i2][i1-1][i0]) +
	    c1[2]*(uc[i2+1][i1][i0] + uc[i2-1][i1][i0]) +
	    c2[0]*(uc[i2][i1][i0+2] + uc[i2][i1][i0-2]) +
	    c2[1]*(uc[i2][i1+2][i0] + uc[i2][i1-2][i0]) +
	    c2[2]*(uc[i2+2][i1][i0] + uc[i2-2][i1][i0]) +
	    c3[0]*(uc[i2][i1][i0+3] + uc[i2][i1][i0-3]) +
	    c3[1]*(uc[i2][i1+3][i0] + uc[i2][i1-3][i0]) +
	    c3[2]*(uc[i2+3][i1][i0] + uc[i2-3][i1][i0]) +
	    c4[0]*(uc[i2][i1][i0+4] + uc[i2][i1][i0-4]) +
	    c4[1]*(uc[i2][i1+4][i0] + uc[i2][i1-4][i0]) +
	    c4[2]*(uc[i2+4][i1][i0] + uc[i2-4][i1][i0])
	    );
      }
    }
  }
  /* boundary conditions - note that uc[-1][i]=0 etc. */
  if (lbc[2]) {
    for (i1=s[1];i1<=e[1];i1++) {
      for (i0=s[0];i0<=e[0];i0++) {
	up[s[2]-2][i1][i0]=-up[s[2]+0][i1][i0];
	up[s[2]-3][i1][i0]=-up[s[2]+1][i1][i0];
	up[s[2]-4][i1][i0]=-up[s[2]+2][i1][i0];
      }
    }
  }
  if (rbc[2]) {
    for (i1=s[1];i1<=e[1];i1++) {
      for (i0=s[0];i0<=e[0];i0++) {
	up[e[2]+2][i1][i0]=-up[e[2]-0][i1][i0];
	up[e[2]+3][i1][i0]=-up[e[2]-1][i1][i0];
	up[e[2]+4][i1][i0]=-up[e[2]-2][i1][i0];
      }
    }
  }
  if (lbc[1]) {
    for (i2=s[2];i2<=e[2];i2++) {
      for (i0=s[0];i0<=e[0];i0++) {
	up[i2][s[1]-2][i0]=-up[i2][s[1]+0][i0];
	up[i2][s[1]-3][i0]=-up[i2][s[1]+1][i0];
	up[i2][s[1]-4][i0]=-up[i2][s[1]+2][i0];
      }
    }
  }
  if (rbc[1]) {
    for (i2=s[2];i2<=e[2];i2++) {
      for (i0=s[0];i0<=e[0];i0++) {
	up[i2][e[1]+2][i0]=-up[i2][e[1]-0][i0];
	up[i2][e[1]+3][i0]=-up[i2][e[1]-1][i0];
	up[i2][e[1]+4][i0]=-up[i2][e[1]-2][i0];
      }
    }
  }
  if (lbc[0]) {
    for (i2=s[2];i2<=e[2];i2++) {
      for (i1=s[1];i1<=e[1];i1++) {
	up[i2][i1][s[0]-2]=-up[i2][i1][s[0]+0];
	up[i2][i1][s[0]-3]=-up[i2][i1][s[0]+1];
	up[i2][i1][s[0]-4]=-up[i2][i1][s[0]+2];
      }
    }
  }
  if (rbc[0]) {
    for (i2=s[2];i2<=e[2];i2++) {
      for (i1=s[1];i1<=e[1];i1++) {
	up[i2][i1][e[0]+2]=-up[i2][i1][e[0]-0];
	up[i2][i1][e[0]+3]=-up[i2][i1][e[0]-1];
	up[i2][i1][e[0]+4]=-up[i2][i1][e[0]-2];
      }
    }

  }
  
}

#include "cstd.h"

/** eta = PML profiles, = 0 in physical region on each axis
    eta[axis][index] - note different lengths on different axes
    ep = 1 - dt^2 eta
    epp = 1/(1 + dt^2 eta) 
    precomputed arrays, stored in fdpars

    maxoff = half-order
    c[axis][offset] = FD coefficient arrays, stored in fdpars
*/

void asg_pstep3d(float *** restrict bulk, 
		 float *** restrict p0, float *** restrict p1, float *** restrict p2,
		 float *** restrict v0, float *** restrict v1, float *** restrict v2,
		 float ** restrict ep, float ** restrict epp,
		 float * restrict sdiv,
		 int * gsc, int * gec, 
		 int * lbc, int * rbc,
		 int maxoff, float ** restrict c) {
  
  int i0, i1, i2;
  int ioff;

  for (i2=gsc[2]; i2 <= gec[2]; i2++) {
    for (i1=gsc[1]; i1 <= gec[1]; i1++) {
      for (i0=gsc[0]; i0 <= gec[0]; i0++) sdiv[i0]=0.0f;
      for (ioff = 0; ioff<maxoff; ioff++) {
	for (i0=gsc[0]; i0 <= gec[0]; i0++) {
	  sdiv[i0] += 
	    (c[0][ioff]*(v0[i2][i1][i0+ioff]-v0[i2][i1][i0-ioff-1]) +
	     c[1][ioff]*(v1[i2][i1+ioff][i0]-v1[i2][i1-ioff-1][i0]) +
	     c[2][ioff]*(v2[i2+ioff][i1][i0]-v2[i2-ioff-1][i1][i0]));
	}
      }
      for (i0=gsc[0]; i0 <= gec[0]; i0++) {
	p0[i2][i1][i0] = (p0[i2][i1][i0]*ep[0][i0] - bulk[i2][i1][i0]*sdiv[i0])*epp[0][i0];
	p1[i2][i1][i0] = (p1[i2][i1][i0]*ep[1][i1] - bulk[i2][i1][i0]*sdiv[i0])*epp[1][i1];
	p2[i2][i1][i0] = (p2[i2][i1][i0]*ep[2][i2] - bulk[i2][i1][i0]*sdiv[i0])*epp[2][i2];
      }
    }
  }

  /* boundary conditions - p is odd about index just before/after comp domain */
  if (lbc[0]) {
    for (i2=gsc[2];i2<=gec[2];i2++) {
      for (i1=gsc[1];i1<=gec[1];i1++) {
	p0[i2][i1][gsc[0]-1]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p0[i2][i1][gsc[0]-ioff-1]=-p0[i2][i1][gsc[0]+ioff-1];
	}
      }
    }
  }
  if (rbc[0]) {
    for (i2=gsc[2];i2<=gec[2];i2++) {
      for (i1=gsc[1];i1<=gec[1];i1++) {
	p0[i2][i1][gec[0]+1]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p0[i2][i1][gec[0]+ioff+1]=-p0[i2][i1][gec[0]-ioff+1];
	}
      }
    }
  }
  if (lbc[1]) {
    for (i2=gsc[2];i2<=gec[2];i2++) {
      for (i0=gsc[0];i0<=gsc[0];i0++) {
	p1[i2][gsc[1]-1][i0]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p1[i2][gsc[1]-ioff-1][i0]=-p1[i2][gsc[1]+ioff-1][i0];
	}
      }
    }
  }
  if (rbc[1]) {
    for (i2=gsc[2];i2<=gec[2];i2++) {
      for (i0=gsc[0];i0<=gsc[0];i0++) {
	p1[i2][gec[1]+1][i0]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p1[i2][gec[1]+ioff+1][i0]=-p1[i2][gec[1]-ioff+1][i0];
	}
      }
    }
  }
  if (lbc[2]) {
    for (i1=gsc[1];i1<=gec[1];i1++) {  
      for (i0=gsc[0];i0<=gsc[0];i0++) {
	p2[gsc[2]-1][i1][i0]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p2[gsc[2]-ioff-1][i1][i0]=-p2[gsc[2]+ioff-1][i1][i0];
	}
      }
    }
  }
  if (rbc[2]) {
    for (i1=gsc[1];i1<=gec[1];i1++) {  
      for (i0=gsc[0];i0<=gsc[0];i0++) {
	p2[gec[2]+1][i1][i0]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p2[gec[2]+ioff+1][i1][i0]=-p2[gec[2]-ioff+1][i1][i0];
	}
      }
    }
  }
}

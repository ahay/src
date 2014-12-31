#include "cstd.h"

void asg_vstep2d(float ** restrict buoy,
		 float ** restrict p0,float ** restrict p1,
		 float ** restrict v0,float ** restrict v1,
		 float ** restrict ev, float ** restrict evp,
		 float ** restrict gradp,
		 int * gsc_v0, int * gec_v0,
		 int * gsc_v1, int * gec_v1,
		 int maxoff,float ** restrict c) {

  int i0, i1;
  int ioff;

  for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {
    for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) gradp[0][i0]=0;      
    for (ioff=0; ioff<maxoff; ioff++) {
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ )       
	gradp[0][i0] += c[0][ioff]*(p0[i1][i0+ioff+1]-p0[i1][i0-ioff]);
    }
    for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) 
      v0[i1][i0] = evp[0][i0]*v0[i1][i0] + ev[0][i0]*buoy[i1][i0]*gradp[0][i0];
  }

  for (i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
    for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) gradp[1][i0]=0;      
    for (ioff=0; ioff<maxoff; ioff++) {
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ )       
	gradp[1][i0] += c[1][ioff]*(p1[i1+ioff+1][i0]-p1[i1-ioff][i0]);
    }
    for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) 
      v1[i1][i0] = evp[1][i1]*v1[i1][i0] + ev[1][i1]*buoy[i1][i0]*gradp[1][i0];
  }

}

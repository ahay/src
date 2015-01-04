#include "cstd.h"

void asg_vstep3d(float *** restrict buoy,
		 float *** restrict p0,float *** restrict p1, float *** restrict p2,
		 float *** restrict v0,float *** restrict v1, float *** restrict v2,
		 float ** restrict ev, float ** restrict evp,
		 float ** restrict gradp,
		 int * gsc_v0, int * gec_v0,
		 int * gsc_v1, int * gec_v1,
		 int * gsc_v2, int * gec_v2,
		 int * lbc, int * rbc,
		 int maxoff,float ** restrict c) {

  int i0, i1, i2;
  int ioff;

  for (i2=gsc_v0[2]; i2 <= gec_v0[2]; i2++) {
    for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) gradp[0][i0]=0;      
      for (ioff=0; ioff<maxoff; ioff++) {
	for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ )       
	  gradp[0][i0] += c[0][ioff]*(p0[i2][i1][i0+ioff+1]-p0[i2][i1][i0-ioff]);
      }
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) 
	v0[i2][i1][i0] = evp[0][i0]*v0[i2][i1][i0] 
	  - ev[0][i0]*0.5*(buoy[i2][i1][i0]+buoy[i2][i1][i0+1])*gradp[0][i0];
    }
  }
  for (i2=gsc_v1[2]; i2 <= gec_v1[2]; i2++) {
    for (i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) gradp[1][i0]=0;      
      for (ioff=0; ioff<maxoff; ioff++) {
	for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ )       
	  gradp[1][i0] += c[1][ioff]*(p1[i2][i1+ioff+1][i0]-p1[i2][i1-ioff][i0]);
      }
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) 
	v1[i2][i1][i0] = evp[1][i1]*v1[i2][i1][i0] 
	  - ev[1][i1]*0.5*(buoy[i2][i1][i0]+buoy[i2][i1+1][i0])*gradp[1][i0];
    }
  }
  for (i2=gsc_v2[2]; i2 <= gec_v2[2]; i2++) {
    for (i1=gsc_v2[1]; i1 <= gec_v2[1]; i1++ ) {
      for (i0=gsc_v2[0]; i0 <= gec_v2[0]; i0++ ) gradp[2][i0]=0;      
      for (ioff=0; ioff<maxoff; ioff++) {
	for (i0=gsc_v2[0]; i0 <= gec_v2[0]; i0++ )       
	  gradp[2][i0] += c[2][ioff]*(p2[i2+ioff+1][i1][i0]-p2[i2-ioff][i1][i0]);
      }
      for (i0=gsc_v2[0]; i0 <= gec_v2[0]; i0++ ) 
	v2[i2][i1][i0] = evp[2][i2]*v2[i2][i1][i0] 
	  - ev[2][i2]*0.5*(buoy[i2][i1][i0]+buoy[i2+1][i1][i0])*gradp[2][i0];
    }
  }

  if (lbc[0]) {
    for (i2=gsc_v0[2]; i2 <= gec_v0[2]; i2++) {
      for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
	for (ioff=1; ioff<maxoff;ioff++) {
	  v0[i2][i1][gsc_v0[0]-ioff]=v0[i2][i1][gsc_v0[0]+ioff-1];
	}
      }
    }
  }
  if (rbc[0]) {
    for (i2=gsc_v0[2]; i2 <= gec_v0[2]; i2++) {
      for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
	for (ioff=1; ioff<maxoff;ioff++) {
	  v0[i2][i1][gec_v0[0]+ioff]=v0[i2][i1][gec_v0[0]-ioff+1];
	}
      }
    }
  }
  if (lbc[1]) {
    for (i2=gsc_v1[2]; i2 <= gec_v1[2]; i2++) {
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	for (ioff=1; ioff<maxoff;ioff++) {
	  v1[i2][gsc_v1[1]-ioff][i0]=v1[i2][gsc_v1[1]+ioff-1][i0];
	}
      }
    }
  }
  if (rbc[1]) {
    for (i2=gsc_v1[2]; i2 <= gec_v1[2]; i2++) {
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	for (ioff=1; ioff<maxoff;ioff++) {
	  v1[i2][gec_v1[1]+ioff][i0]=v1[i2][gec_v1[1]-ioff+1][i0];
	}
      }
    }
  }
  if (lbc[2]) {
    for (i1=gsc_v2[1]; i1 <= gec_v2[1]; i1++) {
      for (i0=gsc_v2[0]; i0 <= gec_v2[0]; i0++) {
	for (ioff=1;ioff<maxoff;ioff++) {
	  v2[gsc_v2[2]-ioff][i1][i0] = v2[gsc_v2[2]+ioff-1][i1][i0];
	}
      }
    }
  }
  if (rbc[2]) {
    for (i1=gsc_v2[1]; i1 <= gec_v2[1]; i1++) {
      for (i0=gsc_v2[0]; i0 <= gec_v2[0]; i0++) {
	for (ioff=1;ioff<maxoff;ioff++) {
	  v2[gec_v2[2]+ioff][i1][i0] = v2[gec_v2[2]-ioff+1][i1][i0];
	}
      }
    }
  }
}

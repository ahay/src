#include <rsf.h>

#include "allp3.h"
#include "apfilt.h"

struct Allpass {
  int nx, ny, nz, nw, nj;
  float*** pp;
};

allpass allpass_init(int nw, int nj, int nx, int ny, int nz, float ***pp)
{
  allpass ap;

  ap = (allpass) sf_alloc(1,sizeof(*ap));

  ap->nw = nw;
  ap->nj = nj;
  ap->nx = nx;
  ap->ny = ny;
  ap->nz = nz;
  ap->pp = pp;

  return ap;
}

void allpass1 (bool der, const allpass ap, float*** xx, float*** yy)
{
  int ix, iy, iz, iw, is;
  float a[7];

  for (iz=0; iz < ap->nz; iz++) {
      for (iy=0; iy < ap->ny; iy++) {
	  for (ix=0; ix < ap->nx; ix++) {
	      yy[iz][iy][ix] = 0.;
	  }
      }
  }
  
  for (iz=0; iz < ap->nz; iz++) {
      for (iy=0; iy < ap->ny-1; iy++) {
	  for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
	      if (der) {
		  aderfilter(ap->nw, ap->pp[iz][iy][ix], a);
	      } else {
		  passfilter(ap->nw, ap->pp[iz][iy][ix], a);
	      }
	      
	      for (iw = 0; iw <= 2*ap->nw; iw++) {
		  is = (iw-ap->nw)*ap->nj;
		  
		  yy[iz][iy][ix] += (xx[iz][iy+1][ix+is] - 
				     xx[iz][iy  ][ix-is]) * a[iw];
	      }
	  }
      }
  }
}

void allpass2 (bool der, const allpass ap, float*** xx, float*** yy)
{
    int ix, iy, iz, iw, is;
    float a[7];
    
    for (iz=0; iz < ap->nz; iz++) {
	for (iy=0; iy < ap->ny; iy++) {
	    for (ix=0; ix < ap->nx; ix++) {
		yy[iz][iy][ix] = 0.;
	    }
	}
    }
    
    for (iz=0; iz < ap->nz-1; iz++) {
	for (iy=0; iy < ap->ny; iy++) {
	    for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
		if (der) {
		    aderfilter(ap->nw, ap->pp[iz][iy][ix], a);
		} else {
		    passfilter(ap->nw, ap->pp[iz][iy][ix], a);
		}
		
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    
		    yy[iz][iy][ix] += (xx[iz+1][iy][ix+is] - 
				       xx[iz  ][iy][ix-is]) * a[iw];
		}
	    }
	}
    }
}

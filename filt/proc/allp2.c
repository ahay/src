#include <rsf.h>

#include "allp2.h"
#include "apfilt.h"

struct Allpass2 {
    int nx, ny, nw, nj;
    float** pp;
};

static allpass2 ap2;

allpass2 allpass2_init(int nw, int nj, int nx, int ny, float **pp)
{
    allpass2 ap;
    
    ap = (allpass2) sf_alloc(1,sizeof(*ap));
    
    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->pp = pp;
    
    return ap;
}

void allpass22_init (allpass2 ap1)
{
    ap2 = ap1;
}

void allpass21_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy)
{
    int i, ix, iy, iw, is, nx, ny;
    float a[7];

    sf_adjnull(adj,add,n1,n2,xx,yy);
  
    nx = ap2->nx;
    ny = ap2->ny;

    for (iy=0; iy < ny-1; iy++) {
	for (ix = ap2->nw*ap2->nj; ix < nx-ap2->nw*ap2->nj; ix++) {
	    passfilter(ap2->nw, ap2->pp[iy][ix], a);
	    i = ix + iy*nx;
	      
	    for (iw = 0; iw <= 2*ap2->nw; iw++) {
		is = (iw-ap2->nw)*ap2->nj;
		  
		if (adj) {
		    xx[i+is+nx] += yy[i]*a[iw];
		    xx[i-is]    -= yy[i]*a[iw];
		} else {
		    yy[i] += (xx[i+is+nx] - xx[i-is]) * a[iw];
		}
	    }
	}
    }
}

void allpass21 (bool der, const allpass2 ap, float** xx, float** yy)
{
    int ix, iy, iw, is;
    float a[7];

    for (iy=0; iy < ap->ny; iy++) {
	for (ix=0; ix < ap->nx; ix++) {
	    yy[iy][ix] = 0.;
	}
    }
  
    for (iy=0; iy < ap->ny-1; iy++) {
	for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
	    if (der) {
		aderfilter(ap->nw, ap->pp[iy][ix], a);
	    } else {
		passfilter(ap->nw, ap->pp[iy][ix], a);
	    }
	      
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = (iw-ap->nw)*ap->nj;
		  
		yy[iy][ix] += (xx[iy+1][ix+is] - 
			       xx[iy  ][ix-is]) * a[iw];
	    }
	}
    }
}

/* 	$Id$	 */

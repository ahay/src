#include <rsf.h>

#include "allp2.h"
#include "apfilt.h"

struct Allpass2 {
    int nx, ny, nw, nj;
    float** pp;
};

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

/* 	$Id: allp2.c,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */

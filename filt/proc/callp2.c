#include <rsf.h>

#include "callp2.h"
#include "apfilt.h"

struct Callpass2 {
    int nx, ny, nw, nj;
    float a[7], d[7];
};

callpass2 callpass2_init(int nw, int nj, int nx, int ny)
{
    callpass2 ap;
    
    ap = (callpass2) sf_alloc(1,sizeof(*ap));
    
    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
   
    return ap;
}

void callpass21_set (callpass2 ap, float p)
{
    aderfilter(ap->nw, p, ap->d);
    passfilter(ap->nw, p, ap->a);
}

void callpass21 (bool der, const callpass2 ap, float** xx, float** yy)
{
    int ix, iy, iw, is;
    float *b;

    for (iy=0; iy < ap->ny; iy++) {
	for (ix=0; ix < ap->nx; ix++) {
	    yy[iy][ix] = 0.;
	}
    }
  
    b = der? ap->d: ap->a;

    for (iy=0; iy < ap->ny-1; iy++) {
	for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = (iw-ap->nw)*ap->nj;
		  
		yy[iy][ix] += (xx[iy+1][ix+is] - 
			       xx[iy  ][ix-is]) * b[iw];
	    }
	}
    }
}

/* 	$Id: callp2.c 704 2004-07-13 18:22:06Z fomels $	 */

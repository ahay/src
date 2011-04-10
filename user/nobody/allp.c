#include <rsf.h>

#include "allp.h"
#include "apfilt.h"

static int n1, n2, nw, nj;
static float** pp;

void allp_init(int nw1, int nj1, int m1, int m2, float **pp1)
{
    nw = nw1;
    nj = nj1;
    n1 = m1;
    n2 = n2;
    pp = pp1;
}

void allp_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
{
    int i1, i2, iw, is, i;
    float a[7];

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (i2=0; i2 < n2-1; i2++) {
	for (i1 = nw*nj; i1 < n1-nw*nj; i1++) {
	    i = i1 + i2*n1;

	    passfilter(nw, pp[i2][i1], a);
	      
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj;
		  
		if (adj) {
		    xx[i+is+n1] += yy[i] * a[iw];
		    xx[i-is]    -= yy[i] * a[iw];
		} else {
		    yy[i] += (xx[i+is+n1] - xx[i-is]) * a[iw];
		}
	    }
	}
    }
}

/* 	$Id$	 */

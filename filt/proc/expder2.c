#include <rsf.h>

#include "expder2.h"

static int n1, n2;
static float *a, *b, *c, *d;

void expder2_init(int m1, int m2, float **aa)
{
    n1 = m1;
    n2 = m2;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];
}

void expder2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    
    int i, j, k;

    if (ny != 4*nx) sf_error("%s: size error: %d != 4*%d",__FILE__,ny,nx);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (j=0; j < n2; j++) {
	for (i=4; i < n1; i++) {
	    k = i + j*n1;
	
	    if (adj) {
		xx[k-1] += b[k]*yy[k] + a[k]*yy[k+nx] + 
		    d[k]*yy[k+2*nx] + c[k]*yy[k+3*nx];
		xx[k-2] -= (2.*b[k]*d[k]*(c[k]*yy[k] + a[k]*yy[k+2*nx]) +
			    (b[k]+2.*a[k]*c[k]*d[k])*yy[k+nx] +
			    (d[k]+2.*a[k]*c[k]*b[k])*yy[k+3*nx]);
		xx[k-3] += b[k]*d[k]*(d[k]*yy[k]+b[k]*yy[k+2*nx]) +
		    d[k]*(a[k]*d[k]+2.*b[k]*c[k])*yy[k+nx] +
		    b[k]*(c[k]*b[k]+2.*d[k]*a[k])*yy[k+3*nx];
		xx[k-4] -= b[k]*d[k]*(d[k]*yy[k+nx] + b[k]*yy[k+3*nx]);
	    } else {
		yy[k] += b[k]*(xx[k-1] + 
			       d[k]*(d[k]*xx[k-3] - 2.*c[k]*xx[k-2])); 
		yy[k+nx] += a[k]*xx[k-1]-(b[k]+2.*a[k]*c[k]*d[k])*xx[k-2] +
		    d[k]*((a[k]*d[k]+2.*b[k]*c[k])*xx[k-3] - 
			  d[k]*b[k]*xx[k-4]);
		yy[k+2*nx] += d[k]*(xx[k-1] + 
				    b[k]*(b[k]*xx[k-3] - 2.*a[k]*xx[k-2]));
		yy[k+3*nx] += c[k]*xx[k-1]-(d[k]+2.*a[k]*c[k]*b[k])*xx[k-2] +
		    b[k]*((c[k]*b[k]+2.*d[k]*a[k])*xx[k-3] - 
			  d[k]*b[k]*xx[k-4]);
	    }
	}
    }
}

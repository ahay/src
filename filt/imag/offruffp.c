#include <rsf.h>

#include "offruffp.h"  

static int nh, nx, num;
static float h0, dh;
static float complex c1, c2;

void offruffp_init (float h0_in, int nh_in, float dh_in,
		    int nx_in, float dx, float w, int num_in)
{
    float w2;

    h0 = h0_in/dx; 
    nh = nh_in; 
    dh = dh_in/dx; 
    nx = nx_in; 
    w2  = w*w; 
    num = num_in;
    c1 = 3.*(9. + w2  - 4.*w*I)/(w2*(3.+w*I));
    c2 = 3.*(w2 - 27. - 8.*w*I)/(w2*(3.+w*I));
}

float complex offruffp_c1 (void)
{
    return c1;
}

float complex offruffp_c2 (void)
{
    return c2;
}

void hderp_lop (bool adj, bool add, int n1, int n2, 
		complex float *x, complex float *y) {
    int ih, ix, i;

    sf_cadjnull (adj,add,n1,n2,x,y);

    for (ih=1; ih < nh; ih++) {
	for (ix=0; ix < nx; ix++) {
	    i = ix + ih*nx;
	    if (adj) {
		x[i]    += y[i];
		x[i-nx] -= y[i];
	    } else {
		y[i] += x[i] - x[i-nx];
	    }
	}
    }
}

void offruffp_lop (bool adj, bool add, int n1, int n2, 
		   complex float *x, complex float *y) {
    int ih, ix, i;
    float h1, h2;
    float complex diag1, diag2, offd1, offd2;

    sf_cadjnull (adj,add,n1,n2,x,y);
    
    for (ih=1; ih < nh; ih++) {
	h1 = h0 + ih*dh; 
	h2 = h1+dh; 
	h1 *= h1; 
	h2 *= h2;

	switch (num) {
	    case 1:
		offd1 = h2; 
		offd2 = h1; 
		diag1 = -2.*h2;
		diag2 = -2.*h1;
		break;
	    case 2:
		offd1 = -h1; 
		offd2 = -h2; 
		diag1 = 2.*h1;
		diag2 = 2.*h2;
		break;
	    default:
		offd1 = 1. - c1*h2 + c2*h1; 
		offd2 = 1. - c1*h1 + c2*h2; 
		diag1 = 12. - 2.*offd1;
		diag2 = 12. - 2.*offd2;
		break;
	}
       
	for (ix=1; ix < nx-1; ix++) {
	    i = ix + ih*nx;

	    if (adj) {
		x[i]    += y[i] * conjf (diag2); 
		x[i-1]  += y[i] * conjf (offd2);
		x[i+1]  += y[i] * conjf (offd2);
		x[i-nx]   -= y[i] * conjf (diag1);
		x[i-nx-1] -= y[i] * conjf (offd1);
		x[i-nx+1] -= y[i] * conjf (offd1);
	    } else {
		y [i] += 
		    x[i]    * diag2 + (x[i-1]    + x[i+1]   ) * offd2 -
		    x[i-nx] * diag1 - (x[i-nx-1] + x[i-nx+1]) * offd1;
	    }
	} /* ix */
    } /* ih */
}


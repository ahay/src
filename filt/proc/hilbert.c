/* From "Closed-form design of maximally flat FIR Hilbert transformers,
 * differentiators, and fractional delayers by power series expansion" by
 * Soo-Chang Pei and Peng-Hua Wang, IEEE Trans. on Circuits and Systems - Part
 * I: Fundamental theory and applications, v. 48, No. 4, 2001, 389-398. */
#include <math.h>

#include <rsf.h>

#include "hilbert.h"

static float c, c2, *h;
static int n, nt;

void hilbert_init(int nt1, int n1, float c1)
{
    n = n1;
    nt = nt1;
    c = 1./(2*sqrtf(c1));
    c2 = c*c;
    h = sf_floatalloc(nt);
}

void hilbert_free(void)
{
    free(h);
}

void hilbert (const float* trace, float* trace2)
{
    int i, it;
    
    for (it=0; it < nt; it++) {
	h[it] = trace[it];
    }

    for (i=n; i >= 1; i--) {
	trace2[0] = h[0] + 4*(h[2]-2*h[1]+h[0])*c2;
	trace2[1] = h[1] + 2*(h[3]-h[2]-h[1]+h[0])*c2;
	for (it=2; it < nt-2; it++) {
	    trace2[it] = h[it]+(h[it+2]-2.*h[it]+h[it-2])*c2;
	}
	trace2[nt-2] = h[nt-2] + 2*(h[nt-1]-h[nt-2]-h[nt-3]+h[nt-4])*c2;
	trace2[nt-1] = h[nt-1] + 4*(h[nt-1]-2*h[nt-2]+h[nt-3])*c2;

	for (it=0; it < nt; it++) {
	    h[it] = trace[it] + trace2[it]*(2*i-1)/(2*i);
	}
    }

    trace2[0] = 2.*(h[1]-h[0])*c;
    for (it=1; it < nt-1; it++) {
	trace2[it] = (h[it+1]-h[it-1])*c;
    }
    trace2[nt-1] = 2.*(h[nt-1]-h[nt-2])*c;
}

    
    
	    

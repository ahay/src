#include <stdlib.h>

#include "ccgstep.h"
#include "alloc.h"

/* precision */
static const float EPSILON=1.e-12;

static float complex* S;  /* model step */
static float complex* Ss; /* residual step */
static bool Allocated = false; /* if S and Ss are allocated */

static double complex dotprod (int n, 
			       const float complex* x, const float complex* y);
static float norm(int n, const float complex* x);

/* ccgstep
   --------
   A step of Claerbout's conjugate-gradient iteration for complex operators.
   nx - model size
   ny - data size
   x[nx] - current model
   g[nx] - gradient
   rr[ny] - residula (d - A x)
   gg[ny] - conjugate gradient */
void sf_ccgstep( bool forget, int nx, int ny, 
		 float complex* x,  const float complex* g, 
		 float complex* rr, const float complex* gg) {
    double complex sds, gdg, gds, sdg, determ, gdr, sdr, alfa, beta;
    int i;
    if (Allocated == false) {
	Allocated = forget = true;
	S  = sf_complexalloc (nx);
	Ss = sf_complexalloc (ny);
    }
    if ( forget == true) {
	for (i = 0; i < nx; i++) {
	    S[i] = 0.;
	}
	for (i = 0; i < ny; i++) {
	    Ss[i] = 0.;
	}    
	beta = 0.0;
	if (norm(ny,gg) <= 0.) return;
	alfa = - dotprod( ny, gg, rr) / dotprod(ny, gg, gg);
    } else {
	/* search plane by solving 2-by-2
	   G . (R - G*alfa - S*beta) = 0
	   S . (R - G*alfa - S*beta) = 0 */
	gdg = dotprod( ny, gg, gg);       
	sds = dotprod( ny, Ss, Ss);       
	gds = dotprod( ny, gg, Ss);    
	sdg = dotprod( ny, Ss, gg);   
	if (gdg == 0. || sds == 0.) return;

	determ = 1.0 - (gds/gdg)*(gds/sds);
	if (creal(determ) > EPSILON) {
	    determ *= gdg * sds;
	} else {
	    determ = gdg * sds * EPSILON;
	}

	gdr = - dotprod( ny, gg, rr);
	sdr = - dotprod( ny, Ss, rr);
	alfa = ( sds * gdr - gds * sdr ) / determ;
	beta = (-sdg * gdr + gdg * sdr ) / determ;
    }
    for (i = 0; i < nx; i++) {
	S[i]  =  alfa * g[i] + beta *  S[i];
	x[i] +=  S[i];
    }
    for (i = 0; i < ny; i++) {
	Ss[i] = alfa * gg[i] + beta * Ss[i];
	rr[i] += Ss[i];
    }
}

void sf_ccgstep_close (void) {
    if (Allocated == true) {
	free (S);
	free (Ss);
	Allocated = false;
    }
}

static double complex dotprod (int n, 
			       const float complex* x, const float complex* y) 
{
    double complex prod;
    float complex xi, yi;
    int i;
    prod = 0.;
    for (i = 0; i < n; i++) {
	xi = x[i];
	yi = y[i];
	prod += creal(xi)*creal(yi) + cimag(xi)*cimag(yi) + 
	    I* (creal(xi)*cimag(yi) - cimag(xi)*creal(yi));
    }
    return prod;
}

static float norm (int n, const float complex* x) {
    int i;
    float xn;

    xn = 0.0;
    for (i=0; i < n; i++) {
	xn += creal(x[i]*conjf(x[i]));
    }
    return xn;
}

/* 	$Id: ccgstep.c,v 1.1 2004/03/13 06:11:03 fomels Exp $	 */


#include <stdlib.h>

#include "cgstep.h"
#include "alloc.h"

/* precision */
static const float EPSILON=1.e-12;

static float* S;  /* model step */
static float* Ss; /* residual step */
static bool Allocated = false; /* if S and Ss are allocated */

static double dotprod (int n, const float* x, const float* y);

/* cgstep
   --------
   A step of Claerbout's conjugate-gradient interation.
   nx - model size
   ny - data size
   x[nx] - current model
   g[nx] - gradient
   rr[ny] - residula (d - A x)
   gg[ny] - conjugate gradient */
void sf_cgstep( bool forget, int nx, int ny, 
		float* x, const float* g, float* rr, const float* gg) {
    double sds, gdg, gds, determ, gdr, sdr, alfa, beta;
    int i;
    if (Allocated == false) {
	Allocated = forget = true;
	S  = sf_floatalloc (nx);
	Ss = sf_floatalloc (ny);
    }
    if ( forget == true) {
	for (i = 0; i < nx; i++) {
	    S[i] = 0.;
	}
	for (i = 0; i < ny; i++) {
	    Ss[i] = 0.;
	}    
	beta = 0.0;
	alfa = dotprod( ny, gg, gg);
	if (alfa <= 0.) return;
	alfa = - dotprod( ny, gg, rr) / alfa;
    } else {
	/* search plane by solving 2-by-2
	   G . (R - G*alfa - S*beta) = 0
	   S . (R - G*alfa - S*beta) = 0 */
	gdg = dotprod( ny, gg, gg);       
	sds = dotprod( ny, Ss, Ss);       
	gds = dotprod( ny, gg, Ss);       
	if (gdg == 0. || sds == 0.) return;

	determ = 1.0 - (gds/gdg)*(gds/sds);
	if (determ > EPSILON) {
	    determ *= gdg * sds;
	} else {
	    determ = gdg * sds * EPSILON;
	}

	gdr = - dotprod( ny, gg, rr);
	sdr = - dotprod( ny, Ss, rr);
	alfa = ( sds * gdr - gds * sdr ) / determ;
	beta = (-gds * gdr + gdg * sdr ) / determ;
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

void sf_cgstep_close (void) {
    if (Allocated == true) {
	free (S);
	free (Ss);
	Allocated = false;
    }
}

static double dotprod (int n, const float* x, const float* y) {
    double prod;
    int i;
    prod = 0.;
    for (i = 0; i < n; i++) {
	prod += ((double) x[i])*y[i];
    }
    return prod;
}

/* 	$Id: cgstep.c,v 1.1 2003/10/21 15:12:39 fomels Exp $	 */


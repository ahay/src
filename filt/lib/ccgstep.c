/* Claerbout's conjugate-gradient iteration for complex numbers. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdlib.h>

#include "ccgstep.h"
#include "alloc.h"

#ifndef __cplusplus
/*^*/

#include "_bool.h"
#include "c99.h"
/*^*/

static const float EPSILON=1.e-12; /* precision */

static float complex* S;  /* model step */
static float complex* Ss; /* residual step */
static bool Allocated = false; /* if S and Ss are allocated */

static double complex dotprod (int n, 
			       const float complex* x, const float complex* y);
/* complex dot product */
static float norm(int n, const float complex* x);
/* complex norm */

void sf_ccgstep( bool forget             /* restart flag */, 
		 int nx                  /* model size */, 
		 int ny                  /* data size */, 
		 float complex* x        /* current model [nx] */,  
		 const float complex* g  /* gradient [nx] */, 
		 float complex* rr       /* data residual [ny] */,
		 const float complex* gg /* conjugate gradient [ny] */) 
/*< Step of Claerbout's conjugate-gradient iteration for complex operators. 
  The data residual is rr = A x - dat
>*/
{
    double complex sds, gdg, gds, sdg, determ, gdr, sdr, alfa, beta;
    int i;
    if (Allocated == false) {
	Allocated = forget = true;
	S  = sf_complexalloc (nx);
	Ss = sf_complexalloc (ny);
    }
    if (forget) {
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

void sf_ccgstep_close (void) 
/*< Free allocated space. >*/ 
{
    if (Allocated) {
	free (S);
	free (Ss);
	Allocated = false;
    }
}

static double complex dotprod (int n, 
			       const float complex* x, const float complex* y)
/* complex dot product */
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

static float norm (int n, const float complex* x) 
/* complex norm */
{
    int i;
    float xn;

    xn = 0.0;
    for (i=0; i < n; i++) {
	xn += creal(x[i]*conjf(x[i]));
    }
    return xn;
}

#endif
/*^*/

/* 	$Id$	 */

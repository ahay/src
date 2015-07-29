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
#include "komplex.h"
#include "blas.h"

#include "_bool.h"
#include "c99.h"
/*^*/

static const double EPSILON=1.e-12; /* precision */

static sf_complex* S;  /* model step */
static sf_complex* Ss; /* residual step */
static bool Allocated = false; /* if S and Ss are allocated */

static sf_double_complex dotprod (int n, 
				  const sf_complex* x, const sf_complex* y);
/* complex dot product */

void sf_ccgstep( bool forget             /* restart flag */, 
		 int nx                  /* model size */, 
		 int ny                  /* data size */, 
		 sf_complex* x        /* current model [nx] */,  
		 const sf_complex* g  /* gradient [nx] */, 
		 sf_complex* rr       /* data residual [ny] */,
		 const sf_complex* gg /* conjugate gradient [ny] */) 
/*< Step of Claerbout's conjugate-gradient iteration for complex operators. 
  The data residual is rr = A x - dat
>*/
{
    sf_double_complex sds, gdg, gds, sdg, determ, gdr, sdr, alfa, beta;
    int i;
    if (Allocated == false) {
	Allocated = forget = true;
	S  = sf_complexalloc (nx);
	Ss = sf_complexalloc (ny);
    }
    if (forget) {
	for (i = 0; i < nx; i++) {
	    S[i] = sf_cmplx(0.0,0.0);
	}
	for (i = 0; i < ny; i++) {
	    Ss[i] = sf_cmplx(0.0,0.0);
	}    

	beta = sf_dcmplx(0.0,0.0);
	if (cblas_scnrm2(ny,gg,1) <= 0.) return;
#ifdef SF_HAS_COMPLEX_H
	alfa = - dotprod( ny, gg, rr) / dotprod(ny, gg, gg);
#else
	alfa = sf_dcneg(sf_dcdiv(dotprod( ny, gg, rr),dotprod(ny, gg, gg)));
#endif
    } else {
	/* search plane by solving 2-by-2
	   G . (R - G*alfa - S*beta) = 0
	   S . (R - G*alfa - S*beta) = 0 */
	gdg = dotprod( ny, gg, gg);       
	sds = dotprod( ny, Ss, Ss);       
	gds = dotprod( ny, gg, Ss);    
	sdg = dotprod( ny, Ss, gg);   
	if (cabs(gdg) == 0. || cabs(sds) == 0.) return;

#ifdef SF_HAS_COMPLEX_H
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
#else
	determ = sf_dcneg(sf_dcmul(sf_dcdiv(gds,gdg),
				   sf_dcdiv(gds,sds)));
	determ.r += 1.0;

	if (creal(determ) > EPSILON) {
	    determ = sf_dcmul(determ,sf_dcmul(gdg,sds));
	} else {
	    determ = sf_dcmul(gdg,sf_dcrmul(sds,EPSILON));
	}

	gdr = sf_dcneg(dotprod( ny, gg, rr));
	sdr = sf_dcneg(dotprod( ny, Ss, rr));
	alfa = sf_dcdiv(sf_dcsub(sf_dcmul(sds,gdr),
				 sf_dcmul(gds,sdr)),determ);
	beta = sf_dcdiv(sf_dcsub(sf_dcmul(gdg,sdr),
				 sf_dcmul(sdg,gdr)),determ);
#endif
    }

    for (i = 0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	S[i]  =  alfa * g[i] + beta *  S[i];
	x[i] +=  S[i];
#else
	S[i]  = sf_cadd(sf_dccmul(alfa,g[i]),
			sf_dccmul(beta,S[i]));
	x[i] = sf_cadd(x[i],S[i]);
#endif
    }
    for (i = 0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
	Ss[i] = alfa * gg[i] + beta * Ss[i];
	rr[i] += Ss[i];
#else
	Ss[i] = sf_cadd(sf_dccmul(alfa,gg[i]),
			sf_dccmul(beta,Ss[i]));
	rr[i] = sf_cadd(rr[i],Ss[i]);
#endif
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

static sf_double_complex dotprod (int n, 
				  const sf_complex* x, const sf_complex* y)
/* complex dot product */
{
    sf_double_complex prod, pi;
    sf_complex xi, yi;
    int i;

    prod = sf_dcmplx(0.,0.);
    for (i = 0; i < n; i++) {
	xi = x[i];
	yi = y[i];
	pi = sf_dcmplx(
	    (double) crealf(xi)*crealf(yi) + 
	    (double) cimagf(xi)*cimagf(yi), 
	    (double) crealf(xi)*cimagf(yi) - 
	    (double) cimagf(xi)*crealf(yi));
#ifdef SF_HAS_COMPLEX_H
	prod += pi;
#else
	prod = sf_dcadd(prod,pi);
#endif
    }
    return prod;
}

/* 	$Id: ccgstep.c 7107 2011-04-10 02:04:14Z ivlad $	 */

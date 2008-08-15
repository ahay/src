/* Bessel function */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
/* 1983 version.  w. fullerton, c3, los alamos scientific lab */

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "bessel.h"

static float xmin, xsml, xmax;

static float cheb_eval(float x, const double *cs, int n)
/* evaluate n-term chebyshev series cs at x 
   r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973) */
{
    int i, ni;
    float c, b0, b1, b2, twox;

    b0 = b1 = b2 = 0.;
    twox = 2.*x;
    
    for (i=0; i < n; i++) {
        b2 = b1;
        b1 = b0;
        ni = n - 1 - i;
        b0 = twox*b1 - b2 + cs[ni];
    }

    c = 0.5 * (b0-b2);
    return c;
}

static int bessel_init(const double *cs, int n)
{
    double err;
    int i, ii;

    err = 0.;
    for (ii=0; ii < n; ii++) {
        i = n - 1 - ii;
	err += fabs(cs[i]);
	if (err > 0.1*FLT_EPSILON) 
	    return (i+1);
    }

    sf_error("%s: precision issue",__FILE__);
    return 0;
}

static float bessel_i1e(float x)
{
    const int nb=11, na1=21, na2=22;
    const double bc[] = {
	-.001971713261099859,
	.407348876675464810,
	.034838994299959456,
	.001545394556300123,
	.000041888521098377,
	.000000764902676483,
	.000000010042493924,
	.000000000099322077,
	.000000000000766380,
	.000000000000004741,
	.000000000000000024
    };
    const double a1c[] = {
        -.02846744181881479,
	-.01922953231443221,
	-.00061151858579437,
	-.00002069971253350,
	+.00000858561914581,
	+.00000104949824671,
	-.00000029183389184,
        -.00000001559378146,
	+.00000001318012367,
	-.00000000144842341,
        -.00000000029085122,
	+.00000000012663889,
	-.00000000001664947,
        -.00000000000166665,
	+.00000000000124260,
	-.00000000000027315,
	+.00000000000002023,
	+.00000000000000730,
	-.00000000000000333,
	+.00000000000000071,
	-.00000000000000006
    };
    const double a2c[] = {
	+.02857623501828014,
	-.00976109749136147,
	-.00011058893876263,
	-.00000388256480887,
	-.00000025122362377,
	-.00000002631468847,
	-.00000000383538039,
	-.00000000055897433,
	-.00000000001897495,
	+.00000000003252602,
	+.00000000001412580,
	+.00000000000203564,
	-.00000000000071985,
	-.00000000000040836,
	-.00000000000002101,
	+.00000000000004273,
	+.00000000000001041,
	-.00000000000000382,
	-.00000000000000186,
	+.00000000000000033,
	+.00000000000000028,
	-.00000000000000003
    };
    static int ntb=0, nta1, nta2;
    float y, b;

    if (0.0==x) return 0.0f;

    if (0==ntb) {
	ntb = bessel_init(bc,nb);
	nta1 = bessel_init(a1c,na1);
	nta2 = bessel_init(a2c,na2);
    }

    y = fabsf(x);

    if (y <= 3.0) {
	if (y < xmin) sf_error("%s: underflow",__FILE__);
	if (y > xsml) {
	    b = x * (.875 + cheb_eval(y*y/4.5-1., bc,ntb));
	} else {
	    b = 0.5*x;
	}
	b *= expf(-y);
    } else {
	if (y <= 8.0) {
	    b = (.375 + cheb_eval ((48./y-11.)/5., a1c, nta1)) / sqrtf(y);
	} else {
	    b = (.375 + cheb_eval (16./y-1.0,      a2c, nta2)) / sqrtf(y);
	}
	b *= SF_SIG(x);
    }

    return b;
}


float bessel_i1 (float x)
/*< Bessel I_1 function >*/
{
    const int ns=11;
    const double cs[] = {
	-.001971713261099859,
	.407348876675464810,
	.034838994299959456,
	.001545394556300123,
	.000041888521098377,
	.000000764902676483,
	.000000010042493924,
	.000000000099322077,
	.000000000000766380,
	.000000000000004741,
	.000000000000000024
    };
    static int nt=0;
    float y, b;

    if (0.0==x) return 0.0f;
    
    if (0==nt) {
	nt = bessel_init(cs,ns);
	xmin = 2.0*FLT_MIN;
	xsml = sqrtf(8.0*FLT_EPSILON);
	xmax = logf (FLT_MAX);
    }

    y = fabsf(x);

    if (y < xmin) sf_error("%s: undeflow",__FILE__);
    if (y > xmax) sf_error("%s: overflow",__FILE__);

    if (y <= 3.0) {
	if (y > xsml) {
	    b = x * (0.875 + cheb_eval(y*y/4.5-1.,cs,nt));
	} else { 
	    b = 0.5*x;
	}
    } else {
	b = exp(y) * bessel_i1e(x);
    }

    return b;
}     





/* Traveltime interpolation interface */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include <rsf.h>
/*^*/

#include "tinterp.h"

static int nn;
static float dx;
static char* type;

float hermite_exp(float t, float dx, float p0, float m0, float p1, float m1);
float hermite_fac(float t, float dx, float p0, float m0, float p1, float m1);
float hermite_ber(float t, float dx, float p0, float m0, float p1, float m1);

void tinterp_init(int nt     /* model dimension */,
		  float ds   /* original source sampling */,
		  char* what /* type of basis function */)
/*< initialization >*/
{
    nn = nt;
    dx = ds;
    type = what;
}

void tinterp_linear(float* p /* interpolated result */,
		    float x  /* position */,
		    float* p0, float* p1 /* values at both ends */)
/*< linear interpolation >*/
{
    int i;
    float t;

    t = x/dx;

    for (i=0; i < nn; i++) {
	p[i] = (1.-t)*p0[i]+t*p1[i];
    }
}

void tinterp_hermite(float* p /* interpolated result */,
		     float x  /* position */,
		     float* p0, float* p1 /* values at both ends */,
		     float* m0, float* m1 /* tangent at both ends */)
/*< cubic Hermite spline interpolation >*/
{
    int i;
    float t;

    t = x/dx;

    switch (type[0]) {
	case 'e': /* expanded */
	    for (i=0; i < nn; i++) {
		p[i] = hermite_exp(t,dx,p0[i],m0[i],p1[i],m1[i]);
	    }
	    break;

	case 'f': /* factorized */
	    for (i=0; i < nn; i++) {
		p[i] = hermite_fac(t,dx,p0[i],m0[i],p1[i],m1[i]);
	    }
	    break;

	case 'b': /* Bernstein */
	    for (i=0; i < nn; i++) {
		p[i] = hermite_ber(t,dx,p0[i],m0[i],p1[i],m1[i]);
	    }
	    break;

	default:
	    sf_error("Basis function not supported");
	    break;
    }
}

void dinterp_hermite(float* p /* interpolated result */,
		     float x  /* position */,
		     float* p0, float* p1 /* values at both ends */,
		     float* m0, float* m1 /* tangent at both ends */)
/*< cubic Hermite spline interpolation (derivative) >*/
{
    int i;
    float t;
    double h00, h10, h01, h11;

    t = x/dx;

    /* all basis functions provide the same coefficients */
    h00 = -6.*t*(1.-t);
    h10 = (1.-t)*(1.-3.*t);
    h01 = 6.*t*(1.-t);
    h11 = 3.*t*t-2.*t;

    for (i=0; i < nn; i++) {
	p[i] = h00*p0[i]+h10*dx*m0[i]+h01*p1[i]+h11*dx*m1[i];
    }
}

void tinterp_close()
/*< close >*/
{
    return;
}

float hermite_exp(float t, float dx, float p0, float m0, float p1, float m1)
/* cubic Hermite with expanded basis */
{
    double h00, h10, h01, h11;

    h00 = 2.*t*t*t-3.*t*t+1.;
    h10 = t*t*t-2.*t*t+t;
    h01 = -2.*t*t*t+3.*t*t;
    h11 = t*t*t-t*t;

    return h00*p0+h10*dx*m0+h01*p1+h11*dx*m1;
}

float hermite_fac(float t, float dx, float p0, float m0, float p1, float m1)
/* cubic Hermite with factorized basis */
{
    double h00, h10, h01, h11;

    h00 = (1.+2.*t)*(1.-t)*(1.-t);
    h10 = t*(1.-t)*(1.-t);
    h01 = t*t*(3.-2.*t);
    h11 = t*t*(t-1.);

    return h00*p0+h10*dx*m0+h01*p1+h11*dx*m1;
}

float hermite_ber(float t, float dx, float p0, float m0, float p1, float m1)
/* cubic Hermite with Bernstein basis */
{
    double h00, h10, h01, h11;

    h00 = (1.-t)*(1.-t)*(1.-t)+3.*t*(1.-t)*(1.-t);
    h10 = t*(1.-t)*(1.-t);
    h01 = 3.*t*t*(1.-t)+t*t*t;
    h11 = -t*t*(1.-t);

    return h00*p0+h10*dx*m0+h01*p1+h11*dx*m1;
}

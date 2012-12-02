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

#include "tinterpsr.h"

static int nn;
static float sdx, rdx;

void tinterp_init(int nt     /* model dimension */,
		  float sds   /* original source sampling */,
		  float rds   /* original source sampling */)
/*< initialization >*/
{
    nn = nt;
    sdx = sds;
    rdx = rds;
}

void tinterp_linear(bool source /* source or receiver */,
		    float* p /* interpolated result */,
		    float x  /* position */,
		    float* p0, float* p1 /* values at both ends */)
/*< linear interpolation >*/
{
    int i;
    float t, dx;

    dx = source?sdx:rdx;
    t = x/dx;

    for (i=0; i < nn; i++) {
	p[i] = (1.-t)*p0[i]+t*p1[i];
    }
}

void dinterp_linear(bool source /* source or receiver */,
		    float* p /* interpolated result */,
		    float x  /* position */,
		    float* p0, float* p1 /* values at both ends */)
/*< linear interpolation (derivative) >*/
{
    int i;
    float dx;

    dx = source?sdx:rdx;

    for (i=0; i < nn; i++) {
	p[i] = (-p0[i]+p1[i])/dx;
    }
}

void tinterp_partial(bool source /* source or receiver */,
		     float* p /* interpolated result */,
		     float x  /* position */,
		     int n1, int n2 /* grid size */,
		     float dd /* grid horizontal sampling */,
		     float* p0, float* p1 /* values at both ends */)
/*< interpolation with fixed relative coordinate >*/
{
    int i, k1, k2;
    float t, dx;

    dx = source?sdx:rdx;
    t = x/dx;
    
    k1 = x/dd+0.5;
    k2 = (dx-x)/dd+0.5;

    for (i=0; i < nn; i++) {
	if (i-k1*n1>=0 && i+k2*n1<nn) {
	    p[i] = (1.-t)*p0[i-k1*n1]+t*p1[i+k2*n1];
	} else {
	    p[i] = (1.-t)*p0[i]+t*p1[i];
	}
    }
}

void dinterp_partial(bool source /* source or receiver */,
		     float* p /* interpolated result */,
		     float x  /* position */,
		     int n1, int n2 /* grid size */,
		     float dd /* grid horizontal sampling */,
		     float* p0, float* p1 /* values at both ends */)
/*< interpolation with fixed relative coordinate (derivative) >*/
{
    int i, k1, k2;
    float t, dx;

    dx = source?sdx:rdx;
    t = x/dx;
    
    k1 = x/dd+0.5;
    k2 = (dx-x)/dd+0.5;

    for (i=0; i < nn; i++) {
	if (i-k1*n1>=0 && i+k2*n1<nn) {
	    p[i] = (-p0[i-k1*n1]+p1[i+k2*n1])/dx;

	    if (i-(k1+1)*n1>=0 && i-(k1-1)*n1<nn)
		p[i] -= (1.-t)*(-p0[i-(k1+1)*n1]+p0[i-(k1-1)*n1])/(2.*dd);
	    else if (i-(k1+1)*n1<0)
		p[i] -= (1.-t)*(-p0[i-k1*n1]+p0[i-(k1-1)*n1])/dd;
	    else
		p[i] -= (1.-t)*(-p0[i-(k1+1)*n1]+p0[i-k1*n1])/dd;
		
	    if (i-(k2+1)*n1>=0 && i-(k2-1)*n1<nn)
		p[i] -= t*(-p1[i-(k2+1)*n1]+p1[i-(k2-1)*n1])/(2.*dd);
	    else if (i-(k2+1)*n1<0)
		p[i] -= t*(-p1[i-k2*n1]+p1[i-(k2-1)*n1])/dd;
	    else
		p[i] -= t*(-p1[i-(k2+1)*n1]+p1[i-k2*n1])/dd;
	} else {
	    p[i] = (-p0[i]+p1[i])/dx;
	}
    }
}

void tinterp_hermite(bool source /* source or receiver */,
		     float* p /* interpolated result */,
		     float x  /* position */,
		     float* p0, float* p1 /* values at both ends */,
		     float* m0, float* m1 /* tangent at both ends */)
/*< cubic Hermite spline interpolation >*/
{
    int i;
    float t, dx;
    double h00, h10, h01, h11;

    dx = source?sdx:rdx;
    t = x/dx;

    h00 = 2.*t*t*t-3.*t*t+1.;
    h10 = t*t*t-2.*t*t+t;
    h01 = -2.*t*t*t+3.*t*t;
    h11 = t*t*t-t*t;

    for (i=0; i < nn; i++) {
	p[i] = h00*p0[i]+h10*dx*m0[i]+h01*p1[i]+h11*dx*m1[i];
    }
}

void dinterp_hermite(bool source /* source or receiver */,
		     float* p /* interpolated result */,
		     float x  /* position */,
		     float* p0, float* p1 /* values at both ends */,
		     float* m0, float* m1 /* tangent at both ends */)
/*< cubic Hermite spline interpolation (derivative) >*/
{
    int i;
    float t, dx;
    double h00, h10, h01, h11;

    dx = source?sdx:rdx;
    t = x/dx;

    h00 = -6.*t*(1.-t);
    h10 = (1.-t)*(1.-3.*t);
    h01 = 6.*t*(1.-t);
    h11 = 3.*t*t-2.*t;

    for (i=0; i < nn; i++) {
	p[i] = (h00*p0[i]+h10*dx*m0[i]+h01*p1[i]+h11*dx*m1[i])/dx;
    }
}

void tinterp_close()
/*< close >*/
{
    return;
}

/* Weighted Sinc interpolation */
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

#include <math.h>
#include <float.h>

#include <rsf.h>

#include "interp_sinc.h"

static float a;
static float (*win)(float,int);

static float besselI(float x)
{
    float a, y;
    int k, kmax;
    
    y = x/2; 
    y = y*y;
    a = 1.;
    kmax = (int) x/(2.*sqrtf(FLT_EPSILON));
    for (k=kmax; k > 0; k--) {
    	a = 1. + a*y/(k*k);
    }

    return a;
}

static float cosine_win (float x, int nw)
{
    return x/cosf(SF_PI * x/nw);
}

static float kaiser_win (float x, int nw)
{
    nw = (nw+1)*0.5; 
    return x*besselI(a)/besselI(a * sqrtf(1-x*x/(nw*nw)));
}

static float lanczos_win (float x, int nw)
{
    nw = (nw+1)*0.5; 
    return x*x/(nw*sinf(SF_PI*x/nw));
}

static float welch_win (float x, int nw)
{
    nw = (nw+1)*0.5; 
    return x/(1.-x*x/(nw*nw));
}

void sinc_init(char wintype /* window type */, 
	       float a_in   /* factor for Kaiser window */)
/*< initialize >*/
{
    switch (wintype) {
	case 'c': 
	    win = cosine_win;
	    break;
	case 'k':
	    win = kaiser_win;
	    a = a_in;
	    break;
	case 'l':
	    win = lanczos_win;
	    break;
	case 'w':
	    win = welch_win;
	    break;
	default:
	    sf_error("%s: window type %c is not implemented",
		     __FILE__,wintype);
	    break;
    }    
}

void sinc_int (float x, int nw, float *w)
/*< interpolation function >*/
{
    int i, nc;
    float y, sx;

    sx = sinf (SF_PI*x)/SF_PI; 
    nc = (nw+1)*0.5; 

    for (i=0; i < nw; i++ ) {
	y = x + nc - i - 1.;
	w[i] = win(y,nw);
    }
    for (i=nc%2; i < nw; i += 2) {
	w[i] = - w[i];
    }
    for (i=0; i < nw; i++ ) {
	if (fabsf(w[i]) > 0.) {
	    w[i] = sx/w[i];
	} else {
	    w[i] = 1.;
	}
    }
}


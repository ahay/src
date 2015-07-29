/* Cubic convolution interpolation. */
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

#include <rsf.h>

#include "interp_cube.h"
  
static void cube2_int (float x, float *w);
static void cube4_int (float x, float *w);
static void cube6_int (float x, float *w);
static void cube8_int (float x, float *w);

void cube_int (float x, int n, float *w)
/*< interpolation function >*/
{
    switch (n) {
	case 8:
	    cube8_int (x, w);
	    break;
	case 6:
	    cube6_int (x, w);
	    break;
	case 4:
	    cube4_int (x, w);
	    break;
	case 2:
	    cube2_int (x, w);
	    break;
	default:
	    sf_error("%s: size %d is not implemented",__FILE__,n);
	    break;
    }
}

static void cube2_int (float x, float *w)
{
    float x2;

    x2 = x*x;
    w[0] = 1. + x2*(2.*x-3.);
    w[1] = (3. - 2.*x)*x2;
}

static void cube4_int (float x, float *w) 
{
    float x2;

    x2 = x*x;
    w[0] = 0.5*x*((2. - x)*x-1.);
    w[1] = 1. + 0.5*x2*(3.*x-5.);
    w[2] = 0.5*x*(1. + (4. - 3.*x)*x);
    w[3] = 0.5*(x-1.)*x2;
}

static void cube6_int (float x, float *w) 
{
    float x2;

    x2 = x*x;
    w[0] = x*(1. + (x-2.)*x)*3./32.;
    w[1] =  x*((41. - 19.*x)*x-22.)/32.;
    w[2] = 1. + x2*(21.*x-37.)/16.;
    w[3] = x*(11. + (26. - 21.*x)*x)/16.;
    w[4] = x*((19.*x-16.)*x-3.)/32.;
    w[5] = (1. - x)*x2*3./32.;
}

static void cube8_int (float x, float *w)
{
    float x2;

    x2 = x*x;
    w[0] = x*((2. - x)*x-1.)*157./8064.;
    w[1] = x*(1384. + x*(1227.*x-2611.))/8064.;
    w[2] = x*((11274. - 4945.*x)*x-6329.)/8064.;
    w[3] = 1. + x2*(9799.*x-17863.)/8064.;
    w[4] = x*(6329. + (11534. - 9799.*x)*x)/8064.;
    w[5] = x*(x*(4945.*x-3561.)-1384.)/8064.;
    w[6] = x*(157. + (1070. - 1227.*x)*x)/8064.;
    w[7] = x2*(x-1.)*157./8064.;
}

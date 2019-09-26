/* B-spline interpolation */
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

#include "spline.h"

#include "error.h"

#include "banded.h"
#include "tridiagonal.h"
/*^*/

static const float s3 = 0.75, s4 = 2./3., s6 = 11./20., s8 = 151./315.;
static const float m4 = 0.625, m6 = 77./144., m8 = 151./320.;  
static const float flt3[] = {0.125}; 
static const float flt4[] = {1./6.};
static const float mom4[] = {0.1875};
static const float mom6[] = {2./9., 1./96.};
static const float flt6[] = {26./120., 1./120.};
static const float flt8[] = {1191./5040., 120./5040., 1./5040.};
static const float mom8[] = {2741./11520., 298./11520., 3./11520.};

sf_bands sf_spline_init (int nw /* interpolator length */, 
			 int nd /* data length */)
/*< initialize a banded matrix >*/
{
    sf_bands slv;
    int   na;
    
    na = (nw>0)?(nw-1)/2:(-nw-1)/2;
    slv = sf_banded_init(nd,na);
    
    switch (nw) {
	case -8:
	    sf_banded_const_define (slv, m8, mom8);
	    break;
	case -6:
	    sf_banded_const_define (slv, m6, mom6);
	    break;
	case -4:
	    sf_banded_const_define (slv, m4, mom4);
	    break;
	case 3:
	    sf_banded_const_define (slv, s3, flt3);
	    break;
	case 4:
	    sf_banded_const_define (slv, s4, flt4);
	    break;
	case 6:
	    sf_banded_const_define (slv, s6, flt6);
	    break;
	case 8:
	    sf_banded_const_define (slv, s8, flt8);
	    break;
	default:
	    sf_error("%s: unsupported spline length  %d",__FILE__,nw);
	    break;
    }
    
    return slv;
}

sf_tris sf_spline4_init (int nd /* data length */)
/*< initialize a tridiagonal matrix for cubic splines >*/
{
    sf_tris slv;
    
    slv = sf_tridiagonal_init(nd);
    sf_tridiagonal_const_define (slv, s4, flt4[0], false);
    
    return slv;
}

void sf_spline4_post (int n            /* total trace length */, 
		      int n1           /* start point */, 
		      int n2           /* end point */, 
		      const float* inp /* spline coefficients */, 
		      float* out       /* function values */)
/*< cubic spline post-filtering >*/
{
    int i;
    float o4;

    o4 = flt4[0];

    for (i=n1; i < n2; i++) {
	if (0==i) {
	    out[i-n1] = s4*inp[i] + o4*inp[i+1];
	} else if (n-1==i) {
	    out[i-n1] = s4*inp[i] + o4*inp[i-1];
	} else {
	    out[i-n1] = s4*inp[i] + o4*(inp[i-1]+inp[i+1]);
	}
    }
}


void sf_spline_post (int nw, int o, int d, int n, 
		     const float *modl, float *datr)
/*< post-filtering to convert spline coefficients to model >*/
{
    const float *flt;
    float t, a0;
    int i, j, iw;

    switch (nw) {
	case -8:
	    a0 = m8;
	    flt = mom8;
	    break;
	case -6:
	    a0 = m6;
	    flt = mom6;
	    break;
	case -4:
	    a0 = m4;
	    flt = mom4;
	    break;
	case 3:
	    a0 = s3;
	    flt = flt3;
	    break;
	case 4:
	    a0 = s4;
	    flt = flt4;
	    break;
	case 6:
	    a0 = s6;
	    flt = flt6;
	    break;
	case 8:
	    a0 = s8;
	    flt = flt8;
	    break;
	default:
	    a0 = 0.;
	    flt = NULL;
	    sf_error("%s: unsupported spline length  %d",__FILE__,nw);
	    break;
    }
    nw = (nw>0)?(nw-1)/2:(-nw-1)/2;

    for (i=0; i < n; i++) {
	j = o + i*d;
	t = a0* modl[j];
	for (iw=0; iw < nw; iw++) {
	    if (i+iw+1 < n)  t += modl[j+(iw+1)*d]*flt[iw];
	    if (i-iw-1 >= 0) t += modl[j-(iw+1)*d]*flt[iw];
	}
	datr[j] = t;
    }
}

void sf_spline2 (sf_bands slv1, sf_bands slv2, 
		 int n1, int n2, float** dat, float* tmp)
/*< 2-D spline pre-filtering >*/
{
    int i1, i2;
    for (i2 = 0; i2 < n2; i2++) {
	sf_banded_solve (slv1, dat[i2]);
    }
    for (i1 = 0; i1 < n1; i1++) {
	for (i2 = 0; i2 < n2; i2++) {
	    tmp[i2] = dat[i2][i1];
	}
	sf_banded_solve (slv2, tmp);
	for (i2 = 0; i2 < n2; i2++) {
	    dat[i2][i1] = tmp[i2];
	}
    }
}

void sf_spline3 (sf_bands slv1, sf_bands slv2, sf_bands slv3, 
		 int n1, int n2, int n3, float*** dat, float* tmp2, float* tmp3)
/*< 3-D spline pre-filtering >*/
{
    int i1, i2, i3;

    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
	    sf_banded_solve (slv1, dat[i3][i2]);
	}
    }

    for (i3 = 0; i3 < n3; i3++) {
	for (i1 = 0; i1 < n1; i1++) {
	    for (i2 = 0; i2 < n2; i2++) {
		tmp2[i2] = dat[i3][i2][i1];
	    }
	    sf_banded_solve (slv2, tmp2);
	    for (i2 = 0; i2 < n2; i2++) {
		dat[i3][i2][i1] = tmp2[i2];
	    }
	}
    }

    for (i2 = 0; i2 < n2; i2++) {
	for (i1 = 0; i1 < n1; i1++) {
	    for (i3 = 0; i3 < n3; i3++) {
		tmp3[i3] = dat[i3][i2][i1];
	    }
	    sf_banded_solve (slv3, tmp3);
	    for (i3 = 0; i3 < n3; i3++) {
		dat[i3][i2][i1] = tmp3[i3];
	    }
	}
    }
}

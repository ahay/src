/* Offset continuation roughening operator. */
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
/*^*/

#include "offruffp.h"  

static int nh, nx, num;
static float h0, dh;
static sf_complex c1, c2;

void offruffp_init (float h0_in, int nh_in, float dh_in /* half-offset axis */,
		    int nx_in, float dx                 /* midpoint axis */, 
		    float w                             /* frequency */, 
		    int num_in                          /* continuation */)
/*< Initialize >*/
{
    float w2;

    h0 = h0_in/dx; 
    nh = nh_in; 
    dh = dh_in/dx; 
    nx = nx_in; 
    w2  = w*w; 
    num = num_in;
#ifdef SF_HAS_COMPLEX_H
    c1 = 3.*sf_cmplx(9. + w2,4.*w)/(w2*sf_cmplx(3.,-w));
    c2 = 3.*sf_cmplx(w2 - 27.,8.*w)/(w2*sf_cmplx(3.,-w));
#else
    c1 = sf_cdiv(sf_cmplx(3.*(9. + w2), 3.*4.*w),sf_cmplx(w2*3.,-w2*w));
    c2 = sf_cdiv(sf_cmplx(3.*(w2 - 27.),3.*8.*w),sf_cmplx(w2*3.,-w2*w));
#endif
}

sf_complex offruffp_c1 (void)
/*< return first coefficient >*/
{
    return c1;
}

sf_complex offruffp_c2 (void)
/*< return second coefficient >*/
{
    return c2;
}

void hderp_lop (bool adj, bool add, int n1, int n2, 
		sf_complex *x, sf_complex *y) 
/*< simple derivative roughening >*/
{
    int ih, ix, i;

    sf_cadjnull (adj,add,n1,n2,x,y);

    for (ih=1; ih < nh; ih++) {
	for (ix=0; ix < nx; ix++) {
	    i = ix + ih*nx;
#ifdef SF_HAS_COMPLEX_H
	    if (adj) {
		x[i]    += y[i];
		x[i-nx] -= y[i];
	    } else {
		y[i] += x[i] - x[i-nx];
	    }
#else
	    if (adj) {
		x[i]    = sf_cadd(x[i],   y[i]);
		x[i-nx] = sf_csub(x[i-nx],y[i]);
	    } else {
		y[i] = sf_cadd(y[i],sf_csub(x[i],x[i-nx]));
	    }
#endif
	}
    }
}

void offruffp_lop (bool adj, bool add, int n1, int n2, 
		   sf_complex *x, sf_complex *y) 
/*< offset continuation roughening >*/
{
    int ih, ix, i;
    float h1, h2;
    sf_complex diag1, diag2, offd1, offd2;

    sf_cadjnull (adj,add,n1,n2,x,y);
    
    for (ih=1; ih < nh; ih++) {
	h1 = h0 + ih*dh; 
	h2 = h1+dh; 
	h1 *= h1; 
	h2 *= h2;

	switch (num) {
	    case 1:
		offd1 = sf_cmplx(h2,0.); 
		offd2 = sf_cmplx(h1,0.); 
		diag1 = sf_cmplx(-2.*h2,0.);
		diag2 = sf_cmplx(-2.*h1,0.);
		break;
	    case 2:
		offd1 = sf_cmplx(-h1,0.); 
		offd2 = sf_cmplx(-h2,0.); 
		diag1 = sf_cmplx(2.*h1,0.);
		diag2 = sf_cmplx(2.*h2,0.);
		break;
	    default:
#ifdef SF_HAS_COMPLEX_H
		offd1 = 1. - c1*h2 + c2*h1; 
		offd2 = 1. - c1*h1 + c2*h2; 
		diag1 = 12. - 2.*offd1;
		diag2 = 12. - 2.*offd2;
#else
		offd1 = sf_cadd(sf_cmplx(1.,0.),
				sf_cadd(sf_crmul(c1,-h2),
					sf_crmul(c2, h1))); 
		offd2 = sf_cadd(sf_cmplx(1.,0.),
				sf_cadd(sf_crmul(c1,-h1),
					sf_crmul(c2, h2))); 
		diag1 = sf_csub(sf_cmplx(12.,0.),
				sf_crmul(offd1,2.));
		diag2 = sf_csub(sf_cmplx(12.,0.),
				sf_crmul(offd2,2.));
#endif
		break;
	}
       
	for (ix=1; ix < nx-1; ix++) {
	    i = ix + ih*nx;

#ifdef SF_HAS_COMPLEX_H
	    if (adj) {
		x[i]    += y[i] * conjf (diag2); 
		x[i-1]  += y[i] * conjf (offd2);
		x[i+1]  += y[i] * conjf (offd2);
		x[i-nx]   -= y[i] * conjf (diag1);
		x[i-nx-1] -= y[i] * conjf (offd1);
		x[i-nx+1] -= y[i] * conjf (offd1);
	    } else {
		y [i] += 
		    x[i]    * diag2 + (x[i-1]    + x[i+1]   ) * offd2 -
		    x[i-nx] * diag1 - (x[i-nx-1] + x[i-nx+1]) * offd1;
	    }
#else
	    if (adj) {
		x[i]    = sf_cadd(x[i],sf_cmul(y[i],conjf (diag2))); 
		x[i-1]  = sf_cadd(x[i-1],sf_cmul(y[i],conjf (offd2)));
		x[i+1]  = sf_cadd(x[i+1],sf_cmul(y[i],conjf (offd2)));
		x[i-nx]   = sf_csub(x[i-nx],sf_cmul(y[i],conjf (diag1)));
		x[i-nx-1] = sf_csub(x[i-nx-1],sf_cmul(y[i],conjf (offd1)));
		x[i-nx+1] = sf_csub(x[i-nx+1],sf_cmul(y[i],conjf (offd1)));
	    } else {
		y [i] = 
		    sf_cadd(y [i],
			    sf_csub(sf_cadd(
					sf_cmul(x[i],diag2),
					sf_cmul(sf_cadd(x[i-1],x[i+1]),
						offd2)),
				    sf_cadd(
					sf_cmul(x[i-nx],diag1),
					sf_cmul(sf_cadd(x[i-nx-1],x[i-nx+1]),
						offd1))));
	    }
#endif
	} /* ix */
    } /* ih */
}


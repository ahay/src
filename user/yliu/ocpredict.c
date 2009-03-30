/* Offset slice prediction by offset continuation */
/*
  Copyright (C) 2009 University of Texas at Austin
    
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

#include "ocpredict.h"
#include "ctridiagonal.h"

void ocpredict_step(bool adj           /* adjoint flag */,
		    bool forw          /* forward or backward */,
		    float dw, int nx,
		    float w            /* log-stretch frequency */,
		    float h, float h2  /* offset position */, 
		    sf_complex *trace  /* common offset slice */)
/*< offset continuation prediction step >*/
{
    int ix, k;
    float w2;
    sf_complex diag, diag2, *in, offd, offd2, c1, c2;
    ctris slv;
 
    in = sf_complexalloc(nx);
    slv = ctridiagonal_init (nx);

    w2 = w*w;

    if (fabsf(w) < dw) {
	return;
/*	for (ix=0; ix < nx; ix++) {
	    trace[ix]=sf_cmplx(0.,0.);
	}
*/   }    
    
#ifdef SF_HAS_COMPLEX_H		
    c1 = 3.*sf_cmplx(9. + w2,4.*w)/(w2*sf_cmplx(3.,- w));
    c2 = 3.*sf_cmplx(w2 - 27.,8.*w)/(w2*sf_cmplx(3.,- w));
#else
    c1 = sf_cdiv(sf_cmplx(3.*(9. + w2),3.*4.*w),
		 sf_cmplx(w2*3.,-w2*w));
    c2 = sf_cdiv(sf_cmplx(3.*(w2 - 27.),3.*8.*w),
		 sf_cmplx(w2*3.,-w2*w));
#endif
    
    for (ix=0; ix < nx; ix++) {
	in[ix] = trace[ix];
    }
    h *= h;
    h2 *= h2; 
    if (forw) {
#ifdef SF_HAS_COMPLEX_H	
	offd  = 1. - c1*h2 + c2*h;
	offd2 = 1. - c1*h  + c2*h2;
	diag  = 12. - 2.*offd;
	diag2 = 12. - 2.*offd2;
#else
	offd  = sf_cadd(sf_cmplx(1.,0.),
			sf_cadd(sf_crmul(c1,-h2),
				sf_crmul(c2,h)));
	offd2 = sf_cadd(sf_cmplx(1.,0.),
			sf_cadd(sf_crmul(c1,-h),
				sf_crmul(c2,h2)));
	diag  = sf_cadd(sf_cmplx(12.,0.),sf_crmul(offd,-2.));
	diag2 = sf_cadd(sf_cmplx(12.,0.),sf_crmul(offd2,-2.));
#endif
    } else {
#ifdef SF_HAS_COMPLEX_H	
	offd  = 1. - c1*h  + c2*h2;
	offd2 = 1. - c1*h2 + c2*h;
	diag  = 12. - 2.*offd;
	diag2 = 12. - 2.*offd2;
#else
	offd  = sf_cadd(sf_cmplx(1.,0.),
			sf_cadd(sf_crmul(c1,-h),
				sf_crmul(c2,h2)));
	offd2 = sf_cadd(sf_cmplx(1.,0.),
			sf_cadd(sf_crmul(c1,-h2),
				sf_crmul(c2,h)));
	diag  = sf_cadd(sf_cmplx(12.,0.),sf_crmul(offd,-2.));
	diag2 = sf_cadd(sf_cmplx(12.,0.),sf_crmul(offd2,-2.));
#endif    
    }

    ctridiagonal_const_define (slv, diag2, offd2);
    
#ifdef SF_HAS_COMPLEX_H	
    trace[0] = diag*in[0] + offd*in[1];
#else
    trace[0] = sf_cadd(sf_cmul(diag,in[0]),
			     sf_cmul(offd,in[1]));
#endif
    for (k = 1; k < nx - 1; k++) {
#ifdef SF_HAS_COMPLEX_H	
	trace[k] = diag*in[k] + offd*(in[k+1]+in[k-1]);
#else
	trace[k] = sf_cadd(sf_cmul(diag,in[k]),
			   sf_cmul(offd,sf_cadd(in[k+1],in[k-1])));
#endif
    }
#ifdef SF_HAS_COMPLEX_H	
    trace[nx-1] = diag*in[nx-1] + offd*in[nx-2];
#else
    trace[nx-1] = sf_cadd(sf_cmul(diag,in[nx-1]),
			  sf_cmul(offd,in[nx-2]));
#endif
    ctridiagonal_solve (slv, trace);

}

/* 	$Id$	 */

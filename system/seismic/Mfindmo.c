/* DMO without stacking by finite-difference offset continuation. */
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
#include <rsf.h>
#include "ctridiagonal.h"

int main(int argc, char* argv[])
{
    ctris slv;
    int nw,nh,nx,iw,ix,ih,n1,n2,ns,ig;
    float w0,dw, h0,dh,dx,w,w2,h,h2;
    sf_complex diag,diag2, c1,c2,offd,offd2;
    sf_complex *in=NULL, *out=NULL, **dat=NULL;
    sf_file cmp=NULL, stk=NULL;

    sf_init (argc,argv);
    cmp = sf_input("in");
    stk = sf_output("out");

    if (SF_COMPLEX != sf_gettype(cmp)) sf_error("Need complex input");
    if (!sf_histint(cmp,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(cmp,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
    dh /= dx;
    h0 /= dx;
    if (!sf_histfloat(cmp,"d3",&dw)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&w0)) sf_error("No o3= in input");

    in = sf_complexalloc(nx);
    out = sf_complexalloc(nx);
    dat = sf_complexalloc2(nx,nh);

    slv = ctridiagonal_init (nx);
    
    for (iw=0; iw < nw; iw++) {
	sf_warning("frequency %d of %d;",iw+1,nw);

	w = 2.*SF_PI*(w0 + iw*dw); 
	w2 = w*w;

	sf_complexread(dat[0],nx*nh,cmp);

	if (fabsf(w) < dw) {
	    for (ix=0; ix < nx; ix++) {
		out[ix] = sf_cmplx(0.,0.);
	    }	 
	    for (ig=0; ig < nh; ig++) {
		sf_complexwrite(out,nx,stk);
	    }
	    continue;
	}

#ifdef SF_HAS_COMPLEX_H
	c1 = 3.*sf_cmplx(9. + w2, 4.*w)/(w2*sf_cmplx(3.,- w));
	c2 = 3.*sf_cmplx(w2 - 27.,8.*w)/(w2*sf_cmplx(3.,- w));
#else
	c1 = sf_cdiv(sf_cmplx(3.*(9. + w2),3.*4.*w),
		     sf_cmplx(w2*3.,-w2*w));
	c2 = sf_cdiv(sf_cmplx(3.*(w2 - 27.),3.*8.*w),
		     sf_cmplx(w2*3.,-w2*w));
#endif

	for (ig=0; ig < nh; ig++) {
	    if ((dh > 0. && h0 >= 0.) || (dh < 0. && h0 <= 0.)) {
		n1 = ig-1; 
		ns = -1;
		n2 = -0.5 - h0/dh;
	    } else {
		n1 = 0;
		ns = 1;
		n2 = ig - h0/dh - 0.5;
	    }
	    
	    for (ix=0; ix < nx; ix++) {
		out[ix] = dat[ig][ix];
	    }

	    for (ih=n1; ih != n2; ih += ns) {
		for (ix=0; ix < nx; ix++) {
		    in[ix] = out[ix];
		}		

		h = h0 + ih*dh; 
		h2 = h+dh*ns; 
		h *= h; 
		h2 *= h2;

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

		ctridiagonal_const_define (slv,diag2,offd2);

#ifdef SF_HAS_COMPLEX_H
		out[0] = diag*in[0] + offd*in[1];
#else
		out[0] = sf_cadd(sf_cmul(diag,in[0]),
				 sf_cmul(offd,in[1]));
#endif
		for (ix=1; ix < nx-1; ix++) {
#ifdef SF_HAS_COMPLEX_H	
		    out[ix] = diag*in[ix] + offd*(in[ix+1]+in[ix-1]);
#else
		    out[ix] = sf_cadd(sf_cmul(diag,in[ix]),
				      sf_cmul(offd,sf_cadd(in[ix+1],in[ix-1])));
#endif
		}
#ifdef SF_HAS_COMPLEX_H	
		out[nx-1] = diag*in[nx-1] + offd*in[nx-2];
#else
		out[nx-1] = sf_cadd(sf_cmul(diag,in[nx-1]),
				    sf_cmul(offd,in[nx-2]));
#endif
		
		ctridiagonal_solve (slv,out);
	    }
	    sf_complexwrite (out,nx,stk);
	} /* ig */

    } /* iw */
    sf_warning(".");
	
    exit(0);
}

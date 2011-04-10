/* Differential shot continuation */
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

#include <float.h>

#include <rsf.h>
/*^*/

#include "shotfill.h"
#include "cbanded.h"

static int nh;
static sf_complex **a, *offd[2], **cc, *s0, b0[3];
static float **c, *diag, s, h1, eps;
static void filt2matrix(const sf_complex *filt,
			float* d, sf_complex *dd);
static void predefine (float w /* log-stretch frequency */);

void shotfill_init (int nh_in, float h0, float dh /* half-offset axis */, 
		    float ds                      /* shot sampling */, 
		    float eps1                    /* regularization */)
/*< initialize >*/
{
    nh = nh_in;

    offd[0] = sf_complexalloc(nh-1);
    offd[1] = sf_complexalloc(nh-2);
    diag = sf_floatalloc(nh);
    /* allocate the pentadiagonal matrix diagonals */

    a = sf_complexalloc2(3,nh);
    /* filter coefficients */

    c = sf_floatalloc2(3,nh);
    cc = sf_complexalloc2(3,nh);

    s0 = sf_complexalloc(nh);

    cbanded_init(nh, 2);
    /* initialize complex banded matrix inversion */

    s = 0.5*ds/dh;
    h1 = h0/dh;
    /* normalize by dh */

    b0[2] = sf_cmplx((s-2.)*(s-1.)/12.,0.);
    b0[1] = sf_cmplx(-(s-2.)*(s+2.)/6.,0.);
    b0[0] = sf_cmplx((s+2.)*(s+1.)/12.,0.);
    /* h-independent part of the filter */

    eps = eps1;
}

void shotfill_close (void)
/*< free allocated storage >*/
{
    free(*a);
    free(a);
    free(offd[0]);
    free(offd[1]);
    free(diag);
    free(*c);
    free(c);
    free(*cc);
    free(cc);

    cbanded_close();
}

static void filt2matrix(const sf_complex *filt,
			float* d, sf_complex *dd)
{
#ifdef SF_HAS_COMPLEX_H
    d[0] = crealf (filt[0] * conjf (filt[0]));
    d[1] = crealf (filt[1] * conjf (filt[1]));
    d[2] = crealf (filt[2] * conjf (filt[2]));
    dd[0] = filt[1] * conjf (filt[0]);
    dd[1] = filt[2] * conjf (filt[1]);
    dd[2] = filt[2] * conjf (filt[0]);
#else
    d[0] = crealf (sf_cmul(filt[0],conjf (filt[0])));
    d[1] = crealf (sf_cmul(filt[1],conjf (filt[1])));
    d[2] = crealf (sf_cmul(filt[2],conjf (filt[2])));
    dd[0] = sf_cmul(filt[1],conjf (filt[0]));
    dd[1] = sf_cmul(filt[2],conjf (filt[1]));
    dd[2] = sf_cmul(filt[2],conjf (filt[0]));
#endif
}

static void predefine (float w /* log-stretch frequency */)
{
    float den, h;
    sf_complex *b, a1, a2, b1, b2;
    int ih;

    for (ih=1; ih < nh-1; ih++) {
	/* loop over offsets, skip first and last */
	h = h1 + ih;
	/* offset normalized by dh */

	den = 12.*h*h + (s-1.)*(s+1.)*w*w;
	
	b = a[ih];

	a1 = sf_cmplx(0.,(s-2.)/w);
	a2 = sf_cmplx(0.,(s+2.)/w);

	b1 = sf_cmplx(2.*h,(s-1.)*w);
	b2 = sf_cmplx(2.*h,(s+1.)*w);

#ifdef SF_HAS_COMPLEX_H
	b[0] = b0[0] - 0.5*h*(a1+b1/den); 
	b[1] = b0[1] + h*(sf_cmplx(0.,s/w)+sf_cmplx(2.*h,s*w)/den);
	b[2] = b0[2] - 0.5*h*(a2+b2/den); 
#else
	b[0] = sf_cadd(b0[0],
		       sf_crmul(sf_cadd(a1,sf_crmul(b1,1.0/den)),-0.5*h)); 
	b[1] = sf_cadd(b0[1],sf_crmul(sf_cadd(sf_cmplx(0.,s/w),
					      sf_crmul(sf_cmplx(2.*h,s*w),
						       1.0/den)),h)); 
	b[2] = sf_cadd(b0[2],
		       sf_crmul(sf_cadd(a2,sf_crmul(b2,1.0/den)),-0.5*h)); 
#endif
	filt2matrix(b,c[ih],cc[ih]);
    }
}

void shotfill_define (float w /* log-stretch frequency */)
/*< Fill the shot continuation matrix >*/
{
    int ih;

    predefine(w);

    diag[0] = c[1][0] + c[1][2] + eps;
    diag[1] = 2.*c[1][1] + c[2][0] + c[2][2] + eps;
    for (ih=2; ih < nh-2; ih++) {
	diag[ih] = c[ih-1][0] + c[ih-1][2] + 2.*c[ih][1]
	    +      c[ih+1][0] + c[ih+1][2] + eps;
    }
    diag[nh-2] = c[nh-3][0] + c[nh-3][2] + 2.*c[nh-2][1] + eps;
    diag[nh-1] = c[nh-2][0] + c[nh-2][2] + eps;

#ifdef SF_HAS_COMPLEX_H
    offd[0][0] = cc[1][0] + cc[1][1];
    offd[1][0] = 2.*cc[1][2];
#else
    offd[0][0] = sf_cadd(cc[1][0],cc[1][1]);
    offd[1][0] = sf_crmul(cc[1][2],2.);
#endif
    for (ih=1; ih < nh-2; ih++) {
#ifdef SF_HAS_COMPLEX_H
	offd[0][ih] = 
	    cc[ih][0] + cc[ih+1][0] +  
	    cc[ih][1] + cc[ih+1][1];
	offd[1][ih] = 2.*cc[ih+1][2];
#else
	offd[0][ih] = 
	    sf_cadd(sf_cadd(cc[ih][0],cc[ih+1][0]),  
		    sf_cadd(cc[ih][1],cc[ih+1][1]));
	offd[1][ih] = sf_crmul(cc[ih+1][2],2.);
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    offd[0][nh-2] = cc[nh-2][0] + cc[nh-2][1];
#else
    offd[0][nh-2] = sf_cadd(cc[nh-2][0],cc[nh-2][1]);
#endif
    
    cbanded_define(diag,offd);
    /* define the pentadiagonal matrix */
}

void shotprop_define (float w /* log-stretch frequency */)
/*< Fill the shot propagation matrix >*/
{
    int ih;

    predefine(w);

    diag[0] = c[1][0] + eps;
    diag[1] = c[1][1] + c[2][0] + eps;
    for (ih=2; ih < nh-2; ih++) {
	diag[ih] = c[ih-1][2] + c[ih][1] + c[ih+1][0] + eps;
    }
    diag[nh-2] = c[nh-3][2] + c[nh-2][1] + eps;
    diag[nh-1] = c[nh-2][2] + eps;

    offd[0][0] = cc[1][0];
    offd[1][0] = cc[1][2];
    for (ih=1; ih < nh-2; ih++) {
#ifdef SF_HAS_COMPLEX_H
	offd[0][ih] = cc[ih+1][0] + cc[ih][1];
#else
	offd[0][ih] = sf_cadd(cc[ih+1][0],cc[ih][1]);
#endif
	offd[1][ih] = cc[ih+1][2];
    }
    offd[0][nh-2] = cc[nh-2][1];
    
    cbanded_define(diag,offd);
    /* define the pentadiagonal matrix */
}

void shotfill_apply (const sf_complex *s1, 
		     const sf_complex *s2 /* input shots [nh] */, 
		     sf_complex *s        /* interpolated shot [nh] */)
/*< interpolate a shot gather >*/
{
    sf_complex *b, sum1, sum2; 
    int ih;

    for (ih=0; ih < nh; ih++) {
	s[ih] = sf_cmplx(0.,0.);
#ifdef SF_HAS_COMPLEX_H
	s0[ih] = 0.5*(s1[ih]+s2[ih]);
#else
	s0[ih] = sf_crmul(sf_cadd(s1[ih],s2[ih]),0.5);
#endif
    }

    /** Compute right-hand side: s=M2'*M1*s1 + M1'*M2*s2 **/

    /*** ih-1 <-> ih+1 ??? ***/
    for (ih=1; ih < nh-1; ih++) {
	b = a[ih];
#ifdef SF_HAS_COMPLEX_H
	sum1 = s1[ih-1]*b[2] + s1[ih]*b[1] + s1[ih+1]*b[0];
	sum2 = 
	    s2[ih-1]*conjf(b[0]) + 
	    s2[ih  ]*conjf(b[1]) +
	    s2[ih+1]*conjf(b[2]);
	s[ih-1] += sum1*b[0] + sum2*conjf(b[2]);
	s[ih  ] += sum1*b[1] + sum2*conjf(b[1]);
	s[ih+1] += sum1*b[2] + sum2*conjf(b[0]);
#else
	sum1 = sf_cadd(sf_cmul(s1[ih-1],b[2]),
		       sf_cadd(sf_cmul(s1[ih],b[1]),
			       sf_cmul(s1[ih+1],b[0])));
	sum2 = sf_cadd(
	    sf_cmul(s2[ih-1],conjf(b[0])),
	    sf_cadd(
		sf_cmul(s2[ih  ],conjf(b[1])),
		sf_cmul(s2[ih+1],conjf(b[2]))));
	s[ih-1] = sf_cadd(s[ih-1],
			  sf_cadd(sf_cmul(sum1,b[0]),
				  sf_cmul(sum2,conjf(b[2]))));
	s[ih  ] = sf_cadd(s[ih],
			  sf_cadd(sf_cmul(sum1,b[1]), 
				  sf_cmul(sum2,conjf(b[1]))));
	s[ih+1] = sf_cadd(s[ih+1],
			  sf_cadd(sf_cmul(sum1,b[2]),
				  sf_cmul(sum2,conjf(b[0]))));
#endif
    }

    cbanded_solve(s);
    /* inversion in place */

}

void shotprop_apply (const sf_complex *s1 /* input shot [nh] */,  
		     sf_complex *s        /* propagated shot [nh] */)
/*< propagate a shot gather >*/
{
    sf_complex *b, sum1; 
    int ih;

    for (ih=0; ih < nh; ih++) {
	s[ih] = sf_cmplx(0.,0.);
	s0[ih] = s1[ih];
    }

    /** Compute right-hand side: s=M2'*M1*s1 **/

    /*** ih-1 <-> ih+1 ??? ***/
    for (ih=1; ih < nh-1; ih++) {
	b = a[ih];
#ifdef SF_HAS_COMPLEX_H
	sum1 = s1[ih-1]*b[2] + s1[ih]*b[1] + s1[ih+1]*b[0];
	s[ih-1] += sum1*b[0];
	s[ih  ] += sum1*b[1];
	s[ih+1] += sum1*b[2];
#else
	sum1 = sf_cadd(sf_cmul(s1[ih-1],b[2]),
		       sf_cadd(sf_cmul(s1[ih],b[1]),
			       sf_cmul(s1[ih+1],b[0])));
	s[ih-1] = sf_cadd(s[ih-1],sf_cmul(sum1,b[0]));
	s[ih  ] = sf_cadd(s[ih  ],sf_cmul(sum1,b[1]));
	s[ih+1] = sf_cadd(s[ih+1],sf_cmul(sum1,b[2]));
#endif
    }

    cbanded_solve(s);
    /* inversion in place */

}

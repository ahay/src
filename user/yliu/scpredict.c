/* Shot gather prediction by shot continuation */
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

#include "scpredict.h"
#include "cbanded.h"

static int nh;
static sf_complex **a, *offd[2], **cc, b0[3];
static float **c, *diag, h1, s;
static void filt2matrix(const sf_complex *filt,
			float* d, sf_complex *dd);
static void predefine (float w /* log-stretch frequency */);

void scpredict_step(bool adj           /* adjoint flag          */,
		    float dw, float w  /* log-stretch frequency */,
		    int nh, float h0, 
		    float dh           /* offset parameters     */,  
		    float ds           /* shot position         */, 
		    float eps          /* regularization        */,
		    sf_complex *trace  /* common shot gather    */)
/*< shot continuation prediction step >*/
{
    sf_complex *b, sum1; 
    int ih;
    sf_complex *temp;

    if (fabsf(w) < dw) {
	return;
    }    
    
    temp = sf_complexalloc(nh);

    offd[0] = sf_complexalloc(nh-1);
    offd[1] = sf_complexalloc(nh-2);
    diag = sf_floatalloc(nh);
    /* allocate the pentadiagonal matrix diagonals */

    a = sf_complexalloc2(3,nh);
    /* filter coefficients */

    c = sf_floatalloc2(3,nh);
    cc = sf_complexalloc2(3,nh);

    cbanded_init(nh, 2);
    /* initialize complex banded matrix inversion */

    s = -0.5*ds/dh;
    h1 = h0/dh;
    /* normalize by dh */

    b0[2] = sf_cmplx((s-2.)*(s-1.)/12.,0.);
    b0[1] = sf_cmplx((s-2.)*(s+2.)/6.,0.);
    b0[0] = sf_cmplx((s+2.)*(s+1.)/12.,0.);
    /* h-independent part of the filter */

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

    for (ih=0; ih < nh; ih++) {
	temp[ih] = sf_cmplx(0.,0.);
    }

    /** Compute right-hand side: s=M2'*M1*s1 **/

    /*** ih-1 <-> ih+1 ??? ***/
    for (ih=1; ih < nh-1; ih++) {
	b = a[ih];
#ifdef SF_HAS_COMPLEX_H
	sum1 = trace[ih-1]*b[2] + trace[ih]*b[1] + trace[ih+1]*b[0];
	temp[ih-1] += sum1*b[0];
	temp[ih  ] += sum1*b[1];
	temp[ih+1] += sum1*b[2];
#else
	sum1 = sf_cadd(sf_cmul(trace[ih-1],b[2]),
		       sf_cadd(sf_cmul(trace[ih],b[1]),
			       sf_cmul(trace[ih+1],b[0])));
	temp[ih-1] = sf_cadd(temp[ih-1],sf_cmul(sum1,b[0]));
	temp[ih  ] = sf_cadd(temp[ih  ],sf_cmul(sum1,b[1]));
	temp[ih+1] = sf_cadd(temp[ih+1],sf_cmul(sum1,b[2]));
#endif
    }

    cbanded_solve(temp);

    for (ih=0; ih < nh; ih++) {
	trace[ih] = temp[ih];
    }
    free(temp);
    scpredict_close();
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
/*    float den, h; */
    sf_complex *b;
    int ih;

    for (ih=1; ih < nh-1; ih++) {
	/* loop over offsets, skip first and last */
	/* h = h1 + ih; */
	/* offset normalized by dh */

	/* den = 12.*h*h + (s-1.)*(s+1.)*w*w; */
	
	b = a[ih];

	b[0] = b0[0]; /* - 0.5*h*(I*(s-2.)/w+(2.*h+I*(s-1.)*w)/den); */
	b[1] = b0[1]; /* + h*(I*s/w+(2.*h+I*s*w)/den); */
	b[2] = b0[2]; /* - 0.5*h*(I*(s+2.)/w+(2.*h+I*(s+1.)*w)/den); */

	filt2matrix(b,c[ih],cc[ih]);
    }
}

void scpredict_close (void)
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


/* 	$Id$	 */

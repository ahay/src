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
static float complex **a, *offd[2], **cc, *s0;
static float **c, *diag, s, h1, b0[3], eps;
static void filt2matrix(const float complex *filt,
			float* d, float complex *dd);

void shotfill_init (int nh_in, float h0, float dh /* half-offset axis */, 
		    float ds                      /* shot sampling */, 
		    float eps1                    /* regularization */)
/*< initialize >*/
{
    nh = nh_in;
    offd[0] = sf_complexalloc(nh-1);
    offd[1] = sf_complexalloc(nh-2);
    diag = sf_floatalloc(nh);

    a = sf_complexalloc2(3,nh);
    c = sf_floatalloc2(3,nh);
    cc = sf_complexalloc2(3,nh);
    s0 = sf_complexalloc(nh);

    cbanded_init(nh, 2);
    
    s = -0.5*ds/dh;
    h1 = h0/dh;

    b0[0] = (2.-s)*(1.-s)/12.;
    b0[1] = (2.-s)*(2.+s)/6.;
    b0[2] = (2.+s)*(1.+s)/12.;

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

static void filt2matrix(const float complex *filt,
			float* d, float complex *dd)
{
    d[0] = crealf (filt[0] * conjf (filt[0]));
    d[1] = crealf (filt[1] * conjf (filt[1]));
    d[2] = crealf (filt[2] * conjf (filt[2]));
    dd[0] = filt[1] * conjf (filt[0]);
    dd[1] = filt[2] * conjf (filt[1]);
    dd[2] = filt[2] * conjf (filt[0]);
}

void shotfill_define (float w /* log-stretch frequency */)
/*< Fill the shot continuation matrix >*/
{
    float den, h;
    float complex *b;
    int ih;

/*
    filt2matrix(pef,d,dd);
*/  

    for (ih=1; ih < nh-1; ih++) {
	h = h1 + ih;

	den = 12.*h*h + (s-1.)*(s+1.)*w*w;
	
	b = a[ih];

	b[0] = b0[0] - 0.5*h*(I*(s-2.)/w+(2.*h+I*(s-1.)*w)/den);
	b[1] = b0[1] + h*(I*s/w+(2.*h+I*s*w)/den);
	b[2] = b0[2] - 0.5*h*(I*(s+2.)/w+(2.*h+I*(s+1.)*w)/den);

	filt2matrix(b,c[ih],cc[ih]);
    }

    diag[0] = c[1][0] + c[1][2] + eps;
    diag[1] = 2.*c[1][1] + c[2][0] + eps;
    for (ih=2; ih < nh-2; ih++) {
	diag[ih] = c[ih-1][0] + c[ih-1][2] + 2.*c[ih][1]
	    +      c[ih+1][0] + c[ih+1][2] + eps;
    }
    diag[nh-2] = c[nh-3][0] + c[nh-3][2] + 2.*c[nh-2][1] + eps;
    diag[nh-1] = c[nh-2][0] + c[nh-2][2] + eps;

    offd[0][0] = cc[1][0] + cc[1][1];
    offd[1][0] = 2.*cc[1][2];
    for (ih=1; ih < nh-2; ih++) {
	offd[0][ih] = 
	    cc[ih][0] + cc[ih+1][0] +  
	    cc[ih][1] + cc[ih+1][1];
	offd[1][ih] = 2.*cc[ih+1][2];
    }
    offd[0][nh-2] = cc[nh-2][0] + cc[nh-2][1];
    
    cbanded_define(diag,offd);
}

void shotfill_apply (const float complex *s1, 
		     const float complex *s2 /* input shots [nh] */, 
		     float complex *s        /* interpolated shot [nh] */)
/*< interpolate a shot gather >*/
{
    float complex *b, sum1, sum2; 
    int ih;

    for (ih=0; ih < nh; ih++) {
	s[ih] = 0.;
	s0[ih] = 0.5*(s1[ih]+s2[ih]);
    }

    for (ih=1; ih < nh-1; ih++) {
	b = a[ih];
	sum1 = s1[ih-1]*b[2] + s1[ih]*b[1] + s1[ih+1]*b[0];
	sum2 = 
	    s2[ih-1]*conjf(b[0]) + 
	    s2[ih  ]*conjf(b[1]) +
	    s2[ih+1]*conjf(b[2]);
	s[ih-1] += sum1*b[0] + sum2*conjf(b[2]);
	s[ih  ] += sum1*b[1] + sum2*conjf(b[1]);
	s[ih+1] += sum1*b[2] + sum2*conjf(b[0]);
    }

    s[0] -= (diag[0]-eps)*s0[0] + 
	conjf(offd[0][0])*s0[1] + conjf(offd[1][0])*s0[2];
    s[1] -= offd[0][0]*s0[0] + (diag[1]-eps)*s0[1] + 
	conjf(offd[0][1])*s0[2] + conjf(offd[1][1])*s0[3];
    for (ih=2; ih < nh-2; ih++) {
	s[ih] -= offd[1][ih-2]*s0[ih-2]+offd[0][ih-1]*s0[ih-1] + 
	    (diag[ih]-eps)*s0[ih] + 
	    conjf(offd[0][ih])*s0[ih+1] + conjf(offd[1][ih])*s0[ih+2];
    }
    s[nh-2] -= offd[1][nh-4]*s0[nh-4]+offd[0][nh-3]*s0[nh-3] + 
	(diag[nh-2]-eps)*s0[nh-2] + conjf(offd[0][nh-2])*s0[nh-1];
    s[nh-1] -= offd[1][nh-3]*s0[nh-3]+offd[0][nh-2]*s0[nh-2] + 
	(diag[nh-1]-eps)*s0[nh-1];

    cbanded_solve(s);

    for (ih=0; ih < nh; ih++) {
	s[ih] += s0[ih];
    }
}


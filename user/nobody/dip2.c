/* 2-D dip estimation */
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

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "dip2.h"
#include "div2.h"
#include "allp2.h"

static float **u1, **u2, **dp;
static int n1, n2, n;
static bool sign;

void dip2_init(int niter          /* number of linear iterations */,
	       int nx, int ny     /* dimensions */, 
	       float fx, float fy /* smoothing */, 
	       bool sign1         /* to keep sign */, 
	       bool gauss         /* to use Gaussian smoothing */,
	       bool verb          /* verbosity flag */)
/*< initialize >*/
{
    n1=nx; n2=ny; n=n1*n2;
    u1 = sf_floatalloc2(n1,n2);
    u2 = sf_floatalloc2(n1,n2);
    dp = sf_floatalloc2(n1,n2);

    div2_init(n1,n2,fx,fy,niter,gauss,verb);
    sign = sign1;
}

void dip2_close(void)
/*< free allocated storage >*/
{
    free (u1[0]); free (u1);
    free (u2[0]); free (u2);
    free (dp[0]); free (dp);
    /*   divide2_close(div0); */
    div2_close();
}

void dip2(int niter   /* number of iterations */, 
	  int nw      /* filter size */, 
	  int nj      /* filter stretch for aliasing */, 
	  bool verb   /* verbosity */, 
	  float **u   /* input data */, 
	  float** p   /* output dip */, 
	  bool **mask /* input mask for known data */)
/*< estimate local dip >*/
{
    int i, iter;
    float mean, usum, psum, ui, dpi, pi;
    allpass2 ap;
 
    ap = allpass2_init (nw,nj,n1,n2,p);

    for (iter =0; iter < niter; iter++) {
	allpass21 (false, ap, u,u2);
	allpass21 (true,  ap, u,u1);
	
	mean = 0.;
	for(i=0; i < n; i++) {
	    ui = u1[0][i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrtf(mean/n);

	usum = 0.;
	psum = 0.;

	for(i=0; i < n; i++) {
	    u1[0][i] /= mean;
	    u2[0][i] /= mean;
	    if (verb) {
		usum += u2[0][i]*u2[0][i];
		psum += p[0][i];
	    }
	}

	if (verb) sf_warning("%d %g %g", iter+1, sqrt(usum/n), psum/n);

	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[0][i]) {
		    u1[0][i] = 0.;
		    u2[0][i] = 0.;
		}
	    }
	}

	div2(u2[0],u1[0],dp[0]);

	for(i=0; i < n; i++) {
	    dpi = dp[0][i];
	    pi = p[0][i];
	    if (sign && 
		fabsf(pi) > FLT_EPSILON && 
		fabsf(pi) > fabsf (dpi)) {
		p[0][i] *= expf(dpi/pi);
	    } else {			
		p[0][i] += dpi;
	    }
	}
    } /* iter */
}

/* 	$Id$	 */

/* Two-slope estimation by plane-wave destruction */
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
/*^*/

#include "twodip2.h"
#include "twodiv2.h"
#include "div2.h"
#include "allp2.h"

static float ***u1, **u2, **u3, ***dp;
static int n1, n2, n, skip;
static const int niter=100;
static bool sign;

static void border(float** u);

void twodip2_init(int nx, int ny     /* data size */, 
		  float fx, float fy /* smoothing radius */, 
		  bool sign1         /* if keep slope signs */, 
		  bool gauss         /* Gaussian versus triangle smoothing */, 
		  bool both          /* both slopes or one */)
/*< initialize >*/
{
    n1=nx; n2=ny; n=n1*n2;
    u2 = sf_floatalloc2(n1,n2);
    u3 = sf_floatalloc2(n1,n2);

    if (both) {
	u1 = sf_floatalloc3(n1,n2,2);
	dp = sf_floatalloc3(n1,n2,2);
	twodiv2_init(2,n1,n2,fx,fy,niter,gauss,u1[0][0]);
    } else {
	u1 = sf_floatalloc3(n1,n2,1);
	dp = sf_floatalloc3(n1,n2,1);
	div2_init(n1,n2,fx,fy,niter,gauss);
    }
    
    sign = sign1;
 }

void twodip2_close(void)
/*< free allocated storage >*/
{
    free (u1[0][0]); free (u1[0]); free (u1);
    free (u2[0]); free (u2);
    free (u3[0]); free (u3);
    free (dp[0][0]); free (dp[0]); free (dp);
    twodiv2_close();
}

void twodip2(int niter        /* number of iterations */, 
	     int nw           /* filter order */, 
	     int nj1, int nj2 /* dealisianing stretch */, 
	     bool verb        /* verbosity flag */, 
	     float **u        /* input data */, 
	     float*** pq      /* output slopes */, 
	     bool **mask      /* mask for missing data */)
/*< estimate slopes >*/
{
    int i, iter;
    float mean, usum, psum, qsum, ui, dpi, pi;
    allpass2 ap, aq;
 
    ap = allpass2_init (nw,nj1,n1,n2,pq[0]);
    aq = allpass2_init (nw,nj2,n1,n2,pq[1]);
    skip = nj1>nj2? 2*nw*nj1 : 2*nw*nj2;

    for (iter =0; iter < niter; iter++) {
	allpass21 (false, ap, u,u3);
	allpass21 (false, aq, u3,u2);
	border(u2);

	allpass21 (true,  aq, u3,u1[1]);
	border(u1[1]);
	
	allpass21 (true,  ap, u,u3);
	allpass21 (false, aq, u3,u1[0]);
	border(u1[0]);

	mean = 0.;
	for(i=0; i < 2*n; i++) {
	    ui = u1[0][0][i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrt (mean/n);

	usum = 0.;
	psum = 0.;
	qsum = 0.;

	for(i=0; i < n; i++) {
	    u1[0][0][i]  /= mean;
	    u1[1][0][i]  /= mean;
	    u2[0][i]     /= mean;
	    if (verb) {
		usum += u2[0][i]*u2[0][i];
		psum += pq[0][0][i];
		qsum += pq[1][0][i];
	    }
	}

	if (verb) sf_warning("%d %g %g %g", iter+1, 
			     sqrt(usum/n), psum/n, qsum/n);

	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[0][i]) {
		    u1[0][0][i] = 0.;
		    u1[1][0][i] = 0.;
		    u2[0][i] = 0.;
		}
	    }
	}

	twodiv2(u2[0],dp[0][0]);
	
	for(i=0; i < 2*n; i++) {
	    dpi = dp[0][0][i];
	    pi = pq[0][0][i];
	    if (sign && 
		fabsf(pi) > FLT_EPSILON && 
		fabsf(pi) > fabsf (dpi)) {
		pq[0][0][i] *= expf(dpi/pi);
	    } else {			
		pq[0][0][i] += dpi;
	    }
	}
    } /* iter */
}


void otherdip2(int niter        /* number of iterations */, 
	       int nw           /* filter order */, 
	       int nj1, int nj2 /* dealising stretch */, 
	       bool verb        /* verbosity flag */, 
	       float **u        /* input data */, 
	       float*** pq      /* output slope */, 
	       bool **mask      /* mask for missing data */)
/*< estimate the second slope only >*/
{
    int i, iter;
    float mean, usum, psum, qsum, ui, dpi, pi;
    allpass2 ap, aq;
 
    ap = allpass2_init (nw,nj1,n1,n2,pq[0]);
    aq = allpass2_init (nw,nj2,n1,n2,pq[1]);
    skip = nj1>nj2? 2*nw*nj1 : 2*nw*nj2;

    for (iter =0; iter < niter; iter++) {
	allpass21 (false, ap, u,u3);
	allpass21 (false, aq, u3,u2);
	border(u2);
	
	allpass21 (true,  ap, u,u3);
	allpass21 (false, aq, u3,u1[0]);
	border(u1[0]);

	mean = 0.;
	for(i=0; i < n; i++) {
	    ui = u1[0][0][i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrt (mean/n);

	usum = 0.;
	psum = 0.;
	qsum = 0.;

	for(i=0; i < n; i++) {
	    u1[0][0][i]  /= mean;
	    u2[0][i]     /= mean;
	    if (verb) {
		usum += u2[0][i]*u2[0][i];
		psum += pq[0][0][i];
	    }
	}

	if (verb) sf_warning("%d %g %g", iter+1, sqrt(usum/n), psum/n);

	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[0][i]) {
		    u1[0][0][i] = 0.;
		    u2[0][i] = 0.;
		}
	    }
	}


	div2(u2[0],u1[0][0],dp[0][0]);

	for(i=0; i < n; i++) {
	    dpi = dp[0][0][i];
	    pi = pq[0][0][i];
	    if (sign && 
		fabsf(pi) > FLT_EPSILON && 
		fabsf(pi) > fabsf (dpi)) {
		pq[0][0][i] *= expf(dpi/pi);
	    } else {			
		pq[0][0][i] += dpi;
	    }
	}
    } /* iter */
}

static void border(float** u)
{
    int i1, i2;

    for (i2=n2-2; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    u[i2][i1]=0.;
	}
    }
    for (i2=0; i2 < n2-2; i2++) {
	for (i1=0; i1 < skip; i1++) {
	    u[i2][i1]=0.;
	}
	for (i1=n1-skip; i1 < n1; i1++) {
	    u[i2][i1]=0.;
	}
    }
}

/* 	$Id$	 */

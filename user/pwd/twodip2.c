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
#include "div2.h"
#include "allp3.h"

static float **u1, *u2, *u3, *dp, *p0;
static int n1, n2, n, skip;
static const int liter=100;
static bool sign;

static void border(float* u);

void twodip2_init(int nx, int ny     /* data size */, 
		  float fx, float fy /* smoothing radius */, 
		  bool sign1         /* if keep slope signs */,
		  bool gauss         /* Gaussian versus triangle smoothing */,
		  bool verb          /* verbosity flag */,
		  bool both          /* both slopes or one */)
/*< initialize >*/
{
    int nd[2], rect[2];

    n1=nx; n2=ny; n=n1*n2;
    u2 = sf_floatalloc(n);
    u3 = sf_floatalloc(n);

    if (both) {
	u1 = sf_floatalloc2(n,2);
	dp = sf_floatalloc(n*2);
	p0 = sf_floatalloc(n*2);
	nd[0] = nx;
	nd[1] = ny;
	rect[0] = fx;
	rect[1] = fy;
	sf_multidivn_init(2,2,n,nd,rect,u1[0],NULL,verb);
    } else {
	u1 = sf_floatalloc2(n,1);
	dp = sf_floatalloc(n);
	p0 = sf_floatalloc(n);
	div2_init(n1,n2,fx,fy,liter,gauss,true);
    }
    
    sign = sign1;
 }

void twodip2_close(void)
/*< free allocated storage >*/
{
    free (*u1); free (u1);
    free (u2);
    free (u3);
    free (dp);
    free (p0);
    sf_multidivn_close();
}

void twodip2(int niter        /* number of iterations */, 
	     int nw           /* filter order */, 
	     int nj1, int nj2 /* dealisianing stretch */,
	     bool drift       /* if shift filter */,
	     bool verb        /* verbosity flag */, 
	     float *u         /* input data */, 
	     float** pq       /* output slopes */, 
	     bool *mask       /* mask for missing data */)
/*< estimate slopes >*/
{
    int i, iter, k;
    float mean, usum, psum, qsum, ui, dpi, pi, lam, usum2;
    allpass ap, aq;
 
    ap = allpass_init (nw,nj1,n1,n2,1,drift,pq[0]);
    aq = allpass_init (nw,nj2,n1,n2,1,drift,pq[1]);
    skip = nj1>nj2? 2*nw*nj1 : 2*nw*nj2;

    allpass1 (false, false, ap, u, u3);
    allpass1 (false, false, aq, u3,u2);
    border(u2);
    
    for (iter =0; iter < niter; iter++) {
	allpass1 (false, true,  aq, u3,u1[1]);
	border(u1[1]);
	
	allpass1 (false, true,  ap, u, u3);
	allpass1 (false, false, aq, u3,u1[0]);
	border(u1[0]);

	mean = 0.;
	for(i=0; i < 2*n; i++) {
	    ui = u1[0][i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrt (mean/n);

	usum = 0.;
	psum = 0.;
	qsum = 0.;

	for(i=0; i < n; i++) {
	    u1[0][i]  /= mean;
	    u1[1][i]  /= mean;
	    u2[i]     /= mean;
	    usum += u2[i]*u2[i];
	    if (verb) {
		psum += pq[0][i];
		qsum += pq[1][i];
	    }
	    p0[i]   = pq[0][i];
	    p0[n+i] = pq[1][i];
	}

	if (verb) sf_warning("%d %g %g %g", iter+1, 
			     sqrt(usum/n), psum/n, qsum/n);

	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[0][i] = 0.;
		    u1[1][i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}

	sf_multidivn(u2,dp,liter);
	
	lam = 1.;
	for (k=0; k < 8; k++) {
	    for(i=0; i < 2*n; i++) {
		dpi = lam*dp[i];
		if (sign) {
		    pi = p0[i];
		    /* dip is nonzero and we keep its sign */ 
		    if (fabsf(pi) > FLT_EPSILON && dpi/pi >= -1.)
			pq[0][i] = p0[i]*(1.+dpi/pi);
		    /* otherwise we don't change it */
		} else {
		    pq[0][i] = p0[i]+dpi;
		}
	    }
	    
	    allpass1 (false, false, ap, u, u3);
	    allpass1 (false, false, aq, u3,u2);
	    border(u2);

	    usum2 = 0.;
	    for(i=0; i < n; i++) {
		usum2 += u2[i]*u2[i];
	    }
	    
	    if (usum2 < usum*mean*mean) break;
	    lam *= 0.5;
	}
    } /* iter */
}

void otherdip2(int niter        /* number of iterations */, 
	       int nw           /* filter order */, 
	       int nj1, int nj2 /* dealising stretch */, 
	       bool drift       /* if shift filter */,
	       bool verb        /* verbosity flag */, 
	       float *u         /* input data */, 
	       float** pq       /* output slope */, 
	       bool *mask       /* mask for missing data */)
/*< estimate the second slope only >*/
{
    int i, iter, k;
    float mean, usum, psum, ui, dpi, pi, lam, usum2;
    allpass ap, aq;
 
    ap = allpass_init (nw,nj1,n1,n2,1,drift,pq[0]);
    aq = allpass_init (nw,nj2,n1,n2,1,drift,pq[1]);
    skip = nj1>nj2? 2*nw*nj1 : 2*nw*nj2;

    allpass1 (false, false, ap, u, u3);
    allpass1 (false, false, aq, u3,u2);
    border(u2);

    for (iter =0; iter < niter; iter++) {
	allpass1 (false, true,  ap, u, u3);
	allpass1 (false, false, aq, u3,u1[0]);
	border(u1[0]);

	mean = 0.;
	for(i=0; i < n; i++) {
	    ui = u1[0][i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrt (mean/n);

	usum = 0.;
	psum = 0.;

	for(i=0; i < n; i++) {
	    u1[0][i]  /= mean;
	    u2[i]     /= mean;
	    usum += u2[i]*u2[i];
	    if (verb) psum += pq[0][i];
	    p0[i] = pq[0][i];
	}

	if (verb) sf_warning("%d %g %g", iter+1, sqrt(usum/n), psum/n);

	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[0][i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}


	div2(u2,u1[0],dp);

	lam = 1.;
	for (k=0; k < 8; k++) {
	    for(i=0; i < n; i++) {
		dpi = lam*dp[i];
		if (sign) {
		    pi = p0[i];
		    /* dip is nonzero and we keep its sign */ 
		    if (fabsf(pi) > FLT_EPSILON && dpi/pi >= -1.)
			pq[0][i] = p0[i]*(1.+dpi/pi);
		    /* otherwise we don't change it */
		} else {
		    pq[0][i] = p0[i]+dpi;
		}
	    }
	    
	    allpass1 (false, false, ap, u, u3);
	    allpass1 (false, false, aq, u3,u2);
	    border(u2);
	    
	    usum2 = 0.;
	    for(i=0; i < n; i++) {
		usum2 += u2[i]*u2[i];
	    }
	    
	    if (usum2 < usum*mean*mean) break;
	    lam *= 0.5;
	}
    } /* iter */
}

static void border(float* u)
{
    int i1, i2;

    for (i2=n2-2; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    u[i2*n1+i1]=0.;
	}
    }
    for (i2=0; i2 < n2-2; i2++) {
	for (i1=0; i1 < skip; i1++) {
	    u[i2*n1+i1]=0.;
	}
	for (i1=n1-skip; i1 < n1; i1++) {
	    u[i2*n1+i1]=0.;
	}
    }
}

/* 	$Id$	 */

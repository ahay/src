/* 2-D omnidirectional dip estimation */
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

#include "odip2.h"
#include "opwd2.h"
#include "matrixdivn.h"

static float *dat, *dp, *p0, **den;
static int n, n1, n2;

void odip2_init(int m1, int m2 /* dimensions */, 
		int* rect      /* smoothing radius [2] */, 
		int niter      /* number of iterations */,    
		bool verb      /* verbosity flag */)
/*< initialize >*/
{
    int nn[2];
    
    n1=m1;
    n2=m2;
    n = n1*n2;

    dat = sf_floatalloc(2*n);
    dp  = sf_floatalloc(2*n);
    p0  = sf_floatalloc(2*n);
    den = sf_floatalloc2(n,4);

    nn[0]=n1;
    nn[1]=n2;

    matrixdivn_init (2, n, nn, rect, den, niter, verb);
}

void odip2_close(void)
/*< free allocated storage >*/
{
    free (dat);
    free (dp);
    free (p0);
    free (*den);
    free (den);
    matrixdivn_close();
}

void odip2(int niter   /* number of nonlinear iterations */, 
	   int nw      /* filter size */, 
	   float *u    /* input data */, 
	   float* p    /* output dips */)
/*< estimate local dip >*/
{
    int i, j, iter, k;
    float usum, usum2, lam, mean, one;
    omni2 ap;
 
    ap = opwd2_init (nw,n1,n2,p,p+n);

    for (i=0; i < n; i++) {
	one = hypotf(p[i],p[i+n])+SF_EPS;
	p[i] /= one;
	p[i+n] /= one;
    } 
    opwd21(false,false,ap,u,dat);
    opwd12(false,false,ap,u,dat+n);
    
    usum = 0.0f;
    for(i=0; i < 2*n; i++) {
	usum += dat[i]*dat[i];
    }
	
    for (iter =0; iter < niter; iter++) {
	opwd21(true,false,ap,u,den[0]);
	opwd21(false,true,ap,u,den[1]);
	opwd12(true,false,ap,u,den[2]);
	opwd12(false,true,ap,u,den[3]);

	for (i=0; i < 2*n; i++) {
	    p0[i] = p[i];
	}

	mean = 0.0f;
	for (j=0; j < 4; j++) {
	    for (i=0; i < n; i++) {
		mean += den[j][i]*den[j][i];
	    }
	}
	mean = sqrtf (mean/(4*n));
	
	for (i=0; i < n; i++) {
	    dat[i]   /= mean;
	    dat[i+n] /= mean;
	}
	for (j=0; j < 4; j++) {
	    for (i=0; i < n; i++) {
		den[j][i] /= mean;
	    }
	}
	
	matrixdivn (dat,dp);

	lam = 1.0f;
	for (k=0; k < 8; k++) {
	    for(i=0; i < 2*n; i++) {
		p[i] = p0[i]+lam*dp[i];
	    }
	    for (i=0; i < n; i++) {
		one = hypotf(p[i],p[i+n])+SF_EPS;
		p[i] /= one;
		p[i+n] /= one;
	    } 

	    opwd21(false,false,ap,u,dat);
	    opwd12(false,false,ap,u,dat+n);

	    usum2 = 0.0f;
	    for(i=0; i < 2*n; i++) {
		usum2 += dat[i]*dat[i];
	    }

	    if (usum2 < usum) break;
	    lam *= 0.5f;
	}

	usum = usum2;
    } /* iter */

    opwd2_close(ap);
}

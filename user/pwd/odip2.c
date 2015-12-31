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

static float *dat, *dp, *p0, **mat;
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
    mat = sf_floatalloc2(n,4);

    nn[0]=n1;
    nn[1]=n2;

    matrixdivn_init (2, n, nn, rect, mat, niter, verb);
}

void odip2_close(void)
/*< free allocated storage >*/
{
    free (dat);
    free (dp);
    free (p0);
    free (*mat);
    free (mat);
    matrixdivn_close();
}

void odip2(int niter   /* number of nonlinear iterations */, 
	   int nw      /* filter size */, 
	   float *u    /* input data */, 
	   float* p    /* output dips */,
	   float eps   /* regularization */)
/*< estimate local dip >*/
{
    int i, iter, k;
    float usum, usum2, lam, mean, one;
    omni2 ap;
 
    ap = opwd2_init (nw,n1,n2,p,p+n);
    opwd21(false,false,ap,u,dat);

    for (i=0; i < n; i++) {
	dat[n+i] = eps*(1.0f-p[i]*p[i]-p[n+i]*p[n+i]);
    }

    for (iter =0; iter < niter; iter++) {
	opwd21(true,false,ap,u,mat[0]);
	opwd21(false,true,ap,u,mat[1]);
	
	for (i=0; i < n; i++) {
	    mat[2][i] = 2*eps*p[i];
	    mat[3][i] = 2*eps*p[n+i];
	}

	usum = 0.0f;
	for(i=0; i < 2*n; i++) {
	    p0[i]   = p[i];
	    usum += dat[i]*dat[i];
	}

	mean = 0.0f;
	for (i=0; i < n; i++) {
	    mean += mat[0][i]*mat[0][i] +  mat[1][i]*mat[1][i];
	}
	mean = sqrtf (mean/(2*n));
	
	for (i=0; i < n; i++) {
	    mat[0][i] /= mean;
	    mat[1][i] /= mean;
	    dat[i] /= mean;
	}

	matrixdivn (dat,dp);

	lam = 1.0f;
	for (k=0; k < 8; k++) {
	    for(i=0; i < 2*n; i++) {
		p[i] = p0[i]+lam*dp[i];
	    }
	    for(i=0; i < n; i++) {
		one = hypotf(p[i],p[n+i]);
		p[i]   /= one;
		p[n+i] /= one;
	    }

	    opwd21(false,false,ap,u,dat);
	    for (i=0; i < n; i++) {
		dat[n+i] = eps*(1.0f-p[i]*p[i]-p[n+i]*p[n+i]);
	    }

	    usum2 = 0.0f;
	    for(i=0; i < 2*n; i++) {
		usum2 += dat[i]*dat[i];
	    }

	    if (usum2 < usum) break;
	    lam *= 0.5f;
	}
    } /* iter */

    opwd2_close(ap);
}

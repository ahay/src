/* Estimate two frequency components */
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

#include "expont2.h"
#include "expder2.h"

static float *u1, *u2, *dp;
static int n1, n2, n;

void twofreq2_init(int nx, int ny     /* data size */, 
		   int fx, int fy /* smoothing radius */)
/*< initialize >*/
{
    int nd[2], rect[2];

    n1=nx; n2=ny; n=n1*n2;
    u1 = sf_floatalloc(n*4);
    u2 = sf_floatalloc(n);
    dp = sf_floatalloc(n*4);

    nd[0] = nx;
    nd[1] = ny;
    rect[0] = fx;
    rect[1] = fy;

    sf_multidivn_init(4,2,n,nd,rect,u1,NULL,true);
}

void twofreq2_close(void)
/*< free allocated storage >*/
{
    free (u1);
    free (u2);
    free (dp);
    sf_multidivn_close();
}

void twofreq2(int niter  /* number of iterations */, 
	      bool verb  /* verbosity flag */, 
	      float *u   /* input data */, 
	      float** pq /* estimated frequencies */,
	      bool *mask /* mask for missing data */)
/*< estimate >*/
{
    int i, iter;
    double mean, usum=0., psum1=0., psum2=0., psum3=0., psum4=0., ui;
 
    expont2_init(n1,n2,pq);
    expder2_init(n1,n2,pq);

    for (iter =0; iter < niter; iter++) {
	expont2_lop (false,false,n,n,u,u2);
	expder2_lop (false,false,n,4*n,u,u1);

	mean = 0.;
	for(i=0; i < 4*n; i++) {
	    ui = u1[i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrt (mean/n);

	if (verb) {
	    usum = 0.;
	    psum1 = psum2 = psum3 = psum4 = 0.;
	}

	for(i=0; i < 4*n; i++) {
	    u1[i] /= mean;
	}

	for(i=0; i < n; i++) {
	    u2[i] /= mean;
	    if (verb) {
		usum += u2[i]*u2[i];
		psum1 += pq[0][i]*pq[0][i];
		psum2 += pq[1][i]*pq[1][i];
		psum3 += pq[2][i]*pq[2][i];
		psum4 += pq[3][i]*pq[3][i];
	    }
	}

	if (verb) sf_warning("%d %g %g %g %g %g", iter+1, 
			     sqrt(usum/n), psum1/n, psum2/n, psum3/n, psum4/n);

	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}

	sf_multidivn(u2,dp,niter);

	for(i=0; i < 4*n; i++) {
	    pq[0][i] += dp[i];
	}
    } /* iter */
}

/* 	$Id: twofreq2.c 714 2004-07-19 20:39:50Z fomels $	 */

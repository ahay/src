#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "twofreq2.h"
#include "twodiv2.h"
#include "expont2.h"
#include "expder2.h"

static float *u1, *u2, *dp;
static int n1, n2, n;
static const int niter=100;

void twofreq2_init(int nx, int ny, float fx, float fy, bool gauss)
{
    n1=nx; n2=ny; n=n1*n2;
    u1 = sf_floatalloc(n*4);
    u2 = sf_floatalloc(n);
    dp = sf_floatalloc(n*4);

    twodiv2_init(4,n1,n2,fx,fy,niter,gauss,u1);
}

void twofreq2_close(void)
{
    free (u1);
    free (u2);
    free (dp);
    twodiv2_close();
}

void twofreq2(int niter, bool verb, float *u, float** pq)
{
    int i, iter;
    float mean, usum, psum1, psum2, psum3, psum4, ui;
 
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
	twodiv2(u2,dp);

	for(i=0; i < 4*n; i++) {
	    pq[0][i] += dp[i];
	}
    } /* iter */
}

/* 	$Id: twofreq2.c 714 2004-07-19 20:39:50Z fomels $	 */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "twodip2.h"
#include "twodiv2.h"
#include "allp2.h"

static float ***u1, **u2, **u3, ***dp;
static int n1, n2, n, skip;
static const int niter=100;
static bool sign;

static void border(float** u);

void twodip2_init(int nx, int ny, float fx, float fy, bool sign1, bool gauss)
{
    n1=nx; n2=ny; n=n1*n2;
    u1 = sf_floatalloc3(n1,n2,2);
    u2 = sf_floatalloc2(n1,n2);
    u3 = sf_floatalloc2(n1,n2);
    dp = sf_floatalloc3(n1,n2,2);

    twodiv2_init(2,n1,n2,fx,fy,niter,gauss,u1[0][0]);
    sign = sign1;
}

void twodip2_close(void)
{
    free (u1[0][0]); free (u1[0]); free (u1);
    free (u2[0]); free (u2);
    free (u3[0]); free (u3);
    free (dp[0][0]); free (dp[0]); free (dp);
    twodiv2_close();
}

void twodip2(int niter, int nw, int nj1, int nj2, 
	     bool verb, float **u, float*** pq)
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

static void border(float** u)
{
    int i1, i2;

    for (i2=n2-2; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    u2[i2][i1]=0.;
	}
    }
    for (i2=0; i2 < n2-2; i2++) {
	for (i1=0; i1 < skip; i1++) {
	    u2[i2][i1]=0.;
	}
	for (i1=n1-skip; i1 < n1; i1++) {
	    u2[i2][i1]=0.;
	}
    }
}

/* 	$Id$	 */

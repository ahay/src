#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "dip2.h"
/* #include "divide2.h" */
#include "divide.h"
#include "allp2.h"

static float **u1, **u2, **dp;
static int n1, n2, n;
static const int niter=100;
static bool sign;
/* static div2 div0; */

void dip2_init(int nx, int ny, float fx, float fy, bool sign1)
{
    n1=nx; n2=ny; n=n1*n2;
    u1 = sf_floatalloc2(n1,n2);
    u2 = sf_floatalloc2(n1,n2);
    dp = sf_floatalloc2(n1,n2);

    /* div0 = divide2_init (n1,n2,eps,lam); */
    divide_init(n1,n2,fx,fy,niter);
    sign = sign1;
}

void dip2_close(void)
{
    free (u1[0]); free (u1);
    free (u2[0]); free (u2);
    free (dp[0]); free (dp);
    /*   divide2_close(div0); */
    divide_close();
}

void dip2(int niter, int nw, int nj, bool verb, float **u, float** p)
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

	mean = sqrt (mean/n);

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

	/* divide2 (div0, u2, u1, dp); */
	divide(u2[0],u1[0],dp[0]);

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

/* 	$Id: dip2.c,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */

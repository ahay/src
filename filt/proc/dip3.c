#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "dip3.h"
#include "divide1.h"
#include "allp3.h"

static float ***u1, ***u2, ***dp;
static int n1, n2, n3, n;
static bool sign;

void dip3_init(int nx, int ny, int nz, float eps, float lam, bool sign1)
{
    n1=nx; n2=ny; n3=nz; n=n1*n2*n3;
    u1 = sf_floatalloc3(n1,n2,n3);
    u2 = sf_floatalloc3(n1,n2,n3);
    dp = sf_floatalloc3(n1,n2,n3);

    divide1_init (n1,eps,lam);
    sign = sign1;
}

void dip3_close(void)
{
    free (u1[0][0]); free (u1[0]); free (u1);
    free (u2[0][0]); free (u2[0]); free (u2);
    free (dp[0][0]); free (dp[0]); free (dp);
    divide1_close();
}

void dip3(int dip, int niter, int nw, int nj, bool verb, 
	  float ***u, float*** p)
{
    int i, iter;
    float mean, usum, psum, ui, dpi, pi;
    allpass ap;
 
    ap = allpass_init (nw,nj,n1,n2,n3,p);

    for (iter =0; iter < niter; iter++) {
	if (dip == 1) {
	    allpass1 (false, ap, u,u2);
	    allpass1 (true,  ap, u,u1);
	} else {
	    allpass2 (false, ap, u,u2);
	    allpass2 (true,  ap, u,u1);
	}
	
	mean = 0.;
	for(i=0; i < n; i++) {
	    ui = u1[0][0][i];
	    mean += ui*ui;
	}
	if (mean == 0.) return;

	mean = sqrt (mean/n);

	usum = 0.;
	psum = 0.;

	for(i=0; i < n; i++) {
	    u1[0][0][i] /= mean;
	    u2[0][0][i] /= mean;
	    if (verb) {
		usum += u2[0][0][i]*u2[0][0][i];
		psum += p[0][0][i];
	    }
	}

	if (verb) sf_warning("%d %g %g", iter+1, sqrt(usum/n), psum/n);

	divide3 (n3, n2, u2, u1, dp);

	for(i=0; i < n; i++) {
	    dpi = dp[0][0][i];
	    pi = p[0][0][i];
	    if (sign && 
		fabsf(pi) > FLT_EPSILON && 
		fabsf(pi) > fabsf (dpi)) {
		p[0][0][i] *= expf(dpi/pi);
	    } else {			
		p[0][0][i] += dpi;
	    }
	}
    } /* iter */
}

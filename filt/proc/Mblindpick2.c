/* Automatic picking from semblance-like panels using shaping regularization.

Takes: < semblance.rsf > pick.rsf
*/
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "divide.h"

int main(int argc, char* argv[])
{
    int nt, ns, nx, n, i, ix, is, it, niter, imax;
    float **slice, *pick0, *pick, *ampl;
    float s0, ds, eps, lam, t, asum, ai, amax, ap, am, num, den, dx;
    bool gauss;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ns)) ns=1;
    nx = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o2",&s0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");
    sf_putint(out,"n2",1);

    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* Vertical smoothness */
    if (!sf_getfloat("lam",&lam)) lam=1.;
    /* Horizontal smoothness */

    eps = sqrtf(12*eps+1.);
    lam = sqrtf(12*lam+1.);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getbool("gauss",&gauss)) gauss = false;
    /* if y, use exact Gaussian for smoothing */

    divide_init(nt,nx,eps,lam,niter,gauss);
    
    slice = sf_floatalloc2(nt,ns);

    n = nt*nx;
    pick0 = sf_floatalloc(n);
    pick = sf_floatalloc(n);
    ampl = sf_floatalloc(n);

    for (ix = 0; ix < nx; ix++) {
	sf_read (slice[0],sizeof(float),nt*ns,in);

	/* pick blind maximum */
	for (it = 0; it < nt; it++) {
	    i = ix*nt+it;

	    imax = 0;
	    amax = 0.;
	    for (is = 0; is < ns; is++) {
		t = slice[is][it];
		if (t > amax) {
		    amax = t;
		    imax = is;
		}
	    }

	    /* quadratic interpolation for sub-pixel accuracy */
	    if (imax > 0 && imax < ns-1) {
		am = slice[imax-1][it];
		ap = slice[imax+1][it];
		num = 0.5*(am-ap);
		den = am+ap-2.*amax;
		dx = num*den/(den*den + FLT_EPSILON);
		if (fabsf(dx) >= 1.) dx=0.;

		ampl[i] = amax - 0.5*dx*dx*den;
		pick0[i] = s0+(imax+dx)*ds;
	    } else {
		ampl[i] = amax;
		pick0[i] = s0+imax*ds;
	    }
	}
    }
    
    /* normalize amplitudes */
    asum = 0.;
    for (i = 0; i < n; i++) {
	ai = ampl[i];
	asum += ai*ai;
    }
    asum = sqrtf (asum/n);

    for(i=0; i < n; i++) {
	ampl[i] /= asum;
	pick0[i] *= ampl[i];
    }

    divide(pick0,ampl,pick);
    sf_write (pick,sizeof(float),n,out);	
    
    exit (0);
}

/* 	$Id: Mblindpick2.c,v 1.2 2004/03/15 06:57:59 fomels Exp $	 */


/* 2-D two dip estimation by plane wave destruction.

Takes: < data.rsf > dip.rsf
*/

#include <math.h>

#include <rsf.h>

#include "twodip2.h"
#include "mask6.h"

int main (int argc, char *argv[])
{
    int n1,n2, n12, niter, nw, nj1, nj2, i;
    float eps, lam, p0, q0, **u, ***p;
    bool verb, sign, gauss, both, **m;
    sf_file in, out, mask, dip1, dip2;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2= in input");
    n12 = n1*n2;
    
    sf_putint(out,"n3",2);

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1; 
    /* vertical smoothness */
    if (!sf_getfloat("lam",&lam)) lam=1; 
    /* horizontal smoothness */

    eps = sqrtf(12*eps+1.);
    lam = sqrtf(12*lam+1.);

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing for first dip */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing for second dip */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getbool("sign",&sign)) sign = false;
    /* if y, keep dip sign constant */
    if (!sf_getbool("gauss",&gauss)) gauss = false;
    /* if y, use exact Gaussian for smoothing */
    if (!sf_getbool("both",&both)) both = true;
    /* if y, estimate both dips */
    
    /* initialize dip estimation */
    twodip2_init(n1, n2, eps, lam, sign, gauss, both);

    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc3(n1,n2,2);

    if (NULL != sf_getstring("dip1")) {
	p0 = 1.;
	dip1 = sf_input("dip1");
    } else {
	if (!sf_getfloat("p0",&p0)) p0=1.;
	/* initial first dip */
	dip1 = NULL;
    }

    if (NULL != sf_getstring("dip2")) {
	q0 = 0.;
	dip2 = sf_input("dip2");
    } else {
	if (!sf_getfloat("q0",&q0)) q0=0.;
	/* initial second dip */
	dip2 = NULL;
    }

    /* initialize dips */
    if (NULL == dip1) {
	for(i=0; i < n12; i++) {
	    p[0][0][i] = p0;
	}
    } else {
	sf_floatread(p[0][0],n12,dip1);
    }
    if (NULL == dip2) {
	for(i=0; i < n12; i++) {
	    p[1][0][i] = q0;
	}
    } else {
	sf_floatread(p[1][0],n12,dip2);
    }
  
    /* initialize mask */
    if (NULL != sf_getstring("mask")) {
	m = sf_boolalloc2(n1,n2);

	mask = sf_input("mask");
	sf_floatread(u[0],n12,mask);
	sf_fileclose(mask);

	mask6 (nw, nj1, nj2, n1, n2, u, m);
    } else {
	m = NULL;
    }

    /* read data */
    sf_floatread(u[0],n12,in);

    /* estimate dip */
    if (both) {
	twodip2(niter, nw, nj1, nj2, verb, u, p, m);
    } else {
	otherdip2(niter, nw, nj1, nj2, verb, u, p, m);
    }

    /* write dips */
    sf_floatwrite(p[0][0],n12*2,out);
     
    exit (0);
}

/* 	$Id$	 */

/* 2-D two dip estimation by plane wave destruction.

Takes: < data.rsf > dip.rsf
*/

#include <math.h>

#include <rsf.h>

#include "twodip2.h"

int main (int argc, char *argv[])
{
    int n1,n2, n12, niter, nw, nj1, nj2, i;
    float eps, lam, p0, q0, **u, ***p;
    bool verb, sign, gauss;
    sf_file in, out;

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

    if (!sf_getfloat("p0",&p0)) p0=1.;
    /* initial first dip */
    if (!sf_getfloat("q0",&q0)) q0=0.;
    /* initial second dip */

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
    if (!sf_getbool("gauss",&gauss)) gauss = true;
    /* if y, use exact Gaussian for smoothing */

    /* initialize dip estimation */
    twodip2_init(n1, n2, eps, lam, sign, gauss);

    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc3(n1,n2,2);

    /* read data */
    sf_floatread(u[0],n12,in);

    /* initialize dip */
    for(i=0; i < n12; i++) {
	p[0][0][i] = p0;
	p[1][0][i] = q0;
    }
  
    /* estimate dip */
    twodip2(niter, nw, nj1, nj2, verb, u, p);

    /* write dips */
    sf_floatwrite(p[0][0],n12*2,out);
     
    sf_close();
    exit (0);
}

/* 	$Id: Mtwodip2.c,v 1.3 2004/04/19 21:51:46 fomels Exp $	 */

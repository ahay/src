/* 2-D two spectral component estimation.

Takes: < data.rsf > dip.rsf
*/

#include <math.h>

#include <rsf.h>

#include "twofreq2.h"
#include "mask4freq.h"

int main (int argc, char *argv[])
{
    int n1,n2, n12, niter, i;
    float eps, lam, p0, q0, p1, q1, *u, **p;
    bool verb, *m;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2= in input");
    n12 = n1*n2;
    
    sf_putint(out,"n3",4);

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1; 
    /* vertical smoothness */
    if (!sf_getfloat("lam",&lam)) lam=1; 
    /* horizontal smoothness */

    eps = sqrtf(12*eps+1.);
    lam = sqrtf(12*lam+1.);

    if (!sf_getfloat("p0",&p0)) p0=1.;
    /* initial first component */
    if (!sf_getfloat("q0",&q0)) q0=1.;
    /* initial first component */
    if (!sf_getfloat("p1",&p1)) p1=-1.;
    /* initial second component */
    if (!sf_getfloat("q1",&q1)) q1=1.;
    /* initial second component */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    /* initialize spectral estimation */
    twofreq2_init(n1, n2, (int) eps, (int) lam);

    u = sf_floatalloc(n12);
    p = sf_floatalloc2(n12,4);

    /* allocate mask */
    if (NULL != sf_getstring("mask")) {
	m = sf_boolalloc(n12);
	mask = sf_input("mask");
    } else {
	mask = NULL;
	m = NULL;
    }

    /* initialize mask */
    if (NULL != mask) {
	sf_floatread(u,n12,mask);
	mask4freq (2, 1, n1, n2, u, m);
    }

    /* read data */
    sf_floatread(u,n12,in);

    /* initialize components */
    for(i=0; i < n12; i++) {
	p[0][i] = p0;
	p[1][i] = q0;
	p[2][i] = p1;
	p[3][i] = q1;
    }
  
    /* estimate components */
    twofreq2(niter, verb, u, p, m);

    /* write components */
    sf_floatwrite(p[0],n12*4,out);

    exit (0);
}

/* 	$Id: Mtwofreq2.c 704 2004-07-13 18:22:06Z fomels $	 */

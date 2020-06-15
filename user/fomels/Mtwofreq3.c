/* 2-D two spectral component estimation.

Takes: < data.rsf > dip.rsf
*/

#include <math.h>

#include <rsf.h>

#include "twofreq3.h"
#include "mask4freq3.h"

int main (int argc, char *argv[])
{
    int n1,n2,n3, n123, niter, i;
    float eps, lam, eps2, p0, q0, p1, q1, p0q, q0q, p1q, q1q, *u, **p;
    bool verb, *m;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("Need n3= in input");
    n123 = n1*n2*n3;
    
    sf_putint(out,"n4",8);

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1; 
    /* vertical smoothness */
    if (!sf_getfloat("lam",&lam)) lam=1; 
    /* horizontal smoothness */
    if (!sf_getfloat("eps2",&eps2)) eps2=1;
    /* third dimension smoothness */

    eps = sqrtf(12*eps+1.);
    eps2 = sqrtf(12*eps2+1.);
    lam = sqrtf(12*lam+1.);

    if (!sf_getfloat("p0",&p0)) p0=1.;
    /* initial first component */
    if (!sf_getfloat("q0",&q0)) q0=1.;
    /* initial first component */
    if (!sf_getfloat("p1",&p1)) p1=-1.;
    /* initial second component */
    if (!sf_getfloat("q1",&q1)) q1=1.;
    /* initial second component */

    if (!sf_getfloat("p0q",&p0q)) p0q=1.;
    /* initial first component */
    if (!sf_getfloat("q0q",&q0q)) q0q=1.;
    /* initial first component */
    if (!sf_getfloat("p1q",&p1q)) p1q=-1.;
    /* initial second component */
    if (!sf_getfloat("q1q",&q1q)) q1q=1.;
    /* initial second component */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    /* initialize spectral estimation */
    twofreq3_init(n1, n2, n3, (int) eps, (int) lam, (int) eps2);

    u = sf_floatalloc(n123);
    p = sf_floatalloc2(n123,8);
    /* q = sf_floatalloc2(n123,4); */
 
    /* allocate mask */
    if (NULL != sf_getstring("mask")) {
	m = sf_boolalloc(n123);
	mask = sf_input("mask");
    } else {
	mask = NULL;
	m = NULL;
    }

    /* initialize mask */
    if (NULL != mask) {
	sf_floatread(u,n123,mask);
	mask4freq3 (2, 1, n1, n2, n3, u, m);
    }

    /* read data */
    sf_floatread(u,n123,in);

    /* initialize components */
    for(i=0; i < n123; i++) {
	p[0][i] = p0;
	p[1][i] = q0;
	p[2][i] = p1;
	p[3][i] = q1;
        p[4][i] = p0q;
        p[5][i] = q0q;
        p[6][i] = p1q;
        p[7][i] = q1q;
    }
  
    /* estimate components */
    twofreq3(niter, verb, u, p, m);

    /* write components */
    sf_floatwrite(p[0],n123*8,out);

    exit (0);
}

/* 	$Id: Mtwofreq3.c 704 2004-07-13 18:22:06Z fomels $	 */

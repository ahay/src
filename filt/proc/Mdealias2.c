/* 2-D trace interpolation to a denser grid using PWD.

Takes: < in.rsf > out.rsf

It may be necessary to bandpass the data before and after dealiasing 
to ensure that the temporal spectrum is banded. Rule of thumb: if 
max(jx,jy)=N, the temporal bandwidth should be 1/N of Nyquist.
*/

#include <math.h>

#include <rsf.h>

#include "dip2.h"
#include "interpd.h"

int main (int argc, char *argv[])
{
    int n1, n2, n12, niter, nf, i, i1, m2, m12, nj;
    bool verb, sign, gauss;
    float eps, lam, d2, p0;
    float **u1, **uu1, **p;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of dip-estimation iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.; 
    /* vertical dip smoothness */
    if (!sf_getfloat("lam",&lam)) lam=1.; 
    /* horizontal dip smoothness */

    eps = sqrtf(12*eps+1.);
    lam = sqrtf(12*lam+1.);

    if (!sf_getint("order",&nf)) nf=1;
    /* [1,2,3] dip filter accuracy */
    if (nf < 1 || nf > 3) sf_error ("accuracy must be between 1 and 3");

    if(SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if(!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if(!sf_getbool("sign",&sign)) sign = false;
    /* if y, keep dip sign constant */
    if (!sf_getbool("gauss",&gauss)) gauss = true;
    /* if y, use exact Gaussian for smoothing */

    if (!sf_getfloat("p0",&p0)) p0=0.;
    /* initial dip */
    if (!sf_getint("nj",&nj)) nj=2;
    /* antialiasing */
    
    m2 = n2*2-1; 
    if (n2 > 1) d2 /= 2;
    
    n12 = n1*n2;
    m12 = n1*m2;
    
    sf_putint(out,"n2",m2); 
    sf_putfloat(out,"d2",d2);

    u1 = sf_floatalloc2(n1,n2);

    p = sf_floatalloc2(n1,n2);

    uu1 = sf_floatalloc2(n1,m2);

    dip2_init(n1, n2, eps, lam, sign, gauss);
    interp_init (n1, 0.0001, verb);

    sf_read(u1[0],sizeof(float),n12,in);

    if (verb) sf_warning("Estimating slopes...");

    for (i1=0; i1 < n2; i1++) {
	for (i=0; i < n1; i++) {
	    p[i1][i]=p0;
	}
    }

    dip2(niter, nf, nj, verb, u1, p);

    if (verb) sf_warning("Interpolating...");

    interp2(n2,u1,uu1,p);
    
    sf_write(uu1[0],sizeof(float),m12,out);
  
    exit (0);
}

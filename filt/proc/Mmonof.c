/* Mono-frequency wavelet estimation.

Takes: < data.rsf ma=ma.rsf > monof.rsf
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "monof.h"

int main(int argc, char* argv[])
{
    int n2, i2, nk, ik, niter, i0;
    float k0, dk, k, a0, f, *data, a;

    bool verb;
    sf_file in, out, ma;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    ma = sf_output("ma");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nk)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");

    if (!sf_getfloat("a0",&a0)) a0=1.;
    /* starting sharpness */
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(ma,"n1",2);
    sf_putint(ma,"nk",nk);
    sf_putfloat(ma,"dk",dk);
    sf_putfloat(ma,"k0",k0);
    sf_fileflush(ma,in);

    data = sf_floatalloc(nk);

    for (i2=0; i2 < n2; i2++) {
	sf_read(data,sizeof(float),nk,in);

	a = 0.;
	i0 = 0;
	for (ik=0; ik < nk; ik++) {
	    f = data[ik];
	    if (f > a) {
		a = f;
		i0 = ik;
	    }
	}

	a = monof(data,i0,niter,a0,nk,2.*SF_PI*dk,verb);

	k = (float) i0;
         
	sf_write(&a,sizeof(float),1,ma);
	sf_write(&k,sizeof(float),1,ma);
        
	sf_write (data,sizeof(float),nk,out);
    }
    
    sf_close();
    exit (0);
}

/* 	$Id: Mmonof.c,v 1.3 2004/04/08 14:03:57 fomels Exp $	 */

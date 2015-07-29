#include <stdio.h>
#include <math.h>

#include <rsf.h>

#include "burg.h"

int main(void)
{
    const int n=10, nc=5, nf=3;
    float **x, a[3];
    sf_complex p, xi;
    int i, ic;

    x = sf_floatalloc2(n,nc);
    p = cexpf(sf_cmplx(0.,0.5));

    for (ic=0; ic < nc; ic++) {
	xi = sf_cmplx(ic+1.,0.);
	x[ic][0] = ic+1.;
	for (i=1; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H	    
	    xi *= p;
#else
	    xi = sf_cmul(xi,p);
#endif
	    x[ic][i] = crealf(xi);
	}
    }

    burg_init (n,nc,nf);
    burg_apply (x[0],a);
    burg_close();

    printf("p=(%g,%g)\nq=%g\na={%g,%g,%g}\n",
	   crealf(p),cimagf(p),
	   2*crealf(p),
	   a[0],a[1],a[2]);
    exit(0);
}
	   
	   

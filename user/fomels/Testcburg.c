#include <stdio.h>
#include <math.h>

#include <rsf.h>

#include "cburg.h"

int main(void)
{
    const int n=10, nc=5, nf=2;
    sf_complex p, **x, a[2];
    int i, ic;

    x = sf_complexalloc2(n,nc);
    p = cexpf(sf_cmplx(0.,0.5));

    for (ic=0; ic < nc; ic++) {
	x[ic][0] = sf_cmplx(ic+1.,0.);
	for (i=1; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[ic][i] = x[ic][i-1]*p;
#else
	    x[ic][i] = sf_cmul(x[ic][i-1],p);
#endif
	}
    }

    cburg_init (n,nc,nf);
    cburg_apply (x[0],a);
    cburg_close();

    printf("p=(%g,%g)\na={(%g,%g),(%g,%g)}\n",
	   crealf(p),cimagf(p),
	   crealf(a[0]),cimagf(a[0]),
	   crealf(a[1]),cimagf(a[1]));

    exit(0);
}
	   
	   

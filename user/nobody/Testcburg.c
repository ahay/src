#include <stdio.h>
#include <math.h>

#include <rsf.h>

#include "cburg.h"

int main(void)
{
    const int n=10, nc=5, nf=2;
    float complex p, **x, a[2];
    int i, ic;

    x = sf_complexalloc2(n,nc);

    p = cexpf(0.5*I);

    for (ic=0; ic < nc; ic++) {
	x[ic][0] = ic+1.;
	for (i=1; i < n; i++) {
	    x[ic][i] = x[ic][i-1]*p;
	}
    }

    cburg_init (n,nc,nf);
    cburg_apply (x,a);
    cburg_close();

    printf("p=(%g,%g)\na={(%g,%g),(%g,%g)}\n",
	   crealf(p),cimagf(p),
	   crealf(a[0]),cimagf(a[0]),
	   crealf(a[1]),cimagf(a[1]));

    exit(0);
}
	   
	   

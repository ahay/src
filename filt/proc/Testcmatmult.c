#include <stdio.h>

#include <rsf.h>

#include "cmatmult.h"

int main (void) {
    const int n = 4;
    float complex **a, b[4], x[4];

    a = sf_complexalloc2(n,n);

    a[0][0] = 10.;
    a[0][1] = 1.+I;
    a[0][2] = 1.+2.*I;
    a[0][3] = 1.+3.*I;

    a[1][0] = 1.-I;
    a[1][1] = 9.;
    a[1][2] = 2.+I;
    a[1][3] = 2.+2.*I;

    a[2][0] = 1.-2.*I;
    a[2][1] = 2.-I;
    a[2][2] = 8.;
    a[2][3] = 3.+I;

    a[3][0] = 1.-3.*I;
    a[3][1] = 2.-2.*I;
    a[3][2] = 3.-I;
    a[3][3] = 7.;

    b[0] = 27.+27.*I;
    b[1] = 39.-9.*I;
    b[2] = 45.+9.*I;
    b[3] = 44.-32.*I;

    cmatmult_init (a);
    
    sf_csolver(cmatmult_lop, sf_ccgstep, n, n, x, b, 2*n, 
	      "verb", true, "end");

    fprintf(stderr,"inverse: (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n", 
	    crealf(x[0]),cimagf(x[0]),
	    crealf(x[1]),cimagf(x[1]),
	    crealf(x[2]),cimagf(x[2]),
	    crealf(x[3]),cimagf(x[3]));

    exit(0);
}

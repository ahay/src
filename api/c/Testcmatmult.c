#include <stdio.h>

#include "cmatmult.h"
#include "alloc.h"
#include "bigsolver.h"
#include "ccgstep.h"
#include "ccdstep.h"

int main (void) {
    const int n = 4;
    sf_complex **a, b[4], x[4];

    a = sf_complexalloc2(n,n);

    a[0][0] = sf_cmplx(10.,0.);
    a[0][1] = sf_cmplx(1.,1.);
    a[0][2] = sf_cmplx(1.,2.);
    a[0][3] = sf_cmplx(1.,3.);

    a[1][0] = sf_cmplx(1.,-1.);
    a[1][1] = sf_cmplx(9.,0.);
    a[1][2] = sf_cmplx(2.,1.);
    a[1][3] = sf_cmplx(2.,2.);

    a[2][0] = sf_cmplx(1.,-2.);
    a[2][1] = sf_cmplx(2.,-1.);
    a[2][2] = sf_cmplx(8.,0.);
    a[2][3] = sf_cmplx(3.,+1.);

    a[3][0] = sf_cmplx(1.,-3.);
    a[3][1] = sf_cmplx(2.,-2.);
    a[3][2] = sf_cmplx(3.,-1.);
    a[3][3] = sf_cmplx(7.,0.);

    b[0] = sf_cmplx(27.,+27.);
    b[1] = sf_cmplx(39.,-9.);
    b[2] = sf_cmplx(45.,+9.);
    b[3] = sf_cmplx(44.,-32.);

    sf_cmatmult_init (a);
    
    sf_csolver(sf_cmatmult_lop, sf_ccgstep, n, n, x, b, 2*n, 
	      "verb", true, "end");

    fprintf(stderr,"inverse: (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n", 
	    crealf(x[0]),cimagf(x[0]),
	    crealf(x[1]),cimagf(x[1]),
	    crealf(x[2]),cimagf(x[2]),
	    crealf(x[3]),cimagf(x[3]));

    sf_ccdstep_init();
    sf_csolver(sf_cmatmult_lop, sf_ccdstep, n, n, x, b, 2*n, 
	      "verb", true, "end");

    fprintf(stderr,"inverse: (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n", 
	    crealf(x[0]),cimagf(x[0]),
	    crealf(x[1]),cimagf(x[1]),
	    crealf(x[2]),cimagf(x[2]),
	    crealf(x[3]),cimagf(x[3]));

    exit(0);
}

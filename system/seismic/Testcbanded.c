#include <stdio.h>

#include "cbanded.h"

int main(void)
{
    float diag[]={4.,2.,3.,4.};
    sf_complex offd1[3];
    sf_complex offd2[2];
    sf_complex *offd[2];

    sf_complex x[4];

    x[0] = sf_cmplx(5.,-2.);
    x[1] = sf_cmplx(-2.,-3.);
    x[2] = sf_cmplx(7.,1.);
    x[3] = sf_cmplx(-3.,9.);

    offd1[0] = sf_cmplx(0.,-1.);
    offd1[1] = sf_cmplx(0.,0.);
    offd1[2] = sf_cmplx(1.,1.);

    offd2[0] = sf_cmplx(1.,0.);
    offd2[1] = sf_cmplx(-1.,-1.);
    
    cbanded_init (4,2);

    offd[0]=offd1;
    offd[1]=offd2;

    cbanded_define(diag,offd);
    cbanded_solve(x);
    cbanded_close();

    printf("(%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",
	   crealf(x[0]),cimagf(x[0]),
	   crealf(x[1]),cimagf(x[1]),
	   crealf(x[2]),cimagf(x[2]),
	   crealf(x[3]),cimagf(x[3]));
    /* (1,0),(-1,1),(2,-1),(-2,2) */

    return 0;
}
    
/* 	$Id$	 */


#include <stdio.h>

#include "ctridiagonal.h"

int main(void)
{
    ctris slv;
    sf_complex diag[4], offd[3], x[4];
    
    diag[0] = sf_cmplx(1.,0.);
    diag[1] = sf_cmplx(2.,-1.);
    diag[2] = sf_cmplx(3.,0.);
    diag[3] = sf_cmplx(4.,0.);

    offd[0] = sf_cmplx(1.,1.);
    offd[1] = sf_cmplx(0.,0.);
    offd[2] = sf_cmplx(1.,-2.);

    x[0] = sf_cmplx(2.,3.);
    x[1] = sf_cmplx(6.,1.);
    x[2] = sf_cmplx(15.,-7.);
    x[3] = sf_cmplx(19.,-2.);

    slv = ctridiagonal_init (4);

    ctridiagonal_define(slv,diag,offd);
    ctridiagonal_solve(slv,x);

    printf("(%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",
	   crealf(x[0]),cimagf(x[0]),
	   crealf(x[1]),cimagf(x[1]),
	   crealf(x[2]),cimagf(x[2]),
	   crealf(x[3]),cimagf(x[3]));

    
    ctridiagonal_const_define(slv,sf_cmplx(1.,2.),sf_cmplx(0.,1.));
    ctridiagonal_solve(slv,x);

    printf("(%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",
	   crealf(x[0]),cimagf(x[0]),
	   crealf(x[1]),cimagf(x[1]),
	   crealf(x[2]),cimagf(x[2]),
	   crealf(x[3]),cimagf(x[3]));

    ctridiagonal_close(slv);

    return 0;
}
    
/* 	$Id$	 */


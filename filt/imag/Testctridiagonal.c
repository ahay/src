#include <stdio.h>

#include "ctridiagonal.h"

int main(void)
{
    ctris slv;
    float complex diag[]={1.,2.-I,3.,4.};
    float complex offd[]={1.+I,0.,1.-2.*I};
    float complex x[]={2.+3.*I,6.+I,15.-7.*I,19-2.*I};
    
    slv = ctridiagonal_init (4);

    ctridiagonal_define(slv,diag,offd);
    ctridiagonal_solve(slv,x);
    ctridiagonal_close(slv);

    printf("(%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",
	   crealf(x[0]),cimagf(x[0]),
	   crealf(x[1]),cimagf(x[1]),
	   crealf(x[2]),cimagf(x[2]),
	   crealf(x[3]),cimagf(x[3]));

    return 0;
}
    
/* 	$Id: Testctridiagonal.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */


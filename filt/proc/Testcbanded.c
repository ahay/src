#include <stdio.h>

#include "cbanded.h"

int main(void)
{
    float diag[]={4.,2.,3.,4.};
    float complex offd1[3];
    float complex offd2[2];
    float complex *offd[2];

    float complex x[]={5.-2.*I,-2.-3*I,7.+I,-3.+9.*I}; 

    offd1[0] = -I;
    offd1[1] = 0.;
    offd1[2] = 1.+I;

    offd2[0] = 1.;
    offd2[1] = -1.-I;
    
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
    
/* 	$Id: Testcbanded.c,v 1.2 2004/04/19 22:03:22 fomels Exp $	 */


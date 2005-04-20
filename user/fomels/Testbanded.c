#include <stdio.h>

#include "banded.h"

int main(void)
{
    bands slv;
    float diag[]={1.,2.,3.,4.};
    float offd1[]={1.,0.,1.};
    float offd2[]={1.,-1.};
    float *offd[2];
    float x[]={2.,1.,5.,-5.};
    
    slv = banded_init (4,2);

    offd[0]=offd1;
    offd[1]=offd2;

    banded_define(slv,diag,offd);
    banded_solve(slv,x);
    banded_close(slv);

    printf("%g %g %g %g\n",x[0],x[1],x[2],x[3]);
    /* 1,-1,2,-2 */

    return 0;
}
    
/* 	$Id$	 */


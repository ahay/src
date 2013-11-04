#include <stdio.h>

#include "gaussel.h"

int main(void)
{
    float a0[4] = {1.0,  2.0, -1.0, -2.0};
    float a1[4] = {1.0,  3.0, -1.0, -2.0};
    float a2[4] = {2.0,  1.0,  1.0,  1.0};
    float a3[4] = {3.0,  1.0,  2.0,  1.0};
    float *a[4];
    float b[4] = {-6.0, -4.0, 11.0, 15.0};
    float x[4];
    
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    a[3] = a3;

    sf_gaussel_init(4);
    sf_gaussel_solve(a,b,x);

    printf("%g %g %g %g\n",x[0],x[1],x[2],x[3]);

    return 0;
}

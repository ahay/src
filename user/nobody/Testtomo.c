#include <stdio.h>

#include <rsf.h>

#include "tomo.h"

int main(void) {
    float dot1[2], dot2[2];
    static int nz=100, nx=100, np=11;
    static float dz=1., dx=2.;

    tomo_init(np, nz, nx, dz, dx);
    sf_dot_test(tomo_lop, nz*nx, np*nx, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

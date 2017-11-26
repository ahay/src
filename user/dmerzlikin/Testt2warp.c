#include <stdio.h>

#include <rsf.h>

#include "t2warp.h"

int main(void) {
    double dot1[2], dot2[2];
    static int n1=100, n2=100, n3=10; 
    static float o1=0, d1=0.004, o2=0, d2=0.001584, eps=0.01;

    t2warp_init(n1,n2,n3,o1,d1,o2,d2,eps);
    sf_dot_test(t2warp, n3*n1, n3*n2, dot1, dot2);
    t2warp_close();

    printf ("%12.3f ? %12.3f\n",dot1[0],dot1[1]);
    printf ("%12.3f ? %12.3f\n",dot2[0],dot2[1]);

    t2warp_init(n1,n2,n3,o1,d1,o2,d2,eps);
    sf_dot_test(t2warp_inv, n3*n2, n3*n1, dot1, dot2);
    t2warp_close();

    printf ("%12.3f ? %12.3f\n",dot1[0],dot1[1]);
    printf ("%12.3f ? %12.3f\n",dot2[0],dot2[1]);

    exit(0);
}

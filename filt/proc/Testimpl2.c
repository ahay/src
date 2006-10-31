#include <stdio.h>
#include <stdlib.h>

#include <rsf.h>

#include "impl2.h"

int main(void) {
    static int n1=100, n2=200;
    static float r1=5., r2=10., tau=1., pclip=50.;
    static bool up=false;
    int i1, i2;
    float dot1[2], dot2[2], **x;

    x = sf_floatalloc2(n1,n2);

    for (i2 = 0; i2 < n2; i2++) {
	for (i1 = 0; i1 < n1; i1++) {
	    x[i2][i1] = ((float) rand())/RAND_MAX;
	}
    }

    impl2_init(r1, r2, n1, n2, tau, pclip, up, false, NULL, 1, NULL);
    impl2_set(x);

    sf_dot_test(impl2_lop, n1*n2, n1*n2, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

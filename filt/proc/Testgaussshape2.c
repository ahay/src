#include <stdio.h>

#include "gaussshape2.h"
#include "freqfilt2.h"

int main(void) {
    float dot1[2], dot2[2], a[3]={4.,4.,1.};
    static int n1=100, n2=100;

    gaussshape2_init(n1, n2);
    gaussshape2_set2(a);
    sf_dot_test(freqfilt2_lop, n1*n2, n1*n2, dot1, dot2);
    gaussshape2_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

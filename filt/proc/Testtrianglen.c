#include <stdio.h>

#include <rsf.h>

#include "trianglen.h"

int main(void) {
    float dot1[2], dot2[2];
    int nbox[]={10,5}, ndat[]={50,100}, n12;

    trianglen_init(2, nbox, ndat);
    n12 = ndat[0]*ndat[1];

    sf_dot_test(trianglen_lop, n12, n12, dot1, dot2);
    trianglen_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

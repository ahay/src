#include <stdio.h>

#include <rsf.h>

#include "smspray.h"

int main(void) {
    double dot1[2], dot2[2];
    int n=20, ns=3;

    smspray_init(n,ns,'t');
    sf_dot_test(smspray_lop, n, n, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

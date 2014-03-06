#include <stdio.h>

#include <rsf.h>

#include "spray.h"

int main(void) {
    double dot1[2], dot2[2];
    int n=20, ns=3, ns2;

    ns2 = spray_init(ns);
    sf_dot_test(spray_lop, n, n*ns2, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

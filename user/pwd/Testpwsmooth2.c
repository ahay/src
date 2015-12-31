#include <stdio.h>

#include <rsf.h>

#include "pwsmooth2.h"

int main(void) {
    double dot1[2], dot2[2];
    int n2=10, ns=3, n1=10;
    int i1, i2;
    float **p1, **p2;

    p1 = sf_floatalloc2(n1,n2);
    p2 = sf_floatalloc2(n1,n2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    p1[i2][i1] = (i1-n1/2)*(i2-n2/2)/((n1-1.)*(n2-1.));
	    p2[i2][i1] = -(i1-n1/4)*(i2-3*n2/4)/((n1-1.)*(n2-1.));
	}
    }

    pwsmooth2_init(ns,n1,n2,1,0.1);
    pwsmooth2_set(p1,p2);
    sf_dot_test(pwsmooth2_lop, n1*n2, n1*n2, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

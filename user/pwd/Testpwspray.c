#include <stdio.h>

#include <rsf.h>

#include "pwspray.h"

int main(void) {
    double dot1[2], dot2[2];
    int n2=10, ns=3, ns2, n1=10;
    int i1, i2;
    float **p;

    p = sf_floatalloc2(n1,n2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    p[i2][i1] = (i1-n1/2)*(i2-n2/2)/((n1-1.)*(n2-1.));
	}
    }

    ns2 = pwspray_init(ns,n1,n2,1,0.1,p);
    sf_dot_test(pwspray_lop, n1*n2, n1*n2*ns2, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

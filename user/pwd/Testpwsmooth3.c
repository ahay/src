#include <stdio.h>

#include <rsf.h>

#include "pwsmooth3.h"

int main(void) {
    double dot1[2], dot2[2];
    int n2=10, ns1=3, ns2=2, n1=10, n3=3;
    int i1, i2, i3;
    float ****p;

    p = sf_floatalloc4(n1,n2,n3,2);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		p[0][i3][i2][i1] = (i1-n1/2)*(i2-n2/2)/((n1-1.)*(n2-1.));
		p[1][i3][i2][i1] = (i1-n1/4)*(i2-n2/3)/((n1-1.)*(n2-1.));
	    }
	}
    }

    pwsmooth3_init(ns1,ns2,n1,n2,n3,1,1,0.1,p);
    sf_dot_test(pwsmooth3_lop, n1*n2*n3, n1*n2*n3, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

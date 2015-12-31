#include <stdio.h>

#include <rsf.h>

#include "pwdsl.h"

int main(void) {
    double dot1[2], dot2[2];
    int n1=10, n2=10, nw=1, rect1=3, rect2=3; 
    int i1, i2;
    float **p;

    p = sf_floatalloc2(n1,n2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    p[i2][i1] = (i1-n1/2)*(i2-n2/2)/((n1-1.)*(n2-1.));
	}
    }

    pwdsl_init(n1, n2, nw, rect1, rect2, 0.1);
    pwdsl_set(p);
    sf_dot_test(pwdsl_lop, n1*n2, n1*n2, dot1, dot2);
    pwdsl_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

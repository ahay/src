#include <stdio.h>

#include <rsf.h>

#include "allp3.h"

int main(void) {
    double dot1[2], dot2[2];
    float *p, *q;
    bool drift = false;
    int n1=10, n2=10, n3=10, nw=1, nj=1; 
    int i, n;

    n = n1*n2*n3;
    p = sf_floatalloc(n);
    q = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	p[i]=0.3*i/n;
	q[i]=0.7;
    }

    allpass3_init (allpass_init(nw,nj,n1,n2,n3,drift,p),
		   allpass_init(nw,nj,n1,n2,n3,drift,q));

    sf_dot_test(allpass3_lop, n, 2*n, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

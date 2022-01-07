#include <stdio.h>

#include <rsf.h>

#include "sfchain.h"

int main(void) {
    double dot1[2], dot2[2];
    int n=16, nf, i;
    float *w, *wf, *x1, *x2, *s;

    w = sf_floatalloc(n);
    
    nf = 1+n/2;
    wf = sf_floatalloc(nf);

    x1 = sf_floatalloc(n);
    x2 = sf_floatalloc(n);
    s = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	w[i] = 1+0.1*i;
	s[i] = i*i;
	x1[i] = 0.01;
	x2[i] = 0.02;
    }

    for (i=0; i < nf; i++) {
	wf[i] = 1.+0.1*i*i;
    }
    
    sfchain_init(n,nf,w,wf,x1,x2,s);

    sf_dot_test(sfchain_lop, 3*n+nf, 3*n, dot1, dot2);
 

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    sf_dot_test(sfchainx_lop, 2*n, 3*n, dot1, dot2);
 
    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

   
    sfchain_close();
    
    exit(0);
}

#include <stdio.h>

#include <rsf.h>

#include "fchain.h"

int main(void) {
    double dot1[2], dot2[2];
    int n=16, nf, i;
    float *w, *wf, *inp;

    w = sf_floatalloc(n);
    
    nf = 1+n/2;
    wf = sf_floatalloc(nf);

    inp = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	w[i] = 1+0.1*i;
	inp[i] = i*i;
    }

    for (i=0; i < nf; i++) {
	wf[i] = 1.+0.1*i*i;
    }
    
    fchain_init(n,nf,w,wf,inp);

    sf_dot_test(fchain_lop, n+nf, n, dot1, dot2);
    
    fchain_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

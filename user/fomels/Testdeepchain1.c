#include <stdio.h>

#include <rsf.h>

#include "deepchain1.h"

int main(void) {

    double dot1[2], dot2[2];
    //int n=16, nk, i;
    int i, nt, nx, nw, ntx, nwx, nf;
    float *w, *wf, *x1, *x2, *x3, *x4, *s;

    nt = 6;
    nx = 4;    
    ntx =nt*nx;
    nw = kiss_fft_next_fast_size((nt+1)/2)+1;
    nwx = nw*nx;
    nf = 2*(nw-1);

    x1 = sf_floatalloc(ntx);
    x2 = sf_floatalloc(ntx);
    x3 = sf_floatalloc(ntx);
    x4 = sf_floatalloc(ntx);
    w = sf_floatalloc(ntx);
    wf = sf_floatalloc(nwx);
    s = sf_floatalloc(ntx);

    for (i=0; i < ntx; i++) {
    x1[i] = 0.01+i;
    x2[i] = 0.02+i;
    x3[i] = 0.03+i;
    x4[i] = 0.04+i;
    w[i] = 1.0;
    s[i] = 0.0+0.1*i;

    }

    for (i=0; i < nwx; i++) {
        wf[i] = 1.0;
     }
    
    printf("parameters: nt %d nw %d nf %d nx %d ntx %d nwx %d || nxx %d nyy %d \n",nt,nw,nf,nx,ntx,nwx,5*ntx+nwx,5*ntx);

    sfdeepchain1_init(nt, nx, nw, w, wf, x1, x2, x3, x4, s);

    sf_dot_test(sfdeepchain1_lop, 5*ntx+nwx, 5*ntx, dot1, dot2);
    
    sfdeepchain1_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

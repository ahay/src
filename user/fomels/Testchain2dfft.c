#include <stdio.h>

#include <rsf.h>

#include "chain2dfft.h"
#include "fft2.h"

int main(void) {
    double dot1[2], dot2[2];
    int n, nt, nx, nt1, nx2, nk, i, nw;
    float *w, *wf, *x1, *x2, *s;
    nt = 7;
    nx = 5;    
    nk = fft2_init(false, 1, nt, nx, &nt1, &nx2);
    n = nt*nx;
    w = sf_floatalloc(n);
    nw = nk/nx2;

    wf = sf_floatalloc(nk);

    x1 = sf_floatalloc(n);
    x2 = sf_floatalloc(n);
    s = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	w[i] = 1.0f;
	s[i] = 3.0f;
	x1[i] = 0.2f;
	x2[i] = 0.6f;

    }

    for (i=0; i < nk; i++) {
     wf[i] = 1.2;
     }
    
    printf("parameters: nt %d nx %d n %d nt1 %d nx2 %d nk %d nw %d || nxx %d nyy %d \n",nt,nx,n,nt1,nx2,nk,nw,3*nt*nx+nk,3*nt*nx);
    sfchain2d_init(nt, nx, nt1, nx2, nk, w, wf, x1, x2, s);
    sf_dot_test(sfchain2d_lop, 3*n+nk, 3*n, dot1, dot2);
    
    sfchain2d_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

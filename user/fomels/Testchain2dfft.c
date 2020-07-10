#include <stdio.h>

#include <rsf.h>

//#include "sfchain.h"
#include "chain2dfft.h"
#include "fft2.h"

int main(void) {
    double dot1[2], dot2[2];
    //int n=16, nk, i;
    int n, nt, nx, nt1, nx2, nk, i;
    float *w, *wf, *x1, *x2, *s;
    nt = 4;
    nx = 4;    
    nk = fft2_init(false, 1, nt, nx, &nt1, &nx2);
    n = nt*nx;
    w = sf_floatalloc(n);
    //sf_warning("nk is %d", nk);

    //nk = 1+n/2;
    wf = sf_floatalloc(nk);

    x1 = sf_floatalloc(n);
    x2 = sf_floatalloc(n);
    s = sf_floatalloc(n);

    for (i=0; i < n; i++) {
/*	w[i] = 1.0f;
	s[i] = 1.0f;
	x1[i] = 1.0f;
	x2[i] = 1.0f;*/
    w[i] = 1+0.1*i;
    s[i] = i*i+1.0;
    x1[i] = 0.01;
    x2[i] = 0.02;
    }

    for (i=0; i < nk; i++) {
/*    wf[i] = 2.0f;
*/    wf[i] = i*i;

     }
    
    printf("***input***\n");
    printf("w is: %g\n",w[0]);
    printf("s is: %g\n",s[0]);
    printf("x1 is: %g\n",x1[0]);
    printf("x2 is: %g\n",x2[0]);
    printf("wf is: %g \n",wf[0]);
    printf("parameters: nt %d nx %d n %d nt1 %d nx2 %d nk %d || nxx %d nyy %d \n",nt,nx,n,nt1,nx2,nk,3*nt*nx+nk,3*nt*nx);

    //sfchain_init(n,nf,w,wf,x1,x2,s);
    sfchain2d_init(nt, nx, nt1, nx2, nk, w, wf, x1, x2, s);
    sf_dot_test(sfchain2d_lop, 3*n+nk, 3*n, dot1, dot2);
    
    sfchain2d_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

#include <stdio.h>

#include <rsf.h>

#include "veltran.h"

int main(void) {
    double dot1[2], dot2[2];
    static bool pull=false;
    static int nx=10, ns=20, nt=100; 
    static float x0=0, dx=0.01, s0=0, ds=0.001, t0=0, dt=0.001, s1=0;

    veltran_init(pull,x0,dx,nx,s0,ds,ns,t0,dt,nt,s1,1.,0,0);
    sf_dot_test(veltran_lop, nt*ns, nt*nx, dot1, dot2);
    veltran_close();

    printf ("%12.3f ? %12.3f\n",dot1[0],dot1[1]);
    printf ("%12.3f ? %12.3f\n",dot2[0],dot2[1]);

    exit(0);
}

#include <stdio.h>

#include <rsf.h>

#include "lr2zortm.h"
#include "fft2.h"

int main(void) {
    double dot1[2], dot2[2];

    int nt, nx, nz, nzpad, nxpad, nk, pad = 1;
    int nzx, ntx, m2;
    int i, im, ix, iz, ik;
    float **lft, **rht;


    nt = 9;
    nz = 5;
    nx = 4;    
    m2 = 3;
    nzx = nz*nx;
    ntx = nt*nx;
    nk = fft2_init(false, pad, nx, nz, &nxpad, &nzpad);

    lft = sf_floatalloc2(nzx,m2);
    rht = sf_floatalloc2(m2,nk);


    for (im = 0; im < m2; im++) {
    for (ix = 0; ix < nx; ix++) {
    for (iz = 0; iz < nz; iz++) {
        i = ix+iz*nx;
        lft[im][i] = (i+1.0)/100.0;
        }
        }
    }

    for (ik = 0; ik < nk; ik++) {
    for (im = 0; im < m2; im++) {
        rht[ik][im] = (i+1.0)/50.0;
        }
    }


    lr2zortm_init(nt, nx, nz, nk, nxpad, nzpad, m2, lft, rht);

    sf_dot_test(lr2zortm_lop, nzx, ntx, dot1, dot2);
    
    printf("finish computing ...\n");

    lr2zortm_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

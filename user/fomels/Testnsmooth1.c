#include <stdio.h>

#include <rsf.h>

#include "nsmooth1.h"

int main(void) {
    double dot1[2], dot2[2];
    int ndat1=50, ndat2=10, ndat, i1, i2;
    float **rect;

    rect = sf_floatalloc2(ndat1,ndat2);

    for (i2=0; i2 < ndat2; i2++) {
	for (i1=0; i1 < ndat1; i1++) {
	    rect[i2][i1] = 1+0.3*((i1+i2)%4);
	}
    }

    nsmooth1_init(ndat1, ndat2, rect);

    ndat = ndat1*ndat2;
    sf_dot_test(nsmooth1_lop, ndat, ndat, dot1, dot2);
    nsmooth1_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

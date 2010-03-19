#include <stdio.h>

#include <rsf.h>

#include "ntriangle1.h"

int main(void) {
    double dot1[2], dot2[2];
    int nbox=5, ndat1=50, ndat2=10, ndat, **rect, **shift, i1, i2;

    rect = sf_intalloc2(ndat1,ndat2);
    shift = sf_intalloc2(ndat1,ndat2);

    for (i2=0; i2 < ndat2; i2++) {
	for (i1=0; i1 < ndat1; i1++) {
	    rect[i2][i1] = 1+(i1+i2)%(nbox-1);
	    shift[i2][i1] = i1%2;
	}
    }

    ntriangle1_init(nbox, ndat1, ndat2, rect, shift);

    ndat = ndat1*ndat2;
    sf_dot_test(ntriangle1_lop, ndat, ndat, dot1, dot2);
    ntriangle1_close();

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

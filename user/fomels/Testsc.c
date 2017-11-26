#include <stdio.h>

#include <rsf.h>

#include "sc.h"

int main(void) {
    double dot1[2], dot2[2];
    int nd=4, nm=2, **indx, *size;

    indx = sf_intalloc2(nm,nd);
    indx[0][0] = 0;
    indx[0][1] = 0;
    indx[1][0] = 0;
    indx[1][1] = 1;
    indx[2][0] = 1;
    indx[2][1] = 0;
    indx[3][0] = 1;
    indx[3][1] = 1;

    size = sf_intalloc(nm);
    size[0] = 2;
    size[1] = 2;
    
    sc_init(nm,indx,size);
    sf_dot_test(sc_lop, 2*nm, nd, dot1, dot2);

    printf ("%12.8f ? %12.8f\n",dot1[0],dot1[1]);
    printf ("%12.8f ? %12.8f\n",dot2[0],dot2[1]);

    exit(0);
}

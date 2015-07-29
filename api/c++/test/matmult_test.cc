#include <stdio.h>

#include <rsf.hh>

#include "matmult_wrap.hh"

int main (void) {
    float** a;
    float x[4], y[5], res[5];
    int i,j, iter;

    a = sf_floatalloc2(4,5);

    a[0][0] = 1.; a[0][1] = 1.;	a[0][2] = 1.; a[0][3] = 0.;
    a[1][0] = 1.; a[1][1] = 2.;	a[1][2] = 0.; a[1][3] = 0.;
    a[2][0] = 1.; a[2][1] = 3.;	a[2][2] = 1.; a[2][3] = 0.;
    a[3][0] = 1.; a[3][1] = 4.;	a[3][2] = 0.; a[3][3] = 1.;
    a[4][0] = 1.; a[4][1] = 5.;	a[4][2] = 1.; a[4][3] = 1.;

    y[0]=3.; y[1]=3.; y[2]=5.; y[3]=7.; y[4]=9.;

    printf ("y = \n");
    for (i=0; i < 5; i ++) {
	printf (" %10.2f",y[i]);
    }
    printf ("\n");
    printf ("a = \n");
    for (j=0; j < 4; j ++) {
	for (i=0; i < 5; i ++) {
	    printf (" %10.2f",a[j][i]);
	}
	printf("\n");
    }
    printf("\n");

    matmult_init(4,5,a); 

    printf ("cgstep\n------\n");
    for (iter =0; iter < 6; iter++) {
	sf_solver( matmult_lop, sf_cgstep, 4, 5, x, y, iter, 
		   "res", res, "end");
	sf_cgstep_close();
	printf ("x = ");
	for (i=0; i < 4; i ++) {
	    printf (" %12.8f",x[i]);
	}
	printf ("\n");
	printf ("res = ");
	for (i=0; i < 5; i ++) {
	    printf (" %12.8f",res[i]);
	}
	printf ("\n");
    }

    exit(0);
}

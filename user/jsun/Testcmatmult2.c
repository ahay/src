#include <stdio.h>

#include <rsf.h>

#include "cmatmult2.h"
#include "cgmres.h"

int main (void) {
    sf_complex** a;
    sf_complex x[4], y[4];
    int i,j, iter;

    a = sf_complexalloc2(4,4);

    /*
    a[0][0] = sf_cmplx(1.,0.); a[0][1] = sf_cmplx(1.,0.); a[0][2] = sf_cmplx(1.,0.); a[0][3] = sf_cmplx(0.,0.);
    a[1][0] = sf_cmplx(1.,0.); a[1][1] = sf_cmplx(2.,0.); a[1][2] = sf_cmplx(0.,0.); a[1][3] = sf_cmplx(0.,0.);
    a[2][0] = sf_cmplx(1.,0.); a[2][1] = sf_cmplx(3.,0.); a[2][2] = sf_cmplx(1.,0.); a[2][3] = sf_cmplx(0.,0.);
    a[3][0] = sf_cmplx(1.,0.); a[3][1] = sf_cmplx(4.,0.); a[3][2] = sf_cmplx(0.,0.); a[3][3] = sf_cmplx(1.,0.);
    x[0]=sf_cmplx(1.,0.); x[1]=sf_cmplx(1.,0.); x[2]=sf_cmplx(1.,0.); x[3]=sf_cmplx(2.,0.);
    y[0]=sf_cmplx(3.,0.); y[1]=sf_cmplx(3.,0.); y[2]=sf_cmplx(5.,0.); y[3]=sf_cmplx(7.,0.);
    */

    a[0][0] = sf_cmplx(1.,2.); a[0][1] = sf_cmplx(1.,0.); a[0][2] = sf_cmplx(1.,0.); a[0][3] = sf_cmplx(0.,0.);
    a[1][0] = sf_cmplx(-1.,3.); a[1][1] = sf_cmplx(2.,0.); a[1][2] = sf_cmplx(2.,-7.); a[1][3] = sf_cmplx(0.,0.);
    a[2][0] = sf_cmplx(1.,0.); a[2][1] = sf_cmplx(3.,-4.); a[2][2] = sf_cmplx(1.,0.); a[2][3] = sf_cmplx(0.,0.);
    a[3][0] = sf_cmplx(5.,2.); a[3][1] = sf_cmplx(4.,-2.); a[3][2] = sf_cmplx(0.,0.); a[3][3] = sf_cmplx(1.,0.);

    x[0]=sf_cmplx(1.,1.); x[1]=sf_cmplx(2.,3.); x[2]=sf_cmplx(2.,1.); x[3]=sf_cmplx(-2.,2.);

    printf ("x = \n");
    for (i=0; i < 4; i ++) {
        printf ("(%6.2f, %6.2f)",crealf(x[i]),cimagf(x[i]));
    }
    printf ("\n");

    printf ("a = \n");
    for (j=0; j < 4; j ++) {
	for (i=0; i < 4; i ++) {
	    printf ("(%6.2f, %6.2f)",crealf(a[j][i]),cimagf(a[j][i]));
	}
	printf("\n");
    }
    printf("\n");

    //    cmatmult2_init(a); 

    cmatmult2(4,x,y,a);

    printf ("y = Ax = \n");
    for (i=0; i < 4; i ++) {
        printf ("(%6.2f, %6.2f)",crealf(y[i]),cimagf(y[i]));
    }
    printf ("\n");

    printf ("cgmres\n------\n");
    cgmres_init(4,4);

    for (iter =0; iter < 5; iter++) {
	for (i=0; i < 4; i ++) {
	    x[i] = sf_cmplx(0.0f,0.0f);
	}
	cgmres( y, x, cmatmult2, a, iter, 0.01*SF_EPS, false);
	printf ("x = ");
	for (i=0; i < 4; i ++) {
	    printf ("(%12.8f, %12.8f)",crealf(x[i]),cimagf(x[i]));
	}
	printf ("\n");
    }    

    exit(0);
}

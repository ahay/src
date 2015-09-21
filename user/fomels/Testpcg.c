#include <stdio.h>
#include <rsf.h>
#include <rsf.h>

#include "pcg.h"

static float tmp[5], **a, *w;

void matmult_lop (bool adj, bool add, 
		  int nx, int ny, float* x, float*y) 
/*< linear operator >*/
{
    int ix, iy;
    sf_adjnull (adj, add, nx, ny, x, y);
    for (ix = 0; ix < nx; ix++) {
	for (iy = 0; iy < ny; iy++) {
	    if (adj) x[ix] += a[iy][ix] * y[iy];
	    else     y[iy] += a[iy][ix] * x[ix];
	}
    }
}


static void normal(int n, const float *inp, float *out)
{
    matmult_lop(false,false,4,5,(float*) inp,tmp);
    matmult_lop(true,false,4,5,out,tmp);
}

static void weight(int n, const float *inp, float *out)
{
    int i;

    for (i=0; i < n; i++) {
	out[i] = w[i]*inp[i];
    }
}

int main (void) 
{
    float x[4], y[5], d[4];
    int i,j;

    a = sf_floatalloc2(4,5);
    w = sf_floatalloc(4);
    sf_randn (4,w);

    printf ("w = ");
    for (i=0; i < 4; i ++) {
	printf (" %12.8f",w[i]);
    }
    printf ("\n");

    a[0][0] = 1.; a[0][1] = 1.;	a[0][2] = 1.; a[0][3] = 0.;
    a[1][0] = 1.; a[1][1] = 2.;	a[1][2] = 0.; a[1][3] = 0.;
    a[2][0] = 1.; a[2][1] = 3.;	a[2][2] = 1.; a[2][3] = 0.;
    a[3][0] = 1.; a[3][1] = 4.;	a[3][2] = 0.; a[3][3] = 1.;
    a[4][0] = 1.; a[4][1] = 5.; a[4][2] = 1.; a[4][3] = 1.;

    y[0]=3.; y[1]=3.; y[2]=5.; y[3]=7.; y[4] = 9.;

    printf ("y = \n");
    for (i=0; i < 5; i ++) {
	printf (" %10.2f",y[i]);
    }
    printf ("\n");
    printf ("a = \n");
    for (j=0; j < 5; j ++) {
	for (i=0; i < 4; i ++) {
	    printf (" %10.2f",a[j][i]);
	}
	printf("\n");
    }
    printf("\n");

    matmult_lop(true,false,4,5,d,y);

    conjgrad(normal, NULL, 4, d, NULL, x, 6, 0.0);

    printf ("x = ");
    for (i=0; i < 4; i ++) {
	printf (" %12.8f",x[i]);
    }
    printf ("\n");

    conjgrad(normal, weight, 4, d, NULL, x, 6, 0.0);

    printf ("x = ");
    for (i=0; i < 4; i ++) {
	printf (" %12.8f",x[i]);
    }
    printf ("\n");

    exit(0);
}

#include "edge.h"

/* a = 0.3*sqrt(3.); b = 0.05/sqrt(3.) */
static const float a=0.51961524227066318806, b=0.028867513459481288225;

void grad2 (int n, const float *x, float *w)
{
    int i;
    float ww;

    w[0] = 0.;
    for (i=1; i < n-1; i++) {
	ww = 0.5*(x[i+1]-x[i-1]);
	w[i] = ww*ww;
    }
}

void grad3 (int n1, int n2, float **x, float **w1, float **w2)
{
    int i1, i2;

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (i2 == 0 || i2 == n2-1 || i1 == 0 || i1 == n1-1) {
		w1[i2][i1] = 0.;
		w2[i2][i1] = 0.;
	    } else {
		w1[i2][i1] =
		    b*(x[i2-1][i1-1] + x[i2+1][i1-1] + 
		       x[i2+1][i1+1] + x[i2-1][i1+1]) +
		    a*(x[i2-1][i1] + x[i2+1][i1] - 2.*x[i2][i1]) +
		    (0.5-2.*b)*x[i2][i1+1] - (0.5+2.*b)*x[i2][i1-1];
		w2[i2][i1] =
		    b*(x[i2-1][i1-1] + x[i2+1][i1-1] +
		       x[i2+1][i1+1] + x[i2-1][i1+1]) +
		    a*(x[i2][i1-1] + x[i2][i1+1] - 2.*x[i2][i1]) +
		    (0.5-2.*b)*x[i2+1][i1] - (0.5+2.*b)*x[i2-1][i1];
	    }
	}
    }
}

void grad9 (int n1, int n2, float **x, float **w)
{
    int i1, i2;
    float w1, w2;

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (i2 == 0 || i2 == n2-1 || i1 == 0 || i1 == n1-1) {
		w[i2][i1] = 0.;
	    } else {
		w1 = b*(x[i2-1][i1-1] + x[i2+1][i1-1] + 
			x[i2+1][i1+1] + x[i2-1][i1+1]) +
		    a*(x[i2-1][i1] + x[i2+1][i1] - 2.*x[i2][i1]) +
		    (0.5-2.*b)*x[i2][i1+1] - (0.5+2.*b)*x[i2][i1-1];
		w2 = b*(x[i2-1][i1-1] + x[i2+1][i1-1] +
			x[i2+1][i1+1] + x[i2-1][i1+1]) +
		    a*(x[i2][i1-1] + x[i2][i1+1] - 2.*x[i2][i1]) +
		    (0.5-2.*b)*x[i2+1][i1] - (0.5+2.*b)*x[i2-1][i1];
		w[i2][i1] = w1*w1 + w2*w2;
	    }
	}
    }
}



  




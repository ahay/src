#include <rsf.h>

#include "cburg.h"

static  int n, nf, nc;
float complex **f, **b;

void cburg_init (int n_in, int nc_in, int nf_in)
{
    n = n_in; 
    nf = nf_in; 
    nc = nc_in;

    f = sf_complexalloc2(n,nc);
    b = sf_complexalloc2(n,nc);
}

void cburg_close(void)
{
    free (*f);
    free (f);
    free (*b);
    free (b);
}

void cburg_apply (float complex **x, float complex *a)
{
    double complex cj, num; 
    double den;
    float complex fi, bi, ai;
    int j, ic, i;

    for (ic=0; ic < nc; ic++) {
	for (i=0; i < n; i++) {
	    b[ic][i] = f[ic][i] = x[ic][i];
	}
    }

    a[0] = 1.;
    for (j=1; j < nf; j++) {
	num = den = 0.;
	for (ic=0; ic < nc; ic++) {
	    for (i=j; i < n; i++) {
		fi = f[ic][i];
		bi = b[ic][i-j];
		num += fi*conj(bi);
		den += creal(fi*conj(fi) + bi*conj(bi));
	    }
	}
	cj = 2.*num/den;
	for (ic=0; ic < nc; ic++) {
	    for (i=j; i < n; i++) {
		fi = f[ic][i];
		bi = b[ic][i-j];
		f[ic][i] -= cj*bi;
		b[ic][i-j] -= conj(cj)*fi;
	    }
	}
	for (i=1; i <= j/2; i++) {
	    ai = a[j-i]-cj*conjf(a[i]);
	    a[i] -= cj*conjf(a[j-i]);
	    a[j-i] = ai;
	}
	a[j] = -cj;
    }
}

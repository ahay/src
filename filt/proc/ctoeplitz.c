#include <rsf.h>

#include "ctoeplitz.h"

static int n;
static float complex *a;

static float complex cdprod (int j, 
			     const float complex* a, const float complex* b) 
{
    int i;
    float complex c;
    c = 0.;
    for (i=1; i <= j; i++) {
	c += a[j-i]*conjf(b[i]);
    }
    return c;
}

void ctoeplitz_init (int n_in)
{
    n = n_in;
    a = sf_complexalloc (n);
    a[0] = 1.;
}

  
/* r - top row of toeplitz matrix */
/* f - right-hand side overwritten with solution */
void ctoeplitz_solve (const float complex *r, float complex *f)
{    
    int i,j;
    float complex e,c,w, bot;
    float v;
    
    v=crealf(r[0]);
    f[0] /= v;
    
    for (j=1; j < n; j++) {
	e = cdprod(j,a,r);
	c = -e/v;

	v += crealf(c)*crealf(e) + cimagf(c)*cimagf(e);
       
	for (i=1; i <= j/2; i++) {
	    bot  = a[j-i] + c*conjf(a[i]);
	    a[i] = a[i] + c*conjf(a[j-i]);
	    a[j-i] = bot;
	}
	a[j] = c;
       
	w = cdprod(j,f,r);
	c = (f[j]-w)/v;
       
	for (i=0; i < j; i++) {
	    f[i] += c*conjf(a[j-i]);
	}
	f[j] = c;
    }
}
  
void ctoeplitz_close()
{
    free (a);
}











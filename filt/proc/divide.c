#include <rsf.h>

#include "divide.h"
#include "gauss2.h"
#include "weight.h"

static int n, niter;
static float *p;

void divide_init(int n1, int n2, float f1, float f2, int niter1)    
{
    n = n1*n2;
    niter = niter1;

    gauss2_init(n1,n2,f1,f2);
    sf_conjgrad_init(n, n, n, 1., 1.e-8, true, false);
    p = sf_floatalloc (n);
}

void divide_close (void)
{
    gauss2_close();
    sf_conjgrad_close();
    free (p);
}

void divide (const float* num, float* den,  float* rat)
{
    weight_init(den);
    sf_conjgrad(weight_lop,gauss2_lop,p,rat,num,niter);
}

/* 	$Id: divide.c,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */

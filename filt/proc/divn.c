#include <rsf.h>

#include "divn.h"
#include "trianglen.h"
#include "weight.h"

static int niter;
static float *p;

void divn_init(int ndim, int nd, int *ndat, int *nbox, int niter1) 
{
    niter = niter1;

    trianglen_init(ndim, nbox, ndat);
    sf_conjgrad_init(nd, nd, nd, nd, 1., 1.e-6, true, false);
    p = sf_floatalloc (nd);
}

void divn_close (void)
{
    trianglen_close();
    sf_conjgrad_close();
    free (p);
}

void divn (float* num, float* den,  float* rat)
{
    weight_init(den);
    sf_conjgrad(NULL, weight_lop,trianglen_lop,p,rat,num,niter); 
}

/* 	$Id: divn.c,v 1.1 2004/04/19 21:55:10 fomels Exp $	 */

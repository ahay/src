#include <rsf.h>

#include "twodiv2.h"
#include "gauss2.h"
#include "freqfilt2.h"
#include "triangle2.h"
#include "repeat.h"
#include "weight2.h"

static int n, niter;
static float *p;
static bool gauss;

void twodiv2_init(int nw, int n1, int n2, float f1, float f2, int niter1, 
		  bool gauss1, float* den) 
{
    n = n1*n2;
    niter = niter1;
    gauss = gauss1;

    if (gauss) {
	gauss2_init(n1,n2,f1,f2);
	repeat_init(n,nw,freqfilt2_lop);
    } else {
	triangle2_init((int) f1, (int) f2, n1, n2);
	repeat_init(n,nw,triangle2_lop);
    }
    sf_conjgrad_init(nw*n, nw*n, n, n, 1., 1.e-6, true, false);
    p = sf_floatalloc (nw*n);
    weight2_init(nw,n,den);
}

void twodiv2_close (void)
{
    if (gauss) {
	gauss2_close();
    } else { 
	triangle2_close();
    }
    sf_conjgrad_close();
    free (p);
    weight2_close();
}

void twodiv2 (float* num, float* rat)
{
    sf_conjgrad(NULL,weight2_lop,repeat_lop,p,rat,num,niter);
}

/* 	$Id$	 */

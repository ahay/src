#include <rsf.h>

#include "div1.h"
#include "gauss.h"
#include "triangle1.h"
#include "weight.h"

static int niter;
static float *p;
static bool gauss;

void div1_init(int n1, float f1, int niter1, bool gauss1) 
{
    niter = niter1;
    gauss = gauss1;

    if (gauss) {
	gauss_init(n1,f1);
    } else {
	triangle1_init((int) f1, n1);
    }
    sf_conjgrad_init(n1, n1, n1, n1, 1., 1.e-6, true, false);
    p = sf_floatalloc (n1);
}

void div1_close (void)
{
    if (gauss) {
	gauss_close();
    } else { 
	triangle1_close();
    }
    sf_conjgrad_close();
    free (p);
}

void div1 (float* num, float* den,  float* rat)
{
    weight_init(den);
    if (gauss) {
	sf_conjgrad(NULL, weight_lop,sf_freqfilt_lop,p,rat,num,niter);
    } else {
	sf_conjgrad(NULL, weight_lop,triangle1_lop,p,rat,num,niter); 
    }
}

/* 	$Id$	 */

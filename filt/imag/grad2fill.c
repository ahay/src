#include <rsf.h>

#include "grad2fill.h"

#include "igrad2.h"

static int n12;
static float* zero;

void grad2fill_init (int m1, int m2)
{
    int i;

    n12 = m1*m2;

    zero = sf_floatalloc (2*n12);

    for (i=0; i < 2*n12; i++) {
	zero[i] = 0.;
    }

    igrad2_init (m1,m2);
}

void grad2fill_close (void)
{
    free (zero);
}

void grad2fill(int niter, float* mm, bool *known)
{
    sf_solver (igrad2_lop, sf_cgstep, n12, 2*n12, mm, zero, niter, 
	       "x0", mm, "known", known, "end");
    sf_cgstep_close ();
}



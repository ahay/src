#include <rsf.h>

#include "lapfill.h"

#include "laplac2.h"
#include "igrad2.h"

static int n12;
static float* zero;
static bool grad;

void lapfill_init (int m1, int m2, bool grad1)
{
    int i;
    int n;

    n12 = m1*m2;
    grad = grad1;

    n = grad? 2*n12: n12;

    zero = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	zero[i] = 0.;
    }

    if (grad) {
	igrad2_init (m1,m2);
    } else {
	laplac2_init (m1,m2);
    }
}

void lapfill_close (void)
{
    free (zero);
}

void lapfill(int niter, float* mm, bool *known)
{
    if (grad) {
	sf_solver (igrad2_lop, sf_cgstep, n12, 2*n12, mm, zero, niter, 
		   "x0", mm, "known", known, "end");
    } else {
	sf_solver (laplac2_lop, sf_cgstep, n12, n12, mm, zero, niter, 
		   "x0", mm, "known", known, "end");
    }
    sf_cgstep_close ();
}

#include <stdlib.h>

#include "dottest.h"
#include "bigsolver.h"
#include "alloc.h"

static float dotprod (int n, const float* x, const float* y);
static void randvec (int n, float* x);

/* dot_test
   --------
   The dot product test to see if the adjoint is coded correctly.
   In the output dot1[0] shpould be equal to dot1[1] 
   (within machine precision),
   and dot2[0] should be equal to dot2[1]. */
void sf_dot_test(sf_operator oper, int nm, int nd, float* dot1, float* dot2) {
    float *mod1, *mod2, *dat1, *dat2;

    mod1 = sf_floatalloc (nm);
    mod2 = sf_floatalloc (nm);
    dat1 = sf_floatalloc (nd);
    dat2 = sf_floatalloc (nd);

    randvec( nm, mod1);
    randvec( nd, dat2);

    oper(false, false, nm, nd, mod1, dat1);
    dot1[0] = dotprod( nd, dat1, dat2);

    oper(true, false, nm, nd, mod2, dat2);
    dot1[1] = dotprod( nm, mod1, mod2);

    oper(false, true, nm, nd, mod1, dat1);
    dot2[0] = dotprod( nd, dat1, dat2);

    oper(true, true, nm, nd, mod2, dat2);
    dot2[1] = dotprod( nm, mod1, mod2);

    free (mod1);
    free (mod2);
    free (dat1);
    free (dat2);
}

/* randvec
   ------
   Fills an array with pseudo-random numbers */
static void randvec (int n, float* x) {
    int i;
    
    for (i = 0; i < n; i++) {
	x[i] = ((float) rand())/RAND_MAX;
    }
}

/* dotprod
   -------
   dot product of two vectors */
static float dotprod (int n, const float* x, const float* y) {
    float prod;
    int i;
    
    prod=0.;
    for (i = 0; i < n; i++) {
	prod += x[i]*y[i];
    }
    return prod;
}


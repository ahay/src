#include <rsf.h>

#include "levint.h"

#include "tcai2.h"
#include "peftc.h"
#include "lint1.h"

void levint1 (int niter, int nm, int nr, int nd, float *coord, 
	      const float *dd, float o1, float d1, float *rr, float eps) 
{
    float *aa;
    bool *ma;
    int ia, ir, na;

    ma = sf_boolalloc(nr);
    na = nr-nm;
    aa = rr+nm;

    lint1_init (o1, d1, coord);
    peftc_init (na, nm, aa, rr);
    tcai2_init (na, nm, aa);

    /* starting guess */
    for (ir=0; ir < nr; ir++) {
	ma[ir] = (ir >= nm);
	rr[ir] = 0.;
    }
    aa[0] = 1.;
    aa[1] = -2.;
    aa[2] = 1.;
    for (ia=3; ia < na; ia++) {
	aa[ia] = 0.;
    }

    sf_solver_reg (lint1_lop, sf_cgstep, peftc_lop, nr, nr, nd, rr, dd, niter,
		   eps, "nlreg", tcai2_lop, "x0", rr, "known", ma, "end");

    /* free filter coefficients */
    for (ir=nm+1; ir < nr; ir++) {
	ma[ir] = false;
    }
    
    sf_cgstep_close();

    sf_solver_reg (lint1_lop, sf_cgstep, peftc_lop, nr, nr, nd, rr, dd, niter,
		   eps, "nlreg", tcai2_lop, "x0", rr, "known", ma, "end");

    sf_cgstep_close();

    free(ma);
}

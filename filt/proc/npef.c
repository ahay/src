#include <rsf.h>

#include "npef.h"
#include "nhconest.h"
#include "npolydiv2.h"

void find_pef(int nd, float *dd, nfilter aa, nfilter rr, 
	      int niter, float eps, int nh) 
{
    int ip, ih, na, np, nr;
    float *flt;

    np = aa->np;
    nr = np*nh;
    flt = sf_floatalloc(nr);

    nhconest_init(dd, aa, nh);
    npolydiv2_init( nr, rr);

    sf_solver_prec(nhconest_lop, sf_cgstep, npolydiv2_lop,
		   nr, nr, nd, flt, dd, niter, eps, "end");
    sf_cgstep_close();
    npolydiv2_close();

    for (ip=0; ip < np; ip++) {
	na = aa->hlx[ip]->nh;
	for (ih=0; ih < na; ih++) {
	    aa->hlx[ip]->flt[ih] = flt[ip*nh + ih];
	}
    }

    free(flt);
}


#include <rsf.h>

#include "signoi.h"

#include "helicon.h"
#include "polydiv.h"

static int niter, nd;
static float eps, *dd;
static filter nn, ss;


void signoi_init(filter nn_in, filter ss_in, 
		 int niter_in, int nd_in, float eps_in)
{
    nn = nn_in;
    ss = ss_in;
    niter = niter_in;
    nd = nd_in;
    eps = eps_in;
    dd = sf_floatalloc(nd);
}

void signoi_lop (bool adj, bool add, int n1, int n2, float *data, float *sign)
{
    helicon_init (nn);
    polydiv_init (nd, ss); 

    if (nd != n1 || nd != n2) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,n1,n2,data,sign);

    helicon_lop (false, false, n1, n1, data, dd);
    sf_solver_prec(helicon_lop, sf_cgstep, polydiv_lop, nd, nd, nd, sign, dd,
		   niter, eps, "end");
    sf_cgstep_close();

    nn++;
    ss++;
}


void signoi_close(void)
{
    free(dd);
    polydiv_close();
}

#include <rsf.h>

#include "expsignoi2.h"
#include "expont.h"

static int n1,n2,n12;
static float eps, *a, *b, *c, *d, *tmp;

void expsignoi2_init (int m1,int m2, float eps1, float **aa)
{
    n1 = m1;
    n2 = m2;
    n12 = n1*n2;
    eps = eps1;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];

    tmp = sf_floatalloc(n12);
}

void expsignoi2_close(void)
{
    free (tmp);
}

void expsignoi2_lop (bool adj, bool add, int ns, int nd, 
		     float *sig, float *dat)
{
    int i;

    if (2*ns != nd) sf_error("%s: wrong size: 2*%d != %d",__FILE__,ns,nd);

    sf_adjnull(adj,add,ns,nd,sig,dat);

    expont_init (n1, n2, a, b);
    expont_lop (adj, true, ns, ns, sig, dat);
    expont_init (n1, n2, c, d);

    for (i=0; i < n12; i++) {
	tmp[i] = adj? dat[i+n12] * eps: sig[i] * eps;
    }

    if (adj) {
	expont_lop (true, true, ns, ns, sig, tmp);
    } else {
	expont_lop (false, true, ns, ns, tmp, dat+n12);
    }
}

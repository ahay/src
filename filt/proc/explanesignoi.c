#include <rsf.h>

#include "explanesignoi.h"
#include "expont.h"
#include "allp2.h"

static int n1,n2,n12;
static float eps, *a, *b, *c, *d, *tmp, *tmp2;
static allpass2 noi, sig;

void explanesignoi_init (int m1,int m2, float eps1, float **aa, 
			 int nw, int nj1, int nj2, float **nn, float **ss)
{
    n1 = m1;
    n2 = m2;
    n12 = n1*n2;
    eps = eps1;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];

    noi = allpass2_init(nw,nj1,n1,n2,nn);
    sig = allpass2_init(nw,nj2,n1,n2,nn);

    tmp = sf_floatalloc(n12);
    tmp2 = sf_floatalloc(n12);
}

void explanesignoi_close(void)
{
    free (tmp);
}

void explanesignoi_lop (bool adj, bool add, int ns, int nd, 
			float *ss, float *dat)
{
    int i;

    if (2*ns != nd) sf_error("%s: wrong size: 2*%d != %d",__FILE__,ns,nd);

    sf_adjnull(adj,add,ns,nd,ss,dat);

    allpass22_init(noi);
    expont_init (n1, n2, a, b);
    sf_chain(allpass21_lop, expont_lop, adj, true, ns, ns, ns, ss, dat, tmp);

    allpass22_init(sig);
    expont_init (n1, n2, c, d);

    for (i=0; i < n12; i++) {
	tmp2[i] = adj? dat[i+n12] * eps: ss[i] * eps;
    }

    if (adj) {
	sf_chain(allpass21_lop, expont_lop, 
		 true, true, ns, ns, ns, ss, tmp2, tmp);
    } else {
	sf_chain(allpass21_lop, expont_lop,
		 false, true, ns, ns, ns, tmp2, dat+n12, tmp);
    }
}

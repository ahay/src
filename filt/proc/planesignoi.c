#include "planesignoi.h"

#include "allp2.h"

static float eps, *tmp;
static allpass2 noi, sig;

void planesignoi_init (int nw, int nj1, int nj2, int nx, int ny, 
		       float **nn, float **ss, float eps1)
{
    eps = eps1;
    noi = allpass2_init(nw,nj1,nx,ny,nn);
    sig = allpass2_init(nw,nj2,nx,ny,ss);
    tmp = sf_floatalloc(nx*ny);
}

void planesignoi_lop (bool adj, bool add, int ns, int nd, float *s, float *d)
{
    int is;

    if (nd != 2*ns) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,ns,nd,s,d);

    allpass22_init(noi);
    allpass21_lop (adj, true, ns, ns, s, d);

    allpass22_init (sig);
    if (adj) {
	for (is=0; is < ns; is++) {
	    tmp[is] = eps*d[ns+is];
	}
	allpass21_lop (true, true, ns, ns, s, tmp);
    } else {
	for (is=0; is < ns; is++) {
	    tmp[is] = s[is];
	}
	allpass21_lop (false, true, ns, ns, tmp, d+ns);
    }
}

void planesignoi_close(void)
{
    free(noi);
    free(sig);
    free(tmp);
}

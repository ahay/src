#include <rsf.h>

#include "peftc.h"

#include "tcai1.h"
#include "tcaf1.h"

static int n;

void peftc_init (int na, int ny, float *aa, float *yy)
{
    n = ny;
    tcai1_init (na, aa);
    tcaf1_init (ny, yy);
}

void peftc_lop (bool adj, bool add, int nx, int nr, float *x, float *r)
{
    sf_adjnull(adj,add,nx,nr,x,r);

    tcai1_lop (adj, true, n,    nr, x,   r);
    tcaf1_lop (adj, true, nx-n, nr, x+n, r);
}

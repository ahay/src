#include <rsf.h>

#include "tcai2.h"
#include "tcai1.h"

static int ny;

void tcai2_init (int na, int nx, float *aa)
{
    ny = nx;
    tcai1_init (na,aa);
}

void tcai2_lop (bool adj, bool add, int nx, int nr, float *x, float *r)
{
    sf_adjnull(adj,add,nx,nr,x,r);
    tcai1_lop (adj,true,ny,nr,x,r);
}








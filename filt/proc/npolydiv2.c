#include <rsf.h>

#include "npolydiv2.h"
#include "npolydiv.h"

static float *tt;

void npolydiv2_init (int nd, nfilter aa)
{
    npolydiv_init (nd, aa);
    tt = sf_floatalloc(nd);
}

void npolydiv2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    if (nx != ny) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
        npolydiv_lop (false,false,nx,ny,yy,tt);
        npolydiv_lop (true,true,nx,ny,xx,tt);
    } else {
        npolydiv_lop (false,false,nx,ny,xx,tt);
        npolydiv_lop (true,true,nx,ny,yy,tt);
    }
}

void npolydiv2_close(void)
{
    npolydiv_close();
    free(tt);
}

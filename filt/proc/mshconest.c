#include <rsf.h>

#include "mshconest.h"
#include "mshelix.h"
#include "hconest.h"

static msfilter msaa;

void mshconest_init(float *x, msfilter msaa_in)
{
    msaa = msaa_in;
    hconest_init(x,msaa->one);
}

void mshconest_lop(bool adj, bool add, int na, int ny, float *a, float *y)
{
    int  is, nx;
    
    if (na != msaa->nh) sf_error("%s: Wrong data dimensions",__FILE__);
    nx = ny/msaa->ns;

    sf_adjnull(adj, add, na, ny, a, y);

    for (is=0; is < msaa->ns; is++) {
	onescale(is,msaa);
	hconest_lop(adj,true,na,nx,a,y+is*nx);
    }
}

/* 	$Id: mshconest.c,v 1.1 2004/06/11 10:51:33 fomels Exp $	 */

#include "nbound.h"
#include "nhelix.h"
#include "bound.h"

void nbound (int ip, int dim, int *nd, int *na, nfilter aa) 
{
    int i, n;
    filter bb;

    n=1;
    for (i=0; i < dim; i++) {
	n *= nd[i];
    }

    aa->mis = sf_boolalloc(n);
    bb = aa->hlx[ip];

    bound (dim, nd, nd, na, bb);

    for (i=0; i < n; i++) {
	aa->mis[i] = bb->mis[i];
    }

    free(bb->mis);
    bb->mis = NULL;
}


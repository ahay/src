#include <rsf.h>

#include "tcaf1.h"

static int nx;
static float *xx;

void tcaf1_init(int ny, float* yy)
{
    nx = ny;
    xx = yy;
}

void tcaf1_lop(bool adj, bool add, int nb, int ny, float *bb, float *yy)
{
    int x, b, y;

    if(ny < nx+nb-1) sf_error("%s: size problem: %d < %d+%d-1",
			      __FILE__,ny,nx,nb);
    sf_adjnull (adj, add, nb, ny, bb, yy);

    for (b=0; b < nb; b++) {
	for (x=0; x < nx; x++) {                  
	    y = x + b;

	    if( adj) bb[b] += yy[y] * xx[x];
	    else     yy[y] += bb[b] * xx[x];
        }
    }
}

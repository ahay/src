#include <assert.h> /* for assert  */
#include <rsf.h>    /* for adjnull */

static int nf;
static float* ff;

void conv_init (int na    /* filter length */, 
		 float* aa /* filter [na] */) 
/*< initialize >*/
{
    nf = na;
    ff = aa;
}

void conv_lop (bool adj, bool add, 
		int nx, int ny, float* xx, float* yy) 
/*< linear operator >*/
{
    int f, x, y;

    assert (ny == nx);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    for (f=0; f < nf; f++) {
	for (y = f-1; y < nx + f-1; y++) {
	    x = y - f + 1;
	    if( adj) {
		/* add code */
	    } else {
		yy[y] += xx[x] * ff[f];
	    }
	}
    }
}

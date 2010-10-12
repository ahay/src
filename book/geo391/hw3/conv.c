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
    int f, x, y, x0, x1;

    assert (ny == nx);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    for (f=0; f < nf; f++) {
	x0 = SF_MAX(0,1-f);
	x1 = SF_MIN(nx,nx+1-f);
	for (x = x0; x < x1; x++) {
	    if( adj) {
		/* add code */
	    } else {
		yy[x+f-1] += xx[x] * ff[f];
	    }
	}
    }
}

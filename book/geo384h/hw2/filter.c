#include <assert.h> /* for assert  */
#include <rsf.h>    /* for adjnull */

static float a, b, c;

void filter_init(float a1, float b1, float c1)
{
    a=a1;
    b=b1;
    c=c1;
}

void filter_lop (bool adj, bool add, 
		int nx, int ny, float* xx, float* yy) 
/*< linear operator >*/
{
    int i;
    float t;

    assert (ny == nx);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    if (adj) {
	/* !!! INSERT CODE !!! */
    } else {
	t = a*xx[0];
	yy[0] += t;
	for (i = 1; i < nx; i++) {
	    t = a*xx[i] + b*xx[i-1] + c*t;
	    yy[i] += t;
	}
    }
}

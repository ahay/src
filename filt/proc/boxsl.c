#include <math.h>

#include <rsf.h>

#include "boxsl.h"
#include "triangle.h"

static int n1, n2, rect;
static float **p, **tmp, *trace1, *trace2, amp;
static triangle tr;

void boxsl_init(int m1, int m2, int rect1, int rect2)
{
    n1 = m1;
    n2 = m2;
    rect = rect2;
    p = sf_floatalloc2(n1,n2+rect);

    tmp = sf_floatalloc2(n1,n2+rect);
    trace1 = sf_floatalloc(n1);
    trace2 = sf_floatalloc(n1);

    tr = triangle_init (rect1,n1);
    amp = 1./rect;
}

void boxsl_set(int m2, float** p1)
{
    int i1, i2;

    if (m2 >= n2+rect) {
	for (i2=0; i2 < n2+rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		p[i2][i1]=p1[i2][i1];
	    }
	}
    } else {
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		p[i2][i1]=p1[i2][i1];
	    }
	}
	for (i2=m2; i2 < n2+rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		p[i2][i1]=p1[m2-1][i1];
	    }
	}
    }
}

void boxsl_close(void)
{
    free(*p);
    free(p);
    free(*tmp);
    free(tmp);
    free(trace1);
    free(trace2);
}

static void roll(bool adj, float** in)
{
    int i1, i2, it;
    float t;

    if (adj) {
	/* roll from the future */
	for (i2=n2+rect-2; i2 >=0; i2--) {	    
	    for (i1=0; i1 < n1; i1++) {
		t = i1 + p[i2][i1];
		if (t < 0.) continue;
		it = t;
		t -= it;
		if (it >=0 && it < n1-1) 
		    in[i2][i1] += in[i2+1][it]*(1.-t) + in[i2+1][it+1]*t;
	    }	    
	}
    } else {
	/* unroll to the future */
	for (i2=0; i2 < n2+rect-1; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		t = i1 + p[i2][i1];
		if (t < 0.) continue;
		it = t;
		t -= it;
		if (it >=0 && it < n1-1) {
		    in[i2+1][it]   += in[i2][i1]*(1.-t);
		    in[i2+1][it+1] += in[i2][i1]*t;
		}
	    }
	}
    }
}

static void shifts(bool adj, float** in, float* out)
{
    int i1, i2, ir, it;
    float t;

    for (i2=0; i2 < n2+rect; i2++) {
	if (adj) {
	    for (i1=0; i1 < n1; i1++) {
		in[i2][i1]  = amp*out[i2*n1+i1];
	    }

	    if (i2 < n2) {
		for (i1=0; i1 < n1; i1++) {
		    trace2[i1] = -amp*out[(i2+rect)*n1+i1];
		}
		
		for (ir=i2+rect-1; ir >= i2; ir--) {
		    for (i1=0; i1 < n1; i1++) {
			trace1[i1] = trace2[i1];
			trace2[i1] = 0.;
		    }
		    for (i1=0; i1 < n1; i1++) {
			t = i1 + p[ir][i1];
			if (t < 0.) continue;
			it = t;
			t -= it;
			if (it >=0 && it < n1-1) 
			    trace2[i1] = trace1[it]*(1.-t) + trace1[it+1]*t;
		    }
		}

		for (i1=0; i1 < n1; i1++) {
		    in[i2][i1] += trace2[i1];
		}
	    }
	} else {
	    if (i2 < n2) {
		for (i1=0; i1 < n1; i1++) {
		    trace2[i1]=in[i2][i1];
		}
		
		for (ir=i2; ir < i2+rect; ir++) {
		    for (i1=0; i1 < n1; i1++) {
			trace1[i1] = trace2[i1];
			trace2[i1]=0.;
		    }
		    for (i1=0; i1 < n1; i1++) {
			t = i1 + p[ir][i1];
			if (t < 0.) continue;
			it = t;
			t -= it;
			if (it >=0 && it < n1-1) {
			    trace2[it]   += trace1[i1]*(1.-t);
			    trace2[it+1] += trace1[i1]*t;
			}
		    }
		}
		
		for (i1=0; i1 < n1; i1++) {
		    out[(i2+rect)*n1+i1] -= amp*trace2[i1];
		}
	    }

	    for (i1=0; i1 < n1; i1++) {
		out[i2*n1+i1] += amp*in[i2][i1];
	    }
	}
    }
}

void boxsl_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
{
    int i1, i2;

    if (nx != n1*n2 || ny != n1*(n2+rect)) 
	sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj, add, nx, ny, x, y);

    if (adj) {
	shifts(adj,tmp,y);
	roll(adj,tmp);
	for (i2=0; i2 < n2; i2++) {	    
	    smooth (tr,0,1,false,tmp[i2]);
	    for (i1=0; i1 < n1; i1++) {
		x[i2*n1+i1] += tmp[i2][i1];
	    }
	}
    } else {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tmp[i2][i1] = x[i2*n1+i1];
	    }
	    smooth (tr,0,1,false,tmp[i2]);
	}
	for (i2=n2; i2 < n2+rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tmp[i2][i1] = 0.;
	    }
	}
	roll(adj,tmp);
	shifts(adj,tmp,y);
    }
}

void trisl_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
{
    if (nx != n1*(n2+rect) || ny != nx) 
	sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj, add, nx, ny, x, y);
    
    boxsl_lop(true, false, n1*n2, ny, tmp[0], adj? y: x);
    boxsl_lop(false, true, n1*n2, ny, tmp[0], adj? x: y);
}


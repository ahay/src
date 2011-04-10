/* Convert warping function to gamma=Vp/Vs using LS and shaping regularization.

Takes: < warp.rsf [warpout=wout.rsf] > gamma.rsf
rect1= rect2= ...

rectN defines the size of the smoothing stencil in N-th dimension.
*/
#include <math.h>

#include <rsf.h>

#include "smoothder.h"

int main(int argc, char* argv[])
{
    int i, niter, nw, nd, dim, n1,n2, i1,i2;
    int n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float **w, **g, t0, dt;
    char key[6];
    sf_file warp, gamma, wout;

    sf_init(argc,argv);
    warp = sf_input("in");
    gamma = sf_output("out");
 
    if (!sf_histfloat(warp,"d1",&dt)) sf_warning("No d1= in input");
    if (!sf_histfloat(warp,"o1",&t0)) sf_warning("No o1= in input");

    dim = sf_filedims (warp,n);

    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
    }

    if (!sf_getint("nw",&nw)) nw=2*rect[0];
    /* trace extension */
    n[0] += 2*nw;

    nd = smoothder_init(dim, rect,n,false,false);
    n1 = n[0];
    n2 = nd/n1;

    w = sf_floatalloc2(n1,n2);
    g = sf_floatalloc2(n1,n2);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    sf_floatread(w[0],nd,warp);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    w[i2][i1] += (t0+i1*dt);
	    
	}
    }

    smoothder(niter, NULL, w[0], g[0]);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    g[i2][i1] = 2.*g[i2][i1]/dt-1.;
	    w[i2][i1] -= (t0+i1*dt);
	}
    }

    sf_floatwrite(g[0],nd,gamma);

    if (NULL != sf_getstring("warpout")) {
	/* optionally, output predicted warp */
	wout = sf_output("warpout");

	sf_floatwrite(w[0],nd,wout);
    }


    exit(0);
}

/* 	$Id$	 */

/* Convert RMS to interval velocity using LS and shaping regularization.

Takes: < vrms.rsf weight=weight.rsf [vrmsout=vout.rsf] > vint.rsf
rect1= rect2= ...

rectN defines the size of the smoothing stencil in N-th dimension.
*/
#include <math.h>

#include <rsf.h>

#include "smoothder.h"

int main(int argc, char* argv[])
{
    int i, niter, nd, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM], n1, n2, i1, i2;
    float **vr, **vi, **wt, wti;
    char key[6];
    sf_file vrms, vint, weight, vout;

    sf_init(argc,argv);
    vrms = sf_input("in");
    vint = sf_output("out");
    weight = sf_input("weight");

    dim = sf_filedims (vrms,n);

    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
    }

    nd = smoothder_init(dim, rect,n);
    n1 = n[0];
    n2 = nd/n1;
    
    vr = sf_floatalloc2(n1,n2);
    vi = sf_floatalloc2(n1,n2);
    wt = sf_floatalloc2(n1,n2);

    sf_read(vr[0],sizeof(float),nd,vrms);
    sf_read(wt[0],sizeof(float),nd,weight);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    wti = 0.;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    wti += wt[i2][i1]*wt[i2][i1];
	}
    }
    if (wti > 0.) wti = sqrtf(n1*n2/wti);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vr[i2][i1] *= vr[i2][i1]*(i1+1.);
	    wt[i2][i1] *= wti/(i1+1.); /* decrease weight with time */	    
	}
    }
    
    smoothder(niter, wt[0], vr[0], vi[0]);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vr[i2][i1] = sqrtf(vr[i2][i1]/(i1+1.));
	    vi[i2][i1] = sqrtf(vi[i2][i1]);
	}
    }

    sf_write(vi[0],sizeof(float),nd,vint);

    if (NULL != sf_getstring("vrmsout")) {
	/* optionally, output predicted vrms */
	vout = sf_output("vrmsout");

	sf_write(vr[0],sizeof(float),nd,vout);
    }

    sf_close();
    exit(0);
}

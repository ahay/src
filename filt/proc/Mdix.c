/* Convert RMS to interval velocity using LS and shaping regularization.

Takes: < vrms.rsf weight=weight.rsf [vrmsout=vout.rsf] > vint.rsf
rect1= rect2= ...

rectN defines the size of the smoothing stencil in N-th dimension.
*/

#include <rsf.h>

#include "dix.h"

int main(int argc, char* argv[])
{
    int i, niter, nd, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *vr, *vi, *wt;
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

    nd = dix_init(dim, rect,n);
    
    vr = sf_floatalloc(nd);
    vi = sf_floatalloc(nd);
    wt = sf_floatalloc(nd);

    sf_read(vr,sizeof(float),nd,vrms);
    sf_read(wt,sizeof(float),nd,weight);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    dix(niter, wt, vr, vi);

    sf_write(vi,sizeof(float),nd,vint);

    if (NULL != sf_getstring("vrmsout")) {
	/* optionally, output predicted vrms */
	vout = sf_output("vrmsout");

	sf_write(vr,sizeof(float),nd,vout);
    }

    sf_close();
    exit(0);
}

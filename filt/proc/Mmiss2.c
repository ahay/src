/* 2-D missing data interpolation.

Takes: < input.rsf > interpolated.rsf
*/

#include <rsf.h>

#include "mask.h"
#include "gaussshape2.h"
#include "freqfilt2.h"

int main(int argc, char* argv[])
{
    int niter, n1, n2, n3, i3, n12, i;
    float *mm, *kk, *pp, filt1, filt2, a[3], eps;
    bool *known, shape, force;
    sf_file in, out, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    mm = sf_floatalloc(n12);
    kk = sf_floatalloc(n12);
    pp = sf_floatalloc(n12);
    known = sf_boolalloc(n12);

    if (NULL != sf_getstring("mask")) {
	/* optional input mask file for known data */
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    if (!sf_getfloat("filt1",&filt1)) filt1=3.;
    if (!sf_getfloat("filt2",&filt2)) filt2=filt1;
    /* smoothing radius */
    
    if (!sf_getfloat("a0",a))   a[0] = (filt1*filt1-1.)/12.;
    if (!sf_getfloat("b0",a+1)) a[1] = 0.;
    if (!sf_getfloat("c0",a+2)) a[2] = (filt2*filt2-1.)/12.;
    /* initial Gaussian shape parameters */
    
    if (!sf_getfloat("eps",&eps)) eps=0.0001;
    /* regularization parameter */

    if (!sf_getbool("shape",&shape)) shape=false;
    /* if y, estimate shaping */

    if (!sf_getbool("force",&force)) force=true;
    /* if y, keep known values */

    mask_init(known);
    gaussshape2_init(n1,n2);
    sf_conjgrad_init(n12, n12, n12, n12, eps, 1.e-9, true, false);

    if (!shape) gaussshape2_set2(a);

    for (i3=0; i3 < n3; i3++) {
	sf_read(mm,sizeof(float),n12,in);

	if (NULL != mask) {
	    sf_read(kk,sizeof(float),n12,mask);
	    for (i=0; i < n12; i++) {
		known[i] = (kk[i] != 0.);
	    }
	} else {
	    for (i=0; i < n12; i++) {
		known[i] = (mm[i] != 0.);
	    }
	}

	if (force) {
	    for (i=0; i < n12; i++) {
		if (known[i]) kk[i]=mm[i];
	    }
	}
 
	if (shape) gaussshape2_set(a, mm, 100);
	sf_conjgrad(NULL, mask_lop, freqfilt2_lop, pp, mm, mm, niter);

	if (force) {
	    for (i=0; i < n12; i++) {
		if (known[i]) mm[i]=kk[i];
	    }
	}

	sf_write(mm,sizeof(float),n12,out);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mmiss2.c,v 1.4 2004/04/05 14:35:11 fomels Exp $	 */

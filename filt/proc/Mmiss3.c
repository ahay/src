/* Missing data interpolation (N-dimensional) using shaping regularization.

Takes: < input.rsf > interpolated.rsf
*/

#include <rsf.h>

#include "mask.h"
#include "trianglen.h"

int main(int argc, char* argv[])
{
    int niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM], n12, i;
    float *mm, *kk, *pp, eps;
    bool *known, force;
    char key[6];
    sf_file in, out, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */

    dim = sf_filedims (in,n);
    n12 = 1;
    for (i=0; i < dim; i++) {
	n12 *= n[i];
	if (n[i] > 1) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	} else {
	    rect[i]=1;
	}
    }

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
    
    if (!sf_getfloat("eps",&eps)) eps=0.0001;
    /* regularization parameter */

    if (!sf_getbool("force",&force)) force=true;
    /* if y, keep known values */

    mask_init(known);
    trianglen_init(dim, rect, n);
    sf_conjgrad_init(n12, n12, n12, n12, eps, 1.e-9, true, false);

    sf_floatread(mm,n12,in);

    if (NULL != mask) {
	sf_floatread(kk,n12,mask);
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
 
    sf_conjgrad(NULL, mask_lop, trianglen_lop, pp, mm, mm, niter);

    if (force) {
	for (i=0; i < n12; i++) {
	    if (known[i]) mm[i]=kk[i];
	}
    }
    
    sf_floatwrite(mm,n12,out);

    sf_close();
    exit (0);
}

/* 	$Id: Mmiss3.c,v 1.1 2004/05/07 03:40:51 fomels Exp $	 */

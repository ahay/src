/* 2-D missing data interpolation by differential offset continuation.

Takes: < input.rsf known=mask.rsf > filled.rsf
*/
#include <math.h>

#include <rsf.h>

#include "offruffp.h"
  
int main(int argc, char* argv[])
{
    int nx, nh, nw, iw, niter, nhx, i, *zero;
    float dx, h0, dh, w0, dw, w;
    float complex *slice, *dat;
    bool simple, *mask;
    sf_file in, out, known;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    known = sf_input("known");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (SF_INT != sf_gettype(known)) sf_error("Need int input in known");

    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");
    nhx = nh*nx;

    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");
    w0 *= 2.*SF_PI;
    dw *= 2.*SF_PI;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    if (!sf_getbool("simple",&simple)) simple=false;
    /* if y, use simple h derivative for regularization */

    slice = sf_complexalloc(nhx);
    mask = sf_boolalloc(nhx);
    zero = sf_intalloc(nhx);
    dat = sf_complexalloc(nhx);

    sf_intread(zero,nhx,known);
    for (i=0; i < nhx; i++) {
	dat[i] = 0.;
	mask[i] = (zero[i] != 0);
    }

    for (iw=0; iw < nw; iw++) {
	sf_warning("frequency %d of %d",iw+1,nw);
     
	w = w0 + iw*dw;
	sf_complexread (slice,nhx,in);

	if (fabsf(w) < dw) {
	    sf_complexwrite(dat,nhx,out);
	} else {
	    offruffp_init (h0, nh, dh, nx, dx, w, 0);
	    sf_csolver (simple? hderp_lop : offruffp_lop, sf_ccgstep, nhx, nhx,
			slice, dat, niter, "x0", slice, "known", mask, "end");
	    sf_ccgstep_close();

	    sf_complexwrite(slice,nhx,out);
	}
    }

    sf_close();
    exit(0);
}



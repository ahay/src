/* 1-D convolution.

Takes: < inp.rsf filtin=filt.rsf > out.rsf
*/
#include <rsf.h>

#include "tcai1.h"
#include "icai1.h"

int main(int argc, char* argv[])
{
    int nx, nf, ny, i2, n2, lag;
    bool each, trans;
    float *xx, *yy, *ff;
    sf_file in, out, filt;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    filt = sf_input("filt");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(filt)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histint(filt,"n1",&nf)) sf_error("No n1= in filtin");

    xx = sf_floatalloc(nx);
    ff = sf_floatalloc(nf);

    if (!sf_getbool("trans",&trans)) trans=false;
    /* if y, transient convolution; if n, internal */
    if (!sf_getbool("each",&each)) each=false;
    /* if y, new filter for each trace */
    if (!sf_getint("lag",&lag)) lag=1;
    /* lag for internal convolution */

    ny = trans? nx+nf-1: nx;
    yy = sf_floatalloc(ny);
    if (trans) sf_putint(out,"n1",ny);

    if (!each) sf_floatread (ff,nf,filt);

    for (i2=0; i2 < n2; i2++) {
        if (each) sf_floatread (ff,nf,filt);
        sf_floatread (xx,nx,in);
        if (trans) {
	    tcai1_init(nf,ff);
	    tcai1_lop (false, false,nx,ny,xx,yy);
	} else {
	    icai1_init(nf,ff,lag);
	    icai1_lop (false, false,nx,ny,xx,yy);
	}
	sf_floatwrite (yy,ny,out);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mconv.c,v 1.2 2004/04/19 21:51:46 fomels Exp $	 */

/* 1-D inverse interpolation with shaping regularization.

Takes: < irregular.rsf head=header.rsf > regular.rsf
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "int1.h"
#include "interp.h"
#include "gauss.h"
#include "freqfilt.h"
#include "triangle1.h"
#include "monofshape.h"

int main (int argc, char* argv[])
{
    bool pef, gauss;
    int id, nd, nt, it, nx, interp, niter;
    float *pp, *mm, *dd, *offset, x0, dx, xmin, xmax, f, filt, eps;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) nd=1;
    if (!sf_histint(in,"n2",&nt)) nt=1;
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    /* create coordinates */
    offset = sf_floatalloc(nd);

    head = sf_input("head");
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float head");
    sf_floatread (offset,nd,head);
    sf_fileclose (head);

    xmin = +FLT_MAX;
    xmax = -FLT_MAX;
    for (id=0; id<nd; id++) {	
	f = offset[id]; 
	if (f < xmin) xmin=f;
	if (f > xmax) xmax=f;
    }
 
    /* create model */
    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    /* number of bins */

    sf_putint(out,"n1",nx);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    /* grid size */
    sf_getfloat ("xmax",&xmax);
    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    /* grid origin */
    sf_putfloat (out,"o1",x0);

    if (!sf_getfloat("dx",&dx)) {
	/* grid sampling */
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }
    sf_putfloat (out,"d1",dx);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* [1,2] forward interpolation method, 1: nearest neighbor, 2: linear */

    switch (interp) {
	case 1:
	    int1_init (offset, x0,dx,nx, bin_int, 1, nd);
	    sf_warning("Using nearest-neighbor interpolation");
	    break;
	case 2:
	    int1_init (offset, x0,dx,nx, lin_int, 2, nd);
	    sf_warning("Using linear interpolation");
	    break;
	default:
	    sf_error("Unsupported interp=%d",interp);
	    break;
    }

    if (!sf_getfloat("filter",&filt)) filt=3.;
    /* smoothing length */

    if (filt < 1.) sf_error("filt must be >= 1, got filt=%g",filt);

    pp = sf_floatalloc(nx);
    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nx;
    /* number of conjugate-gradient iterations */
    if (!sf_getfloat("eps",&eps)) eps=1./nd;
    /* regularization parameter */
    if (!sf_getbool("pef",&pef)) pef=false;
    /* if y, use monofrequency regularization */
    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* if y, use Gaussian shaping */

    if (gauss) {
	gauss_init (nx, filt);
    } else {
	triangle1_init ((int) filt, nx);
    }
    sf_conjgrad_init(nx, nx, nd, nd, eps, 1.e-9, true, false);

    if (pef) monofshape_init(nx);

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_floatread (dd,nd,in);

	if (gauss) {
	    sf_conjgrad(NULL, int1_lop, freqfilt_lop, pp, mm, dd, niter);
	} else {
	    sf_conjgrad(NULL, int1_lop, triangle1_lop, pp, mm, dd, niter);
	}

	if (pef) {
	    monofshape_set(0.1,nx,mm,100);
	    sf_conjgrad(NULL, int1_lop, monofshape_lop, pp, mm, dd, niter);
	}
	
    	sf_floatwrite (mm,nx,out);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mshapebin1.c,v 1.6 2004/04/19 21:51:46 fomels Exp $	 */

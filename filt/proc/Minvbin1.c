#include <float.h>
#include <math.h>

#include <rsf.h>

#include "int1.h"
#include "interp.h"
#include "tcai1.h"
#include "cgstep.h"
#include "bigsolver.h"

int main (int argc, char* argv[])
{
    int id, nd, nt, it, nx, interp, filt, meth, niter;
    float *mm, *dd, *offset, *aa, x0, dx, xmin, xmax, f, eps;
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
    sf_read (offset,sizeof(float),nd,head);
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

    sf_putint(out,"n1",nx);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    sf_getfloat ("xmax",&xmax);
    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    sf_putfloat (out,"o1",x0);

    if (!sf_getfloat("dx",&dx)) {
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }
    sf_putfloat (out,"d1",dx);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=1;

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

    if (!sf_getint("filter",&filt)) filt=1;
    if (!sf_getint("meth",&meth)) meth=1;

    switch (meth) {
	case 1:
	    filt++;
	    sf_warning("Using differential regularization");	    
	    aa = sf_floatalloc(filt);
	    tcai1_init(filt,aa);
	    switch (filt) {
		case 2:
		    aa[0] = 1.;
		    aa[1] = -1.;
		    break;
		default:
		    aa[0] = 1.;
		    aa[1] = -2.;
		    aa[2] = 1.;
		    break;
	    }
	    break;
    }   

    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nx;
    if (!sf_getfloat("eps",&eps)) eps=0.2;

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_read (dd,sizeof(float),nd,in);
	switch (meth) {
	    case 1:
		solver_reg(int1_lop, cgstep, tcai1_lop, nx+filt, nx, nd, 
			   mm, dd, niter, eps, "end");
		break;
	}
	cgstep_close();
	sf_write (mm,sizeof(float),nx,out);
    }
    
    exit(0);
}



/* Data binning in 1-D slices.

Takes: < input.rsf head=header.rsf > binned.rsf
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "int1.h"
#include "interp.h"

int main (int argc, char* argv[])
{
    int id, nd, nt, it, ix, nx, interp;
    float *mm, *count, *dd, *offset, x0, dx, xmin, xmax, f, clip;
    sf_file in, out, head, fold;

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
    /* Number of bins */

    sf_putint(out,"n1",nx);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    /* grid dimensions */
    sf_getfloat ("xmax",&xmax);

    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    /* grid origin */
    sf_putfloat (out,"o1",x0);

    if (!sf_getfloat("dx",&dx)) {
	/* grid spacing */
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }
    sf_putfloat (out,"d1",dx);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=1;
    /* [1,2] interpolation method, 1: nearest neighbor, 2: linear */

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

    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);
    count  = sf_floatalloc(nx);

    /* compute fold */
    for (id=0; id<nd; id++) {
	dd[id]=1.;
    }

    int1_lop (true, false,nx,nd,count,dd);
 
    if (NULL != sf_getstring("fold")) {
	/* output fold file (optional) */
	fold = sf_output("fold");
	sf_putint(fold,"n1",nx);
	sf_putint(fold,"n2",1);
	sf_putfloat(fold,"o1",x0);
	sf_putfloat(fold,"d1",dx);
	sf_write (count,sizeof(float),nx,fold);
	sf_fileclose (fold);
    }

    if (!sf_getfloat("clip",&clip)) clip = FLT_EPSILON;
    /* clip for fold normalization */

    for (ix=0; ix<nx; ix++) {
	if (clip < count[ix]) count[ix]=1./fabsf(count[ix]);
	else                  count[ix]=0.;
    }

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_read (dd,sizeof(float),nd,in);
	int1_lop (true,false,nx,nd,mm,dd);
	for (ix=0; ix<nx; ix++) {
	    mm[ix] *= count[ix];
	}
	sf_write (mm,sizeof(float),nx,out);
    }
    
    exit(0);
}

/* 	$Id: Mbin1.c,v 1.4 2003/10/01 14:38:31 fomels Exp $	 */


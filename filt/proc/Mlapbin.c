/* Data binning in 2-D slices by inverse interpolation.

Takes: < input.rsf head=header.rsf > binned.rsf
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "int2.h"
#include "interp.h"
#include "laplac2.h"
#include "gauss2.h"

int main (int argc, char* argv[])
{
    bool gauss;
    int id, nk, nd, nm, nt, it, nx, ny, n2, xkey, ykey, interp, niter;
    float *pp, *mm, *dd, **xy, *hdr, filt1, filt2;
    float x0, y0, dx, dy, xmin, xmax, ymin, ymax, f, dt, t0, eps;
    char *xk, *yk;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nt)) sf_error("Need n2= in in");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (NULL != (xk = sf_getstring("xk"))) {
	/* x key name */
	xkey = sf_segykey(xk);
    }  else if (!sf_getint("xkey",&xkey)) {
	/* x key number (if no xk), default is sx */
	xkey = sf_segykey("sx");
    }
    if (NULL != (yk = sf_getstring("yk"))) {
	/* y key name */
	ykey = sf_segykey(yk);
    }  else if (!sf_getint("ykey",&ykey)) {
	/* y key number (if no yk), default is sy */
	ykey = sf_segykey("sy");
    }

    /* create coordinates */
    xy = sf_floatalloc2(2,nd);
    head = sf_input("head");

    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

    hdr = sf_floatalloc(nk);

    ymin = xmin = +FLT_MAX;
    ymax = xmax = -FLT_MAX;
    for (id=0; id<nd; id++) {	
	sf_read (hdr,sizeof(float),nk,head);
	f = hdr[xkey]; 
	if (f < xmin) xmin=f;
	if (f > xmax) xmax=f;
	xy[id][0] = f;
	f = hdr[ykey]; 
	if (f < ymin) ymin=f;
	if (f > ymax) ymax=f;
	xy[id][1] = f;
    }

    sf_fileclose (head);

    /* create model */
    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    /* Number of bins in x */
    if (!sf_getint ("ny",&ny)) sf_error("Need ny=");
    /* Number of bins in y */

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_putint(out,"n3",nt);

    if (sf_histfloat(in,"o2",&t0)) sf_putfloat(out,"o3",t0);
    if (sf_histfloat(in,"d2",&dt)) sf_putfloat(out,"d3",dt);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    sf_getfloat ("xmax",&xmax);
    sf_getfloat ("ymin",&ymin);
    /* Grid dimensions */
    sf_getfloat ("ymax",&ymax);

    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);
    if (ymax <= ymin) sf_error ("ymax=%f <= ymin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    if (!sf_getfloat("y0",&y0)) y0=ymin; 
    /* grid origin */

    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"o2",y0);

    if (!sf_getfloat("dx",&dx)) {
	/* bin size in x */
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }

    if (!sf_getfloat("dy",&dy)) {
	/* bin size in y */
	if (1 >= nx) {
	    dy = dx;
	} else {
	    dy = (ymax-ymin)/(ny-1);
	}
    }

    sf_putfloat (out,"d1",dx);
    sf_putfloat (out,"d2",dy);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* [1,2] interpolation method, 1: nearest neighbor, 2: bi-linear */

    switch (interp) {
	case 1:
	    int2_init (xy, x0,y0,dx,dy,nx,ny, bin_int, 1, nd);
	    sf_warning("Using nearest-neighbor interpolation");
	    break;
	case 2:
	    int2_init (xy, x0,y0,dx,dy,nx,ny, lin_int, 2, nd);
	    sf_warning("Using linear interpolation");
	    break;
	case 3:
	    sf_error("Unsupported interp=%d",interp);
	    break;
    }

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nm;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1./nd;
    /* regularization parameter */
    if (!sf_getbool("gauss",&gauss)) gauss=true;
    /* if y, use gaussian shaping */

    if (gauss) {
	if (!sf_getfloat("filt1",&filt1)) filt1=3.;
	if (!sf_getfloat("filt2",&filt2)) filt2=filt1;
	/* smoothing length for shaping */

	if (filt1 < 1. || filt2 < 1.) 
	    sf_error("wrong filt1=%g or filt2=%g",filt1,filt2);

	gauss2_init(nx,ny,filt1,filt2);
	pp = sf_floatalloc(nm);
    } else {
	laplac2_init(nx,ny);
	pp = NULL;
    }

    sf_conjgrad_init(nm, nm, nd, eps, 1.e-9, true, false);

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_read (dd,sizeof(float),nd,in);

	if (gauss) {
	    sf_conjgrad(int2_lop, gauss2_lop, pp, mm, dd, niter);
	} else {
	    sf_solver_reg(int2_lop,sf_cgstep,laplac2_lop,
			  nm,nm,nd,mm,dd,niter,eps,"end");
	    sf_cgstep_close();
	}

	sf_write (mm,sizeof(float),nm,out);
    }

    exit(0);
}

/* 	$Id: Mlapbin.c,v 1.1 2004/02/14 07:01:42 fomels Exp $	 */

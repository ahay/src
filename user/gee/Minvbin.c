/* Data interpolation in 2-D slices using helix preconditioning. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "npolydiv.h"
#include "regrid.h"
#include "gradint2.h"

int main (int argc, char* argv[])
{
    int id, nk, nd, nm, nt, it, nx, ny, n2, xkey, ykey, na, ia, niter;
    int n[2], m[2], interp, np, *nh, *ph, ip;
    float *mm, *dd, **xy, *hdr;
    float x0, y0, dx, dy, xmin, xmax, ymin, ymax, f, dt, t0, a0, eps;
    char *lagfile, *nhfile, *header;
    bool stat, der;
    sf_filter aa;
    nfilter naa=NULL;
    sf_file in, out, head, flt, lag, nhh, pch=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nt)) sf_error("Need n2= in in");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    /* create coordinates */
    xy = sf_floatalloc2(2,nd);

    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);

    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

    if (!sf_getint("xkey",&xkey)) sf_error("Need xkey=");
    /* x key number */
    if (!sf_getint("ykey",&ykey)) sf_error("Need ykey=");
    /* y key number */

    if (xkey < 0 || xkey >= nk) 
	sf_error("xkey=%d is out of the range [0,%d]",xkey,nk-1);
    if (ykey < 0 || ykey >= nk) 
	sf_error("ykey=%d is out of the range [0,%d]",ykey,nk-1);

    hdr = sf_floatalloc(nk);

    ymin = xmin = +FLT_MAX;
    ymax = xmax = -FLT_MAX;
    for (id=0; id<nd; id++) {	
	sf_floatread (hdr,nk,head);
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

    if (sf_histfloat(in,"o2",&t0)) sf_putfloat(out,"o3",t0);
    if (sf_histfloat(in,"d2",&dt)) sf_putfloat(out,"d3",dt);

    /* let user overwrite */
    sf_getfloat ("xmax",&xmax);
    sf_getfloat ("xmin",&xmin);
    sf_getfloat ("ymax",&ymax);
    sf_getfloat ("ymin",&ymin);
    /* Grid dimensions */
    
    if (xmax < xmin) sf_error ("xmax=%f < xmin=%f",xmax,xmin);
    if (ymax < ymin) sf_error ("ymax=%f < ymin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    if (!sf_getfloat("y0",&y0)) y0=ymin; 
    /* grid origin */

    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"o2",y0);

    /* create model */
    if (!sf_getint ("nx",&nx)) nx = (int) (xmax - xmin + 1.);
    /* Number of bins in x */
    if (!sf_getint ("ny",&ny)) ny = (int) (ymax - ymin + 1.);
    /* Number of bins in y */

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_putint(out,"n3",nt);

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

    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation length */

    if (!sf_getbool("der",&der)) der=false;
    /* if y, apply derivative filter on the residual */

    if (der) {
	gradint2_init (xy, x0,y0, dx,dy, nx,ny, sf_spline_int, interp, nd);
    } else {
	sf_int2_init (xy, x0,y0, dx,dy, nx,ny, sf_spline_int, interp, nd);
    }

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);

    if (!sf_getbool("stat",&stat)) stat=true;
    /* stationary or nonstationary filter */

    flt = sf_input("filt");

    if (NULL != (lagfile = sf_getstring("lag")) /* file with filter lags */
	|| 
	NULL != (lagfile = sf_histstring(flt,"lag"))) {
	lag = sf_input(lagfile);
    } else {
	lag = NULL;
    }

    n[0] = nx;
    n[1] = ny;
    if (!sf_getints ("n",m,2) && 
	(NULL == lag || !sf_histints (lag,"n",m,2))) {
	m[0] = nx;
	m[1] = ny;
    }

    if (stat) {
	if (!sf_histint(flt,"n1",&na)) sf_error("No n1= in filt");
	aa = sf_allocatehelix (na);

	if (!sf_histfloat(flt,"a0",&a0)) a0=1.;
	sf_floatread (aa->flt,na,flt);

	for( ia=0; ia < na; ia++) {
	    aa->flt[ia] /= a0;
	}
	
	if (NULL != lag) {	    
	    sf_intread(aa->lag,na,lag);
	    sf_fileclose(lag);
	} else {
	    for( ia=0; ia < na; ia++) {
		aa->lag[ia] = ia+1;
	    }
	}
	
	sf_fileclose(flt);

	regrid (2, m, n, aa);
	sf_polydiv_init (nm, aa);
    } else {
	if (NULL == (nhfile = sf_getstring("nh")) /* file with filter sizes */
	    &&
	    NULL == (nhfile = sf_histstring(flt,"nh")))
	    sf_error("Need nh=");

	nhh = sf_input(nhfile);

	if (SF_INT != sf_gettype(nhh)) sf_error("Need int nh");
	if (!sf_histint(nhh,"n1",&np)) sf_error("No n1= in nh");

	ph = sf_intalloc(nm);
	nh = sf_intalloc(np);

	sf_intread(nh,np,nhh);
	sf_fileclose(nhh);

	naa = nallocate(np,nm,nh,ph);
	free(ph);

	pch = sf_input("pch");

	if (NULL == lag || SF_INT != sf_gettype(lag)) 
	    sf_error("Need int lag=");

	for (ip=0; ip < np; ip++) {
	    na = nh[ip];
	    aa = naa->hlx[ip];
	    sf_intread(aa->lag,na,lag);
	    sf_floatread(aa->flt,na,flt);
	    regrid (2, m, n, aa);
	}
	free (nh);
	sf_fileclose(lag);
	sf_fileclose(flt);

	npolydiv_init(nm, naa);
    }

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=0.01; 
    /* regularization parameter */

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_warning("%d of %d;",it+1,nt);
	
	sf_floatread (dd,nd,in);

	if (stat) {
	    if (der) {
		sf_solver_prec (gradint2_lop, sf_cgstep, sf_polydiv_lop,
				nm, nm, nd, mm, dd, niter, eps, "end");
	    } else {
		sf_solver_prec (sf_int2_lop, sf_cgstep, sf_polydiv_lop,
				nm, nm, nd, mm, dd, niter, eps, "end");
	    }
	} else {
	    sf_intread(naa->pch,nm,pch);

	    sf_solver_prec (sf_int2_lop, sf_cgstep, npolydiv_lop,
			    nm, nm, nd, mm, dd, niter, eps, "end");
	}

	sf_cgstep_close();
	sf_floatwrite (mm,nm,out);
    }
    sf_warning(".");


    exit(0);
}

/* 	$Id: Mbin.c 991 2005-02-10 11:42:19Z fomels $	 */

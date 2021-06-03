/*  Velocity transform for generating velocity spectra and its corresponding hyperbolic modeling. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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
#include <rsf.h>

#include "velxf.h"
#include "l1.h"
#include "l1step.h"

int main(int argc, char* argv[])
{
    bool *mask, adj, robust;
    int it,ix,is, nt,nx,ns, ntx, nts, niter, nliter;
    float *data, *modl, *x2, *z2, *s;
    float dt,dx,ds, x, z, ot,ox,os, perc,fact, eps;
    char *type;
    sf_file cmp, vtr;

    sf_init(argc,argv);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adj = 0: from velocity-domain(t,s) to cmp-gather(t,x)
       adj = 1: from cmp-gather(t,x) to velocity-domain(t,s) */

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations (invoked if adj=y) */


    if (adj) {
	cmp = sf_input("in");
	vtr = sf_output("out");

	if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(cmp,"n2",&nx)) sf_error("No n2= in input");

	if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
   
	if (!sf_histfloat(cmp,"o1",&ot)) sf_error("No o1= in input");
	if (!sf_histfloat(cmp,"o2",&ox)) sf_error("No o2= in input");

	if (!sf_getint("ns",&ns)) ns=nx;
	/* number of slowness values */
	if (!sf_getfloat("ds",&ds)) ds=0.001;
	/* slowness sampling */
	if (!sf_getfloat("os",&os)) os=0.00000001;
	/* slowness origin */

	sf_putint(vtr,"n2",ns);
	sf_putfloat(vtr,"d2",ds);
	sf_putfloat(vtr,"o2",os);
    } else {
	vtr = sf_input("in");
	cmp = sf_output("out");

	if (!sf_histint(vtr,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(vtr,"n2",&ns)) sf_error("No n2= in input");

	if (!sf_histfloat(vtr,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(vtr,"d2",&ds)) sf_error("No d2= in input");
   
	if (!sf_histfloat(vtr,"o1",&ot)) sf_error("No o1= in input");
	if (!sf_histfloat(vtr,"o2",&os)) sf_error("No o2= in input");

	if (!sf_getint("nx",&nx)) nx=ns;
	/* number of offset values */
	if (!sf_getfloat("dx",&dx)) dx=10.;
	/* offset sampling */
	if (!sf_getfloat("ox",&ox)) ox=0.;
	/* offset origin */

	sf_putint(cmp,"n2",nx);
	sf_putfloat(cmp,"d2",dx);
	sf_putfloat(cmp,"o2",ox);
    }

    ntx = nt*nx;
    nts = nt*ns;

    data = sf_floatalloc(ntx);
    modl = sf_floatalloc(nts);
    mask = sf_boolalloc(nx);
    x2 = sf_floatalloc(nx);
    z2 = sf_floatalloc(nt);
    s  = sf_floatalloc(ns);

    if (adj) {
	sf_floatread(data, ntx, cmp);
    } else {
	sf_floatread(modl, nts, vtr);
    }

    for (ix=0; ix < nx; ix++) {
	mask[ix] = true;
	x = ox + dx*ix;
	x2[ix] = x*x;
    }

    for (it=0; it < nt; it++) {
	z = ot + dt*it;
	z2[it] = z*z;
    }

    for (is=0; is < ns; is++) {
	s[is] = os + ds*is;
    }

    velxf_init(nt, nx, ns,
	       ot, dt,
	       x2, z2, s, mask);

    if (adj && niter > 0) {
	    if (!sf_getbool("robust",&robust)) robust=false;
	    /* robust inversion */
	    
	    if (robust) {
		if (!sf_getfloat("perc",&perc)) perc=90.;
		/* threshold percentage for robust inversion */

		if (!sf_getfloat("fact",&fact)) fact=1.5;
		/* threshold factor for robust inversion */
		
		if (!sf_getint("nliter",&nliter)) nliter=10;
		/* number of POCS iterations for robust inversion */
		
		if (!sf_getfloat("eps",&eps)) eps=1.;
		/* regularization parameter for robust inversion */

		if (NULL == (type = sf_getstring("type"))) type="threshold";
		/* thresholding type */

		l1_init(ntx,nliter,perc,fact,type,true);
		sf_solver_reg(velxf,l1step,sf_copy_lop,nts,nts,ntx,modl,data,niter,eps,"verb",true,"end");
	    } else {		
		sf_solver(velxf,sf_cgstep,nts,ntx,modl,data,niter,"verb",true,"end");
	    }
    } else {
	velxf(adj,false,nts,ntx,modl,data);
    }

    if (adj) {
	sf_floatwrite(modl, nts, vtr);
    } else {
	sf_floatwrite(data, ntx, cmp);
    }

    exit(0);
}



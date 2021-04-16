/* Hyperbolic Radon transform with conjugate-directions inversion */
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

#include <rsf.h>
#include "veltran.h"

int main(int argc, char* argv[])
{
    int nt, n3, nx, nv, i3, ntx, ntv, ic, nc;
    int niter, miter, psun1, psun2;
    bool adj;
    float o1,d1, x0,dx, v0,dv, anti, s1,s0,ds, perc,fact;
    float *cmp=NULL, *vscan=NULL, *error=NULL, *mask=NULL;
    sf_file in=NULL, out=NULL, err=NULL, msk=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    n3 = sf_leftsize(in,2);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */
    if (!sf_getint("miter",&miter)) miter=2;
    /* conjugate-direction memory */
    if (!sf_getint("psun1",&psun1)) psun1=1;
    /* amplitude type for adjoint */
    if (!sf_getint("psun2",&psun2)) psun2=1;
    /* amplitude type for forward */
    if (!sf_getfloat("anti",&anti)) anti=1.;
    /* antialiasing */

    if(niter > 0)  adj=false;

    if (adj) {
	if (!sf_histint(in,"n2",&nv))   sf_error("No n2= in input");
	if (!sf_histfloat(in,"v0",&v0)) sf_error("No v0= in input");
	if (!sf_histfloat(in,"dv",&dv)) sf_error("No dv= in input");

	if (!sf_getint("nx",&nx) && !sf_histint(in,"nx",&nx)) nx=nv;
	if (!sf_getfloat("x0",&x0) && !sf_histfloat(in,"x0",&x0))
	    sf_error("need x0=");
	if (!sf_getfloat("dx",&dx) && !sf_histfloat(in,"dx",&dx))
	    sf_error("need dx=");
    } else {
	if (!sf_histint(in,"n2",&nx))   sf_error("No n2= in input");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	
	if (!sf_getint("nv",&nv) && !sf_histint(in,"nv",&nv)) nv=nx;
	if (!sf_getfloat("v0",&v0) && !sf_histfloat(in,"v0",&v0))
	    sf_error("need v0=");
	if (!sf_getfloat("dv",&dv) && !sf_histfloat(in,"dv",&dv))
	    sf_error("need dv=");

	sf_putfloat(out,"x0",x0);
	sf_putfloat(out,"dx",dx);
	sf_putint(out,"nx",nx);
		
	sf_putfloat(out,"v0",v0);
	sf_putfloat(out,"dv",dv);
    }

    if (!sf_getfloat("s02",&s1)) s1=0.;
    /* reference slowness squared (for antialiasing) */

    s0 = 1./(v0 +(nv-1)*dv); s0 *= s0; ds = (1./(v0*v0) - s0)/(nv-1);
/*    s0 = ds;		               ds = (1./(v0*v0) - s0)/(nv-1); */

    if (adj) {
	sf_putstring(out,"label2","Offset");
	sf_putstring(out,"unit2","km");
	sf_putfloat(out,"o2",x0);
	sf_putfloat(out,"d2",dx);
	sf_putint(out,"n2",nx);
    } else {
	sf_putstring(out,"label2","Slowness Squared");
	sf_putstring(out,"unit2","s\\^2\\_/km\\^2");
	sf_putfloat(out,"o2",s0);
	sf_putfloat(out,"d2",ds);
	sf_putint(out,"n2",nv);
    }

    if(niter > 0 && NULL != sf_getstring("error")) {
	err = sf_output("error");
	
	sf_putint(err,"n1",niter);
	sf_putfloat(err,"o1",1.);
	sf_putfloat(err,"d1",1.);
	sf_putint(err,"n2",1);
	sf_putstring(err,"label1","Iteration Number");
	sf_putstring(err,"label2","Relative Squared Error");

	error = sf_floatalloc(niter);
    } else {
	err = NULL;
	error = NULL;
    }

    ntx = nt*nx;
    ntv = nt*nv;

    if (NULL != sf_getstring("mask")) {
	msk = sf_input("mask");
	mask = sf_floatalloc(ntv);
    } else {
	msk = NULL;
	mask = NULL;
    }

    if (!sf_getint("ncycle",&nc)) nc=0;
    /* number of sharpening cycles */

    if (nc > 0) {
	if (!sf_getfloat("perc",&perc)) perc=50.0;
	/* percentage for sharpening */

	if (!sf_getfloat("fact",&fact)) fact=0.5;
	/* factor for sharpening */

	sf_sharpen_init(ntv,perc,fact);
	if (NULL == mask) mask = sf_floatalloc(ntv);
	sf_conjgrad_init(ntv,ntv,ntx,ntx,1.0,1.e-6,true,false);
    }

    cmp = sf_floatalloc(ntx);
    vscan = sf_floatalloc(ntv);

    veltran_init (true, x0, dx, nx, s0, ds, nv, o1, d1, nt, 
		  s1, anti, psun1, psun2);

    for (i3=0; i3 < n3; i3++) {
	if( adj) {
	    sf_floatread(vscan,ntv,in);
	} else {		
	    sf_floatread(cmp,ntx,in);
	}
	
	if(0==niter) {
	    veltran_lop((bool) !adj,false,ntv,ntx,vscan,cmp);
	} else { 
	    if (nc > 0) {
		for (ic=0; ic < nc; ic++) {
		    sf_conjgrad(NULL,veltran_lop,sf_weight_lop,
				mask,vscan,cmp,niter);
		    sf_sharpen(vscan);
		}
	    } else {
		if (NULL != msk) sf_floatread(mask,ntv,msk);
		
		sf_cdstep_init();
		sf_solver (veltran_lop, sf_cdstep, 
			   ntv, ntx, vscan, cmp, niter, 
			   "err", error, "nmem", 0, "nfreq", miter, 
			   "mwt", mask, "end");
		sf_cdstep_close();
		
		if(NULL != err) sf_floatwrite(error,niter,err);
	    }
	}

	if( adj) {
	    sf_floatwrite(cmp,ntx,out);
	} else {		
	    sf_floatwrite(vscan,ntv,out);
	}
    }

    exit(0);
}

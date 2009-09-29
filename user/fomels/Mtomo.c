/* Simple tomography test. */
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

#include "tomo.h"

int main(int argc, char* argv[]) {
    bool adj;
    int np, nz, nx, nt, ns, nq, dq, niter, rect1, rect2;
    float dz, dx, eps, tol;
    float *s, *t, *p;
    sf_file slow, time;

    /* initialize parameters */
    sf_init(argc,argv);

    /* read parameter from the command line adj=y or adj=n */
    if (!sf_getbool("adj",&adj)) adj=false;
    /* if n, generate traveltime from slowness; 
       if y, invert slowness from traveltime */

    if (adj) {
	time = sf_input("in");
	slow = sf_output("out");

	if (!sf_getint("niter",&niter)) niter=100;
	/* maximum number of iterations */

	if (!sf_getfloat("eps",&eps)) eps=1.;
	/* scaling parameter */

	if (!sf_getfloat("tol",&tol)) tol=1.e-7;
	/* tolerance for stopping iterations */

	if (!sf_getint("rect1",&rect1)) rect1=1;
	if (!sf_getint("rect2",&rect2)) rect2=1;
	/* smoothing length in z and x */

	if (!sf_histint(time,"n1",&np)) sf_error("No n1= in input");
	if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(time,"d2",&dx)) sf_error("No d2= in input");
    
	if (!sf_getint("nz",&nz) && !sf_histint(time,"nz",&nz)) 
	    sf_error("Need nz=");
	if (!sf_getfloat("dz",&dz) && !sf_histfloat(time,"dz",&dz)) 
	    sf_error("Need dz=");

	sf_putint(slow,"n1",nz);
	sf_putfloat(slow,"d1",dz);
	sf_putfloat(slow,"o1",0.);

	if (!sf_getint("ns",&nq) && !sf_histint(time,"n3",&nq)) 
	    sf_error("Need ns=");
	if (!sf_getint("ds",&dq) && !sf_histint(time,"d3",&dq)) 
	    sf_error("Need ds=");

	sf_putint(slow,"n3",1);
    } else {
	slow = sf_input("in");
	time = sf_output("out");

	if (!sf_histint(slow,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histint(slow,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(slow,"d1",&dz)) sf_error("No d1= in input");
	if (!sf_histfloat(slow,"d2",&dx)) sf_error("No d2= in input");
    
	if (!sf_getint("np",&np)) np=11;

	sf_putint(time,"n1",np);
	sf_putfloat(time,"d1",2./(np-1));
	sf_putfloat(time,"o1",-1.);

	sf_putint(time,"nz",nz);
	sf_putfloat(time,"dz",dz);    
	
	if (!sf_getint("ns",&nq)) nq=1;
	/* number of depth steps */
	if (!sf_getint("ds",&dq)) dq=nz;
	/* step size */
	
	sf_putint(time,"n3",nq);
	sf_putfloat(time,"d3",dq*dz);
	sf_putfloat(time,"o3",0.);

	sf_putint(time,"ds",dq);
    }

    ns = nz*nx;
    nt = nq*np*nx;

    /* allocate space */
    s = sf_floatalloc(ns);
    t = sf_floatalloc(nt);
    p = adj? sf_floatalloc(ns): NULL;

    tomo_init(np, nq, dq, nz, nx, dz, dx);

    if (adj) {
	sf_floatread(t, nt, time);

	/* initializtion for shaping (smoothing) */ 
	sf_triangle2_init(rect1,rect2,nz,nx,1);

	/* intitialize conjugate gradients */
	sf_conjgrad_init(ns, ns, nt, nt, eps, tol, true, false);
    
	/* invert */
	sf_conjgrad(NULL, tomo_lop, sf_triangle2_lop, p, s, t, niter);
    
	sf_floatwrite(s, ns, slow);

    } else {
	sf_floatread(s, ns, slow);

	tomo_lop(false, false, ns, nt, s, t);

	sf_floatwrite(t, nt, time);
    }
    

    exit(0);
}

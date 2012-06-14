/* 2-D Post-stack Kirchhoff redatuming. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kirdat.h"

int main(int argc, char* argv[])
{
    bool verb;
    int it, nt, ix, nx, left, right, ic, aper, shift, tap;
    float datum, length, t0, dt, x0, dx, dist, tau, delta, coef;
    float **tr_in, **tr_out, **table;
    sf_file in, out, green;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getfloat("datum",&datum)) sf_error("Need datum=");
    /* datum depth */

    if (!sf_getint("aperture",&aper)) aper=50;
    /* aperture (number of traces) */

    if (!sf_getint("taper",&tap)) tap=10;
    /* taper (number of traces) */

    if (!sf_getfloat("length",&length)) length=0.025;
    /* filter length (in seconds) */

    /* read input */
    if (!sf_histint(in,"n1",&nt)) sf_error("No nt=");
    if (!sf_histint(in,"n2",&nx)) sf_error("No nx=");
    
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No t0=");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No dt=");
    
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No x0=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No dx=");

    tr_in = sf_floatalloc2(nt,nx);
    sf_floatread(tr_in[0],nt*nx,in);

    /* allocate memory for output */
    tr_out = sf_floatalloc2(nt,nx);

    /* read Green's function (traveltime table) */
    green = sf_input("green");

    table = sf_floatalloc2(nx,nx);
    sf_floatread(table[0],nx*nx,green);
    sf_fileclose(green);

    /* initialize */
    filt_init(dt,length);

#ifdef _OPENMP
#pragma omp parallel for private(left,right,ic,coef,tau,dist,shift,it,delta)
#endif
    for (ix=0; ix < nx; ix++) {
	if (verb) sf_warning("Processing grid %d of %d.",ix+1,nx);	
	
	/* aperture */
	left  = (ix-aper < 0)?    0:    ix-aper;
	right = (ix+aper > nx-1)? nx-1: ix+aper;
	
	for (ic=left; ic <= right; ic++) {
	    	    
	    /* taper coefficient */
	    coef = 1.;
	    coef *= (ic-left  >= tap)? 1.: (ic-left)/tap;
	    coef *= (right-ic >= tap)? 1.: (right-ic)/tap;

	    /* time delay */
	    tau = 2.*table[ic][ix];
	    
	    /* distance */
	    dist = datum*datum+(ic-ix)*dx*(ic-ix)*dx;
	    
	    /* filter (tau dependent) */
	    filt_set(tau);
	    
	    shift = 0;
	    delta = 0.;
	    for (it=0; it < nt; it++) {
		if (((float)it)*dt < tau) 
		    continue;
		else if (shift == 0)
		    delta = (((float)it*dt)-tau)/dt;
		
		tr_out[ix][it] += coef/SF_PI
		    *dx*datum*tau/dist
		    *pick(delta,tr_in[ic],shift);
		shift++;
	    }
	}
    }
    
    /* write output */
    sf_floatwrite(tr_out[0],nt*nx,out);

    exit(0);
}

/* Post-stack Stolt modeling/migration. 

Requires the input to be cosine-transformed over the lateral axes.
*/
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

#include <math.h>
#include <rsf.h>
#include "shprefilter.h"

int main(int argc, char* argv[])
{
    int nt,nx,ny, iw,ix,iy, nf, nw;
    float dw, dt, dx,dy, t0, y, w,st,sq, *str=NULL, *trace2=NULL, *trace=NULL, vel, a, b, x;
    char *intp;
    sf_bands spl;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) nx=1;
    if (!sf_histint(in,"n3",&ny)) ny=1;

    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");
    /* Constant velocity (use negative velocity for modeling) */
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    
    intp = sf_getstring("interp");
     if (NULL == intp) intp[0]='s' ;

    if (!sf_getint ("pad",&nw)) nw=nt;
    /* padding on the time axis */
    nw=2*kiss_fft_next_fast_size(nw-1);

    sf_cosft_init(nw/2+1);
    dw = 2 * SF_PI/(nw*dt);

    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dy)) dy=dx;
    dx *= SF_PI * fabsf (vel);
    dy *= SF_PI * fabsf (vel);

    if (!sf_getfloat("stretch", &st) && !sf_histfloat(in,"stretch",&st)) st=1.;
    /* Stolt stretch parameter */
    if (vel < 0) st = 2.-st;
    a = (1.-1./st);
    b = 1./st;

    trace = sf_floatalloc(nw);
    trace2 = sf_floatalloc(nw);
    str = sf_floatalloc(nw);
    
    if (intp[0] == 's') {
	    if (!sf_getint("nf",&nf)) nf=2;
	    /* Interpolation accuracy */
	    spl = (nf > 2)? sf_spline_init (nf, nw): NULL;
    }
    for (iy = 0; iy < ny; iy++) {
	sf_warning("%d of %d;",iy+1,ny);
	y = iy*dy;
	y *= y;
	for (ix = 0; ix < nx; ix++) {
	    x = ix*dx;
	    x = st*(x*x + y);  

	    sf_floatread(trace,nt,in);
	    for (iw = nt; iw < nw; iw++) { /* pad */
		trace[iw]=0.;
	    }
	    sf_cosft_frw (trace,0,1);

	    for (iw = 0; iw < nw; iw++) {
		w = iw*dw;
		sq = (vel < 0)? w*w - x: w*w + x;
		if (sq > 0.) {
		    sq = sqrtf(sq);
		    str[iw] = a*w + b*sq;
		    /* trace[iw] *= (a + b*w/sq); Jacobian */
		} else {
		    str[iw] = - 2.*dw;
		    trace[iw] = 0.;
		}
	    }

	    if (intp[0] == 's') {
		    sf_int1_init (str, 0., dw, nw, sf_spline_int, nf, nw);	    
		    if (nf > 2) sf_banded_solve (spl, trace);
	    } else {
		    sf_int1sh_init (str, 0., dw, nw, sf_lin_int, nf, nw);
		    shprefilter(nw,trace); 
	    }
	    sf_int1_lop (false,false,nw,nw,trace,trace2);
	    sf_cosft_inv (trace2,0,1);
	    sf_floatwrite(trace2,nt,out);
	}
    }
    sf_warning(".");

    exit (0);
}

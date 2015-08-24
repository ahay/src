/* Post-stack Stolt modeling/migration. 

Requires the input to be cosine-transformed over the lateral axes.

August 2014 program of the month:
http://ahay.org/blog/2014/08/03/program-of-the-month-sfstolt/
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
#include "fint1.h"

static float a, b, x, vel;

static float stolt(float w, int iw) {
    float sq;

    sq = (vel < 0)? w*w - x: w*w + x;
    if (sq > 0.) {
	sq = sqrtf(sq);
	sq = a*w + b*sq;
    }
    return sq;
}

int main(int argc, char* argv[])
{
    fint1 map;
    int nt,nx,ny, iw,ix,iy, nf, nw, mute;
    float dw, dt, dx,dy, x0,y0, t0, y, st, *trace=NULL, minstr;
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

    if (!sf_getint ("pad",&nw)) nw=nt;
    /* padding on the time axis */
    nw=2*kiss_fft_next_fast_size(nw-1);

    sf_cosft_init(nw/2+1);
    dw = 2 * SF_PI/(nw*dt);

    if (!sf_histfloat(in,"o2",&x0)) x0=0.0; 
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_histfloat(in,"o3",&y0)) y0=0.0; 
    if (!sf_histfloat(in,"d3",&dy)) dy=dx;

    x0 *= SF_PI * fabsf (vel);
    y0 *= SF_PI * fabsf (vel);

    dx *= SF_PI * fabsf (vel);
    dy *= SF_PI * fabsf (vel);	

    if (!sf_getfloat("stretch", &st) && !sf_histfloat(in,"stretch",&st)) st=1.;
    /*( stretch=1 Stolt stretch parameter )*/
    if (1. != st) sf_warning("stretch=%g",st);

    if (vel > 0) st = 2.-st;
    a = (1.-1./st);
    b = 1./st;

    if (!sf_getint("extend",&nf)) nf=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("minstr",&minstr)) minstr=0.0;
    /* minimum stretch allowed */

    trace = sf_floatalloc(nw);
    map = fint1_init (nf, nw, mute);

    for (iy = 0; iy < ny; iy++) {
	sf_warning("%d of %d;",iy+1,ny);
	y = y0+iy*dy;
	y *= y;
	for (ix = 0; ix < nx; ix++) {
	    x = x0+ix*dx;
	    x = st*(x*x + y);  

	    sf_floatread(trace,nt,in);
	    for (iw = nt; iw < nw; iw++) { /* pad */
		trace[iw]=0.;
	    }
	    sf_cosft_frw (trace,0,1);

	    fint1_set(map,trace);
	    stretch(map,stolt,nw,dw,0.,nw,dw,0.,trace,minstr);

	    sf_cosft_inv (trace,0,1);
	    sf_floatwrite(trace,nt,out);
	}
    }
    sf_warning(".");

    exit (0);
}

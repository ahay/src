/* Prestack Stolt modeling/migration. */
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

int main(int argc, char* argv[])
{
    fint1 str;
    bool inv, stack, depth;
    int nt,nw,nx,ny,nh,nf, it,iw,ix,iy,ih, iw2;
    float dw,dx,dy,dh, x,y,h,xh, vel, w0, wh, w2, sq;
    float *trace=NULL, *keep=NULL;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* y: modeling, n: migration */

    if (!sf_getbool("depth",&depth)) depth=false;
    /* y: depth migration, n: time migration */

    if (!sf_getbool ("stack",&stack)) stack=true;
    /* if y: stack migrated image */

    if (inv) { /* modelling */
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nx)) nx=1;
	if (!sf_histint(in,"n3",&ny)) ny=1;
	
	if (!sf_getint ("nh",&nh)) sf_error("Need nh=");
	/* number of offsets */

	if (!sf_getfloat ("dh",&dh)) sf_error("Need dh=");
	/* offset sampling */

	dh = 1./(2*(nh-1)*dh);

	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"d3",&dy)) dy=dx;

	sf_putint(out,"n2",nh); sf_putfloat(out,"d2",dh);
	sf_putint(out,"n3",nx); sf_putfloat(out,"d3",dx);
	sf_putint(out,"n4",ny); sf_putfloat(out,"d4",dy);
	sf_putfloat(out,"o4",0.);
    } else { /* migration */
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nh)) nh=1;
	if (!sf_histint(in,"n3",&nx)) nx=1;
	if (!sf_histint(in,"n4",&ny)) ny=1;

	if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"d3",&dx)) sf_error("No d3= in input");
	if (!sf_histfloat(in,"d4",&dy)) dy=dx;

	if (stack) {
	    sf_putint(out,"n2",nx); sf_putfloat(out,"d2",dx);
	    sf_putint(out,"n3",ny); sf_putfloat(out,"d3",dy);
	    sf_putint(out,"n4",1);
	}
    }

    if (!sf_getfloat ("vel",&vel)) sf_error("Need vel=");
    /* constant velocity */

    if (!sf_histfloat(in,"o1",&w0)) w0=0.;
    if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");

    if (!sf_getint ("pad",&nw)) nw=nt;
    /* padding on the time axis */
    nw=2*(nw-1);

    sf_cosft_init(nw);
    dw = 2.*SF_PI/(nw*dw);
    dh *= 2.*SF_PI;
    dx *= 2.*SF_PI;
    dy *= 2.*SF_PI;

    if (depth) {
	if (inv) {
	    dw *= vel;
	} else {
	    dw *= 1./vel;
	}
    } else {
	dh *= vel;
	dx *= vel;
	dy *= vel;
    }

    trace = sf_floatalloc(nw);

    if (stack) keep = sf_floatalloc(nw);

    if (!sf_getint("extend",&nf)) nf=4;
    /* trace extension */

    str = fint1_init(nf,nw,0);

    for (iy = 0; iy < ny; iy++) {
	y = iy*dy;
	y *= y;
	for (ix = 0; ix < nx; ix++) {
	    x = ix*dx;
	    x *= x;

	    if (inv) {
		sf_floatread(trace,nt,in);
		for (it=nt; it < nw; it++) { /* pad */
		    trace[it]=0.;
		}
		
		sf_cosft_frw (trace,0,1);
		fint1_set(str,trace);
	    } else if (stack) {
		for (it=0; it < nw; it++) {
		    keep[it] = 0.;
		}
	    }

	    for (ih = 0; ih < nh; ih++) {
		h = ih * dh;
		h *= h;
		xh = x*h;
		h += x + y;

		if (!inv) {
		    sf_floatread(trace,nt,in);
		    for (it=nt; it < nw; it++) { /* pad */
			trace[it]=0.;
		    }

		    sf_cosft_frw (trace,0,1);
		    fint1_set(str,trace);
		}

		for (iw = 0; iw < nw; iw++) {
		    w2 = iw*dw;
		    w2 *= w2;

		    if (inv) { /* modeling */
			wh = w2-h;
			sq = wh*wh - 4.*xh;

			if (wh > 0. && sq > 0.) {
			    w2 = sqrtf(0.5*(wh + sqrtf (sq)))/dw;
			    iw2 = w2;
			    trace[iw] = (iw2 < nw)? 
				fint1_apply(str,iw2,w2-iw2,false):0.;
			} else {
			    trace[iw] = 0.;
			}
		    } else { /* migration */
			if (w2 == 0.) {
			    trace[iw] = 0.;
			} else {
			    w2 = sqrtf (w2 + h + xh/w2)/dw;
			    iw2 = w2;
			    trace[iw] = (iw2 < nw)? 
				fint1_apply(str,iw2,w2-iw2,false):0.;
			}
		    }
		}

		if (inv || !stack) {
		    sf_cosft_inv (trace,0,1);
		    sf_floatwrite(trace,nt,out);
		} else {
		    for (iw=0; iw < nw; iw++) {
			keep[iw] += trace[iw];
		    }
		}
	    } /* h */
	    if (!inv && stack) {
		sf_cosft_inv (keep,0,1);
		sf_floatwrite(keep,nt,out);
	    }
	} /* x */
    } /* y */

    exit (0);
}

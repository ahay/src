/* Local radial median filtering. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <stdio.h>
#include <math.h>

#include "median.h"

int main (int argc, char* argv[]) 
{
    int n3, i3, nfw, nt, i1, nx, i2, iw, m, tempn1;

    float vmin, vmax, dx, x0, dt, t0, tmax, tmin, tp, xp, vi, xi, ti, tempxi, tempti, nti, s;

    bool boundary;
    
    float *dat, *window, *outp;

    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");


    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input"); 
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input"); 
    
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input"); 
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input"); 

    if (!sf_getint("nfw",&nfw)) sf_error("Need nfw=");
    /* filter window of median filter */

    if (nfw < 1)  sf_error("Need positive integer input"); 
    if (nfw%2 == 0)  nfw = (nfw+1);
    m=(nfw-1)/2;
    
    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/

    if (!sf_getfloat("vmin",&vmin)) sf_error("Need vmin=");
    /* minimum velocity */
    
    if (!sf_getfloat("vmax",&vmax)) sf_error("Need vmax=");
    /* maximum velocity */

    if (!sf_getfloat("tmin",&tmin)) sf_error("Need tmin=");
    /* zero-offset time for maximum velocity */
    
    if (!sf_getfloat("tmax",&tmax)) sf_error("Need tmax=");
    /* zero-offset time for mimumum velocity */

    xp = 0.;
    tp = 0.;

    if (vmin > 0. && vmax > 0. && vmax!=vmin) {
	xp = (tmax - tmin)/(1./vmax - 1./vmin);
	tp = tmax + xp/vmin;
    } else {
	sf_error("Need vmin and vmax > 0.");
    }

    n3 = sf_leftsize(in,2);

    dat = sf_floatalloc(nt*nx);
    window = sf_floatalloc(nfw);
    outp = sf_floatalloc(nt*nx);

    for (i3 = 0; i3 < n3; i3 ++) {
	sf_floatread(dat,nt*nx,in);

	for (i2 = 0; i2 < nx; i2 ++) {
	    for (i1 = 0; i1 < nt; i1 ++) {
		xi = i2*dx + x0;
		ti = i1*dt + t0;
		if (ti != tp) {
		    vi = (xi-xp)/(ti-tp);
		    if (vi ==0.) vi +=1.e-10;
		} else {
		    vi = (xi-xp)/(ti-tp+1.e-20);
		    if (vi ==0.) vi +=1.e-10;
		}
		if (vi <= vmin || vi >= vmax) {
		    outp[nt*i2+i1] = 0.;
		} else {
		    window[m] = dat[nt*i2+i1];
		    for (iw = 1; iw <= m; iw ++) {
			if ((i2+iw) > (nx-1) ) { /* boundary condition */
			    if (boundary) {
				window[m+iw] = window[m+iw-1];
			    } else {
				window[m+iw] = 0.;
			    }
			} else {
			    tempxi = (i2+iw)*dx+x0;
			    tempti = ti+(tempxi-xi)/vi;
			    nti = (tempti-t0)/dt;
			    tempn1 = (int)(nti);
			    s = nti-tempn1;
			    if ((tempn1+1) > (nt-1)) { /* boundary condition */
				if (boundary) {
				    window[m+iw] = window[m+iw-1];
				} else {
				    window[m+iw] = 0.;
				}
			    } else {
				window[m+iw]=s*dat[nt*(i2+iw)+tempn1+1]-(1.-s)*dat[nt*(i2+iw)+tempn1];
			    }
			}
		    }
		    for (iw = 1; iw <= m; iw ++) {
			if ((i2-iw) < 0 ) { /* boundary condition */
			    if (boundary) {
				window[m-iw] = window[m-iw+1];
			    } else {
				window[m-iw] = 0.;
			    }
			} else {
			    tempxi = (i2-iw)*dx+x0;
			    tempti = ti+(tempxi-xi)/vi;
			    nti = (tempti-t0)/dt;
			    tempn1 = (int)(nti);
			    s = nti-tempn1;
			    if (tempn1 < 0) { /* boundary condition */
				if (boundary) {
				    window[m-iw] = window[m-iw+1];
				} else {
				    window[m-iw] = 0.;
				}
			    } else {
				window[m-iw]=s*dat[nt*(i2-iw)+tempn1+1]-(1.-s)*dat[nt*(i2-iw)+tempn1];
			    }
			}
		    }
		    outp[nt*i2+i1] = medianfilter(window,nfw);
		}
	    }
	}
	sf_floatwrite(outp,nt*nx,out);
    }

    exit (0);
}

/* 	$Id$	 */

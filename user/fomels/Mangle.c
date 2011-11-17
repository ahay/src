/* Illustration of angle gathers.
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

int main(int argc, char* argv[])
{
    int iw, ih, im, nw, nh, nm;
    float w0, dw, dh, dm, w, m, h, vel, *trace;
    sf_file angle;
    
    sf_init(argc,argv);
    angle=sf_output("out");
    sf_setformat(angle,"native_float");

    if (!sf_getint("nw",&nw)) nw=513;
    if (!sf_getint("nm",&nm)) nm=257;
    if (!sf_getint("nh",&nh)) nh=257;

    if (!sf_getfloat("dw",&dw)) dw=1./(2*(nw-1)*0.004);
    if (!sf_getfloat("dm",&dm)) dm=1./(2*(nm-1)*0.01);
    if (!sf_getfloat("dh",&dh)) dh=1./(2*(nh-1)*0.01);

    if (!sf_getfloat("w0",&w0)) w0=dw; 

    if (!sf_getfloat("vel",&vel)) vel=2.;

    sf_putint(angle,"n1",nh);
    sf_putfloat(angle,"o1",0.);
    sf_putfloat(angle,"d1",dh);
    sf_putstring(angle,"label1","Offset Wavenumber");
    sf_putstring(angle,"unit1","1/km");

    sf_putint(angle,"n2",nm);
    sf_putfloat(angle,"o2",0.);
    sf_putfloat(angle,"d2",dm);
    sf_putstring(angle,"label2","Midpoint Wavenumber");
    sf_putstring(angle,"unit2","1/km");

    sf_putint(angle,"n3",nw);
    sf_putfloat(angle,"o3",w0);
    sf_putfloat(angle,"d3",dw);
    sf_putstring(angle,"label3","Frequency");
    sf_putstring(angle,"unit3","1/s");

    w0 *= 2.*SF_PI/vel;
    dw *= 2.*SF_PI/vel;
    dm *= 2.*SF_PI;
    dh *= 2.*SF_PI;

    trace = sf_floatalloc(nh);

    for (iw=0; iw < nw; iw++) {
	sf_warning("frequency %d of %d;",iw,nw-1); 
	w = w0+iw*dw;
	w *= 4.*w;
	for (im=0; im < nm; im++) {
	    m = im*dm;
	    m *= m;
	    for (ih=0; ih < nh; ih++) {
		h = ih*dh;
		h *= h;
		h = w + m - h;
		if (h*h >= 4.*m*w) {
		    h = (h + sqrtf(h*h - 4.*m*w))/(2.*w);
		    if (h >= 0. && h <= 1.) {
			h = acosf(sqrtf(h))*180./SF_PI;		
			trace[ih] = h;
		    } else {
			trace[ih] = -dh;
		    }
		} else {
		    trace[ih] = -dh;
		}
	    }
	    sf_floatwrite(trace,nh,angle);
	}
    }
    sf_warning(".");

    exit(0);
}

/* 	$Id$	 */


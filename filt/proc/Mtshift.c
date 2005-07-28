/* Compute shift from pseudo-v to pseudo-tan(theta) */
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

int main (int argc, char* argv[])
{
    fint1 sft;
    int it,ix,iz,ih,ia, nt,nx,nw, nh,na;
    float h0, f, dh, v, da, a0;
    float **gath, **gath2, *vel, *trace;
    sf_file in, out, velocity;

    sf_init (argc,argv);
    in       = sf_input("in");
    velocity = sf_input("velocity");
    out      = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint  (in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint  (in,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");

    if (!sf_getint("na",&na)) na=nh;
    /* tangent samples */
    if (!sf_getfloat("da",&da)) da=1./(nh-1);
    /* tangent sampling */
    if (!sf_getfloat("a0",&a0)) a0=0.;
    /* tangent origin */

/*    sf_putfloat(out,"o2",0.);*/

    sf_putfloat(out,"o2",a0);
    sf_putfloat(out,"d2",da);
    sf_putint  (out,"n2",na);

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    gath = sf_floatalloc2(nt,nh);
    gath2 = sf_floatalloc2(nt,na);
    trace = sf_floatalloc(nh);
    vel   = sf_floatalloc(nt);

    sft = fint1_init (nw, nh);
    
    for (ix = 0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);	
	sf_floatread (gath[0],nt*nh,in);

	for (it = 0; it < nt; it++) {
	    for (ih = 0; ih < nh; ih++) {
		trace[ih] = gath[ih][it];
	    }
	    fint1_set(sft,trace);
	    v = vel[it];
	    
	    for (ia=0; ia < na; ia++) {
		f = (v*hypotf(a0+ia*da,1.) - h0)/dh;
		iz = f;
		if (iz >= 0 && iz < nh) {
		    gath2[ia][it] = fint1_apply(sft,iz,f-iz,false);
		} else {
		    gath2[ia][it] = 0.;
		}
	    }
	}
	    
	sf_floatwrite (gath2[0],nt*na,out);
    }
	
    exit (0);
}

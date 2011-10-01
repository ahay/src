/* Normal moveout in TXP domain. */
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

#include "fint1.h"

float x, p, *vel;

static float nmo_map(float t, int it) {
    float v;
    
    v = vel[it];
    v *= v;

    t = t+p*x*x/(x+t*p*v+SF_EPS);
    
    return t;
}

int main (int argc, char* argv[])
{
    fint1 nmo;
    int nt, ix, nx, ip, np, mute, nw;
    float t0, dt, x0, dx, p0, dp, str, *trace;
    sf_file cmp, nmod, velocity;

    sf_init (argc,argv);
    cmp      = sf_input("in");
    velocity = sf_input("velocity");
    nmod     = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint  (cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint  (cmp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_histint  (cmp,"n3",&np)) sf_error("No n3= in input");
    if (!sf_histfloat(cmp,"d3",&dp)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&p0)) sf_error("No o3= in input");

    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    nmo = fint1_init (nw, nt, mute);

    trace = sf_floatalloc(nt);
    vel   = sf_floatalloc(nt);

    sf_floatread (vel,nt,velocity);

    for (ip = 0; ip < np; ip++) {
	sf_warning("slope %d of %d;",ip+1,np);

	p = p0+ip*dp;

	for (ix = 0; ix < nx; ix++) {
	    x = x0+ix*dx;

	    sf_floatread (trace,nt,cmp);

	    fint1_set(nmo,trace);

	    stretch(nmo,nmo_map,nt,dt,t0,nt,dt,t0,trace,str);
	    
	    sf_floatwrite (trace,nt,nmod);
	}
    }
    
    sf_warning(".");

    exit (0);
}

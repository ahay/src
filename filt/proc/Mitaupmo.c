/* Inverse normal moveout in tau-p domain. */
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

#include "stretch4.h"

int main (int argc, char* argv[])
{
    map4 nmo; /* using cubic spline interpolation */
    int it,ix,ip, nt,nx, np, nw;
    float dt, t0, p, p0, f, ft, dp, eps;
    float *trace, *vel, *str, *out;
    sf_file cmp, nmod, velocity;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&p0)) sf_error("No o2= in input");

    nx = sf_leftsize(cmp,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = stretch4_init (nt, t0, dt, nt, eps);
    
    for (ix = 0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);	

	for (ip = 0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    p *= p;

	    sf_floatread (trace,nt,cmp);
	    
	    f = 0.;
	    for (it=0; it < nt; it++) {
		str[it] = t0+f*dt;

		ft = vel[it];
		ft = 1.-p*ft*ft;

		if (ft < 0.) {
		    for (; it < nt; it++) {
			str[it]=t0-10.*dt;			
		    }
		    break;
		}

		f += sqrtf(ft);
	    }

	    stretch4_define (nmo,str);
	    stretch4_apply (nmo,trace,out);
	    
	    sf_floatwrite (out,nt,nmod);
	}
    }
	
    exit (0);
}



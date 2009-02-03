/* Slope-based tau-p moveout. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <float.h>

#include <rsf.h>

#include "stretch4.h"

int main (int argc, char* argv[])
{
    map4 nmo;
    int it,ix,ip, nt,nx, np;
    float dt, t0, p, p0, t, f, dp, eps;
    float *trace, *slope, *dsldt, *str, *cos;
    sf_file inp, nmod, cos2, dip, dipt;

    sf_init (argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    dipt = sf_input("dipt");
    nmod = sf_output("out");
    cos2 = sf_output("cos2");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(inp,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&p0)) sf_error("No o2= in input");
    p0 /= dp;

    nx = sf_leftsize(inp,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    slope = sf_floatalloc(nt);
    dsldt = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    cos = sf_floatalloc(nt);

    nmo = stretch4_init (nt, t0, dt, nt, eps);

    eps = 100.*FLT_EPSILON;
    
    for (ix = 0; ix < nx; ix++) { /* midpoints */
	for (ip = 0; ip < np; ip++) { /* offset */
	    p = p0+ip;

	    sf_floatread (trace,nt,inp);  
	    sf_floatread (slope,nt,dip);  
	    sf_floatread (dsldt,nt,dipt); 

	    for (it=0; it < nt; it++) { /* time */
		t = t0 + it*dt;

		f = t - p*slope[it]*dt; 
		    
		if (f < 0. || f < t) {
		    str[it] = t0-10.*dt;
		    cos[it] = 0.;
		} else {
		    str[it] = sqrtf(t*f); /* t -> tau */
		    cos[it] = 4*t/(t+f-p*t*dsldt[it]+eps)-t/(f+eps);
		}
	    }

	    stretch4_define (nmo,str);
	    
	    stretch4_apply (nmo,trace,trace);	    
	    sf_floatwrite (trace,nt,nmod);
	    
	    stretch4_apply (nmo,cos,cos);	    
	    sf_floatwrite (cos,nt,cos2);
	}
    }
    
    exit (0);
}

/* 	$Id$	 */

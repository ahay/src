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
#include <rsf.h>

int main (int argc, char* argv[])
{
    sf_map4 nmo;
    char *type;
    int it,ix,ip, nt,nx, np;
    float dt, t0, p, p0, t, f, dp, eps, v0, c0, cos2, dp1;
    float *trace=NULL, *slope=NULL, *dsldt=NULL, *str=NULL, *vel=NULL;
    sf_file inp=NULL, nmod=NULL, vel2=NULL, dip=NULL, dipt=NULL;

    sf_init (argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    dipt = sf_input("dipt");
    nmod = sf_output("out");
    vel2 = sf_output("vel2");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(inp,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&p0)) sf_error("No o2= in input");
    p0 /= dp;
    dp1 = 1./(dp*dp);

    nx = sf_leftsize(inp,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    if (!sf_getfloat("v0",&v0)) v0=0.;
    /* initial velocity */
    v0 *= dp;

    if (NULL == (type = sf_getstring("type"))) type="dix";
    /* transform type */

    trace = sf_floatalloc(nt);
    slope = sf_floatalloc(nt);
    dsldt = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);

    nmo = sf_stretch4_init (nt, t0, dt, nt, eps);

    eps = 100.*SF_EPS;

    for (ix = 0; ix < nx; ix++) { /* midpoints */
	for (ip = 0; ip < np; ip++) { /* offset */
	    p = p0+ip;
	    c0 = p*v0;
	    c0 = 1.-c0*c0;

	    sf_floatread (trace,nt,inp);
	    sf_floatread (slope,nt,dip);
	    sf_floatread (dsldt,nt,dipt);

	    for (it=0; it < nt; it++) { /* time */
		t = t0 + it*dt;

		f = t - p*slope[it]*dt*c0;

		if (f < 0. || f < t) {
		    str[it] = t0-10.*dt;
		    vel[it] = 0.;
		} else {
		    str[it] = sqrtf(t*f); /* t -> tau */
		    switch (type[0]) {
			case 'd':
			    cos2 = t*c0*(4.0/(t+f-c0*p*t*dsldt[it]+eps)-1.0/
					 (f+eps));
			    vel[it] = (1-cos2)/(p*p+eps);
			    break;
			case 'l':
			    vel[it] = dsldt[it]/(p*(p*dsldt[it]-1.0)+eps);
			    break;
			default:
			    sf_error("Unknown type \"%s\"",type);
			    break;
		    }
		    vel[it] *= dp1;
		}
	    }

	    sf_stretch4_define (nmo,str);

	    sf_stretch4_apply (false,nmo,trace,trace);
	    sf_floatwrite (trace,nt,nmod);

	    sf_stretch4_apply (false,nmo,vel,vel);
	    sf_floatwrite (vel,nt,vel2);
	}
    }

    exit (0);
}

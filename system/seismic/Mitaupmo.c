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
    bool interval;
    map4 nmo; /* using cubic spline interpolation */
    int it,ix,ip, nt,nx, np;
    float dt, t0, p, p0, f, ft, dp, eps, at;
    float *trace=NULL, *vel=NULL, *str=NULL, *out=NULL, *n=NULL;
    sf_file cmp=NULL, nmod=NULL, velocity=NULL, eta=NULL;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (!sf_getbool("interval",&interval)) interval=true;
    /* use interval velocity */

    if (NULL == sf_getstring("eta")) {
	eta = NULL;
    } else {
	eta = sf_input("eta");
    }

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

    n = (NULL == eta)? NULL: sf_floatalloc(nt);

    nmo = stretch4_init (nt, t0, dt, nt, eps);

    for (ix = 0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);	
	if (NULL != eta) sf_floatread (n,nt,eta);	

	for (ip = 0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    p *= p;

	    sf_floatread (trace,nt,cmp);
	    
	    f = 0.;
	    for (it=0; it < nt; it++) {
		ft = vel[it];
		ft *= ft;
		
		if (NULL == n) {
		    ft = 1.-p*ft;
		} else {
		    at = n[it];
		    ft = (1.-(1.+2.*at)*p*ft)/(1.-2.*at*p*ft);
		}

		if (ft < 0.) {
		    for (; it < nt; it++) {
			str[it]=t0-10.*dt;
		    }
		    break;
		}

		if (interval) {
		    if (it==0) f = sqrtf(ft)*t0/dt;
		    str[it] = f*dt;
		    f += sqrtf(ft);
		} else {		    
		    f = sqrtf(ft);
		    str[it] = f*(t0+it*dt);
		}
	    }

	    stretch4_define (nmo,str);
	    stretch4_apply (nmo,trace,out);

	    sf_floatwrite (out,nt,nmod);
	}
    }

    exit (0);
}

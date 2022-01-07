/* Slope-based tau-p 3D moveout. */
/*
  Copyright (C) 2010 Politecnico di Milano

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
    int it,ix,ip1, ip2,nt,nx, np1, np2;
    float dt, t0, t, eps, f, p1, p10, dp1, p2, p20, dp2;
    float *trace=NULL, *slope1=NULL, *slope2=NULL, *str=NULL;
    sf_file inp=NULL, nmod=NULL, dip1=NULL, dip2=NULL ,tau0=NULL;

    sf_init (argc,argv);
    inp = sf_input("in");
    if (NULL != sf_getstring("dip1")) dip1 = sf_input("dip1"); /*slope field mesaure along dimension 2*/
    else sf_error("Need dip1 input");
    if (NULL != sf_getstring("dip2")) dip2 = sf_input("dip2"); /*slope field mesaure along dimension 3*/
    else sf_error("Need dip2 input");
    nmod = sf_output("out");
    if (NULL != sf_getstring("tau0")) tau0 = sf_output("tau0"); /*tau0(tau,p) */
    else tau0 = NULL ;
    

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(inp,"n2",&np1)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dp1)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&p10)) sf_error("No o2= in input");

    if (!sf_histint(inp,"n3",&np2)) sf_error("No n3= in input");
    if (!sf_histfloat(inp,"d3",&dp2)) sf_error("No d3= in input");
    if (!sf_histfloat(inp,"o3",&p20)) sf_error("No o3= in input");

    p10 /= dp1; /* slope axis is normalized: no need for slope normalization dp1=1 and dp2=1*/
    p20 /= dp2;

    nx = sf_leftsize(inp,3);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    slope1 = sf_floatalloc(nt);
    slope2 = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    /*cos = sf_floatalloc(nt);*/

    nmo = sf_stretch4_init (nt, t0, dt, nt, eps);

    eps = SF_EPS;

    for (ix = 0; ix < nx; ix++) { /* midpoints */

	for (ip2 = 0; ip2 < np2; ip2++) { /* slope 2 (dimension 3)*/
	    p2 = p20+ip2;

	    for (ip1 = 0; ip1 < np1; ip1++) { /* slope 1 (dimension 2)*/
		p1 = p10+ip1;

		sf_floatread (trace,nt,inp);
		sf_floatread (slope1,nt,dip1);
		sf_floatread (slope2,nt,dip2);

		for (it=0; it < nt; it++) { /* time */
		    t = t0 + it*dt;
		
		    f = t - p1*slope1[it]*dt - p2*slope2[it]*dt;

		    if (f < 0. || f < t) {
			str[it] = t0-10.*dt;
		    } else {
			str[it] = sqrtf(t*f); /* t -> tau */
		    }

		}

		sf_stretch4_define (nmo,str,false);

		sf_stretch4_apply (false,nmo,trace,trace);
		sf_floatwrite (trace,nt,nmod);
		sf_floatwrite (str,nt,tau0);

	    }
	}
    }

    exit (0);
}

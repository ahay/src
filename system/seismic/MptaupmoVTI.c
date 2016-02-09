/* Slope-based tau-p moveout in VTI. */
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
    int it,ix,ip, nt,nx, np;
    float dt, t0, p, p0, t, dp, eps, f; 
    float N,D;
    float *trace=NULL, *slope=NULL, *curv=NULL, *str=NULL;
    sf_file inp=NULL, nmod=NULL, dip=NULL, ddip=NULL ,tau0=NULL;

    sf_init (argc,argv);
    inp = sf_input("in");
    if (NULL != sf_getstring("dip")) dip = sf_input("dip"); /*slope field */
	else dip = NULL ;
	if (NULL != sf_getstring("ddip")) ddip = sf_input("ddip"); /*curvature field */
	else ddip = NULL ;
    nmod = sf_output("out");
	if (NULL != sf_getstring("tau0")) tau0 = sf_output("tau0"); /*tau0(tau,p) */
	else tau0 = NULL ;
    
    /*cos2 = sf_output("cos2");*/

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(inp,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&p0)) sf_error("No o2= in input");
    p0 /= dp; /* slope axis is normalized: no need for slope normalization dp=1*/

    nx = sf_leftsize(inp,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    slope = sf_floatalloc(nt);
    curv = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    /*cos = sf_floatalloc(nt);*/

    nmo = sf_stretch4_init (nt, t0, dt, nt, eps);

    eps = SF_EPS;

    for (ix = 0; ix < nx; ix++) { /* midpoints */
	for (ip = 0; ip < np; ip++) { /* slope */
	    p = p0+ip;

	    sf_floatread (trace,nt,inp);
	    sf_floatread (slope,nt,dip);
	    sf_floatread (curv,nt,ddip);

	    for (it=0; it < nt; it++) { /* time */
		t = t0 + it*dt;
		
	/*N = 3*TAU*R(irow,icol)+TAU*P*Q(irow,icol)-3*R(irow,icol)^2*P;
        D = 3*TAU*R(irow,icol)+TAU*P*Q(irow,icol)+R(irow,icol)^2*P;*/
		N=3*t*slope[it]*(dt)+t*p*curv[it]*(dt)-3*slope[it]*slope[it]*(dt*dt)*p;
		D=3*t*slope[it]*(dt)+t*p*curv[it]*(dt)+slope[it]*slope[it]*(dt*dt)*p;
		/*f = t - p*slope[it]*dt;*/
		f = N/D; 
		
		if (f < 0.) {
		    str[it] = t0-10.*dt;
		} else { 
		    str[it] = t*sqrtf(f); /* t -> tau */
		} 
	    }

	    sf_stretch4_define (nmo,str);

	    sf_stretch4_apply (false,nmo,trace,trace);
	    sf_floatwrite (trace,nt,nmod);
	    sf_floatwrite (str,nt,tau0);	
	    /*sf_stretch4_apply (nmo,cos,cos);
	    sf_floatwrite (cos,nt,cos2);*/
	}
    }

    exit (0);
}

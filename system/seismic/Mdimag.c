/* Diffraction imaging in the plane-wave domain. */
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

static float p, s, v, x0;

static float diffrac(float t, int it) 
{
    return (hypotf(0.5*t,(s-x0)/v)+0.5*t*sqrtf(1-p*p*v*v)+p*(s-x0));
}

int main(int argc, char* argv[])
{
    fint1 nmo;
    bool sembl;
    int it,ip,is,iv, nt,np,ns,nv, ib,ie,nb,i, nw, mute;
    float amp, dt, dp, t0, p0, v0, dv, num, den, str, s0, ds;
    float *trace=NULL, **stack=NULL, **stack2=NULL;
    sf_file cmp=NULL, scan=NULL;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&np)) sf_error("No n2= in input");

    ns = sf_leftsize(cmp,2);
    if (!sf_histfloat(cmp,"o3",&s0)) sf_error("No o3= in input");
    if (!sf_histfloat(cmp,"d3",&ds)) sf_error("No d3= in input");

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
     if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");

    if (!sf_histfloat(cmp,"o2",&p0)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d2",&dp)) sf_error("No d2= in input");

    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    if (!sf_getint("nv",&nv))   sf_error("Need nv=");

    sf_putfloat(scan,"o2",v0);
    sf_putfloat(scan,"d2",dv);
    sf_putint(scan,"n2",nv);

    stack =  sf_floatalloc2(nt,nv);
    stack2 = sembl? sf_floatalloc2(nt,nv) : NULL;

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("str",&str)) str=0.;
    /* maximum stretch allowed */

    if (!sf_getfloat("x0",&x0)) sf_error("Need x0=");

    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt,mute);

    for (is=0; is < ns; is++) {
	sf_warning("shot %d of %d",is+1,ns);
	s = s0 + is * ds;

	for (it=0; it < nt*nv; it++) {
	    stack[0][it] = 0.;
	    if (sembl) stack2[0][it] = 0.;
	}

	for (ip=0; ip < np; ip++) {
	    sf_floatread(trace,nt,cmp); 
	    p = p0 + ip * dp;

	    for (it=0; it < nt; it++) {
		trace[it] /= nt*np;
	    }
	    fint1_set(nmo,trace);

	    for (iv=0; iv < nv; iv++) {
		v = v0 + iv * dv;

		stretch(nmo,diffrac,nt,dt,t0,nt,dt,t0,trace,str);

		for (it=0; it < nt; it++) {
		    amp = trace[it];
		    if (sembl) stack2[iv][it] += amp*amp;
		    stack[iv][it] += amp;
		}
	    } /* v */
	} /* P */
	
	if (sembl) {
	    for (iv=0; iv < nv; iv++) {
		for (it=0; it < nt; it++) {
		    ib = it-nb;
		    ie = it+nb+1;
		    if (ib < 0) ib=0;
		    if (ie > nt) ie=nt;
		    num = 0.;
		    den = 0.;
		    for (i=ib; i < ie; i++) {
			num += stack[iv][i]*stack[iv][i];
			den += stack2[iv][i];
		    }
		    den *= np;
		    trace[it] = (den > 0.)? num/den: 0.;
		}
		sf_floatwrite(trace,nt,scan);
	    }
	} else {
	    sf_floatwrite (stack[0],nt*nv,scan);
	}
    } /* s */

    exit(0);
}

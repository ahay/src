/* Velocity analysis.

Inverse of sfvelmod
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

#include "fint1.h"

static float v;

static float hyperb(float t) 
{ 
    return hypotf(t,v); 
} 

int main(int argc, char* argv[])
{
    fint1 nmo;
    bool sembl, half, slow;
    int it,ih,ix,iv, nt,nh,nx,nv, ib,ie,nb,i, nw, CDPtype;
    float amp, dt, dh, t0, h0, v0, dv, h, num, den, dy;
    float *trace, **stack, **stack2;
    sf_file cmp, scan;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
    if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");

    if (!sf_getfloat("v0",&v0) && !sf_histfloat(cmp,"v0",&v0)) 
	sf_error("Need v0=");
    if (!sf_getfloat("dv",&dv) && !sf_histfloat(cmp,"dv",&dv)) 
	sf_error("Need dv=");
    if (!sf_getint("nv",&nv) && !sf_histint(cmp,"nv",&nv)) 
	sf_error("Need nv=");

    sf_putfloat(scan,"o2",v0);
    sf_putfloat(scan,"d2",dv);
    sf_putint(scan,"n2",nv);

    sf_putfloat(scan,"h0",h0);
    sf_putfloat(scan,"dh",dh);
    sf_putint(scan,"nh",nh);

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    if (half) {
	dh *= 2.;
	h0 *= 2.;
    }

    CDPtype=1;
    if (sf_histfloat(cmp,"d3",&dy)) {
	CDPtype=0.5+0.5*dh/dy;
	if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
    } 	    
    sf_warning("CDPtype=%d",CDPtype);

    sf_putstring(scan,"label2",slow? "slowness": "velocity");

    trace = sf_floatalloc(nt);
    stack =  sf_floatalloc2(nt,nv);
    stack2 = sembl? sf_floatalloc2(nt,nv) : NULL;

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt);

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */

    for (ix=0; ix < nx; ix++) {
	sf_warning("cmp %d of %d",ix+1,nx);

	for (it=0; it < nt*nv; it++) {
	    stack[0][it] = 0.;
	    if (sembl) stack2[0][it] = 0.;
	}

	for (ih=0; ih < nh; ih++) {
	    h = h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
	    sf_floatread(trace,nt,cmp); 

	    for (it=0; it < nt; it++) {
		trace[it] /= nt*nh;
	    }
	    fint1_set(nmo,trace);

	    for (iv=0; iv < nv; iv++) {
		v = v0 + iv * dv;
		v = slow? h*v: h/v;

		stretch(nmo,hyperb,nt,dt,t0,nt,dt,t0,trace);

		for (it=0; it < nt; it++) {
		    amp = fabsf(v)*trace[it];
		    stack[iv][it] += amp;
		    if (sembl) stack2[iv][it] += amp*amp;
		}
	    } /* v */
	} /* h */
	
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
		    trace[it] = (den > 0.)? num/den: 0.;
		}
		sf_floatwrite(trace,nt,scan);
	    }
	} else {
	    sf_floatwrite (stack[0],nt*nv,scan);
	}
    } /* x */

    sf_close();
    exit(0);
}

/* 	$Id: Mvscan.c,v 1.6 2004/06/25 18:08:43 fomels Exp $	 */

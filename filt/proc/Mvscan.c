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

static float hyperb(float t, int it) 
{ 
    return hypotf(t,v); 
} 

int main(int argc, char* argv[])
{
    fint1 nmo;
    bool sembl, half, slow, dsembl, asembl, weight;
    int it,ih,ix,iv, nt,nh,nx,nv, ib,ie,nb,i, nw, CDPtype, mute, *mask;
    float amp, amp2, dt, dh, t0, h0, v0, dv, h, num, den, dy, str, sh=0., sh2=0.;
    float *trace, **stack, **stack2, **stackh, *hh;
    char *time, *space, *unit;
    size_t len;
    sf_file cmp, scan, offset, msk;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
    if (sembl || !sf_getbool("diffsemblance",&dsembl)) dsembl=false;
    /* if y, compute differential semblance */
    if (sembl || dsembl || !sf_getbool("avosemblance",&asembl)) asembl=false;
    /* if y, compute AVO-friendly semblance */
    if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */
    if (!sf_getbool("weight",&weight)) weight=true;
    /* if y, apply pseudo-unitary weighting */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	hh = sf_floatalloc(nh);

	h0 = dh = 0.;
    } else {
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	
	sf_putfloat(scan,"h0",h0);
	sf_putfloat(scan,"dh",dh);
	sf_putint(scan,"nh",nh);

	if (sf_histfloat(cmp,"d3",&dy)) {
	    CDPtype=half? 0.5+dh/dy: 0.5+0.5*dh/dy;
	    if (0 == CDPtype) CDPtype=1;
	    if (1 != CDPtype) {
		sf_histint(cmp,"CDPtype",&CDPtype);
		sf_warning("CDPtype=%d",CDPtype);
	    }
	}

	offset = NULL;
	hh = NULL;
    }

    if (NULL != sf_getstring("mask")) {
	msk = sf_input("mask");
	mask = sf_intalloc(nh);
    } else {
	msk = NULL;
	mask = NULL;
    }

    if (!sf_getfloat("v0",&v0) && !sf_histfloat(cmp,"v0",&v0)) 
	sf_error("Need v0=");
    if (!sf_getfloat("dv",&dv) && !sf_histfloat(cmp,"dv",&dv)) 
	sf_error("Need dv=");
    if (!sf_getint("nv",&nv) && !sf_histint(cmp,"nv",&nv)) 
	sf_error("Need nv=");

    sf_putfloat(scan,"o2",v0);
    sf_putfloat(scan,"d2",dv);
    sf_putint(scan,"n2",nv);

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */
    sf_putstring(scan,"label2",slow? "slowness": "velocity");

    if (NULL != (time = sf_histstring(cmp,"unit1")) &&
	NULL != (space = sf_histstring(cmp,"unit2"))) {
	len = strlen(time)+strlen(space)+2;
	unit = sf_charalloc(len);
	if (slow) {
	    snprintf(unit,len,"%s/%s",time,space);
	} else {
	    snprintf(unit,len,"%s/%s",space,time);
	}
	sf_putstring(scan,"unit2",unit);
    }

    stack =  sf_floatalloc2(nt,nv);
    stack2 = (sembl || dsembl || asembl)? sf_floatalloc2(nt,nv) : NULL;
    stackh = asembl? sf_floatalloc2(nt,nv) : NULL;

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */

    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt,mute);

    for (ix=0; ix < nx; ix++) {
	sf_warning("cmp %d of %d",ix+1,nx);

	for (it=0; it < nt*nv; it++) {
	    stack[0][it] = 0.;
	    if (sembl || asembl) stack2[0][it] = 0.;
	    if (asembl) stackh[0][it] = 0.;
	}

	if (NULL != offset) sf_floatread(hh,nh,offset);
	if (NULL != msk) sf_intread(mask,nh,msk);

	if (asembl) sh = sh2 = 0.;

	for (ih=0; ih < nh; ih++) {
	    sf_floatread(trace,nt,cmp); 
	    if (NULL != msk && 0==mask[ih]) continue;

	    h = (NULL != offset)? hh[ih]: 
		h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
	    if (half) h *= 2.;

	    if (asembl) {
		sh  += h;
		sh2 += h*h;
	    }

	    for (it=0; it < nt; it++) {
		trace[it] /= nt*nh;
	    }
	    fint1_set(nmo,trace);
	    
	    for (iv=0; iv < nv; iv++) {
		v = v0 + iv * dv;
		v = slow? h*v: h/v;

		stretch(nmo,hyperb,nt,dt,t0,nt,dt,t0,trace,str);

		for (it=0; it < nt; it++) {
		    amp = weight? fabsf(v)*trace[it]: trace[it];
		    if (dsembl) {
			if (ih > 0) {
			    amp2 = amp - stack2[iv][it];
			    stack[iv][it] += amp2*amp2;
			}
			stack2[iv][it] = amp;
		    } else {
			if (sembl || asembl) stack2[iv][it] += amp*amp;
			if (asembl) stackh[iv][it] += amp*h;

			stack[iv][it] += amp;
		    }
		}
	    } /* v */
	} /* h */
	
	if (sembl || asembl) {
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
			
			/* (h2 s^2 - 2 h s sh + N sh^2)/((-h^2 + h2 N) s2) */
			if (asembl) 
			    num += (nh*stackh[iv][i]*stackh[iv][i] - 
				    2.*sh*stack[iv][i]*stackh[iv][i])/sh2;
		    }
		    if (asembl) den *= (nh-sh*sh/sh2);
		    trace[it] = (den > 0.)? num/den: 0.;
		}
		sf_floatwrite(trace,nt,scan);
	    }
	} else {
	    sf_floatwrite (stack[0],nt*nv,scan);
	}
    } /* x */

    exit(0);
}

/* 	$Id$	 */

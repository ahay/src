/* Stolt stretch. */
/*
  Copyright (C) 2007 University of Texas at Austin

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
#include "vt2w.h"

static int ix;
static float **str;
static float map(float t, int it) { return str[ix][it]; }

int main(int argc, char* argv[])
{
    int nt, ns, nx, n2, it, nstr, ifix;
    float dt, t0, eps, v0, *v=NULL, *trace=NULL, *out=NULL, *ww=NULL, wsum;
    char buffer[20];
    sf_map4 stolt;
    fint1 istolt;
    bool inv;
    sf_file in=NULL, st=NULL, vel=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    st = sf_output("out");
    vel = sf_input("velocity");

    if (!sf_getbool ("inv", &inv)) inv=false;
    /* if y, inverse stretch */

    if (!sf_getint("nstretch",&nstr)) nstr=1;
    /* number of steps */

    if (inv) {
		if (!sf_histint(in,"n1",&ns)) sf_error("No n1= in input");
		if (!sf_histint(vel,"n1",&nt)) sf_error("No n1= in velocity");
		sf_putint(st,"n1",nt);
    } else {
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_getint("pad",&ns)) ns=nt;
	/* time axis padding */
	sf_putint(st,"n1",ns);
	sf_putint(st,"nstretch", nstr);
    }
    nx = sf_leftsize(in,1);
    n2 = sf_leftsize(vel,1);

    if (n2 != 1 && n2 != nx) sf_error("Wrong number of traces in velocity");

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0) || 0. != t0) sf_error("Need o1=0 in input");
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */
    if (!sf_getfloat("vel",&v0)) sf_error("Need vel=");
    /* reference velocity */

    trace = sf_floatalloc(nt);
    v = sf_floatalloc(nt);
    str = sf_floatalloc2(nt,nx);
    out = sf_floatalloc(ns);
    ww = sf_floatalloc(nx);

    if (inv) {
	stolt = NULL;
	istolt = fint1_init(4,ns,0);
    } else {
	stolt =  sf_stretch4_init (ns,t0,dt,nt,eps);
	istolt = NULL;
    }

    wsum = 0.;
    for (ix=0; ix < nx; ix++) {
	if (0==ix || n2 > 1) {
	    sf_floatread(v,nt,vel);
	    ww[ix] = vt2w (nt,v,str[ix]);
	} else if (1 == n2) {
	    ww[ix] = ww[0];
	    for (it=0; it < nt; it++) {
		str[ix][it] = str[0][it];
	    } 
	}
	wsum += ww[ix];
    }

    sf_warning("%d traces processed",nx);

    if (1==nstr) {
	sf_putfloat(st,"stretch", wsum/nx);
    } else {
	for (ix=0; ix < nstr; ix++) {
	    ifix=SF_MAX(1,(int) (nx*ix/(nstr-1)));
	    sf_warning("%d %d",ix,ifix);

	    snprintf(buffer,20,"stretch%d",ix);
	    sf_putfloat (st,buffer,ww[ifix]);

	    snprintf(buffer,20,"node%d",ix);
	    sf_putint (st,buffer,ifix);
	}
    }

    for (ix=0; ix < nx; ix++) {
	for (it=0; it < nt; it++) {
	    str[ix][it] *= dt/v0;
	}
	if (inv) {
	    sf_floatread(out,ns,in);
	    fint1_set(istolt,out);	    
	    stretch(istolt, map, ns, dt, t0, nt, dt, t0, trace, 0.);
	    sf_floatwrite (trace,nt,st);
	} else {
	    sf_floatread (trace,nt,in);
	    sf_stretch4_define (stolt, str[ix],false);
	    sf_stretch4_apply (false,stolt, trace, out);
	    sf_floatwrite (out,ns,st);
	}
    }

    exit(0);
}

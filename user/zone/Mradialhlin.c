/* Radial transform.with shifted-linear interpolation*/
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

#include <rsf.h>
#include "shprefilter.h"

int main(int argc, char* argv[])
{
    bool inv;
    int nt, nx, nv, it, ix, iv, nw, n3, i3, ntr, ntm, im;
    float *trace=NULL, *modl=NULL, *r=NULL;
    float vmin, vmax, dv, dx, x0, t0, t, dt, tp, xp;
    char *unit=NULL, *space=NULL, *time=NULL, *intp=NULL;
    size_t len;
    sf_bands spl=NULL;
    sf_file in, out;
    const float tau=0.21;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getint("nw",&nw)) nw=2;
    /* accuracy level */

    if (!sf_histint(in,"n2",&nt)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dt)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&t0)) sf_error("No o2= in input");
 
    if (!sf_getfloat("tp",&tp)) tp=t0;
    if (!sf_getfloat("xp",&xp)) xp=0.;
     
     intp = sf_getstring("interp");
     if (NULL == intp) intp[0]='s' ;

    if (inv) {
	if (!sf_histint(in,"n1",&nv)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dv)) sf_error("No d1= in input"); 
	if (!sf_histfloat(in,"o1",&vmin)) sf_error("No o1= in input"); 

	if (!sf_getint("nx",&nx) && !sf_histint(in,"nx",&nx)) 
	    sf_error("Need nx=");
	/* number of offsets (if inv=y) */

	if (!sf_getfloat("x0",&x0) && !sf_histfloat(in,"x0",&x0)) 
	    sf_error("Need x0=");
	/* offset origin (if inv=y) */

	if (!sf_getfloat("dx",&dx) && !sf_histfloat(in,"dx",&dx)) 
	    sf_error("Need dx=");
	/* offset sampling (if inv=y) */

	sf_putint(out,"n1",nx);
	sf_putfloat(out,"d1",dx);
	sf_putfloat(out,"o1",x0);
	sf_putstring(out,"label1","Offset");

	if (NULL != (unit=sf_histstring(in,"unit1"))) {
	    space=strchr(unit,'/');
	    if (*space == '\0') sf_putstring(out,"unit1",unit);
	}

	if (intp[0] == 's') {
		if (nw > 2) spl = sf_spline_init (nw, nv);
	}
	
	ntr = nv;
	ntm = nx;
    } else {
	if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input"); 
	if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input"); 

	if (!sf_getint("nv",&nv)) sf_error("Need nv=");
	/* number of velocities (if inv=n) */

	if (!sf_getfloat("vmin",&vmin)) sf_error("Need vmin=");
	/* minimum velocity (if inv=n) */

	if (!sf_getfloat("vmax",&vmax)) sf_error("Need vmax=");
	/* maximum velocity (if inv=n) */

	dv = (vmax - vmin)/(nv-1.);

	sf_putint(out,"n1",nv);
	sf_putfloat(out,"d1",dv);
	sf_putfloat(out,"o1",vmin);

	sf_putstring(out,"label1","Velocity");

	if (NULL != (time = sf_histstring(in,"unit2")) &&
	    NULL != (space = sf_histstring(in,"unit1"))) {
	    len = strlen(time)+strlen(space)+2;
	    unit = sf_charalloc(len);
	    snprintf(unit,len,"%s/%s",space,time);
	    sf_putstring(out,"unit1",unit);
	    free(time);
	    free(space);
	}

	sf_putint(out,"nx",nx);
	sf_putfloat(out,"x0",x0);
	sf_putfloat(out,"dx",dx);

	if (intp[0] == 's') {
		if (nw > 2) spl = sf_spline_init (nw, nx);
	}

	ntr = nx;
	ntm = nv;
    }
    n3 = sf_leftsize(in,2);

    trace = sf_floatalloc(ntr);
    modl =  sf_floatalloc(ntm);
    r = sf_floatalloc(ntm);

    for (i3=0; i3 < n3; i3++) {
	for (it=0; it < nt; it++) {
	    t = t0 + it*dt;

	    sf_floatread (trace,ntr,in);
	    if (intp[0] == 's') {
	    	if (nw > 2) sf_banded_solve (spl,trace);
	    }

	    if (t > tp) {
		if (inv) {
		    for (ix=0; ix < nx; ix++) {
			r[ix] = (x0-xp+ix*dx)/(t-tp);
		    }
		    if (intp[0] == 's') {
			sf_int1_init (r, vmin, dv, nv, sf_spline_int, nw, nx, 0.0);
		    } else {
			sf_int1_init (r, vmin, dv, nv, sf_lin_int, 2, nx, tau);
			shprefilter(ntr,trace); 
		    }
		} else {
		    for (iv=0; iv < nv; iv++) {
			r[iv] = xp+(vmin+iv*dv)*(t-tp);
		    }
			if (intp[0] == 's') {
			    sf_int1_init (r, x0,   dx, nx, sf_spline_int, nw, nv, 0.0);
			} else {
			    sf_int1_init (r, x0,   dx, nx, sf_lin_int, 2, nv, tau);
			    shprefilter(ntr,trace); 
			}
		}

		sf_int1_lop (false,false,ntr,ntm,trace,modl);
	    } else {
		for (im=0; im < ntm; im++) {
		    modl[im] = 0.;
		}
	    }

	    sf_floatwrite (modl,ntm,out);
	}
    }

    exit(0);
}

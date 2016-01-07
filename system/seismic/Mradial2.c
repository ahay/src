/* Another version of radial transform. */
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
#include "stretch4.h"

int main(int argc, char* argv[])
{
    map4 mo;
    bool inv;
    int nt, nx, nv, it, ix, iv, n3, i3;
    float *data, *modl, *r;
    float v0, vmax, dv, dx, x0, t0, t, dt, tp, eps;
    char *unit=NULL, *space=NULL, *time=NULL;
    size_t len;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    if (!sf_histint(in,"n2",&nt)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dt)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&t0)) sf_error("No o2= in input");
 
    if (!sf_getfloat("tp",&tp)) tp=t0;

    if (inv) {
	if (!sf_histint(in,"n1",&nv)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dv)) sf_error("No d1= in input"); 
	if (!sf_histfloat(in,"o1",&v0)) sf_error("No o1= in input"); 

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
    } else {
	if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input"); 
	if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input"); 

	if (!sf_getint("nv",&nv)) sf_error("Need nv=");
	/* number of velocities (if inv=n) */

	if (!sf_getfloat("vmin",&v0)) sf_error("Need vmin=");
	/* minimum velocity (if inv=n) */

	if (!sf_getfloat("vmax",&vmax)) sf_error("Need vmax=");
	/* maximum velocity (if inv=n) */

	dv = (vmax - v0)/(nv-1.);

	sf_putint(out,"n1",nv);
	sf_putfloat(out,"d1",dv);
	sf_putfloat(out,"o1",v0);

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
    }
    n3 = sf_leftsize(in,2);

    mo = stretch4_init (nx, x0, dx, nv, eps);

    data = sf_floatalloc(nx);
    modl = sf_floatalloc(nv);
    r = sf_floatalloc(nv);

    for (i3=0; i3 < n3; i3++) {
	for (it=0; it < nt; it++) {
	    t = t0 + it*dt;

	    if (inv) {
		sf_floatread(modl,nv,in);
	    } else {
		sf_floatread(data,nx,in);
	    }

	    if (t > tp) {
		for (iv=0; iv < nv; iv++) {
		    r[iv] = (v0+iv*dv)*(t-tp);
		}
		stretch4_define (mo,r);

		if (inv) {
		    stretch4_apply (false,mo,modl,data);
		} else {
		    stretch4_invert (false,mo,modl,data);
		}
	    } else {
		if (inv) {
		    for (ix=0; ix < nx; ix++) {
			data[ix] = 0.;
		    }
		} else {
		    for (iv=0; iv < nv; iv++) {
			modl[iv] = 0.;
		    }
		}
	    }

	    if (inv) {
		sf_floatwrite(data,nx,out);
	    } else {
		sf_floatwrite (modl,nv,out);
	    }
	}
    }

    exit(0);
}

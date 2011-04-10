/* Normal moveout in tau-p domain. */
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

static bool interval;
static float p, f, t0, dt, *vel, *velx;

static float nmo_map(float t, int it) {
    float ft, fi;

    ft = vel[it];
    ft = 1.-p*ft*ft;

    if (interval) {
	fi = f*dt;
	if (ft > 0.)
	    f += sqrtf(ft);
    } else {
	f = sqrtf(ft);
	fi = f*(t0+it*dt);
    }

    return fi;
}

static float anmo_map(float t, int it) {
    float ft,fi,gt;

    ft = vel[it];
    gt = velx[it];
    gt = 1.-p*gt*gt;
    gt = gt/(gt+p*ft*ft);
    fi = f*dt;
    if (gt > 0.)
	f += sqrtf(gt);
    return fi;
}

int main (int argc, char* argv[])
{
    fint1 nmo;
    int ip, ix, nt, np, nx, nw, mute, np2;
    float dp, p0, str;
    float *trace, *paxis;
    mapfunc map;
    sf_file taup, nmod, velocity, velocityx, slope;

    sf_init (argc,argv);
    taup = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(taup)) sf_error("Need float input");

    if (!sf_histint(taup,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(taup,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(taup,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(taup,"n2",&np)) sf_error("No n2= in input");


    nx = sf_leftsize(taup,2);

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */
    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch */

    if (!sf_getint("extend",&nw)) nw=4;
    /* interpolation accuracy */

    if (!sf_getbool("interval",&interval)) interval=true;
    /* use interval velocity */

    if (NULL != sf_getstring("slope")) {
	slope = sf_input("slope");
	np2 = sf_filesize(slope);
	if (np2 != np && np2 != np*nx) sf_error("Wrong dimensions in slope");

	paxis = sf_floatalloc(np2);
	sf_floatread (paxis,np2,slope);
	sf_fileclose(slope);
    } else {
        if (!sf_histfloat(taup,"d2",&dp)) sf_error("No d2= in input");
        if (!sf_histfloat(taup,"o2",&p0)) sf_error("No o2= in input");

        np2 = np;
        paxis = sf_floatalloc(np2);
        for (ip=0; ip < np; ip++)
        	paxis[ip] = p0 + ip*dp;
    }

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);

    if (NULL != sf_getstring("velx")) {
	velocityx = sf_input("velx");
	velx = sf_floatalloc(nt);
	map = anmo_map;
    } else {
	velocityx = NULL;
	velx = NULL;
	map = nmo_map;
    }

    nmo = fint1_init (nw, nt, mute);

    for (ix=0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);
	if (NULL != velx) sf_floatread (velx,nt,velocityx);

	for (ip=0; ip < np; ip++) {

		p = (np2 == np) ?  paxis[ip] : paxis[ix*np+ip];

	    p *= p;

	    sf_floatread (trace,nt,taup);
	    fint1_set(nmo,trace);

	    f = 0.;
	    stretch(nmo,map,nt,dt,t0,nt,dt,t0,trace,str);
	    sf_floatwrite (trace,nt,nmod);
	}
    }

    exit (0);
}

/* Kirchhoff 2.5-D modeling with analytical Green's functions. */
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
#include <float.h>
#include <math.h>

#include <rsf.h>

#include "kirmod.h"
#include "ricker.h"
#include "aastretch.h"

int main(int argc, char* argv[]) 
{
    int nx, nt, ns, nh, is, ih, ix, it;
    float *rfl, *crv, *shot, *trace, *ts, *tg, vel[5], slow;
    float dx, x0, dt, t0, ds, s0, dh, h0, r0, *time, *ampl, *delt, freq;
    maptype type = CONST;
    aamap map;
    sf_file refl, curv, modl, shots;

    sf_init(argc,argv);
    curv = sf_input("in");
    modl = sf_output("out");

    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint(curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");

    /*** Initialize trace ***/

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time sampling */
    if (!sf_getfloat("t0",&t0)) t0=0.;
    /* time origin */

    trace = sf_floatalloc(nt+1);

    sf_putint(modl,"n1",nt);
    sf_putfloat(modl,"d1",dt);
    sf_putfloat(modl,"o1",t0);

    /*** Initialize shots ***/

    if (NULL != sf_getstring("shots")) {
	shots = sf_input("shots");
	
	if (!sf_histint(shots,"n1",&ns)) sf_error("No n1= in shots");
    } else {
	shots = NULL;

	if (!sf_getint("ns",&ns)) ns=nx;
	/* number of shots */
	if (!sf_getfloat("s0",&s0)) s0=x0;
	/* first shot */
	if (!sf_getfloat("ds",&ds)) ds=dx;
	/* shot increment */

	sf_putfloat(modl,"o3",s0);
	sf_putfloat(modl,"d3",ds);
    }

    shot = sf_floatalloc(ns);

    if (NULL != shots) {
	sf_floatread(shot,ns,shots);
	sf_fileclose(shots);
    } else {
	for (is=0; is < ns; is++) {
	    shot[is] = s0+is*ds;
	}
    }

    sf_putint(modl,"n3",ns);

    /*** Initialize offsets ***/

    if (!sf_getint("nh",&nh)) nh=nx;
    /* number of offsets */
    if (!sf_getfloat("h0",&h0)) h0=0.;
    /* first offset */
    if (!sf_getfloat("dh",&dh)) dh=dx;
    /* offset increment */

    sf_putint(modl,"n2",nh);
    sf_putfloat(modl,"o2",h0);
    sf_putfloat(modl,"d2",dh);

    /*** Allocate space ***/
    
    kirmod_init(ns, s0, ds, nh, h0, dh);

    /*** Initialize reflector ***/

    crv = sf_floatalloc(nx);
    rfl = sf_floatalloc(nx);

    sf_floatread(crv,nx,curv);

    if (NULL != sf_getstring("refl")) {
	refl = sf_input("refl");
	sf_floatread(rfl,nx,refl);
	sf_fileclose(refl);
    } else {
	if (!sf_getfloat("r0",&r0)) r0=1.;
	/* reflectivity */
	for (ix=0; ix < nx; ix++) {
	    rfl[ix] = r0;
	}
    }

    /*** Initialize velocity ***/

    if (!sf_getfloat("vel",vel)) sf_error("Need vel=");
    /* velocity */
    
    if (!sf_getfloat("gradx",&vel[2])) vel[2]=0.;
    if (!sf_getfloat("gradz",&vel[1])) vel[1]=0.;
    /* velocity gradient */

    if (vel[2]==0. && vel[1]==0.) {
	type = CONST; 
        /* constant velocity */
    } else {
	type = S2LIN; 
        /* linear slowness squared */

	slow = 1./(vel[0]*vel[0]);
	/* slowness squared */
	vel[2] *= -slow/vel[0];
	vel[1] *= -slow/vel[0];
	vel[0] = slow;
    }

    if (!sf_getfloat("refx",&vel[4])) vel[4]=x0;
    if (!sf_getfloat("refz",&vel[3])) vel[3]=0.;
    /* reference coordinates for velocity */

    /*** Initialize stretch ***/
    map = aastretch_init (nt, t0, dt, nx);

    time = sf_floatalloc(nx);
    ampl = sf_floatalloc(nx);
    delt = sf_floatalloc(nx);

    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);

    /*** Compute traveltime table ***/

    kirmod_table (type, nx, x0, dx, crv, vel);

    /*** Main loop ***/
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < nt; it++) {
		trace[it] = 0.;
	    }

	    for (ix=0; ix < nx; ix++) {
		ts = kirmod_map(is,nh,ix);
		tg = kirmod_map(is,ih,ix);
		time[ix] = ts[0] + tg[0];
		ampl[ix] = 1./sqrt(ts[1]*tg[1]*(ts[1]+tg[1])+0.001*dt);
		delt[ix] = tg[2]*dh; 
		/* 2.5-D amplitude? */
	    }

	    aastretch_define (map,time,delt,ampl);
	    aastretch_apply (map,rfl,trace);

	    /* convolve with Ricker wavelet */
	    sf_freqfilt(nt,trace);

	    sf_floatwrite(trace,nt,modl);
	}
    }
  
    exit(0);
}
    

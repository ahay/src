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
#include "kirmod2.h"
#include "ricker.h"
#include "aastretch.h"

int main(int argc, char* argv[]) 
{
    int nx, nt, ns, nh, nc, nxc, is, ih, ix, ic;
    float **rfl, **rgd, **crv, **dip, *shot, *trace, vel[5], vel2[5];
    float slow, dx, x0, dt, t0, ds, s0, dh, h0, r0, **time, **ampl, **delt, freq;
    float theta, ava, amp;
    char *type, *type2;
    bool twod;
    surface inc, ref;
    ktable ts, tg;
    sf_file refl, curv, modl, shots;

    sf_init(argc,argv);
    curv = sf_input("in");
    modl = sf_output("out");

    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint  (curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histint  (curv,"n2",&nc)) nc=1; /* number of reflectors */
    nxc = nx*nc;

    /*** Initialize trace ***/

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time sampling */
    if (!sf_getfloat("t0",&t0)) t0=0.;
    /* time origin */

    trace = sf_floatalloc(nt+1);

    sf_putint  (modl,"n1",nt);
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

    if (!sf_getint  ("nh",&nh)) nh=nx;
    /* number of offsets */
    if (!sf_getfloat("h0",&h0)) h0=0.;
    /* first offset */
    if (!sf_getfloat("dh",&dh)) dh=dx;
    /* offset increment */

    sf_putint  (modl,"n2",nh);
    sf_putfloat(modl,"o2",h0);
    sf_putfloat(modl,"d2",dh);

    /*** Initialize reflector ***/

    crv = sf_floatalloc2(nx,nc);
    rfl = sf_floatalloc2(nx,nc);
    rgd = sf_floatalloc2(nx,nc);
    dip = sf_floatalloc2(nx,nc);

    sf_floatread(crv[0],nxc,curv);

    /* reflectivity (A) */
    if (NULL != sf_getstring("refl")) {
	refl = sf_input("refl");
	sf_floatread(rfl[0],nxc,refl);
	sf_fileclose(refl);
    } else {
	if (!sf_getfloat("r0",&r0)) r0=1.;
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		rfl[ic][ix] = r0;
	    }
	}
    }

    /* AVO gradient (B/A) */
    if (NULL != sf_getstring("rgrad")) {
	refl = sf_input("rgrad");
	sf_floatread(rgd[0],nxc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		rgd[ic][ix] = 0.;
	    }
	}
    }

    /* reflector dip */
    if (NULL != sf_getstring("dip")) {
	refl = sf_input("dip");
	sf_floatread(dip[0],nxc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		dip[ic][ix] = 0.;
	    }
	}
    }

    /*** Initialize velocity ***/

    if (!sf_getfloat("vel",&vel[0])) sf_error("Need vel=");
    /* velocity */
    
    if (!sf_getfloat("gradx",&vel[2])) vel[2]=0.;
    if (!sf_getfloat("gradz",&vel[1])) vel[1]=0.;
    /* velocity gradient */

    type = sf_getstring("type");
    /* type of velocity ('c': constant, 's': linear sloth, 'v': linear velocity) */
    if (NULL==type) {
	type= (vel[2]==0. && vel[1]==0.)?"const":"veloc";
    } else if (vel[2]==0. && vel[1]==0.) {
	free(type);
	type = "const"; 
    } else if ('s'==type[0]) {
	/* linear slowness squared */
	
	slow = 1./(vel[0]*vel[0]);
	/* slowness squared */
	vel[2] *= -2.*slow/vel[0];
	vel[1] *= -2.*slow/vel[0];
	vel[0] = slow;     
    } else if ('v' != type[0]) {
	sf_error("Unknown type=%s",type);
    }

    if (!sf_getbool("twod",&twod)) twod=false;
    /* 2-D or 2.5-D */
	
    if (!sf_getfloat("refx",&vel[4])) vel[4]=x0;
    if (!sf_getfloat("refz",&vel[3])) vel[3]=0.;
    /* reference coordinates for velocity */

    if (!sf_getfloat("vel2",&vel2[0])) vel2[0]=vel[0];
    /* converted velocity */
    
    if (!sf_getfloat("gradx2",&vel2[2])) vel2[2]=vel[2];
    if (!sf_getfloat("gradz2",&vel2[1])) vel2[1]=vel[1];
    /* converted velocity gradient */

    type2 = sf_getstring("type2");
    if (NULL==type2) {	
	type2=type;
    } else if (vel2[2]==0. && vel2[1]==0.) {
	free(type2);
	type2 = "const"; 
    } else if ('s'==type2[0]) {
	/* linear slowness squared */
	
	slow = 1./(vel2[0]*vel2[0]);
	/* slowness squared */
	vel2[2] *= -slow/vel2[0];
	vel2[1] *= -slow/vel2[0];
	vel2[0] = slow;     
    } else if ('v' != type2[0]) {
	sf_error("Unknown type=%s",type2);
    }

    if (!sf_getfloat("refx2",&vel2[4])) vel2[4]=vel[4];
    if (!sf_getfloat("refz2",&vel2[3])) vel2[3]=vel[3];
    /* reference coordinates for converted velocity */
    
    /*** Allocate space ***/
    
    inc = kirmod2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc);
    ref = (strcmp(type,type2) ||
	   vel2[0] != vel[0] || 
	   vel2[1] != vel[1] ||
	   vel2[2] != vel[2] ||
	   vel2[3] != vel[3] ||
	   vel2[4] != vel[4])?
	kirmod2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc):
	inc;
    
    /*** Initialize stretch ***/
    aastretch_init (nt, t0, dt, nxc);

    time = sf_floatalloc2(nx,nc);
    ampl = sf_floatalloc2(nx,nc);
    delt = sf_floatalloc2(nx,nc);

    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);

    /*** Compute traveltime table ***/

    kirmod2_table (inc, type[0], twod, crv, dip, vel);
    if (ref != inc) kirmod2_table (ref, type2[0], twod, crv, dip, vel2);

    /*** Main loop ***/
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (ix=0; ix < nx; ix++) {
		for (ic=0; ic < nc; ic++) {
		    ts = kirmod2_map(inc,is,nh,ix,ic);
		    tg = kirmod2_map(ref,is,ih,ix,ic);

		    time[ic][ix] = ts->t + tg->t;

		    theta = sinf(0.5*(tg->an - ts->an));
		    ava = 1.+rgd[ic][ix]*theta*theta;
		    if (ref != inc) ava *= theta;
		    
		    amp = 0.5*(ts->tn + tg->tn)/(ts->a * tg->a * sqrtf(ts->ar + tg->ar) + FLT_EPSILON);

		    ampl[ic][ix] = ava*amp*dx; 
		    delt[ic][ix] = fabsf(ts->tx+tg->tx)*dx; 
		}
	    }

	    aastretch_define (time[0],delt[0],ampl[0]);
	    aastretch_lop (false,false,nxc,nt,rfl[0],trace);

	    /* convolve with Ricker wavelet */
	    sf_freqfilt(nt,trace);
	    
	    sf_floatwrite(trace,nt,modl);
	}
    }
  
    exit(0);
}
    

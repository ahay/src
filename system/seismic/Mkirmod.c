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
    float **rfl, **rgd, **crv, **dip, *trace, *trace2;
    float **time, **ampl, **delt, freq, theta, ava, amp, obl;
    float slow, dx, x0, dt, t0, ds, s0, dh, h0, r0;
    char *type, *type2;
    bool twod, verb, adj, lin;
    surface inc, ref;
    velocity vel, vel2;
    ktable ts, tg, **tss, **tgs;
    sf_file data, refl, curv, modl;
    
    sf_init(argc,argv);

    if (!sf_getbool("lin",&lin)) lin=false;
    /* if linear operator */

    if (lin) {
	if (!sf_getbool("adj",&adj)) adj=false;
	/* adjoint flag */

	if (adj) {
	    modl = sf_input("in");
	    data = sf_output("out");
	} else {
	    data = sf_input("in");
	    modl = sf_output("out");
	}
	curv = sf_input("curv");
    } else {
	adj = false;
	curv = sf_input("in");
	modl = sf_output("out");
	data = NULL;
    }
    
    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint  (curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histint  (curv,"n2",&nc)) nc=1; /* number of reflectors */
    nxc = nx*nc;
    
    /*** Initialize trace ***/
    
    if (adj) {
	if (!sf_histint(modl,"n1",&nt)) sf_error("Need nt=");
	if (!sf_histfloat(modl,"d1",&dt)) dt=0.004;
	if (!sf_histfloat(modl,"o1",&t0)) t0=0.;


	/*** Initialize shots ***/
    
	if (!sf_histint(modl,"n3",&ns)) ns=nx;
	if (!sf_histfloat(modl,"o3",&s0)) s0=x0;
	if (!sf_histfloat(modl,"d3",&ds)) ds=dx;

	/*** Initialize offsets ***/

	if (!sf_histint  (modl,"n2",&nh)) nh=nx;
	if (!sf_histfloat(modl,"o2",&h0)) h0=0.;
	if (!sf_histfloat(modl,"d2",&dh)) dh=dx;


	sf_putint(data,"n1",nx);
	sf_putfloat(data,"d1",dx);
	sf_putfloat(data,"o1",x0);
	sf_putint(data,"n2",nc);
	sf_putint(data,"n3",1);
    } else {
	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples */
	if (!sf_getfloat("dt",&dt)) dt=0.004;
	/* time sampling */
	if (!sf_getfloat("t0",&t0)) t0=0.;
	/* time origin */
    
	sf_putint  (modl,"n1",nt);
	sf_putfloat(modl,"d1",dt);
	sf_putfloat(modl,"o1",t0);
	sf_putstring(modl,"label1","Time");
	sf_putstring(modl,"unit1","s");

	/*** Initialize shots ***/
	
	if (!sf_getint("ns",&ns)) ns=nx;
	/* number of shots */
	if (!sf_getfloat("s0",&s0)) s0=x0;
	/* first shot */
	if (!sf_getfloat("ds",&ds)) ds=dx;
	/* shot increment */
	
	sf_putfloat(modl,"o3",s0);
	sf_putfloat(modl,"d3",ds);
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
    }

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    trace = sf_floatalloc(nt);
    trace2 = sf_floatalloc(nt);

    /*** Initialize reflector ***/

    crv = sf_floatalloc2(nx,nc);
    rfl = sf_floatalloc2(nx,nc);
    rgd = sf_floatalloc2(nx,nc);
    dip = sf_floatalloc2(nx,nc);
    
    sf_floatread(crv[0],nxc,curv);
    
    if (!lin) {
	/* reflectivity (A) */
	if (NULL != sf_getstring("refl")) {
	    refl = sf_input("refl");
	    sf_floatread(rfl[0],nxc,refl);
	    sf_fileclose(refl);
	} else {
	    if (!sf_getfloat("r0",&r0)) r0=1.;
	    /* normal reflectivity (if constant) */
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][ix] = r0;
		}
	    }
	}
    }

    if (NULL != sf_getstring("rgrad")) {
	/* AVO gradient file (B/A) */
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

    if (NULL != sf_getstring("dip")) {
	/* reflector dip file */
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
    vel  = (velocity) sf_alloc(1,sizeof(*vel));
    vel2 = (velocity) sf_alloc(1,sizeof(*vel2));

    if (!sf_getfloat("vel",&(vel->v0))) sf_error("Need vel=");
    /* velocity */
    
    if (!sf_getfloat("gradx",&(vel->gx))) (vel->gx)=0.;
    if (!sf_getfloat("gradz",&(vel->gz))) (vel->gz)=0.;
    /* velocity gradient */

    type = sf_getstring("type");
    /* type of velocity ('c': constant, 's': linear sloth, 'v': linear velocity) */
    if (NULL==type) {
	type= ((vel->gx)==0. && (vel->gz)==0.)?"const":"veloc";
    } else if ((vel->gx)==0. && (vel->gz)==0.) {
	free(type);
	type = "const"; 
    } else if ('s'==type[0]) {
	/* linear slowness squared */
	
	slow = 1./((vel->v0)*(vel->v0));
	/* slowness squared */
	(vel->gx) *= -2.*slow/(vel->v0);
	(vel->gz) *= -2.*slow/(vel->v0);
	(vel->v0) = slow;     
    } else if ('v' != type[0]) {
	sf_error("Unknown type=%s",type);
    }

    if (!sf_getbool("twod",&twod)) twod=false;
    /* 2-D or 2.5-D */
	
    if (!sf_getfloat("refx",&(vel->x0))) (vel->x0)=x0;
    if (!sf_getfloat("refz",&(vel->z0))) (vel->z0)=0.;
    /* reference coordinates for velocity */

    if (!sf_getfloat("vel2",&(vel2->v0))) (vel2->v0)=(vel->v0);
    /* converted velocity */
    
    if (!sf_getfloat("gradx2",&(vel2->gx))) (vel2->gx)=(vel->gx);
    if (!sf_getfloat("gradz2",&(vel2->gz))) (vel2->gz)=(vel->gz);
    /* converted velocity gradient */

    type2 = sf_getstring("type2");
    /* type of the receiver side branch */
    if (NULL==type2) {	
	type2=type;
    } else if ((vel2->gx)==0. && (vel2->gz)==0.) {
	free(type2);
	type2 = "const"; 
    } else if ('s'==type2[0]) {
	/* linear slowness squared */
	
	slow = 1./((vel2->v0)*(vel2->v0));
	/* slowness squared */
	(vel2->gx) *= -slow/(vel2->v0);
	(vel2->gz) *= -slow/(vel2->v0);
	(vel2->v0) = slow;     
    } else if ('v' != type2[0]) {
	sf_error("Unknown type=%s",type2);
    }

    if (!sf_getfloat("refx2",&(vel2->x0))) (vel2->x0)=(vel->x0);
    if (!sf_getfloat("refz2",&(vel2->z0))) (vel2->z0)=(vel->z0);
    /* reference coordinates for converted velocity */
    
    /*** Allocate space ***/
    
    inc = kirmod2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc);
    if (strcmp(type,type2) ||
	(vel2->v0) != (vel->v0) || 
	(vel2->gz) != (vel->gz) ||
	(vel2->gx) != (vel->gx) ||
	(vel2->z0) != (vel->z0) ||
	(vel2->x0) != (vel->x0)) {
	ref = kirmod2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc);
    } else {
	ref = inc;
    }
    
    /*** Initialize stretch ***/
    aastretch_init (nt, t0, dt, nxc);

    time = sf_floatalloc2(nx,nc);
    ampl = sf_floatalloc2(nx,nc);
    delt = sf_floatalloc2(nx,nc);

    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);

    /*** Compute traveltime table ***/

    kirmod2_table (inc, vel, type[0], twod, crv, dip);
    if (ref != inc) kirmod2_table (ref, vel2, type2[0], twod, crv, dip);

    if (lin) {
	if (adj) {
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][ix] = 0.0;
		}
	    }
	} else {
	    sf_floatread(rfl[0],nxc,data);
	}
    }

    tss = (ktable**) sf_alloc(nc,sizeof(ktable*));
    tgs = (ktable**) sf_alloc(nc,sizeof(ktable*));
    for (ic=0; ic < nc; ic++) {
	tss[ic] = (ktable*) sf_alloc(nx,sizeof(ktable));
	tgs[ic] = (ktable*) sf_alloc(nx,sizeof(ktable));
    }

    /*** Main loop ***/
    for (is=0; is < ns; is++) {
	if (verb) sf_warning("shot %d of %d",is+1,ns);

	for (ih=0; ih < nh; ih++) {

	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    ts = kirmod2_map(inc,is,nh,ix,ic);
		    tg = kirmod2_map(ref,is,ih,ix,ic);

		    time[ic][ix] = ts->t + tg->t;
		    delt[ic][ix] = fabsf(ts->tx+tg->tx)*dx;

		    tss[ic][ix] = ts;
		    tgs[ic][ix] = tg;
		}
	    }

	    aastretch_define (time[0],delt[0],NULL);

	    if (adj) {
		sf_floatread(trace2,nt,modl);

		/* convolve with Ricker wavelet */
		sf_freqfilt_lop(true,false,nt,nt,trace,trace2);
	    
		aastretch_lop (true,false,nxc,nt,ampl[0],trace);
	    }

	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    ts = tss[ic][ix];
		    tg = tgs[ic][ix];

		    obl = 0.5*(ts->tn + tg->tn);
		    amp = ts->a * tg->a * sqrtf(ts->ar + tg->ar) + FLT_EPSILON;

		    if (lin) {
			if (adj) {
			    rfl[ic][ix] += ampl[ic][ix]*obl*dx/amp;
			} else {
			    ampl[ic][ix] = rfl[ic][ix]*obl*dx/amp;
			}
		    } else {
			theta = 0.5*(SF_SIG(tg->tx)*tg->an - 
				     SF_SIG(ts->tx)*ts->an);
			theta = sinf(theta);
			
			ava = rfl[ic][ix]+rgd[ic][ix]*theta*theta;
			if (ref != inc) ava *= theta;
			
			obl = 0.5*(ts->tn + tg->tn);
			amp = ts->a * tg->a * sqrtf(ts->ar + tg->ar) + FLT_EPSILON;
			
			ampl[ic][ix] = ava*obl*dx/amp;
		    }
		}
	    }

	    if (!adj) {
		aastretch_lop (false,false,nxc,nt,ampl[0],trace);

		/* convolve with Ricker wavelet */
		sf_freqfilt_lop(false,false,nt,nt,trace,trace2);
	    
		sf_floatwrite(trace2,nt,modl);
	    }
	}
    }

    if (lin && adj) sf_floatwrite(rfl[0],nxc,data);
  
    exit(0);
}
    

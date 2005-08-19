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
#include "kirmod3.h"
#include "ricker.h"
#include "aastretch.h"

int main(int argc, char* argv[]) 
{
    int nx, ny, nt, ns, nsx, nsy, nh, nc, nxyc, is, isx, isy, ih, ix, iy, ic;
    float ***rfl, ***rgd, ***crv, ***dipx, ***dipy, **shot, *trace;
    float slow, dx, x0, dt, t0, dsx, dsy, sx0, sy0, dh, h0, r0, **time, **ampl, **delt, freq;
    float theta, ava, amp, obl;
    char *type, *type2;
    bool twod;
    surface inc, ref;
    velocity vel, vel2;
    ktable ts, tg;
    sf_file refl, curv, modl, shots;
    
    sf_init(argc,argv);
    curv = sf_input("in");
    modl = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint  (curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histint  (curv,"n2",&ny)) sf_error("No n2= in input");
    if (!sf_histfloat(curv,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(curv,"o2",&y0)) sf_error("No o2= in input");
    if (!sf_histint  (curv,"n3",&nc)) nc=1; /* number of reflectors */
    nxyc = nx*ny*nc;
    
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
	
	if (!sf_histint(shots,"n1",&two) || 2 != two) sf_error("Need n1=2 in shots");
	if (!sf_histint(shots,"n2",&ns)) sf_error("No n2= in shots");

	sf_putint(modl,"n4",ns);
    } else {
	shots = NULL;
	
	if (!sf_getint("nsx",&nsx)) nsx=nx;
	/* number of inline shots */
	if (!sf_getfloat("sx0",&sx0)) sx0=x0;
	/* first inline shot */
	if (!sf_getfloat("dsx",&dsx)) dsx=dx;
	/* inline shot increment */

	sf_putfloat(modl,"o4",sx0);
	sf_putfloat(modl,"d4",dsx);
	sf_putint(modl,"n4",nsx);

	if (!sf_getint("nsy",&nsy)) nsy=ny;
	/* number of crossline shots */
	if (!sf_getfloat("sy0",&sy0)) sy0=y0;
	/* first crossline shot */
	if (!sf_getfloat("dsy",&dsy)) dsy=dy;
	/* crossline shot increment */

	sf_putfloat(modl,"o5",sy0);
	sf_putfloat(modl,"d5",dsy);
	sf_putint(modl,"n5",nsy);

	ns = nsx*nsy;
    }
    
    shot = sf_floatalloc2(2,ns);

    if (NULL != shots) {
	sf_floatread(shot[0],2*ns,shots);
	sf_fileclose(shots);
    } else {	
	for (isy=0; isy < nsy; isy++) {
	    for (isx=0; isx < nsx; isx++) {
		is = isx + isy*nsx;
		shot[is][0] = sx0+isx*dsx;
		shot[is][1] = sy0+isy*dsy;
	    }
	}
    }
    
    /*** Initialize offsets ***/

    if (!sf_getint  ("nhx",&nhx)) nhx=nx;
    /* number of inline offsets */
    if (!sf_getfloat("hx0",&hx0)) hx0=0.;
    /* first inline offset */
    if (!sf_getfloat("dhx",&dhx)) dhx=dx;
    /* inline offset increment */

    sf_putint  (modl,"n2",nhx);
    sf_putfloat(modl,"o2",hx0);
    sf_putfloat(modl,"d2",dhx);
    
    if (!sf_getint  ("nhy",&nhy)) nh=ny;
    /* number of crossline offsets */
    if (!sf_getfloat("hy0",&hy0)) hy0=0.;
    /* first crossline offset */
    if (!sf_getfloat("dhy",&dhy)) dhy=dy;
    /* crossline offset increment */

    sf_putint  (modl,"n3",nhy);
    sf_putfloat(modl,"o3",hy0);
    sf_putfloat(modl,"d3",dhy);

    /*** Initialize reflector ***/

    crv = sf_floatalloc3(nx,ny,nc);
    rfl = sf_floatalloc3(nx,ny,nc);
    rgd = sf_floatalloc3(nx,ny,nc);
    dip = sf_floatalloc3(nx,ny,nc);
    
    sf_floatread(crv[0][0],nxyc,curv);
    
    /* reflectivity (A) */
    if (NULL != sf_getstring("refl")) {
	refl = sf_input("refl");
	sf_floatread(rfl[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	if (!sf_getfloat("r0",&r0)) r0=1.;
	/* constant reflectivity */
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][iy][ix] = r0;
		}
	    }
	}
    }

    /* AVO gradient (B/A) */
    if (NULL != sf_getstring("rgrad")) {
	refl = sf_input("rgrad");
	sf_floatread(rgd[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    rgd[ic][iy][ix] = 0.;
		}
	    }
	}
    }

    /* reflector inline dip */
    if (NULL != sf_getstring("dipx")) {
	refl = sf_input("dipx");
	sf_floatread(dipx[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    dipx[ic][iy][ix] = 0.;
		}
	    }
	}
    }

    /* reflector crossline dip */
    if (NULL != sf_getstring("dipy")) {
	refl = sf_input("dipy");
	sf_floatread(dipy[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    dipy[ic][iy][ix] = 0.;
		}
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

    /*** Main loop ***/
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (ix=0; ix < nx; ix++) {
		for (ic=0; ic < nc; ic++) {
		    ts = kirmod2_map(inc,is,nh,ix,ic);
		    tg = kirmod2_map(ref,is,ih,ix,ic);

		    time[ic][ix] = ts->t + tg->t;

		    theta = sinf(0.5*(SF_SIG(tg->tx)*tg->an - SF_SIG(ts->tx)*ts->an));
		    ava = 1.+rgd[ic][ix]*theta*theta;
		    if (ref != inc) ava *= theta;

		    obl = 0.5*(ts->tn + tg->tn);
		    
		    amp = ts->a * tg->a * sqrtf(ts->ar + tg->ar) + FLT_EPSILON;

		    ampl[ic][ix] = ava*obl*dx/amp; 
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
    

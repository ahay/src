/* Kirchhoff 2-D modeling for linear slowness squared. */
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

#include "kirmod.h"

int main(int argc, char* argv[]) 
{
    int nx, nt, ns, nh, is, ih, ix, it;
    float *rfl, *crv, *recv, *shot, *trace, *ts, *tg, slow2[5];
    float dx, x0, dt, t0, ds, s0, dh, h0, r0;
    maptype type = CONST;
    sf_file refl, curv, modl, shots, recvs;

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

    trace = sf_floatalloc(nt);

    sf_putint(modl,"n1",nt);
    
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

    if (NULL != sf_getstring("recvs")) {
	recvs = sf_input("recvs");
	
	if (!sf_histint(recvs,"n1",&nh)) sf_error("No n1= in recvs");
    } else {
	recvs = NULL;

	if (!sf_getint("nh",&nh)) nh=nx;
	/* number of recvs */
	if (!sf_getfloat("h0",&h0)) h0=0.;
	/* first recv */
	if (!sf_getfloat("dh",&dh)) dh=dx;
	/* recv increment */
    }

    recv = sf_floatalloc(nh);

    if (NULL != recvs) {
	sf_floatread(recv,nh,recvs);
	sf_fileclose(recvs);
    } else {
	for (ih=0; ih < nh; ih++) {
	    recv[ih] = h0+ih*dh;
	}
    }
    
    sf_putint(modl,"n2",nh);

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

    /*** Initialize slowness ***/

    if (!sf_getfloat("slow2",slow2)) {
	/* slowness squared */
	if (!sf_getfloat("vel",slow2)) sf_error("Need slow2= or vel=");
	/* velocity (if no slow2=) */
	slow2[0] = 1./(slow2[0]*slow2[0]);
    } 

    if (!sf_getfloat("refx",&slow2[4])) slow2[4]=x0;
    if (!sf_getfloat("refz",&slow2[3])) slow2[3]=0.;
    /* reference coordinates for slowness */
    
    if (!sf_getfloat("gradx",&slow2[2])) slow2[2]=0.;
    if (!sf_getfloat("gradz",&slow2[1])) slow2[1]=0.;
    /* gradient of slowness squared */

    slow2[2] *= 0.5;
    slow2[1] *= 0.5;
    
    /*** Compute traveltime table ***/

    kirmod_table (type, nx, x0, dx, crv, slow2);

    /*** Main loop ***/
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < nt; it++) {
		trace[it] = 0.;
	    }

	    for (ix=0; ix < nx; ix++) {
		ts = kirmod_map(is,nh,ix);
		tg = kirmod_map(is,ih,ix);
	    }

	    sf_floatwrite(trace,nt,modl);
	}
    }
  
    exit(0);
}
    

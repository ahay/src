/* Post-stack 2-D/3-D v(z) time modeling/migration with Gazdag phase-shift. */
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

#include "gazdag.h"

int main (int argc, char *argv[])
{
    bool verb;
    int nt, nt2;	/* number of time samples */
    int nz;		/* number of migrated time samples */
    int nx, ny;		/* number of wavenumbers */
    int it,ix,iy,iz;          /* loop counters 	*/
    
    float dt;		/* time sampling interval 	*/
    float dz;		/* migrated time sampling interval */
    float dx,dy;	/* wave number sampling interval */
    float x,y;          /* wave number */
    float x0,y0;        /* wave number origin */
    float *vt, v0;	/* velocity v(t)		*/
    float *vz, vz0;	/* vertical velocity v(t)		*/
    float *n, n0;	/* eta		*/

    float *p,*q;	/* input, output data		*/

    bool inv;           /* modeling or migration        */
    bool depth;         /* time or depth migration      */
    float eps;          /* dip filter constant          */   

    char *rule;         /* phase-shuft interpolation rule */

    sf_file in, out, vel, velz, eta;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* if y, modeling; if n, migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* stabilization parameter */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool("depth",&depth)) depth = false;
    /* if true, depth migration */

    if (!sf_histint(in,"n2",&nx)) nx = 1;
    if (!sf_histfloat(in,"d2",&dx)) 
	sf_error ("No d2= in input");
    if (!sf_histfloat(in,"o2",&x0)) x0=0.;

    if (!sf_histint(in,"n3",&ny)) ny = 1;
    if (!sf_histfloat(in,"d3",&dy)) dy=dx;
    if (!sf_histfloat(in,"o3",&y0)) y0=0.;

    dx *= 2.*SF_PI; x0 *= 2.*SF_PI;
    dy *= 2.*SF_PI; y0 *= 2.*SF_PI;

    if (NULL == sf_getstring("velocity")) {
	vel = NULL;
	velz = NULL;
	eta = NULL;
    } else {
	vel = sf_input("velocity");
	if (!sf_histint(vel,"n1",&nz)) 
	    sf_error ("No n1= in velocity");
	if (!sf_histfloat(vel,"d1",&dz)) 
	    sf_error ("No d1= in velocity");
	
	if (NULL == sf_getstring("velz")) {
	    velz = NULL;
	    eta = NULL;
	} else {
	    velz = sf_input("velz");
	    if (NULL == sf_getstring("eta")) sf_error("Need eta=");
	    eta = sf_input("eta");
	}
    }

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");

	if (!sf_getint("nt",&nt)) {
	    /* Length of time axis (for modeling) */
	    if (depth) {
		sf_error ("nt= must be supplied");
	    } else {
		nt=nz;
	    }
	} else {
	    sf_putint(out,"n1",nt);
	}

	if (!sf_getfloat("dt",&dt)) {
	    /* Sampling of time axis (for modeling) */
	    if (depth) {
		sf_error ("dt= must be supplied");
	    } else {
		dt=dz;
	    }
	} else {
	    sf_putfloat(out,"d1",dt);
	}

	sf_putfloat(out,"o1",0.);
    } else { /* migration */
	if (!sf_histint(in,"n1",&nt)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error ("No d1= in input");

	if (NULL == vel) {
	    if (!sf_getint("nz",&nz)) {
		/* Length of depth axis (for migration, if no velocity file) */
		if (depth) {
		    sf_error("Need nz=");
		} else {
		    nz = nt;
		}
	    }
	    if (!sf_getfloat("dz",&dz)) {
		/* Sampling of depth axis (for migration, if no velocity file) */
		if (depth) {
		    sf_error("Need dz=");
		} else {
		    dz = dt;
		}
	    }
	}

	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"o1",0.);
    }
    
    vt = sf_floatalloc(nz);
    vz = sf_floatalloc(nz);
    n = sf_floatalloc(nz);

    if (NULL == (rule = sf_getstring("rule"))) rule="simple";
    /* phase-shift interpolation rule (simple, midpoint, linear) */

    if (NULL == vel) {
	/* file with velocity */
	if (!sf_getfloat("vel",&v0)) sf_error ("vel= must be supplied");
	/* Constant velocity (if no velocity file) */

	if (!sf_getfloat("vz",&vz0)) vz0=v0;
	/* Constant vertical velocity (if no velocity file) */

	if (!sf_getfloat("n",&n0)) n0=0.0;
	/* Constant eta (if no velocity file) */	

	for (iz=0; iz < nz; iz++) {
	    vt[iz] = v0;
	    vz[iz] = vz0;
	    n[iz] = n0;
	}
    } else {
	sf_floatread(vt,nz,vel);
	sf_fileclose(vel);

	if ('a' == rule[0]) {
	    if (NULL == velz || NULL == eta) sf_error("Need velz= and eta=");
	    sf_floatread(vz,nz,velz);
	    sf_floatread(n,nz,eta);

	    sf_fileclose(velz);
	    sf_fileclose(eta);
	}  else {
	    for (iz=0; iz < nz; iz++) {
		vz[iz] = vt[iz];
		n[iz] = 0.0;
	    }
	}	
    }

    /* vt -> 1/4 vt^2 */ 
    for (iz=0; iz < nz; iz++) {
	vt[iz] *= 0.25*vt[iz];
	vz[iz] *= 0.25*vz[iz];

	if (depth) {
	    vt[iz] = 1./vt[iz];
	    vz[iz] = 1./vz[iz];
	}
    }

    /* determine frequency sampling */    
    if (!sf_getint("pad",&nt2)) nt2 = 2*kiss_fft_next_fast_size((nt+1)/2);

    p = sf_floatalloc(nt2);
    q = sf_floatalloc(nz);

    gazdag_init (eps, nt2, dt, nz, dz, vt, vz, n, depth, rule[0]);

    for (iy=0; iy < ny; iy++) {
	y = y0+iy*dy;
	y *= y;

	for (ix=0; ix < nx; ix++) {
	    x = x0+ix*dx;
	    x = x*x+y;
	
	    if (verb) sf_warning("wavenumber %d of %d;",iy*nx+ix+1,nx*ny);
	
	    if (inv) {
		sf_floatread(q,nz,in);
	    } else {
		sf_floatread(p,nt,in);
		if (nt != nt2) {
		    for (it=nt; it < nt2; it++) {
			p[it]=0.;
		    }
		}
	    }

	    gazdag(inv,x,p,q);

	    if (inv) {
		sf_floatwrite(p,nt,out);
	    } else {
		sf_floatwrite(q,nz,out);
	    }
	}
    } 
    if (verb) sf_warning(".");

    exit (0);
}

/* 	$Id$	 */

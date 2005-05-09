/* Post-stack 2-D v(z) time modeling/migration with Gazdag phase-shift. */
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
    int nt, nt2;	/* number of time samples */
    int nz;		/* number of migrated time samples */
    int nk;		/* number of wavenumbers */
    int ik,iz;          /* loop counters 	*/
    
    float dt;		/* time sampling interval 	*/
    float dz;		/* migrated time sampling interval */
    float dk;	        /* wave number sampling interval */
    float k;            /* wave number */
    float k0;           /* wave number origin */
    float *vt, v0;	/* velocity v(t)		*/
    float *p,*q;	/* input, output data		*/

    bool inv;           /* modeling or migration        */
    bool depth;         /* time or depth migration      */
    float eps;          /* dip filter constant          */   

    char *rule;         /* phase-shuft interpolation rule */

    sf_file in, out, vel;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; if n, migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* Stabilization parameter */

    if (!sf_getbool("depth",&depth)) depth = false;
    /* if true, depth migration */

    if (!sf_histint(in,"n2",&nk)) nk = 1;
    if (!sf_histfloat(in,"d2",&dk)) 
	sf_error ("No d2= in input");
    if (!sf_histfloat(in,"o2",&k0)) k0=0.;
    dk *= 2.*SF_PI;
    k0 *= 2.*SF_PI;

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
	if (NULL == sf_getstring("velocity")) {
	    if (!sf_getint("nz",&nz)) nz = nt;
	    /* Length of depth axis (for migration, if no velocity file) */
	    if (!sf_getfloat("dz",&dz)) dz = dt;
	    /* Sampling of depth axis (for migration, if no velocity file) */
	} else {
	    vel = sf_input("velocity");
	    if (!sf_histint(vel,"n1",&nz)) 
		sf_error ("No n1= in velocity");
	    if (!sf_histfloat(vel,"d1",&dz)) 
		sf_error ("No d1= in velocity");
	}
	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"o1",0.);
    }
    
    vt = sf_floatalloc(nz);
    if (NULL == sf_getstring("velocity")) {
	/* file with velocity */
	if (!sf_getfloat("vel",&v0)) sf_error ("vel= must be supplied");
	/* Constant velocity (if no velocity file) */
	for (iz=0; iz < nz; iz++) {
	    vt[iz] = v0;
	}
    } else {
	vel = sf_input("velocity");
	sf_floatread(vt,nz,vel);
	sf_fileclose(vel);
    }

    /* vt -> 1/4 vt^2 */ 
    for (iz=0; iz < nz; iz++) {
	vt[iz] *= 0.25*vt[iz];
	if (depth) vt[iz] = 1./vt[iz];
    }

    /* determine frequency sampling */    
    nt2 = nt;
    if (nt%2) nt2++;

    p = sf_floatalloc(nt2);
    q = sf_floatalloc(nz);

    if (NULL == (rule = sf_getstring("rule"))) rule="simple";
    /* phase-shift interpolation rule (simple, midpoint, linear) */

    gazdag_init (eps, nt2, dt, nz, dz, vt, depth, rule[0]);

    for (ik=0; ik < nk; ik++) {
	k = k0+ik*dk;
	k *= k;
	
	if (inv) {
	    sf_floatread(q,nz,in);
	} else {
	    sf_floatread(p,nt,in);
	    if (nt != nt2) p[nt]=0.;
	}

	gazdag(inv,k,p,q);

	if (inv) {
	    sf_floatwrite(p,nt,out);
	} else {
	    sf_floatwrite(q,nz,out);
	}
    } 
    
    exit (0);
}

/* 	$Id$	 */

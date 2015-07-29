/* Time-to-depth conversion in V(z). 

July 2013 program of the month:
http://www.ahay.org/rsflog/index.php?/archives/345-Program-of-the-month-sftime2depth.html
*/
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

#include "fint1.h"

int main (int argc, char *argv[])
{
    int nt, nz, nx, iz, it, ix, nw, n2;
    bool slow, twoway, intime;
    float t0, dt, z0, dz, t=0., z=0., eps;
    float *time, *depth, *vel, *str;
    fint1 fnt;
    sf_map mp;
    sf_file in, out, velocity;

    sf_init(argc, argv);
    in       = sf_input("in");
    velocity = sf_input("velocity");
    out      = sf_output("out");

    if (!sf_histint (in,"n1",&nt)) sf_error ("No n1= in input");
    if (!sf_histfloat (in,"d1",&dt)) sf_error ("No d1= in input");
    if (!sf_histfloat (in,"o1",&t0)) t0 = 0.;

    nx = sf_leftsize(in,1);

    if (!sf_getbool("intime",&intime)) intime=false;
    /* y if velocity is in time rather than depth */

    if ((intime || !sf_histint(velocity,"n1",&nz)) && !sf_getint ("nz",&nz)) {
	/* Number of depth samples (default: n1) */
	nz = nt; 
    } else {
	sf_putfloat(out,"n1",nz);
    }
    if ((intime || !sf_histfloat(velocity,"d1",&dz)) && !sf_getfloat ("dz",&dz)) {	
	/* Depth sampling (default: d1) */
	dz = dt; 
    } else {
	sf_putfloat(out,"d1",dz);
    }
    if ((intime || !sf_histfloat(velocity,"o1",&z0)) && !sf_getfloat ("z0",&z0)) z0 = 0.; 
    /* Depth origin */
    if (!sf_histint(velocity,"n2",&n2)) n2=1;

    sf_putfloat(out,"o1",z0);

    if (!sf_getint ("extend",&nw)) nw = 4;
    /* Interpolation accuracy */
    if (!sf_getbool ("slow",&slow)) slow = false;
    /* If y, input slowness; if n, velocity */

    time = sf_floatalloc (nt);
    depth = sf_floatalloc (nz);
    vel = sf_floatalloc (intime? nt: nz);

    if (!sf_getbool ("twoway",&twoway)) twoway = true;
    /* if y, two-way traveltime */
    if (twoway) {
	dz *= 2.;
	z0 *= 2.;
    }

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    if (intime) {
	fnt = NULL;
	mp = sf_stretch_init (nz, z0, dz, nt, eps, false);
	str = sf_floatalloc(nt);
    } else {
	fnt = fint1_init (nw, nt, 0);
	mp = NULL;
	str = sf_floatalloc(nz);
    }

    for (ix = 0; ix < nx; ix++) {
	sf_floatread (time,nt,in);
	if (0 == ix || n2 > 1) {
	    sf_floatread (vel,intime? nt: nz,velocity);
	    if (intime) {
		for (it=0; it < nt; it++) {
		    if (it == 0) {
			t =  slow? t0/(vel[0]*dt): t0*vel[0]/dt;
		    } else {
			t += slow? 1./vel[it-1]: vel[it-1];
		    }		 
		    str[it] = t*dt;		 
		}
		sf_stretch_define (mp,str);
	    } else {
		for (iz = 0; iz < nz; iz++) {
		    if (iz == 0) {
			z =  slow? z0*vel[0]/dz: z0/(dz*vel[0]);
		    } else {
			z += slow? vel[iz-1]: 1./vel[iz-1];
		    }
		
		    str[iz] = (z*dz-t0)/dt;
		}
	    }
	}

	if (intime) {
	    sf_stretch_apply (mp,time,depth);
	} else {
	    fint1_set (fnt, time);
	    
	    for (iz = 0; iz < nz; iz++) {
		t = str[iz];
		it = t;
		t = t - it;
		
		if (it < -nw/2) {
		    depth[iz] = 0.;
		} else if (it > nt + nw/2) {
		    depth[iz] = time[nt-1];
		} else {
		    depth[iz] = fint1_apply (fnt, it, t, false);
		}		
	    }
	}

	sf_floatwrite (depth,nz,out);
    }


    exit (0);
}

/* 	$Id: Mtime2depth.c 10817 2013-08-03 16:36:21Z sfomel $	 */

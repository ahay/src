/* Second-order cell ray tracing with locally parabolic rays in 3-D.

Takes: > rays.rsf

Rays and wavefronts can be displayed with sfplotrays.
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

#include <math.h>

#include <rsf.h>

#include "celltrace3.h"

int main(int argc, char* argv[])
{
    bool velocity;
    int is, nz, ny, nx, im, nm, order, nshot, ndim, nsr, two;
    int nt, nr, ir, it, na, nb, ia, ib;
    float da=0., a0, amax, db=0., b0, bmax, *t;
    float x[3], p[3], dz, dy, dx, z0, y0, x0, **trj, *slow, **s, *a, *b;
    celltrace3 ct;
    char *trajname;
    FILE *traj;
    sf_file shots, vel, angles, time;

    sf_init (argc,argv);
    vel = sf_input("in");
    time = sf_output("out");

    trajname = sf_getstring("traj");

    if (NULL != trajname) {
	traj = fopen(trajname,"wb");
	if (NULL==traj) sf_error("Cannot open file \"%s\" for writing",trajname);
    } else {
	traj = NULL;
    }

    /* get 2-D grid parameters */
    if (!sf_histint(vel,"n1",&nz))   sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&ny))   sf_error("No n2= in input");
    if (!sf_histint(vel,"n3",&nx))   sf_error("No n3= in input");

    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"d3",&dx)) sf_error("No d3= in input");

    if (!sf_histfloat(vel,"o1",&z0)) z0=0.;
    if (!sf_histfloat(vel,"o2",&y0)) y0=0.;
    if (!sf_histfloat(vel,"o3",&x0)) x0=0.;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    /* If y, the input is velocity; if n, slowness */
    if(!sf_getint("order",&order)) order=4;
    /* Interpolation accuracy */

    if (!sf_getint("nt",&nt)) nt=nx*nz;
    /* number of time steps */

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
	/* file with shot locations */
	shots = sf_input("shotfile");
	if (!sf_histint(shots,"n1",&ndim) || 3 != ndim) 
	    sf_error("Must have n1=3 in shotfile");
	if (!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
  
	s = sf_floatalloc2 (ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose (shots);
    } else {
	nshot = 1;
	ndim = 3;

	s = sf_floatalloc2 (ndim,nshot);

	if (!sf_getfloat("zshot",s[0]))   s[0][0]=0.;
	/* shot location in depth (if shotfile is not specified) */
	if (!sf_getfloat("yshot",s[0]+1)) s[0][1]=y0 + 0.5*(ny-1)*dy;
	/* shot location inline (if shotfile is not specified) */
	if (!sf_getfloat("xshot",s[0]+2)) s[0][2]=x0 + 0.5*(nx-1)*dx;
	/* shot location crossline (if shotfile is not specified) */
	
	sf_warning("Shooting from z=%g, y=%g, x=%g",s[0][0],s[0][1],s[0][2]);
    }

    if (NULL != sf_getstring("anglefile")) {
	/* file with initial angles */
	angles = sf_input("anglefile");

	if (!sf_histint(angles,"n1",&nr)) sf_error("No n1= in anglefile");
	if (!sf_histint(angles,"n2",&two) || 2 != two) sf_error("Need n2=2 in anglefile");

	a = sf_floatalloc(nr);
	b = sf_floatalloc(nr);
    } else {
	angles = NULL;

	if (!sf_getint("na",&na)) sf_error("Need na=");
	/* Number of azimuths (if anglefile is not specified) */
	if (!sf_getint("nb",&nb)) sf_error("Need nb=");
	/* Number of inclinations (if anglefile is not specified) */

	if (!sf_getfloat("a0",&a0)) a0 = 0.; 
	/* First azimuth angle in degrees (if anglefile is not specified) */
	if (!sf_getfloat("amax",&amax)) amax=360.;
	/* Maximum azimuth angle in degrees (if anglefile is not specified) */

	if (!sf_getfloat("b0",&b0)) b0 = 0.; 
	/* First inclination angle in degrees (if anglefile is not specified) */
	if (!sf_getfloat("bmax",&bmax)) bmax=180.;
	/* Maximum inclination angle in degrees (if anglefile is not specified) */
	
	/* convert degrees to radians */
	a0 *= SF_PI/180.;
	amax *= SF_PI/180.;
	b0 *= SF_PI/180.;
	bmax *= SF_PI/180.;

	/* figure out angle spacing */
	da = (na > 1)? (amax - a0)/(na-1) : 0.;
	db = (nb > 1)? (bmax - b0)/(nb-1) : 0.;

	nr = na*nb;

	a = sf_floatalloc(nr);
	b = sf_floatalloc(nr);

	for (ir=ib=0; ib < nb; ib++) {
	    for (ia=0; ia < na; ia++, ir++) {
		b[ir] = b0 + ib*db;
		a[ir] = a0 + ia*da;
	    }
	}
    }
 
    /* specify output dimensions */
    nsr = nr*nshot;
    if (NULL != traj) fwrite(&nsr,sizeof(int),1,traj);
	    
    sf_putint(time,"n1",nr);
    sf_putint(time,"n2",nshot);

    /* get slowness */
    nm = nz*nx;
    slow = sf_floatalloc(nm);

    sf_floatread(slow,nm,vel);

    if (vel) { /* convert to slowness */
	for(im = 0; im < nm; im++){
	    slow[im] = 1./slow[im];
	}
    }

    /* initialize ray tracing object */
    ct = celltrace3_init (order, nt, nz, ny, nx, dz, dy, dx, z0, y0, x0, slow);
    
    free (slow);

    trj = sf_floatalloc2 (ndim,nt);
    t = sf_floatalloc(nr);

    for( is = 0; is < nshot; is++) { /* loop over shots */
	/* initialize angles */
	if (NULL != angles) { /* read angles from file */
	    sf_floatread(a,nr,angles);
	    sf_floatread(b,nr,angles);
	}

	for (ir = 0; ir < nr; ir++) { /* loop over rays */
	    /* initialize position */
	    x[0] = s[is][0]; 
	    x[1] = s[is][1];
	    x[2] = s[is][2];
	    
	    /* initialize direction */
	    
	    p[2] = +sinf(a[ir])*sinf(b[ir]);
	    p[1] = -cosf(a[ir])*sinf(b[ir]);
	    p[0] = cosf(b[ir]);
	    
	    t[ir] = cell_trace3 (ct, x, p, &it, trj);
	    
	    if (NULL != traj) {
		if (it < 0) it = -it; /* keep side-exiting rays */
		fwrite(&it,sizeof(int),1,traj);
		fwrite(trj[0],sizeof(float),(it+1)*ndim,traj);
	    }
	}
	
	sf_floatwrite(t,nr,time);
    }


    exit (0);
}

/* 	$Id$	 */

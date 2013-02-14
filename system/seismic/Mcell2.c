/* Second-order cell ray tracing with locally parabolic rays.

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

int main(int argc, char* argv[])
{
    bool velocity, lsint;
    int is, nz, nx, im, nm, order, nshot, ndim, nsr;
    int nt, nr, ir, it;
    float da=0., a0, amax, *t;
    float x[2], p[2], dz, dx, z0, x0, **trj, *slow, **s, *a;
    sf_celltrace ct;
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
    if (!sf_histint(vel,"n2",&nx))   sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&z0)) z0=0.;
    if (!sf_histfloat(vel,"o2",&x0)) x0=0.;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    /* If y, the input is velocity; if n, slowness */
    if(!sf_getint("order",&order)) order=4;
    /* Interpolation accuracy */
    if (!sf_getbool("lsint",&lsint)) lsint=false;
    /* if use least-squares interpolation */

    if (!sf_getint("nt",&nt)) nt=nx*nz;
    /* number of time steps */

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
	/* file with shot locations */
	shots = sf_input("shotfile");
	if (!sf_histint(shots,"n1",&ndim) || 2 != ndim) 
	    sf_error("Must have n1=2 in shotfile");
	if (!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
  
	s = sf_floatalloc2 (ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose (shots);
    } else {
	nshot = 1;
	ndim = 2;

	s = sf_floatalloc2 (ndim,nshot);

	if (!sf_getfloat("zshot",s[0]))   s[0][0]=0.;
	/* shot location in depth (if shotfile is not specified) */
	if (!sf_getfloat("yshot",s[0]+1)) s[0][1]=x0 + 0.5*(nx-1)*dx;
	/* shot location in lateral (if shotfile is not specified) */
	
	sf_warning("Shooting from z=%f, x=%f",s[0][0],s[0][1]);
    }

    if (NULL != sf_getstring("anglefile")) {
	/* file with initial angles */
	angles = sf_input("anglefile");

	if (!sf_histint(angles,"n1",&nr)) sf_error("No n1= in anglefile");
    } else {
	angles = NULL;

	if (!sf_getint("nr",&nr)) sf_error("Need nr=");
	/* Number of angles (if anglefile is not specified) */
	if (!sf_getfloat("a0",&a0)) a0 = 0.; 
	/* First angle in degrees (if anglefile is not specified) */
	if (!sf_getfloat("amax",&amax)) amax=360.;
	/* Maximum angle in degrees (if anglefile is not specified) */

	/* convert degrees to radians */
	a0 = a0*SF_PI/180.;
	amax = amax*SF_PI/180.;

	/* figure out angle spacing */
	da = (nr > 1)? (amax - a0)/(nr-1) : 0.;
    }

    a = sf_floatalloc(nr);
 
    /* specify output dimensions */
    nsr = nr*nshot;
    if (NULL != traj) fwrite(&nsr,sizeof(int),1,traj);
	    
    sf_putint(time,"n1",nr);
    sf_putint(time,"n2",nshot);

    /* get slowness */
    nm = nz*nx;
    slow = sf_floatalloc(nm);

    sf_floatread(slow,nm,vel);

    if (velocity) { /* convert to slowness */
	for(im = 0; im < nm; im++){
	    slow[im] = 1./slow[im];
	}
    }

    /* initialize ray tracing object */
    ct = sf_celltrace_init (lsint, order, nt, nz, nx, dz, dx, z0, x0, slow);
    
    free (slow);

    trj = sf_floatalloc2 (ndim,nt);
    t = sf_floatalloc(nr);

    for( is = 0; is < nshot; is++) { /* loop over shots */
	/* initialize angles */
	if (NULL != angles) {
	    sf_floatread(a,nr,angles);
	} else {
	    for (ir = 0; ir < nr; ir++) {
		a[ir] = a0+da*ir;
	    }
	}
	for (ir = 0; ir < nr; ir++) { /* loop over rays */
	    /* initialize position */
	    x[0] = s[is][0]; 
	    x[1] = s[is][1];

	    /* initialize direction */
	    p[0] = -cosf(a[ir]);
	    p[1] = sinf(a[ir]);

	    t[ir] = sf_cell_trace (ct, x, p, &it, trj);

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

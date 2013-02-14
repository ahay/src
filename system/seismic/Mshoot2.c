/* 2-D ray shooting. */
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
#include <float.h>

#include <rsf.h>

static sf_celltrace ct;
static float xs[2], xr, t;

static float shooting(float a);

int main(int argc, char* argv[])
{
    bool velocity, lsint;
    int is, nz, nx, im, nm, order, nshot, ndim, nt, nr, ir, i1, i2, ia;
    float a, b, fb, r0, dr, rmax, xmax, a0, amax, da, r1, r2, tol;
    float dz, dx, z0, x0, *angle, *slow, **shot, *time;
    sf_file shots, vel, out;

    sf_init (argc,argv);
    vel = sf_input("in");
    out = sf_output("out");

    /* get 2-D grid parameters */
    if (!sf_histint(vel,"n1",&nz))   sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&nx))   sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&z0)) z0=0.;
    if (!sf_histfloat(vel,"o2",&x0)) x0=0.;
    xmax = x0 + (nx-1)*dx;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    /* If y, the input is velocity; if n, slowness */
    if(!sf_getint("order",&order)) order=4;
    /* Interpolation order */
    if (!sf_getbool("lsint",&lsint)) lsint=false;
    /* if use least-squares interpolation */

    if (!sf_getint("nt",&nt)) nt=nx*nz;
    /* Maximum number of time steps */

    if (!sf_getint("nr",&nr)) nr=1;
    /* number of recievers */
    if (!sf_getfloat("r0",&r0)) r0=x0;
    /* first receiver */
    if (!sf_getfloat("dr",&dr)) dr=dx;
    /* receiver increment */
    rmax = r0 + (nr-1)*dr;

    /* sanity check */
    if (dr <= 0.) sf_error ("Need positive dr");
    if (r0 < x0 || r0 > xmax || rmax < x0 || rmax > xmax)
	sf_error ("Range of receivers [%f,%f] "
		  "is outside of the model range [%f,%f]",
		  r0,rmax,x0,xmax);
 
    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
	/* file with shot locations */
	shots = sf_input("shotfile");
	if (!sf_histint(shots,"n1",&ndim) || 2 != ndim) 
	    sf_error("Must have n1=2 in shotfile");
	if (!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
  
	shot = sf_floatalloc2 (ndim,nshot);
	sf_floatread(shot[0],ndim*nshot,shots);
	sf_fileclose (shots);
    } else {
	nshot = 1;
	ndim = 2;

	shot = sf_floatalloc2 (ndim,nshot);

	if (!sf_getfloat("zshot",shot[0]))   shot[0][0]=0.;
	/* shot coordinates (if no shotfile) */
	if (!sf_getfloat("yshot",shot[0]+1)) shot[0][1]=x0 + 0.5*(nx-1)*dx;
	
	sf_warning("Shooting from z=%f, x=%f",shot[0][0],shot[0][1]);
    }

    /* specify output dimensions */
    sf_putint (out,"n1",nr);
    sf_putfloat (out,"d1",dr);
    sf_putfloat (out,"o1",r0);
    sf_putint (out,"n2",nshot);
	    
    /* get velocity */
    nm = nz*nx;
    slow = sf_floatalloc(nm);

    sf_floatread(slow,nm,vel);

    if (vel) { /* convert to slowness */
	for(im = 0; im < nm; im++){
	    slow[im] = 1./slow[im];
	}
    }

    /* initialize ray tracing object */
    ct = sf_celltrace_init (lsint, order, nt, nz, nx, dz, dx, z0, x0, slow);
    
    free (slow);

    angle = sf_floatalloc (nr);
    time  = sf_floatalloc (nr);

    if (!sf_getfloat("tol",&tol)) tol=0.01;
    /* Shooting tolerance (in degrees) */
    tol *= SF_PI/180.; /* 1/100 degree */

    for( is = 0; is < nshot; is++) { /* loop over shots */
	/* initialize position */
	xs[0] = shot[is][0]; 
	xs[1] = shot[is][1];

	/* initialize directions */
	da = 0.01*SF_PI;
	a0   = atan2f(r0  -xs[1],xs[0]-z0);
	amax = atan2f(rmax-xs[1],xs[0]-z0);

	xr = r0;
	r1 = 0.;
	for (a0 -= da; a0 > - SF_PI; a0 -= da) {
	    r1 = shooting(a0);
	    if (r1 <= 0.) break;
	}
	r1 += xr;

	xr = rmax;
	for (amax += da; amax < SF_PI; amax += da) {
	    if (shooting(amax) >= 0.) break;
	}
	
	for (ir=0; ir < nr; ir++) {
	    time[ir] = angle[ir] = FLT_MAX;
	}

	da = (amax-a0)/nr;
	for (ia=0; ia < nr; ia++) {
	    a = a0+(ia+1)*da;
	    r2 = xr+shooting(a);
	    if (r1 < r2) {
		i1 = ceilf((r1-r0)/dr);
		i2 = floorf((r2-r0)/dr);
	    } else {
		i1 = ceilf((r2-r0)/dr);
		i2 = floorf((r1-r0)/dr);
	    }
	    for (ir=i1; ir <= i2; ir++) {
		if (ir >= nr) break;
		if (ir < 0) continue;
		xr = r0+ir*dr; /* target */

		b = sf_zero(shooting,a-da,a,r1-xr,r2-xr,tol,false);
		fb = shooting(b);

		if (fabsf(fb) > 0.5*dr) 
		    sf_warning("insufficient accuracy=%f "
			       "on shot %d, receiver %d",fb,is+1,ir+1);
 
		if (t < time[ir]) { /* select first arrivals */
		    angle[ir] = b;
		    time[ir] = t;
		}
	    }
	    r1 = r2;
      	}
	sf_floatwrite(angle,nr,out);
    }


    exit (0);
}

static float shooting(float a)
{
    int it;
    float x[2], p[2];

    x[0]=xs[0];
    x[1]=xs[1];
    
    p[0] = -cosf(a);
    p[1] = sinf(a);

    t = sf_cell_trace (ct, x, p, &it, NULL);

    return (x[1]-xr);
}

/* 	$Id$	 */

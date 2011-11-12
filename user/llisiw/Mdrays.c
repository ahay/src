/* 2D dynamic ray tracing by a Runge-Kutta integrator. 
 Angle is 90 deg along y-axis and 180 deg along z-axis, clockwise.
*/
/*
  Copyright (C) 2011 University of Texas at Austin

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

#include "draytrace.h"

int main(int argc, char* argv[])
{
    bool velocity, verb;
    int is, n[2], im, order, nshot, ndim;
    int nt, nt1, nr, ir, it, i, iz, iy, x0, *m;
    float dt, da=0., a0, amax, v0, deg2rad, shift, real, imag;
    float x[2], p[2], d[2], o[2], **traj, **dire, *slow, **s, *a, tempx[2];
    sf_complex ***dynaM, ***dynaN, ***dynaK, **ctime;
    raytrace rt;
    sf_file shots, vel, rays, angles, dmat, gbeam, proj, mask;

    sf_init (argc,argv);
    vel = sf_input("in");
    gbeam = sf_output("out");

    /* get 2-D grid parameters */
    if (!sf_histint(vel,"n1",n))     sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",n+1))   sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",d))   sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",o))   o[0]=0.;
    if (!sf_histfloat(vel,"o2",o+1)) o[1]=0.;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    /* If y, input is velocity; if n, slowness */
    if(!sf_getint("order",&order)) order=4;
    /* Interpolation order */

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* Number of time steps */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* Sampling in time */

    if(!sf_getbool("verb",&verb)) verb=true;
    /* Verbosity flag */

    if (!sf_getfloat("shift",&shift)) shift=0.5;
    /* Complex source shift */

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
        /* file with shot locations [zshot,yshot,nshot] */
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

	/* shot coordinates (if no shotfile) */
	if (!sf_getfloat("zshot",&s[0][0])) s[0][0]=0.;
        if (!sf_getfloat("yshot",&s[0][1])) s[0][1]=o[1]+0.5*(n[1]-1)*d[1];
    }

    deg2rad = SF_PI/180.;

    if (NULL != sf_getstring("anglefile")) {
        /* file with initial angles [nr,nshot] */
        angles = sf_input("anglefile");

        if (!sf_histint(angles,"n1",&nr)) sf_error("No n1= in anglefile");
    } else {
        angles = NULL;

        if (!sf_getint("nr",&nr)) sf_error("Need nr=");
        /* number of angles (if no anglefile) */
        if (!sf_getfloat("a0",&a0)) a0 = 0.;
        /* minimum angle (if no anglefile) */
        if (!sf_getfloat("amax",&amax)) amax=360.;
        /* maximum angle (if no anglefile) */

        /* convert degrees to radians */
        a0   *= deg2rad;
        amax *= deg2rad;

        /* figure out angle spacing */
        da = (nr>1)? (amax-a0)/(nr-1): 0.;
    }

    a = sf_floatalloc(nr);

    /* specify output dimensions */
    nt1 = nt+1;

    /* output rays contains only trajectory of central rays */
    if (NULL != sf_getstring("rays")) {
	    rays = sf_output("rays");
	    sf_putint(rays,"n1",nt1);
	    sf_putint(rays,"n2",nr);
	    if (nshot>1) sf_putint(rays,"n3",nshot);
	    sf_putfloat(rays,"o1",0.);
	    sf_putfloat(rays,"d1",dt);
	    if (NULL == angles) {
		sf_putfloat(rays,"o2",a0/deg2rad);
		sf_putfloat(rays,"d2",da/deg2rad);
	    }
	    sf_putstring(rays,"label1","time");
	    sf_putstring(rays,"label2","angle");
	    sf_putstring(rays,"unit2","deg");
	    sf_putstring(rays,"label3","shot");
	    sf_settype(rays,SF_COMPLEX);
	    sf_fileflush(rays,NULL);
	    sf_settype(rays,SF_FLOAT);
    } else {
	rays = NULL;
    }

    /* output dynaK contains complex matrix along each central ray */
    if (NULL != sf_getstring("dmat")) {
	dmat = sf_output("dmat");
	sf_putint(dmat,"n1",2);
	sf_putint(dmat,"n2",2);
	sf_putint(dmat,"n3",nt1);
	sf_putint(dmat,"n4",nr);
	if (nshot>1) sf_putint(dmat,"n5",nshot);
	sf_putfloat(dmat,"o3",0.);
	sf_putfloat(dmat,"d3",dt);
	if (NULL == angles) {
	    sf_putfloat(dmat,"o4",a0/deg2rad);
	    sf_putfloat(dmat,"d4",da/deg2rad);
	}
	sf_putstring(dmat,"label3","time");
	sf_putstring(dmat,"label4","angle");
	sf_putstring(dmat,"unit4","deg");
	sf_putstring(dmat,"label5","shot");
	sf_settype(dmat,SF_COMPLEX);
    } else {
	dmat = NULL;
    }

    /* output proj contains projection of grid points on central ray */
    if (NULL != sf_getstring("proj")) {
	proj = sf_output("proj");
	sf_putint(proj,"n3",nr);
	if (nshot>1) sf_putint(proj,"n4",nshot);
	if (NULL == angles) {
	    sf_putfloat(proj,"o3",a0/deg2rad);
	    sf_putfloat(proj,"d3",da/deg2rad);
	}
	sf_putstring(proj,"label3","angle");
	sf_putstring(proj,"unit3","deg");
	sf_putstring(proj,"label4","shot");
	sf_settype(proj,SF_INT);
    } else {
	proj = NULL;
    }

    /* output mask contains projection of central ray onto grid */
    if (NULL != sf_getstring("mask")) {
	mask = sf_output("mask");
	sf_putint(mask,"n3",nr);
	if (nshot>1) sf_putint(mask,"n4",nshot);
	if (NULL == angles) {
	    sf_putfloat(mask,"o3",a0/deg2rad);
	    sf_putfloat(mask,"d3",da/deg2rad);
	}
	sf_putstring(mask,"label3","angle");
	sf_putstring(mask,"unit3","deg");
	sf_putstring(mask,"label4","shot");
	sf_settype(mask,SF_INT);
    } else {
	mask = NULL;
    }

    /* output gbeam contains complex traveltime for each central ray */
    sf_putint(gbeam,"n3",nr);
    if (nshot>1) sf_putint(gbeam,"n4",nshot);
    if (NULL == angles) {
        sf_putfloat(gbeam,"o3",a0/deg2rad);
        sf_putfloat(gbeam,"d3",da/deg2rad);
    }
    sf_putstring(gbeam,"label3","angle");
    sf_putstring(gbeam,"unit3","deg");
    sf_putstring(gbeam,"label4","shot");
    sf_settype(gbeam,SF_COMPLEX);

    /* get slowness squared */
    slow = sf_floatalloc(n[0]*n[1]);

    sf_floatread(slow,n[0]*n[1],vel);

    for(im = 0; im < n[0]*n[1]; im++){
        v0 = slow[im];
        slow[im] = velocity? 1./(v0*v0): v0*v0;
    }

    /* initialize ray tracing object */
    rt = trace_init (2,nt,dt,n,o,d,slow,order);

    /* safe to free after trace_init */
    free(slow);
    
    /* allocate temporary memory */
    traj = sf_floatalloc2(ndim,nt1);
    dire = sf_floatalloc2(ndim,nt1);
    dynaM = sf_complexalloc3(ndim,ndim,nt1);
    dynaN = sf_complexalloc3(ndim,ndim,nt1);
    dynaK = sf_complexalloc3(ndim,ndim,nt1);
    ctime = sf_complexalloc2(n[0],n[1]);

    for( is = 0; is < nshot; is++) { /* loop over shots */
        /* initialize angles */
        if (NULL != angles) {
            sf_floatread(a,nr,angles);
        } else {
            for (ir = 0; ir < nr; ir++) {
                a[ir] = a0+da*ir;
            }
        }

        if (verb) sf_warning("Shooting from (z=%g, x=%g):",s[is][0],s[is][1]);

        for (ir = 0; ir < nr; ir++) { /* loop over rays */
            /* initialize position */
            x[0] = s[is][0]; 
            x[1] = s[is][1];

            /* initialize direction */
            p[0] = -cosf(a[ir]);
            p[1] = sinf(a[ir]);

            it = trace_ray (rt,x,p,shift,traj,dire,dynaM,dynaN);
            if (it < 0) it = -it; /* keep side/buttom exiting rays */

	    if (verb) sf_warning("Ray angle=%g exit at t=%g",a[ir],it==0?nt*dt:it*dt);
            
	    /* write central ray trajectory and/or associated complex matrix */
	    for (i=0; i < nt1; i++) {
		if (0==it || i < it) {
		    dray_assemble (dynaM[i],dynaN[i],dynaK[i]);

		    if (NULL != rays) sf_floatwrite(traj[i],ndim,rays);
		    if (NULL != dmat) sf_complexwrite(dynaK[i][0],ndim*ndim,dmat);
		} else {
		    dray_assemble (dynaM[it],dynaN[it],dynaK[i]);

		    if (NULL != rays) sf_floatwrite(traj[it],ndim,rays);
		    if (NULL != dmat) sf_complexwrite(dynaK[i][0],ndim*ndim,dmat);
		}
	    }
	    
	    /* write central ray mask */
	    if (NULL != mask) {
		m = sf_intalloc(n[0]*n[1]);
		
		for (i=0; i < n[0]*n[1]; i++) m[i] = 0;
		dray_central (traj,it==0?nt1:it,o,d,n,m);

		sf_intwrite (m,n[0]*n[1],mask);
	    }

	    /* construct Gaussian beam */
	    for (iy=0; iy < n[1]; iy++) {
		tempx[1] = o[1]+iy*d[1];

		for (iz=0; iz < n[0]; iz++) {
		    tempx[0] = o[0]+iz*d[0];

		    /* search for projection on central ray, return projection */
		    x0 = dray_search (traj,it==0?nt1:it,tempx);
		    
		    if (NULL != proj) sf_intwrite(&x0,1,proj);

		    /* quadratic function away from central ray */
		    /* T(x) = T(x0)+p(x0)*(x-x0)+1/2*(x-x0)*K(x0)*(x-x0) */
		    real = x0*dt + dire[x0][0]*(tempx[0]-traj[x0][0])+dire[x0][1]*(tempx[1]-traj[x0][1])+
			0.5*crealf((tempx[0]-traj[x0][0])*(tempx[0]-traj[x0][0])*dynaK[x0][0][0]+
				   (tempx[0]-traj[x0][0])*(tempx[1]-traj[x0][1])*dynaK[x0][0][1]+
				   (tempx[1]-traj[x0][1])*(tempx[0]-traj[x0][0])*dynaK[x0][1][0]+
				   (tempx[1]-traj[x0][1])*(tempx[1]-traj[x0][1])*dynaK[x0][1][1]);
		    
		    imag = 0.+
			0.5*cimagf((tempx[0]-traj[x0][0])*(tempx[0]-traj[x0][0])*dynaK[x0][0][0]+
				   (tempx[0]-traj[x0][0])*(tempx[1]-traj[x0][1])*dynaK[x0][0][1]+
				   (tempx[1]-traj[x0][1])*(tempx[0]-traj[x0][0])*dynaK[x0][1][0]+
				   (tempx[1]-traj[x0][1])*(tempx[1]-traj[x0][1])*dynaK[x0][1][1]);
		    
		    ctime[iy][iz] = sf_cmplx(real,imag);
		}
	    }

	    sf_complexwrite (ctime[0],n[0]*n[1],gbeam);
	}
    }

    if (verb) sf_warning(".");

    exit (0);
}

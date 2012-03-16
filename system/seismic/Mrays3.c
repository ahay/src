/* Ray tracing by a Runge-Kutta integrator in 3-D.

Takes: > rays.rsf

Rays can be plotted with sfplotrays.
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

#include "raytrace.h"

int main(int argc, char* argv[])
{
    bool velocity, sym, escvar;
    int is, n[3], im, nm, order, nshot, ndim, two, na, nb;
    int nt, nt1, nr, ir, it, i, ia, ib;
    float t, dt, da=0., a0, amax, db=0., b0, bmax, v0;
    float x[3], p[3], d[3], o[3], **traj, *slow, **s, *a, *b;
    raytrace rt;
    sf_file shots, vel, rays, angles;

    sf_init (argc,argv);
    vel = sf_input("in");
    rays = sf_output("out");

    /* get 2-D grid parameters */
    if (!sf_histint(vel,"n1",n))     sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",n+1))   sf_error("No n2= in input");
    if (!sf_histint(vel,"n3",n+2))   sf_error("No n3= in input");
    if (!sf_histfloat(vel,"d1",d))   sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"d3",d+2)) sf_error("No d3= in input");
    if (!sf_histfloat(vel,"o1",o))   o[0]=0.;
    if (!sf_histfloat(vel,"o2",o+1)) o[1]=0.;
    if (!sf_histfloat(vel,"o3",o+2)) o[2]=0.;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    /* If y, input is velocity; if n, slowness */
    if(!sf_getint("order",&order)) order=4;
    /* Interpolation order */

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* Number of time steps */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* Sampling in time */

    if (!sf_getbool("sym",&sym)) sym=true;
    /* if y, use symplectic integrator */

    if(!sf_getbool("escvar",&escvar)) escvar=false;
    /* If y - output escape values, n - trajectories */

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
        /* file with shot locations */
        shots = sf_input("shotfile");
        if (!sf_histint(shots,"n1",&ndim) || 3 != ndim) 
            sf_error("Must have n1=2 in shotfile");
        if (!sf_histint(shots,"n2",&nshot)) 
            sf_error("No n2= in shotfile");
  
        s = sf_floatalloc2 (ndim,nshot);
        sf_floatread(s[0],ndim*nshot,shots);
        sf_fileclose (shots);
    } else {
        nshot = 1;
        ndim = 3;

        s = sf_floatalloc2 (ndim,nshot);

        if (!sf_getfloat("zshot",&s[0][0])) s[0][0]=o[0];
        /* shot location in depth (if shotfile is not specified) */
        if (!sf_getfloat("yshot",&s[0][1])) s[0][1]=o[1] + 0.5*(n[1]-1)*d[1];
        /* shot location inline (if shotfile is not specified) */
        if (!sf_getfloat("xshot",&s[0][2])) s[0][2]=o[2] + 0.5*(n[2]-1)*d[2];
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
    nt1 = nt+1;
    if (escvar) {
        sf_putint (rays, "n1", 4);
        sf_putfloat (rays, "o1", 0.0);
        sf_putfloat (rays, "d1", 1.0);
        sf_putstring (rays, "label1", "Escape variable");
        sf_putstring (rays, "unit1", "");
        if (NULL == angles) {
            sf_putint (rays, "n2", na);
            sf_putfloat (rays, "d2", da*180.0/SF_PI);
            sf_putfloat (rays, "o2", a0*180.0/SF_PI);
            sf_putstring (rays, "label2", "Azimuth");
            sf_putstring (rays, "unit2", "Degrees");
            sf_putint (rays, "n3", nb);
            sf_putfloat (rays, "d3", db*180.0/SF_PI);
            sf_putfloat (rays, "o3", b0*180.0/SF_PI);
            sf_putstring (rays, "label3", "Inclination");
            sf_putstring (rays, "unit3", "Degrees");
            sf_putint (rays,"n4",nshot);
            sf_putfloat (rays, "d4", 1.0);
            sf_putfloat (rays, "o4", 0.0);
            sf_putstring (rays, "label4", "Shots");
            sf_putstring (rays, "unit4", "");
        } else {
            sf_putint (rays,"n2",nr);
            sf_putfloat (rays, "d2", 1.0);
            sf_putfloat (rays, "o2", 0.0);
            sf_putstring (rays, "label2", "Angles");
            sf_putstring (rays, "unit2", "");
            sf_putint (rays,"n3",nshot);
            sf_putfloat (rays, "d3", 1.0);
            sf_putfloat (rays, "o3", 0.0);
            sf_putstring (rays, "label3", "Shots");
            sf_putstring (rays, "unit3", "");
        }
    } else {
        sf_putint (rays,"n1",ndim);
        sf_putint (rays,"n2",nt1);
        sf_putint (rays,"n3",nr);
        sf_putint (rays,"n4",nshot);
    }
            
    /* get slowness squared */
    nm = n[0]*n[1]*n[2];
    slow = sf_floatalloc(nm);

    sf_floatread(slow,nm,vel);

    for(im = 0; im < nm; im++){
        v0 = slow[im];
        slow[im] = velocity? 1./(v0*v0): v0*v0;
    }

    /* initialize ray tracing object */
    rt = raytrace_init (3, sym, nt, dt, n, o, d, slow, order);
    
    free (slow);

    traj = sf_floatalloc2 (sym? ndim: 2*ndim,nt1);

    for( is = 0; is < nshot; is++) { /* loop over shots */
        /* initialize angles */
        if (NULL != angles) {
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

            it = trace_ray (rt, x, p, traj);
            if (it < 0) it = -it; /* keep side-exiting rays */

            if (escvar) {
                /* Write escape variables only */
                sf_floatwrite (traj[it],ndim,rays); /* z, x, y */
                t = it*dt; /* t */
                sf_floatwrite (&t, 1, rays);
            } else {
                for (i=0; i < nt1; i++) {
                    if (0==it || it > i) {
                        sf_floatwrite (traj[i],ndim,rays);
                    } else {
                        sf_floatwrite (traj[it],ndim,rays);
                    }
                }
            }
        }
    }

    exit (0);
}

/*         $Id$         */

/* Dynamic ray tracing by a Runge-Kutta integrator. */
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
    bool velocity, sym, verb, escvar;
    int is, n[2], im, nm, order, nshot, ndim;
    int nt, nt1, nr, ir, it, i;
    float t, dt, da=0., a0, amax, v0, deg2rad, shift;
    float x[2], p[2], d[2], o[2], **traj, *slow, **s, *a;
    sf_complex ***dynaM, ***dynaN;
    raytrace rt;
    sf_file shots, vel, rays, angles, dynaP, dynaX;

    sf_init (argc,argv);
    vel = sf_input("in");
    rays = sf_output("out");
    dynaP = sf_output("dynaP");
    dynaX = sf_output("dynaX");

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

    if (!sf_getbool("sym",&sym)) sym=true;
    /* if y, use symplectic integrator */

    if(!sf_getbool("verb",&verb)) verb=true;
    /* Verbosity flag */
    if(!sf_getbool("escvar",&escvar)) escvar=false;
    /* If y - output escape values, n - trajectories */

    if (!sf_getfloat("shift",&shift)) shift=0.5;
    /* Complex source shift */

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
        /* file with shot locations */
        shots = sf_input("shotfile");
        if (!sf_histint(shots,"n1",&ndim) || 2 != ndim) 
            sf_error("Must have n1=2 in shotfile");
        if (!sf_histint(shots,"n2",&nshot))
            sf_error("No n2= in shotfile");

        if (sf_histfloat(shots,"o2",&t)) sf_putfloat(rays,"o3",t);
        if (sf_histfloat(shots,"d2",&t)) sf_putfloat(rays,"d3",t);

        s = sf_floatalloc2 (ndim,nshot);
        sf_floatread(s[0],ndim*nshot,shots);
        sf_fileclose (shots);
    } else {
        nshot = 1;
        ndim = 2;

        s = sf_floatalloc2 (ndim,nshot);

        if (!sf_getfloat("zshot",&s[0][0]))   s[0][0]=0.;
        /* shot coordinates (if no shotfile) */
        if (!sf_getfloat("yshot",&s[0][1])) s[0][1]=o[1] + 0.5*(n[1]-1)*d[1];
    }

    deg2rad = SF_PI/180.;

    if (NULL != sf_getstring("anglefile")) {
        /* file with initial angles */
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
        da = (nr > 1)? (amax - a0)/(nr-1) : 0.;
    }

    a = sf_floatalloc(nr);

    /* specify output dimensions */
    nt1 = nt+1;
    sf_putint (rays,"n1",escvar ? 4 : nt1);
    sf_putint (rays,"n2",nr);
    if( nshot > 1 ) sf_putint (rays,"n3",nshot);
    sf_putfloat(rays,"o1",0.);
    sf_putfloat(rays,"d1",dt);
    if (NULL == angles) {
        sf_putfloat(rays,"o2",a0/deg2rad);
        sf_putfloat(rays,"d2",da/deg2rad);
    }
    sf_putstring(rays,"label1","time");
    sf_putstring(rays,"label2","angle");
    sf_putstring(rays,"unit2","deg");
    if (!escvar)
        sf_settype (rays,SF_COMPLEX);
    sf_fileflush (rays,NULL);
    sf_settype (rays,SF_FLOAT);

    sf_putint (dynaP,"n1",2);
    sf_putint (dynaP,"n2",2);
    sf_putint (dynaP,"n3",escvar ? 4 : nt1);
    sf_putint (dynaP,"n4",nr);
    if( nshot > 1 ) sf_putint (dynaP,"n5",nshot);
    sf_putfloat(dynaP,"o3",0.);
    sf_putfloat(dynaP,"d3",dt);
    if (NULL == angles) {
        sf_putfloat(dynaP,"o4",a0/deg2rad);
        sf_putfloat(dynaP,"d4",da/deg2rad);
    }
    sf_putstring(dynaP,"label3","time");
    sf_putstring(dynaP,"label4","angle");
    sf_putstring(dynaP,"unit4","deg");
    sf_settype (dynaP,SF_COMPLEX);

    sf_putint (dynaX,"n1",2);
    sf_putint (dynaX,"n2",2);
    sf_putint (dynaX,"n3",escvar ? 4 : nt1);
    sf_putint (dynaX,"n4",nr);
    if( nshot > 1 ) sf_putint (dynaX,"n5",nshot);
    sf_putfloat(dynaX,"o3",0.);
    sf_putfloat(dynaX,"d3",dt);
    if (NULL == angles) {
        sf_putfloat(dynaX,"o4",a0/deg2rad);
        sf_putfloat(dynaX,"d4",da/deg2rad);
    }
    sf_putstring(dynaX,"label3","time");
    sf_putstring(dynaX,"label4","angle");
    sf_putstring(dynaX,"unit4","deg");
    sf_settype (dynaX,SF_COMPLEX);

    /* get slowness squared */
    nm = n[0]*n[1];
    slow = sf_floatalloc(nm);

    sf_floatread(slow,nm,vel);

    for(im = 0; im < nm; im++){
        v0 = slow[im];
        slow[im] = velocity? 1./(v0*v0): v0*v0;
    }

    /* initialize ray tracing object */
    rt = raytrace_init (2, sym, nt, dt, n, o, d, slow, order);

    free (slow);
    
    traj = sf_floatalloc2 (sym? ndim: 2*ndim,nt1);
    dynaM = sf_complexalloc3 (sym? ndim: 2*ndim,sym? ndim: 2*ndim,nt1);
    dynaN = sf_complexalloc3 (sym? ndim: 2*ndim,sym? ndim: 2*ndim,nt1);

    for( is = 0; is < nshot; is++) { /* loop over shots */
        /* initialize angles */
        if (NULL != angles) {
            sf_floatread(a,nr,angles);
        } else {
            for (ir = 0; ir < nr; ir++) {
                a[ir] = a0+da*ir;
            }
        }
        if (verb)
            sf_warning ("Shooting from z=%g, x=%g;", s[is][0], s[is][1]);

        for (ir = 0; ir < nr; ir++) { /* loop over rays */
            /* initialize position */
            x[0] = s[is][0]; 
            x[1] = s[is][1];

            /* initialize direction */
            p[0] = -cosf(a[ir]);
            p[1] = sinf(a[ir]);

            it = trace_ray (rt, x, p, shift, traj, dynaM, dynaN);
            if (it < 0) it = -it; /* keep side-exiting rays */

            if (escvar) {
                /* Write escape variables only */
                sf_floatwrite (traj[it],ndim,rays); /* z, x */
                t = it*dt; /* t */
                sf_floatwrite (&t,1,rays);
                /* a */
                if (it > 0) {
                    i = it >= 2 ? it - 2 : it - 1;
                    /* Escape vector */
                    x[0] = traj[it][0];
                    x[1] = traj[it][1];
                    x[0] -= traj[i][0];
                    x[1] -= traj[i][1];
                    /* Dot product with unit vector pointing upward */
                    t = sqrt(x[0]*x[0] + x[1]*x[1]); /* Length */
                    t = acos(-x[0]/t);
                    if (x[1] < 0) t = -t;
                } else
                    t = a[ir];
                t /= deg2rad;
                sf_floatwrite (&t,1,rays);
            } else {
                /* Write full trajectory */
                for (i=0; i < nt1; i++) {
                    if (0==it || it > i) {
                        sf_floatwrite (traj[i],ndim,rays);
			sf_complexwrite (dynaM[i][0],2*2,dynaP);
			sf_complexwrite (dynaN[i][0],2*2,dynaX);
                    } else {
                        sf_floatwrite (traj[it],ndim,rays);
			sf_complexwrite (dynaM[it][0],2*2,dynaP);
			sf_complexwrite (dynaN[it][0],2*2,dynaX);
                    }
                }
            }
        }
    }
    if (verb)
        sf_warning (".");

    exit (0);
}

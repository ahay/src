/* Escape tables by ray tracing with escape equations in 2-D. */
/*
  Copyright (C) 2012 University of Texas at Austin

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <rsf.h>

#include "esc_point2.h"
#include "esc_tracer2.h"

#define TRAJ2_COMPS 4
#define TRAJ2_Z 0
#define TRAJ2_X 1
#define TRAJ2_A 2
#define TRAJ2_T 3

typedef struct {
    int       it, nt;
    float     dt;
    float   **pnts;
    sf_file   traj;
} sf_escrt2_traj_cbud;

/* Callback to output points along a ray trajectory; this one is
   called by sf_esc_tracer2 at every step during ray tracing in phase space;
   each next call to this routine can be more separated by more than one
   dt step in time with the previous one; we use interpolation in here
   to fill in the intermediate trajectory points */
void sf_escrt2_traj (float z, float x, float a, int it, void *ud) {
    int i, j, n;
    float z0, x0, a0, dz, dx, da, a1;
    sf_escrt2_traj_cbud *cbud = (sf_escrt2_traj_cbud*)ud;

    if (it < cbud->nt) {
        n = it - cbud->it;
        if (it != 0 && n > 1) {
            /* Interpolate values inbetween the last one
               and the new one */
            z0 = cbud->pnts[cbud->it][TRAJ2_Z];
            x0 = cbud->pnts[cbud->it][TRAJ2_X];
            a0 = cbud->pnts[cbud->it][TRAJ2_A];
            /* Make sure that there is no -pi->pi
               jump between the two points in the angle dimension */
            if (a0 > SF_PI/2.0 && a < -SF_PI/2.0)
                a0 -= 2.0*SF_PI;
            else if (a0 < -SF_PI/2.0 && a > SF_PI/2.0)
                a0 += 2.0*SF_PI;
            dz = (z - z0)/(float)n;
            dx = (x - x0)/(float)n;
            da = (a - a0)/(float)n;
            for (i = cbud->it + 1; i < it; i++) {
                j = i - cbud->it;
                cbud->pnts[i][TRAJ2_Z] = z0 + j*dz;
                cbud->pnts[i][TRAJ2_X] = x0 + j*dx;
                a1 = a0 + j*da;
                /* Keep angle in [-pi;pi] range */
                if (a1 < -SF_PI)
                    a1 += 2.0*SF_PI;
                else if (a1 > SF_PI)
                    a1 -= 2.0*SF_PI;
                cbud->pnts[i][TRAJ2_A] = a1;
            }
        }
        cbud->pnts[it][TRAJ2_Z] = z;
        cbud->pnts[it][TRAJ2_X] = x;
        cbud->pnts[it][TRAJ2_A] = a;
        cbud->it = it;
    }
}

static char* sf_escrt3_warnext (sf_file input) {
    int icpu = 0, ncpu = 1, maxlen;
    char *host = sf_gethost ();
    char *ext = NULL;

    if (input) {
        if (!sf_histint (input, "icpu", &icpu)) icpu = 0;
        if (!sf_histint (input, "ncpu", &ncpu)) ncpu = 1;
    }
    maxlen = strlen (host) + 30;
    ext = (char*)malloc (maxlen*sizeof(char));
    snprintf (ext, maxlen, "[%s:%d/%d]", host, icpu + 1, ncpu);
    return ext;
}

int main (int argc, char* argv[]) {
    int nz, nx, na, ia, ix, iz, i, it, nt, ic, nc = 1, fz, lz, itr;
    float dz, oz, dx, ox, da, oa, z, x, dt, df, md, aper;
    float ***e;
    sf_file spdom, vspline = NULL, out, traj = NULL;
    sf_escrt2_traj_cbud *tdata = NULL; 
    char *ext = NULL;
    bool verb, parab;
    sf_esc_slowness2 esc_slow;
    sf_esc_tracer2 *esc_tracers;
    sf_esc_point2 *esc_points;
    sf_timer timer;

    sf_init (argc, argv);

    if (!sf_stdin ()) {
        spdom = NULL;
    } else {
        spdom = sf_input ("in");
        /* Spatial (z,x) domain */
    }

    out = sf_output ("out");
    /* Escape values */

    /* Spatial dimensions */
    if (spdom) {
        if (!sf_histint (spdom, "n1", &nz)) sf_error ("No n1= in input");
        if (!sf_histint (spdom, "n2", &nx)) sf_error ("No n2= in input");
        if (!sf_histfloat (spdom, "d1", &dz)) sf_error ("No d1= in input");
        if (!sf_histfloat (spdom, "o1", &oz)) sf_error ("No o1= in input");
        if (!sf_histfloat (spdom, "d2", &dx)) sf_error ("No d2= in input");
        if (!sf_histfloat (spdom, "o2", &ox)) sf_error ("No o2= in input");
    }
    ext = sf_escrt3_warnext (spdom);
    if (!sf_getint ("nz", &nz) && !spdom) sf_error ("Need nz=");
    /* Number of samples in z axis */
    if (!sf_getfloat ("oz", &oz) && !spdom) sf_error ("Need oz=");
    /* Beginning of z axis */
    if (!sf_getfloat ("dz", &dz) && !spdom) sf_error ("Need dz=");
    /* Sampling of z axis */
    if (!sf_getint ("nx", &nx) && !spdom) sf_error ("Need nx=");
    /* Number of samples in x axis */
    if (!sf_getfloat ("ox", &ox) && !spdom) sf_error ("Need ox=");
    /* Beginning of x axis */
    if (!sf_getfloat ("dx", &dx) && !spdom) sf_error ("Need dx=");
    /* Sampling of x axis */

    if (!sf_getint ("na", &na)) na = 360;
    /* Number of phase angles */
    da = 2.0*SF_PI/(float)na;
    oa = -SF_PI + 0.5*da;

    if (!sf_getfloat ("df", &df)) df = 0.25;
    /*< Maximum distance to travel per step (fraction of the cell size) >*/

    if (!sf_getfloat ("md", &md)) md = SF_HUGE;
    /* Maximum distance for a ray to travel (default - up to model boundaries) */
    if (md != SF_HUGE)
        md = fabsf (md);

    if (!sf_getfloat ("aper", &aper)) aper = SF_HUGE;
    /* Maximum aperture in x and y directions from current point (default - up to model boundaries) */
    if (aper != SF_HUGE)
        aper = fabsf (aper);

#ifdef _OPENMP
    if (!sf_getint ("nc", &nc)) nc = 0;
    /* Number of threads to use for ray tracing (OMP_NUM_THREADS by default) */
    if (nc)
        omp_set_num_threads (nc); /* User override */
    else
        nc = omp_get_max_threads (); /* Current default */
    sf_warning ("%s Using %d threads", ext, omp_get_max_threads ());
#endif

    if (!sf_getbool ("parab", &parab)) parab = true;
    /* y - use parabolic approximation of trajectories, n - straight line */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (sf_getstring ("traj")) {
        /* Trajectory output */
        traj = sf_output ("traj");
        if (!sf_getint ("nt", &nt)) nt = 1001;
        /* Number of time samples for each trajectory */
        if (!sf_getfloat ("dt", &dt)) dt = 0.001;
        /* Time sampling */
        tdata = (sf_escrt2_traj_cbud*)sf_alloc (nc*na, sizeof(sf_escrt2_traj_cbud));
        /* Time sampling */
        for (itr = 0; itr < nc*na; itr++) {
            tdata[itr].it = 0;
            tdata[itr].nt = nt;
            tdata[itr].dt = dt;
            tdata[itr].pnts = sf_floatalloc2 (TRAJ2_COMPS - 1, nt);
        }
    }

    e = sf_floatalloc3 (ESC2_NUM, na, nc);

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness2_init (vspline, verb);

    /* Make room for escape variables in output */
    if (spdom)
        sf_shiftdim2 (spdom, out, 1);

    sf_putint (out, "n1", ESC2_NUM);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Escape variable");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", na);
    sf_putfloat (out, "d2", da*180.0/SF_PI);
    sf_putfloat (out, "o2", oa*180.0/SF_PI);
    sf_putstring (out, "label2", "Angle");
    sf_putstring (out, "unit2", "Degrees");
    sf_putint (out, "n3", nz);
    sf_putfloat (out, "o3", oz);
    sf_putfloat (out, "d3", dz);
    if (!spdom) {
        sf_putstring (out, "label3", "Depth");
        sf_putstring (out, "unit3", "");
    }
    sf_putint (out, "n4", nx);
    sf_putfloat (out, "o4", ox);
    sf_putfloat (out, "d4", dx);
    if (!spdom) {
        sf_putstring (out, "label4", "Lateral");
        sf_putstring (out, "unit4", "");
    }

    if (traj) {
        if (spdom)
            sf_shiftdim2 (spdom, traj, 1);
        sf_putint (traj, "n1", TRAJ2_COMPS - 1);
        sf_putfloat (traj, "o1", 0.0);
        sf_putfloat (traj, "d1", 1.0);
        sf_putstring (traj, "label1", "Trajectory components");
        sf_putstring (traj, "unit1", "");
        sf_putint (traj, "n2", nt);
        sf_putfloat (traj, "o2", 0.0);
        sf_putfloat (traj, "d2", dt);
        sf_putstring (traj, "label2", "Time");
        sf_putstring (traj, "unit2", "s");
        sf_putint (traj, "n3", na);
        sf_putfloat (traj, "d3", da*180.0/SF_PI);
        sf_putfloat (traj, "o3", oa*180.0/SF_PI);
        sf_putstring (traj, "label3", "Angle");
        sf_putstring (traj, "unit3", "Degrees");
        sf_putint (traj, "n4", nz);
        sf_putfloat (traj, "o4", oz);
        sf_putfloat (traj, "d4", dz);
        if (!spdom) {
            sf_putstring (out, "label4", "Depth");
            sf_putstring (out, "unit4", "");
        }
        sf_putint (traj, "n5", nx);
        sf_putfloat (traj, "o5", ox);
        sf_putfloat (traj, "d5", dx);
        if (!spdom) {
            sf_putstring (out, "label5", "Lateral");
            sf_putstring (out, "unit5", "");
        }
    }

    esc_tracers = (sf_esc_tracer2*)sf_alloc (nc, sizeof(sf_esc_tracer2));
    esc_points = (sf_esc_point2*)sf_alloc (nc, sizeof(sf_esc_point2));
    for (ic = 0; ic < nc; ic++) {
        esc_tracers[ic] = sf_esc_tracer2_init (esc_slow);
        sf_esc_tracer2_set_parab (esc_tracers[ic], parab);
        if (md != SF_HUGE)
            sf_esc_tracer2_set_mdist (esc_tracers[ic], md);
        sf_esc_tracer2_set_df (esc_tracers[ic], df);
        esc_points[ic] = sf_esc_point2_init ();
    }

    timer = sf_timer_init ();
    /* Ray tracing loop */
    for (ix = 0; ix < nx; ix++) {
        x = ox + ix*dx;
        /* Set aperture */
        if (aper != SF_HUGE) {
            for (ic = 0; ic < nc; ic++) {
                sf_esc_tracer2_set_xmin (esc_tracers[ic], x - aper);
                sf_esc_tracer2_set_xmax (esc_tracers[ic], x + aper);
            }
        }
        if (verb)
            sf_warning ("%s Shooting from lateral location %d of %d at %g;",
                        ext, ix + 1, nx, x);
        /* Loop over chunks */
        for (ic = 0; ic < (nz/nc + ((nz % nc) != 0)); ic++) {
            fz = ic*nc;
            lz = (ic + 1)*nc - 1;
            if (lz >= nz)
                lz = nz - 1;
            sf_timer_start (timer);
#ifdef _OPENMP
#pragma omp parallel for                   \
            schedule(static,1)             \
            private(iz,ia,z,it,i,itr) \
            shared(fz,lz,ix,na,nz,nx,oa,oz,ox,da,dz,dx,x,tdata,esc_tracers,esc_points,e,out,traj)
#endif
            for (iz = fz; iz <= lz; iz++) {
                z = oz + iz*dz;
                for (ia = 0; ia < na; ia++) {
                    if (traj) {
                        itr = (iz - fz)*na + ia;
                        sf_esc_tracer2_set_trajcb (esc_tracers[iz - fz], sf_escrt2_traj, dt,
                                                   (void*)&tdata[itr]);
                    }
                    sf_esc_tracer2_compute (esc_tracers[iz - fz], z, x, oa + ia*da,
                                            0.0, 0.0, esc_points[iz - fz], NULL);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC2_NUM; i++)
                        e[iz - fz][ia][i] = sf_esc_point2_get_esc_var (esc_points[iz - fz], i);
                    if (traj) {
                        /* Fill the rest of the trajectory with the last point */
			itr = (iz - fz)*na + ia;
                        for (it = tdata[itr].it + 1; it < tdata[itr].nt; it++) {
                            for (i = 0; i < TRAJ2_COMPS - 1; i++)
                                tdata[itr].pnts[it][i] = tdata[itr].pnts[tdata[itr].it][i];
                        }
                    }
                } /* Loop over a */
            } /* Loop over z */
            sf_timer_stop (timer);
            sf_floatwrite (e[0][0], (size_t)(lz - fz + 1)*(size_t)na*(size_t)ESC2_NUM, out);
            if (tdata) {
                for (itr = 0; itr < (lz - fz + 1)*na; itr++) {
                    sf_floatwrite (tdata[itr].pnts[0],
                                   (size_t)tdata[itr].nt*(size_t)(TRAJ2_COMPS - 1), traj);
                }
            }
        } /* Loop over z chunks */
    } /* Loop over x */
    if (verb) {
        sf_warning (".");
        sf_warning ("%s Total kernel time: %g s, per depth point: %g s",
                    ext, sf_timer_get_total_time (timer)/1000.0,
                    (sf_timer_get_total_time (timer)/(float)((size_t)nx*(size_t)nz))/1000.0);
    }
    sf_timer_close (timer);

    for (ic = 0; ic < nc; ic++) {
        sf_esc_point2_close (esc_points[ic]);
        sf_esc_tracer2_close (esc_tracers[ic]);
    }
    free (esc_points);
    free (esc_tracers);
    if (traj) {
        for (itr = 0; itr < nc*na; itr++) {
            free (tdata[itr].pnts[0]);
            free (tdata[itr].pnts);
        }
        free (tdata);
    }
    sf_esc_slowness2_close (esc_slow);

    free (e[0][0]);
    free (e[0]);
    free (e);

    sf_fileclose (vspline);
    if (traj)
        sf_fileclose (traj);

    return 0;
}


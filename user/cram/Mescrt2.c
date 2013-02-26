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

int main (int argc, char* argv[]) {
    int nz, nx, na, ia, ix, iz, i, it;
    float dz, oz, dx, ox, da, oa, z, x;
    float **e;
    sf_file spdom, vspline = NULL, out;
    sf_escrt2_traj_cbud tdata; 

    bool verb, parab;
    sf_esc_slowness2 esc_slow;
    sf_esc_tracer2 esc_tracer;
    sf_esc_point2 esc_point;

    sf_init (argc, argv);

    if (!sf_stdin ()) {
        spdom = NULL;
    } else {
        spdom = sf_input ("in");
        /* Spatial (z,x) domain */
    }

    out = sf_output ("out");
    /* Escape values */

    tdata.traj = NULL;

    /* Spatial dimensions */
    if (spdom) {
        if (!sf_histint (spdom, "n1", &nz)) sf_error ("No n1= in input");
        if (!sf_histint (spdom, "n2", &nx)) sf_error ("No n2= in input");
        if (!sf_histfloat (spdom, "d1", &dz)) sf_error ("No d1= in input");
        if (!sf_histfloat (spdom, "o1", &oz)) sf_error ("No o1= in input");
        if (!sf_histfloat (spdom, "d2", &dx)) sf_error ("No d2= in input");
        if (!sf_histfloat (spdom, "o2", &ox)) sf_error ("No o2= in input");
    }
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

    if (!sf_getbool ("parab", &parab)) parab = true;
    /* y - use parabolic approximation of trajectories, n - straight line */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (sf_getstring ("traj")) {
        /* Trajectory output */
        tdata.traj = sf_output ("traj");
        if (!sf_getint ("nt", &tdata.nt)) tdata.nt = 1001;
        /* Number of time samples for each trajectory */
        if (!sf_getfloat ("dt", &tdata.dt)) tdata.dt = 0.001;
        /* Time sampling */
        tdata.it = 0;
        tdata.pnts = sf_floatalloc2 (TRAJ2_COMPS - 1, tdata.nt);
    }

    e = sf_floatalloc2 (ESC2_NUM, na);

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

    if (tdata.traj) {
        if (spdom)
            sf_shiftdim2 (spdom, tdata.traj, 1);
        sf_putint (tdata.traj, "n1", TRAJ2_COMPS - 1);
        sf_putfloat (tdata.traj, "o1", 0.0);
        sf_putfloat (tdata.traj, "d1", 1.0);
        sf_putstring (tdata.traj, "label1", "Trajectory components");
        sf_putstring (tdata.traj, "unit1", "");
        sf_putint (tdata.traj, "n2", tdata.nt);
        sf_putfloat (tdata.traj, "o2", 0.0);
        sf_putfloat (tdata.traj, "d2", tdata.dt);
        sf_putstring (tdata.traj, "label2", "Time");
        sf_putstring (tdata.traj, "unit2", "s");
        sf_putint (tdata.traj, "n3", na);
        sf_putfloat (tdata.traj, "d3", da*180.0/SF_PI);
        sf_putfloat (tdata.traj, "o3", oa*180.0/SF_PI);
        sf_putstring (tdata.traj, "label3", "Angle");
        sf_putstring (tdata.traj, "unit3", "Degrees");
        sf_putint (tdata.traj, "n4", nz);
        sf_putfloat (tdata.traj, "o4", oz);
        sf_putfloat (tdata.traj, "d4", dz);
        if (!spdom) {
            sf_putstring (out, "label4", "Depth");
            sf_putstring (out, "unit4", "");
        }
        sf_putint (tdata.traj, "n5", nx);
        sf_putfloat (tdata.traj, "o5", ox);
        sf_putfloat (tdata.traj, "d5", dx);
        if (!spdom) {
            sf_putstring (out, "label5", "Lateral");
            sf_putstring (out, "unit5", "");
        }
    }

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness2_init (vspline, verb);

    if (tdata.traj)
        esc_tracer = sf_esc_tracer2_init (esc_slow,
                                          sf_escrt2_traj, tdata.dt, (void*)&tdata);
    else
        esc_tracer = sf_esc_tracer2_init (esc_slow,
                                          NULL, 0.0, NULL);
    sf_esc_tracer2_set_parab (esc_tracer, parab);

    esc_point = sf_esc_point2_init ();

    for (ix = 0; ix < nx; ix++) {
        x = ox + ix*dx;
        if (verb)
            sf_warning ("Shooting from lateral location %d of %d at %g;", ix + 1, nx, x);
        for (iz = 0; iz < nz; iz++) {
            z = oz + iz*dz;
            for (ia = 0; ia < na; ia++) {
                sf_esc_tracer2_compute (esc_tracer, z, x, oa + ia*da,
                                        0.0, 0.0, esc_point);
                /* Copy escape values to the output buffer */
                for (i = 0; i < ESC2_NUM; i++)
                    e[ia][i] = sf_esc_point2_get_esc_var (esc_point, i);
                if (tdata.traj) {
                    /* Fill the rest of the trajectory with the last point */
                    for (it = tdata.it + 1; it < tdata.nt; it++) {
                        for (i = 0; i < TRAJ2_COMPS - 1; i++)
                            tdata.pnts[it][i] = tdata.pnts[tdata.it][i];
                    }
                    sf_floatwrite (tdata.pnts[0], (size_t)tdata.nt*
                                                  (size_t)(TRAJ2_COMPS - 1), tdata.traj);
                }
            } /* Loop over a */
            sf_floatwrite (e[0], (size_t)na*(size_t)ESC2_NUM, out);
        } /* Loop over z */
    } /* Loop over x */
    if (verb)
        sf_warning (".");

    sf_esc_point2_close (esc_point);
    sf_esc_tracer2_close (esc_tracer);
    sf_esc_slowness2_close (esc_slow);

    if (tdata.traj) {
        free (tdata.pnts[0]);
        free (tdata.pnts);
        sf_fileclose (tdata.traj);
    }
    free (e[0]);
    free (e);

    sf_fileclose (vspline);

    return 0;
}


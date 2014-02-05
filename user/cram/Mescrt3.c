/* Escape tables by ray tracing with escape equations in 3-D. */
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

#include "esc_tracer3.h"

#define TRAJ3_COMPS 6
#define TRAJ3_Z 0
#define TRAJ3_X 1
#define TRAJ3_Y 2
#define TRAJ3_B 3
#define TRAJ3_A 4
#define TRAJ3_T 5

typedef struct {
    int       it, nt;
    float     dt;
    float   **pnts;
} sf_escrt3_traj_cbud;

/* Callback to output points along a ray trajectory; this one is
   called by sf_esc_tracer3 at every step during ray tracing in phase space;
   each next call to this routine can be more separated by more than one
   dt step in time with the previous one; we use interpolation in here
   to fill in the intermediate trajectory points */
void sf_escrt3_traj (float z, float x, float y, float b,
                     float a, int it, void *ud) {
    int i, j, n;
    float z0, x0, y0, b0, a0, dz, dx, dy, db, da, b1, a1;
    sf_escrt3_traj_cbud *cbud = (sf_escrt3_traj_cbud*)ud;

    if (it < cbud->nt) {
        n = it - cbud->it;
        if (it != 0 && n > 1) {
            /* Interpolate values inbetween the last one
               and the new one */
            z0 = cbud->pnts[cbud->it][TRAJ3_Z];
            x0 = cbud->pnts[cbud->it][TRAJ3_X];
            y0 = cbud->pnts[cbud->it][TRAJ3_Y];
            b0 = cbud->pnts[cbud->it][TRAJ3_B];
            a0 = cbud->pnts[cbud->it][TRAJ3_A];

            /* Make sure that angles do not jump
               across 2.0*pi->0 boundary in azimuth */
            if (a0 > 3.0*SF_PI/2.0 && a < SF_PI/2.0)
                a0 -= 2.0*SF_PI;
            else if (a0 < SF_PI/2.0 && a > 3.0*SF_PI/2.0)
                a0 += 2.0*SF_PI;
            else if ((a0 - a) > SF_PI) {
                a0 -= SF_PI;
                b0 = b < SF_PI/2.0
                     ? 2.0*SF_PI - b0 : -b0;
            } else if ((a - a0) > SF_PI) {
                a0 += SF_PI;
                b0 = b < SF_PI/2.0
                     ? 2.0*SF_PI - b0 : -b0;
            }
            /* Steps */
            dz = (z - z0)/(float)n;
            dx = (x - x0)/(float)n;
            dy = (y - y0)/(float)n;
            db = (b - b0)/(float)n;
            da = (a - a0)/(float)n;
            /* Fill the trajectory in */
            for (i = cbud->it + 1; i < it; i++) {
                j = i - cbud->it;
                a1 = a0 + j*da;
                b1 = b0 + j*db;
                /* Keep b in [0; pi] and a in [0; 2*pi] range */
                if (a1 < 0.0)
                    a1 += 2.0*SF_PI;
                else if (a > 2.0*SF_PI)
                    a1 -= 2.0*SF_PI;
                if (b1 < 0.0) {
                    a1 += SF_PI;
                    b1 = fabs (b);
                } else if (b1 > SF_PI) {
                    a1 -= SF_PI;
                    b1 = 2.0*SF_PI - b;
                }
                if (a1 < 0.0)
                    a1 += 2.0*SF_PI;
                else if (a > 2.0*SF_PI)
                    a1 -= 2.0*SF_PI;
                cbud->pnts[i][TRAJ3_B] = b1;
                cbud->pnts[i][TRAJ3_A] = a1;
                cbud->pnts[i][TRAJ3_Z] = z0 + j*dz;
                cbud->pnts[i][TRAJ3_X] = x0 + j*dx;
                cbud->pnts[i][TRAJ3_Y] = y0 + j*dy;
            }
        }
        cbud->pnts[it][TRAJ3_Z] = z;
        cbud->pnts[it][TRAJ3_X] = x;
        cbud->pnts[it][TRAJ3_Y] = y;
        cbud->pnts[it][TRAJ3_B] = b;
        cbud->pnts[it][TRAJ3_A] = a;
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
    int nz, nx, ny, nb, na, ib, ia, iz, ix, iy, i, it, nt, ic, nc = 1, fz, lz, itr = 0;
    float dz, oz, dx, ox, dy, oy, db, ob, da, oa, z, x, y, a, dt, df, md, aper;
    float ****e;
    sf_file spdom, vspline = NULL, out, traj = NULL;
    sf_escrt3_traj_cbud *tdata = NULL; 
    char *ext = NULL;
    bool verb, parab;
    sf_esc_slowness3 esc_slow;
    sf_esc_tracer3 *esc_tracers;
    sf_esc_point3 *esc_points;
    sf_timer timer;

    sf_init (argc, argv);

    if (!sf_stdin ()) {
        spdom = NULL;
    } else {
        spdom = sf_input ("in");
        /* Spatial (z,x,y) domain */
    }

    out = sf_output ("out");
    /* Escape values */

    /* Spatial dimensions */
    if (spdom) {
        if (!sf_histint (spdom, "n1", &nz)) sf_error ("No n1= in input");
        if (!sf_histint (spdom, "n2", &nx)) sf_error ("No n2= in input");
        if (!sf_histint (spdom, "n3", &ny)) sf_error ("No n3= in input");
        if (!sf_histfloat (spdom, "d1", &dz)) sf_error ("No d1= in input");
        if (!sf_histfloat (spdom, "o1", &oz)) sf_error ("No o1= in input");
        if (!sf_histfloat (spdom, "d2", &dx)) sf_error ("No d2= in input");
        if (!sf_histfloat (spdom, "o2", &ox)) sf_error ("No o2= in input");
        if (!sf_histfloat (spdom, "d3", &dy)) sf_error ("No d3= in input");
        if (!sf_histfloat (spdom, "o3", &oy)) sf_error ("No o3= in input");
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
    if (!sf_getint ("ny", &ny) && !spdom) sf_error ("Need ny=");
    /* Number of samples in y axis */
    if (!sf_getfloat ("oy", &oy) && !spdom) sf_error ("Need oy=");
    /* Beginning of y axis */
    if (!sf_getfloat ("dy", &dy) && !spdom) sf_error ("Need dy=");
    /* Sampling of y axis */

    if (!sf_getint ("na", &na)) na = 360;
    /* Number of azimuth phase angles */
    da = 2.0*SF_PI/(float)na;
    oa = 0.5*da;

    if (!sf_getint ("nb", &nb)) nb = 180;
    /* Number of inclination phase angles */
    db = SF_PI/(float)nb;
    ob = 0.5*db;

    if (!sf_getfloat ("df", &df)) df = 0.1;
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
        tdata = (sf_escrt3_traj_cbud*)sf_alloc (nc*na*nb, sizeof(sf_escrt3_traj_cbud));
        for (itr = 0; itr < nc*na*nb; itr++) {
            tdata[itr].it = 0;
            tdata[itr].nt = nt;
            tdata[itr].dt = dt;
            tdata[itr].pnts = sf_floatalloc2 (TRAJ3_COMPS - 1, nt);
        }
    }

    e = sf_floatalloc4 (ESC3_NUM, nb, na, nc);

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness3_init (vspline, verb);

    /* Make room for escape variables in output */
    if (spdom)
        sf_shiftdimn (spdom, out, 1, 3);

    sf_putint (out, "n1", ESC3_NUM);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Escape variable");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", nb);
    sf_putfloat (out, "d2", db*180.0/SF_PI);
    sf_putfloat (out, "o2", ob*180.0/SF_PI);
    sf_putstring (out, "label2", "Inclination");
    sf_putstring (out, "unit2", "Degrees");
    sf_putint (out, "n3", na);
    sf_putfloat (out, "d3", da*180.0/SF_PI);
    sf_putfloat (out, "o3", oa*180.0/SF_PI);
    sf_putstring (out, "label3", "Azimuth");
    sf_putstring (out, "unit3", "Degrees");

    sf_putint (out, "n4", nz);
    sf_putfloat (out, "o4", oz);
    sf_putfloat (out, "d4", dz);
    if (!spdom) {
        sf_putstring (out, "label4", "Depth");
        sf_putstring (out, "unit4", "");
    }
    sf_putint (out, "n5", nx);
    sf_putfloat (out, "o5", ox);
    sf_putfloat (out, "d5", dx);
    if (!spdom) {
        sf_putstring (out, "label5", "X");
        sf_putstring (out, "unit5", "");
    }
    sf_putint (out, "n6", ny);
    sf_putfloat (out, "o6", oy);
    sf_putfloat (out, "d6", dy);
    if (!spdom) {
        sf_putstring (out, "label6", "Y");
        sf_putstring (out, "unit6", "");
    }
    /* Save min/max possible escape values */
    sf_putfloat (out, "Zmin", sf_esc_slowness3_oz (esc_slow));
    sf_putfloat (out, "Zmax", sf_esc_slowness3_oz (esc_slow) +
                              (sf_esc_slowness3_nz (esc_slow) - 1)*
                              sf_esc_slowness3_dz (esc_slow));
    sf_putfloat (out, "Xmin", sf_esc_slowness3_ox (esc_slow));
    sf_putfloat (out, "Xmax", sf_esc_slowness3_ox (esc_slow) +
                              (sf_esc_slowness3_nx (esc_slow) - 1)*
                              sf_esc_slowness3_dx (esc_slow));
    sf_putfloat (out, "Ymin", sf_esc_slowness3_oy (esc_slow));
    sf_putfloat (out, "Ymax", sf_esc_slowness3_oy (esc_slow) +
                              (sf_esc_slowness3_ny (esc_slow) - 1)*
                              sf_esc_slowness3_dy (esc_slow));

    if (traj) {
        if (spdom)
            sf_shiftdimn (spdom, traj, 1, 4);
        sf_putint (traj, "n1", TRAJ3_COMPS - 1);
        sf_putfloat (traj, "o1", 0.0);
        sf_putfloat (traj, "d1", 1.0);
        sf_putstring (traj, "label1", "Escape variable");
        sf_putstring (traj, "unit1", "");
        sf_putint (traj, "n2", nt);
        sf_putfloat (traj, "o2", 0.0);
        sf_putfloat (traj, "d2", dt);
        sf_putstring (traj, "label2", "Time");
        sf_putstring (traj, "unit2", "s");
        sf_putint (traj, "n3", nb);
        sf_putfloat (traj, "d3", db*180.0/SF_PI);
        sf_putfloat (traj, "o3", ob*180.0/SF_PI);
        sf_putstring (traj, "label3", "Inclination");
        sf_putstring (traj, "unit3", "Degrees");
        sf_putint (traj, "n4", na);
        sf_putfloat (traj, "d4", da*180.0/SF_PI);
        sf_putfloat (traj, "o4", oa*180.0/SF_PI);
        sf_putstring (traj, "label4", "Azimuth");
        sf_putstring (traj, "unit4", "Degrees");
        
        sf_putint (traj, "n5", nz);
        sf_putfloat (traj, "o5", oz);
        sf_putfloat (traj, "d5", dz);
        if (!spdom) {
            sf_putstring (traj, "label5", "Depth");
            sf_putstring (traj, "unit5", "");
        }
        sf_putint (traj, "n6", nx);
        sf_putfloat (traj, "o6", ox);
        sf_putfloat (traj, "d6", dx);
        if (!spdom) {
            sf_putstring (traj, "label6", "X");
            sf_putstring (traj, "unit6", "");
        }
        sf_putint (traj, "n7", ny);
        sf_putfloat (traj, "o7", oy);
        sf_putfloat (traj, "d7", dy);
        if (!spdom) {
            sf_putstring (traj, "label7", "Y");
            sf_putstring (traj, "unit7", "");
        }
    }

    esc_tracers = (sf_esc_tracer3*)sf_alloc (nc, sizeof(sf_esc_tracer3));
    esc_points = (sf_esc_point3*)sf_alloc (nc, sizeof(sf_esc_point3));
    for (ic = 0; ic < nc; ic++) {
        esc_tracers[ic] = sf_esc_tracer3_init (esc_slow);
        sf_esc_tracer3_set_parab (esc_tracers[ic], parab);
        if (md != SF_HUGE)
            sf_esc_tracer3_set_mdist (esc_tracers[ic], md);
        sf_esc_tracer3_set_df (esc_tracers[ic], df);
        esc_points[ic] = sf_esc_point3_init ();
    }

    timer = sf_timer_init ();
    /* Ray tracing loop */
    for (iy = 0; iy < ny; iy++) {
        y = oy + iy*dy;
        /* Set aperture */
        if (aper != SF_HUGE) {
            for (ic = 0; ic < nc; ic++) {
                sf_esc_tracer3_set_ymin (esc_tracers[ic], y - aper);
                sf_esc_tracer3_set_ymax (esc_tracers[ic], y + aper);
            }
        }
        for (ix = 0; ix < nx; ix++) {
            x = ox + ix*dx;
            /* Set aperture */
            if (aper != SF_HUGE) {
                for (ic = 0; ic < nc; ic++) {
                    sf_esc_tracer3_set_xmin (esc_tracers[ic], x - aper);
                    sf_esc_tracer3_set_xmax (esc_tracers[ic], x + aper);
                }
            }
            if (verb)
                sf_warning ("%s Shooting from lateral location %d of %d at y=%g, x=%g;",
                            ext, iy*nx + ix + 1, ny*nx, y, x);
            /* Loop over chunks */
            for (ic = 0; ic < (nz/nc + ((nz % nc) != 0)); ic++) {
                fz = ic*nc;
                lz = (ic + 1)*nc - 1;
                if (lz >= nz)
                    lz = nz - 1;
                sf_timer_start (timer);
#ifdef _OPENMP
#pragma omp parallel for                       \
                schedule(static,1)             \
                private(iz,ia,ib,a,z,it,i,itr) \
                shared(fz,lz,iy,ix,nb,na,nz,nx,ny,ob,oa,oz,ox,oy,db,da,dz,dx,dy,x,y,tdata,esc_tracers,esc_points,e,out,traj)
#endif
                for (iz = fz; iz <= lz; iz++) {
                    z = oz + iz*dz;
                    for (ia = 0; ia < na; ia++) {
                        a = oa + ia*da;
                        for (ib = 0; ib < nb; ib++) {
                            if (traj) {
                                itr = (iz - fz)*na*nb + ia*nb + ib;
                                sf_esc_tracer3_set_trajcb (esc_tracers[iz - fz], sf_escrt3_traj, dt,
                                                           (void*)&tdata[itr]);
                            }
                            sf_esc_tracer3_compute (esc_tracers[iz - fz], z, x, y, ob + ib*db, a,
                                                    0.0, 0.0, esc_points[iz - fz], NULL, NULL);
                            /* Copy escape values to the output buffer */
                            for (i = 0; i < ESC3_NUM; i++)
                                e[iz - fz][ia][ib][i] = sf_esc_point3_get_esc_var (esc_points[iz - fz], i);
                            if (traj) {
                                /* Fill the rest of the trajectory with the last point */
                                for (it = tdata[itr].it + 1; it < tdata[itr].nt; it++) {
                                    for (i = 0; i < TRAJ3_COMPS - 1; i++)
                                        tdata[itr].pnts[it][i] = tdata[itr].pnts[tdata[itr].it][i];
                                }
                            }
                        } /* Loop over b */
                    } /* Loop over a */
                } /* Loop over z */
                sf_timer_stop (timer);
                sf_floatwrite (e[0][0][0], (size_t)(lz - fz + 1)*(size_t)nb*(size_t)na*(size_t)ESC3_NUM,
                               out);
                if (tdata) {
                    for (itr = 0; itr < (lz - fz + 1)*na*nb; itr++) {
                        sf_floatwrite (tdata[itr].pnts[0],
                                       (size_t)tdata[itr].nt*(size_t)(TRAJ3_COMPS - 1), traj);
                    }
                }
            } /* Loop over z chunks */
        } /* Loop over x */
    } /* Loop over y */
    if (verb) {
        sf_warning (".");
        sf_warning ("%s Total kernel time: %g s, per depth point: %g s",
                    ext, sf_timer_get_total_time (timer)/1000.0,
                    (sf_timer_get_total_time (timer)/(float)((size_t)nx*(size_t)ny*(size_t)nz))/1000.0);
    }
    sf_timer_close (timer);

    for (ic = 0; ic < nc; ic++) {
        sf_esc_point3_close (esc_points[ic]);
        sf_esc_tracer3_close (esc_tracers[ic]);
    }
    free (esc_points);
    free (esc_tracers);
    if (traj) {
        for (itr = 0; itr < nc*na*nb; itr++) {
            free (tdata[itr].pnts[0]);
            free (tdata[itr].pnts);
        }
        free (tdata);
    }
    sf_esc_slowness3_close (esc_slow);

    free (e[0][0][0]);
    free (e[0][0]);
    free (e[0]);
    free (e);
    free (ext);

    sf_fileclose (vspline);
    if (traj)
        sf_fileclose (traj);

    return 0;
}


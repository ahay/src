/* Escape tables by stitching of escape solutions in supercells in 3-D. */
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

#include "esc_scgrid3.h"

int main (int argc, char* argv[]) {
    int nz, nx, ny, nb, na, iz, ix, iy, ia, ib;
    int icpu = 0, ncpu = 1, morder = 2, nc = 1;
    float dz, oz, dx, ox, dy, oy, db, ob, da, oa;
    float z, x, y;
    float ***e;
    sf_file spdom, vspline = NULL, scgrid = NULL, scdaemon = NULL, 
            out;

    bool verb, parab;
    sf_esc_slowness3 esc_slow;
    sf_esc_tracer3 esc_tracer;
    sf_esc_scgrid3 esc_scgrid;

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
        if (!sf_histint (spdom, "icpu", &icpu)) icpu = 0;
        /* Current CPU number */
        if (!sf_histint (spdom, "ncpu", &ncpu)) ncpu = 1;
        /* Total number of CPUs */
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

#ifdef _OPENMP
    if (!sf_getint ("nc", &nc)) nc = 1;
    /* Number of threads to use for interpolation */
    omp_set_num_threads (nc);
    omp_set_num_threads (nc);
    sf_warning ("Using %d threads, omp_get_max_threads()=%d",
                nc, omp_get_max_threads ());
#endif

    if (!sf_getbool ("parab", &parab)) parab = true;
    /* y - use parabolic approximation of trajectories, n - straight line */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    e = sf_floatalloc3 (ESC3_NUM, nb, na);

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    if (!sf_getstring ("scgrid")) sf_error ("Need scgrid=");
    /* Grid of supercells of local escape solutions */
    scgrid = sf_input ("scgrid");

    if (sf_getstring ("scdaemon")) {
        /* Daemon for distributed computation */
        scdaemon = sf_input ("scdaemon");
    }

    if (!sf_getint ("morder", &morder)) morder = 1;
    /* Order of interpolation accuracy in the angular domain (1-3) */
   
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

    esc_tracer = sf_esc_tracer3_init (esc_slow, NULL, 0.0, NULL);
    sf_esc_tracer3_set_parab (esc_tracer, parab);

    esc_scgrid = sf_esc_scgrid3_init (scgrid, scdaemon, esc_tracer, morder,
                                      (float)icpu/(float)ncpu, verb);

    if (scdaemon)
        sf_fileclose (scdaemon);

    for (iy = 0; iy < ny; iy++) {
        y = oy + iy*dy;
        for (ix = 0; ix < nx; ix++) {
            x = ox + ix*dx;
            if (verb)
                sf_warning ("Projecting from lateral location %d of %d at y=%g, x=%g;",
                            iy*nx + ix + 1, ny*nx, y, x);
            for (iz = 0; iz < nz; iz++) {
                z = oz + iz*dz;
                if (sf_esc_tracer3_inside (esc_tracer, &z, &x, &y, false)) {
                    sf_esc_scgrid3_compute (esc_scgrid, z, x, y, oa, da, ob, db, na, nb, e[0][0]);
                } else {
                    for (ia = 0; ia < na; ia++) {
                        for (ib = 0; ib < nb; ib++) {
                            e[ia][ib][ESC3_Z] = z;
                            e[ia][ib][ESC3_X] = x;
                            e[ia][ib][ESC3_Y] = y;
                            e[ia][ib][ESC3_T] = 0.0;
#ifdef ESC_EQ_WITH_L
                            e[ia][ib][ESC3_L] = 0.0;                            
#endif
                        }
                    }
                }
                sf_floatwrite (e[0][0], (size_t)nb*(size_t)na*(size_t)ESC3_NUM,
                               out);
            } /* Loop over z */
        } /* Loop over x */
    } /* Loop over y */
    if (verb)
        sf_warning (".");

    sf_esc_scgrid3_close (esc_scgrid, verb);
    sf_esc_tracer3_close (esc_tracer);
    sf_esc_slowness3_close (esc_slow);

    free (e[0][0]);
    free (e[0]);
    free (e);

    sf_fileclose (scgrid);
    sf_fileclose (vspline);

    return 0;
}


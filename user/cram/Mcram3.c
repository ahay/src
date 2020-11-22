/* 3-D angle-domain Kirchhoff migration based on escape tables. */
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

#include "cram_gather3.h"
#include "cram_point3.h"
#include "esc_point3.h"

static char* sf_cram3_warnext (sf_file input) {
    int icpu, ncpu, maxlen;
    char *host = sf_gethost ();
    char *ext = NULL;

    if (!sf_histint (input, "icpu", &icpu)) icpu = 0;
    if (!sf_histint (input, "ncpu", &ncpu)) ncpu = 1;
    maxlen = strlen (host) + 30;
    ext = (char*)malloc (maxlen*sizeof(char));
    snprintf (ext, maxlen, "[%s:%d/%d]", host, icpu + 1, ncpu);
    return ext;
}

int main (int argc, char* argv[]) {
    size_t i, j, is, ih, nh;
    float sx, sy, gx, gy;
    float gxmin, gxmax, gymin, gymax;
    int iz, ix, iy, nz, nx, ny, nb, na, nt, np;
    float dt, db, da, t0, b0, a0, dbx, dby, dxm, dym;
    float dz, z0, dx, x0, dy, y0, zd, z, x, y, zf, vconst = 1.5;
    float oazmin = 180.0, oazmax = 180.0, dazmin = 180.0, dazmax = 180.0;
    float xbmin, xbmax, ybmin, ybmax, zbmin, zbmax;
    float oaz, daz, armin, armax;
    float ***esc, *s;
    sf_file esct, data = NULL, survey = NULL, vz = NULL, ddaemon = NULL;
    sf_file imag = NULL, hits = NULL, oimag = NULL, dimag = NULL,
            osmap = NULL, dsmap = NULL, oimap = NULL, dimap = NULL;
    char *ext = NULL;
    bool amp, mute, outaz, extrap, inorm;
    sf_cram_data2 cram_data;
    sf_cram_survey3 cram_survey;
    sf_cram_slowness3 cram_slowness;
    sf_cram_gather3 cram_gather;
    sf_cram_rbranch3 *cram_rbranches;
    sf_cram_point3 *cram_points;

    sf_init (argc, argv);

    esct = sf_input ("in");
    /* Escape tables */

    /* Phase space dimensions (also, imaging dimensions) */
    if (!sf_histint (esct, "n1", &na)) sf_error ("No n1= in input");
    if (na != ESC3_NUM) sf_error ("Need n1=%d in input", ESC3_NUM);
    if (!sf_histint (esct, "n2", &nb)) sf_error ("No n2= in input");
    if (!sf_histint (esct, "n3", &na)) sf_error ("No n3= in input");
    if (!sf_histint (esct, "n4", &nz)) sf_error ("No n4= in input");
    if (!sf_histint (esct, "n5", &nx)) sf_error ("No n5= in input");
    if (!sf_histint (esct, "n6", &ny)) sf_error ("No n6= in input");
    if (!sf_histfloat (esct, "d2", &db)) sf_error ("No d2= in input");
    if (!sf_histfloat (esct, "o2", &b0)) sf_error ("No o2= in input");
    if (!sf_histfloat (esct, "d3", &da)) sf_error ("No d3= in input");
    if (!sf_histfloat (esct, "o3", &a0)) sf_error ("No o3= in input");
    if (!sf_histfloat (esct, "d4", &dz)) sf_error ("No d4= in input");
    if (!sf_histfloat (esct, "o4", &z0)) sf_error ("No o4= in input");
    if (!sf_histfloat (esct, "d5", &dx)) sf_error ("No d5= in input");
    if (!sf_histfloat (esct, "o5", &x0)) sf_error ("No o5= in input");
    if (!sf_histfloat (esct, "d6", &dy)) sf_error ("No d6= in input");
    if (!sf_histfloat (esct, "o6", &y0)) sf_error ("No o6= in input");

    if (!sf_histfloat (esct, "Xmin", &xbmin)) sf_error ("No Xmin= in input");
    if (!sf_histfloat (esct, "Xmax", &xbmax)) sf_error ("No Xmax= in input");
    if (!sf_histfloat (esct, "Ymin", &ybmin)) sf_error ("No Ymin= in input");
    if (!sf_histfloat (esct, "Ymax", &ybmax)) sf_error ("No Ymax= in input");
    if (!sf_histfloat (esct, "Zmin", &zbmin)) sf_error ("No Zmin= in input");
    if (!sf_histfloat (esct, "Zmax", &zbmax)) sf_error ("No Zmax= in input");

    ext = sf_cram3_warnext (esct);

    /* Surface depth */
    zd = zbmin + 0.25*dz;

    if (!sf_getbool ("amp", &amp)) amp = true;
    /* n - do not apply amplitude correction weights */
    if (!sf_getbool ("extrap", &extrap)) extrap = false;
    /* y - extrapolate migrated samples in gathers */
    if (!sf_getbool ("mute", &mute)) mute = false;
    /* y - mute signal in constant z plane before stacking */
    if (!sf_getbool ("outaz", &outaz)) outaz = true;
    /* n - stack azimuth direction before output */
    if (!sf_getbool ("inorm", &inorm)) inorm = false;
    /* y - normalize gathers for illumination */

    if (mute) {
        if (!sf_getfloat ("oazmin", &oazmin)) oazmin = 180.0;
        /* Maximum allowed scattering angle at z min */
        if (!sf_getfloat ("oazmax", &oazmax)) oazmax = 180.0;
        /* Maximum allowed scattering angle at z max */
        if (!sf_getfloat ("dazmin", &dazmin)) dazmin = 180.0;
        /* Maximum allowed dip angle at z min */
        if (!sf_getfloat ("dazmax", &dazmax)) dazmax = 180.0;
        /* Maximum allowed dip angle at z max */
        if (oazmin < 0.0) oazmin = 0.0;
        if (oazmax < 0.0) oazmax = 0.0;
        if (dazmin < 0.0) dazmin = 0.0;
        if (dazmax < 0.0) dazmax = 0.0;
    }

    if (!sf_getfloat ("dbx", &dbx)) dbx = 10.0*dx;
    /* Size of search bins in x */
    if (!sf_getfloat ("dby", &dby)) dby = 10.0*dy;
    /* Size of search bins in y */

    if (!sf_getfloat ("dxm", &dxm)) dxm = 5.0*dx;
    /* Taper length in x */
    if (!sf_getfloat ("dym", &dym)) dym = 5.0*dy;
    /* Taper length in y */

    if (!sf_getfloat ("armin", &armin)) armin = 0.01*dy*dx;
    /* Minimum allowed area for an exit ray branch */
    if (!sf_getfloat ("armax", &armax)) armax = 100.0*dy*dx;
    /* Maximum allowed area for an exit ray branch */

    if (!sf_getint ("np", &np)) np = 1;
    /* number of image points to buffer before accessing data */
    if (np > nz*nx*ny)
        np = nz*nx*ny;

    esc = sf_floatalloc3 (ESC3_NUM, nb, na);

    if (sf_getstring ("data")) {
        /* Processed prestack data */
        data = sf_input ("data");
    } else {
        sf_error ("Need data=");
    }

    if (sf_getstring ("ddaemon")) {
        /* Daemon for distributed data storage */
        ddaemon = sf_input ("ddaemon");
    }

    if (sf_getstring ("survey")) {
        /* Survey info for input data */
        survey = sf_input ("survey");
    } else {
        sf_error ("Need survey=");
    }

    /* Data dimensions */
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in data");
    if (!sf_histfloat (data, "d1", &dt)) sf_error ("No d1= in data");
    if (!sf_histfloat (data, "o1", &t0)) sf_error ("No o1= in data");

    if (sf_getstring ("vz")) {
        /* Velocity model for amplitude weights */
        vz = sf_input ("vz");
    } else {
        if (!sf_getfloat ("vconst", &vconst)) vconst = 1.5;
        /* Constant velocity, if vz= is not used */
    }

#ifdef _OPENMP
    int nc = 1;
    if (!sf_getint("nc", &nc))
        nc = 0;
    /* Number of threads to use for ray tracing (OMP_NUM_THREADS by default) */
    if (nc)
        omp_set_num_threads (nc); /* User override */
    else
        nc = omp_get_max_threads (); /* Current default */
    sf_warning ("%s Using %d threads", ext, omp_get_max_threads ());
#endif

    /* Data object */
    cram_data = sf_cram_data2_init (data, ddaemon);
    if (ddaemon)
        sf_fileclose (ddaemon);
    /* Survey object */
    cram_survey = sf_cram_survey3_init (survey, true);
    sf_fileclose (survey);
    /* Slowness object */
    cram_slowness = sf_cram_slowness3_init (vz, vconst);

    /* Image and gathers accumulator object */
    cram_gather = sf_cram_gather3_init (nb, na, nz, b0, db, a0, da,
                                        oazmax, dazmax, outaz, inorm);

    imag = sf_output ("out");
    /* Image (z, x, y) */
    sf_cram_gather3_setup_image_output (cram_gather, esct, imag);

    if (sf_getstring ("hits")) {
        /* Image illumination (z, x, y) */
        hits = sf_output ("hits");
        sf_cram_gather3_setup_image_output (cram_gather, esct, hits);
    }

    if (sf_getstring ("agath")) {
        /* Scattering angle gathers (angle, azimuth, z, x, y) */
        oimag = sf_output ("agath");
        sf_cram_gather3_setup_oangle_output (cram_gather, esct, oimag);
        if (sf_getstring ("imap")) {
            /* SCattering gathers illumination (angle, azimuth, z, x, y) */
            oimap = sf_output ("imap");
            sf_cram_gather3_setup_oangle_output (cram_gather, esct, oimap);
        }
        if (sf_getstring ("smap")) {
            /* Scattering gathers energy (angle, azimuth, z, x, y) */
            osmap = sf_output ("smap");
            sf_cram_gather3_setup_oangle_output (cram_gather, esct, osmap);
        }
    }

    if (sf_getstring ("dipagath")) {
        /* Dip angle gathers (angle, azimuth, z, x, y) */
        dimag = sf_output ("dipagath");
        sf_cram_gather3_setup_dangle_output (cram_gather, esct, dimag);
        if (sf_getstring ("dipimap")) {
            /* Dip gathers illumination (angle, azimuth, z, x, y) */
            dimap = sf_output ("dipimap");
            sf_cram_gather3_setup_dangle_output (cram_gather, esct, dimap);
        }
        if (sf_getstring ("dipsmap")) {
            /* Dip gathers energy (angle, azimuth, z, x, y) */
            dsmap = sf_output ("dipsmap");
            sf_cram_gather3_setup_dangle_output (cram_gather, esct, dsmap);
        }
    }

    cram_points = (sf_cram_point3*)sf_alloc (np, sizeof (sf_cram_point3));
    cram_rbranches = (sf_cram_rbranch3*)sf_alloc (np, sizeof (sf_cram_rbranch3));
    s = sf_floatalloc (np);

    for (j = 0; j < np; j++) {
        /* Exit ray branches object */
        cram_rbranches[j] = sf_cram_rbranch3_init (nb, na, zd, t0 + (nt - 1)*dt,
                                                   xbmin, xbmax, ybmin, ybmax,
                                                   dbx, dby, cram_slowness);
        sf_cram_rbranch3_set_arminmax (cram_rbranches[j], armin, armax);
        /* Subsurface point image object */
        cram_points[j] = sf_cram_point3_init (nb, b0, db, na, a0, da,
                                              oimag != NULL, dimag != NULL,
                                              cram_data, cram_survey, cram_slowness, cram_rbranches[j]);
        sf_cram_point3_set_amp (cram_points[j], amp);
        sf_cram_point3_set_taper (cram_points[j], dxm, dym);
        sf_cram_point3_set_extrap (cram_points[j], extrap);
    }

    if (np > 1)
        sf_warning ("%s Buffering escape tables for %d points", ext, np);
    j = 0;
    for (iy = 0; iy < ny; iy++) { /* Loop over image y */
        y = y0 + iy*dy;
        for (ix = 0; ix < nx; ix++) { /* Loop over image x */
            x = x0 + ix*dx;
            sf_warning ("%s Lateral %d of %d (x=%g, y=%g)", ext, iy*nx + ix + 1, nx*ny, x, y);
            for (iz = 0; iz < nz; iz++) { /* Loop over image z */
                z = z0 + iz*dz;
                /* Escape variables for this (z, x, y) location */
                sf_floatread (esc[0][0], ESC3_NUM*na*nb, esct);
                if (mute) {
                    zf = (z - zbmin)/(zbmax - zbmin);
                    /* Maximum scattering angle for this z */
                    oaz = oazmax - zf*(oazmax - oazmin);
                    /* Maximum dip angle for this z */
                    daz = dazmax - zf*(dazmax - dazmin);
                    sf_cram_point3_set_mute (cram_points[j], oaz, daz);
                }
                if (z > zd) {
                    s[j] = sf_cram_slowness3_get_value (cram_slowness, z, x, y);
                    sf_cram_rbranch3_set_escapes (cram_rbranches[j], esc);
                } else
                    s[j] = -1.0;
                j++;
                if (np == j || (iy == (ny - 1) && ix == (nx - 1) && iz == (nz - 1))) {
                    /* Buffer is full or last set of escape tables is here - access data now */
                    if (np != 1)
                        sf_warning ("%s Migrating data, z=%g, x=%g, y=%g", ext, z, x, y);
                    /* Loop over known sources */
                    is = sf_cram_survey3_get_first_source (cram_survey, &sx, &sy, &nh,
                                                           &gxmin, &gxmax, &gymin, &gymax);
                    while (is != (size_t)-1) {
                        /* Loop over known receivers */
                        ih = sf_cram_survey3_get_first_receiver (cram_survey, is,
                                                                 &i, &gx, &gy);
                        while (ih != (size_t)-1) {
#ifdef _OPENMP
#pragma omp parallel for                        \
                            schedule(dynamic,1) \
                            private(j)          \
                            shared(np,sx,sy,gx,gy,gxmin,gxmax,gymin,gymax,nh,i,s,cram_points)
#endif
                            for (j = 0; j < np; j++) {
                                if (s[j] > 0.0)
                                    sf_cram_point3_compute_one_trace (cram_points[j], sx, sy,
                                                                      gx, gy, gxmin, gxmax,
                                                                      gymin, gymax, s[j],
                                                                      i, nh);
                            }
                            ih = sf_cram_survey3_get_next_receiver (cram_survey, is, ih,
                                                                    &i, &gx, &gy);
                        } /* Loop over known receivers */
                        is = sf_cram_survey3_get_next_source (cram_survey, is, &sx, &sy, &nh,
                                                              &gxmin, &gxmax, &gymin, &gymax);
                    } /* Loop over known sources */
                    /* Write image */
                    if (iy == (ny - 1) && ix == (nx - 1) && iz == (nz - 1))
                        np = j;
                    for (j = 0; j < np; j++) {
                        sf_cram_gather3_image_output (cram_gather, cram_points[j], imag, hits);
                        /* Add to the angle gathers */
                        if (oimag)
                            sf_cram_gather3_oangle_output (cram_gather, cram_points[j],
                                                           oimag, osmap, oimap);
                        /* Dip gathers */
                        if (dimag)
                            sf_cram_gather3_dangle_output (cram_gather, cram_points[j],
                                                           dimag, dsmap, dimap);
                        sf_cram_point3_reset (cram_points[j]);
                        s[j] = -1.0;
                    }
                    j = 0;
                } /* np == j */
            } /* iz */
        } /* ix */
    } /* iy */
    /* In trace-by-trace mode, the actual migration takes place now,
       once all the escape functions are known */

    for (i = 0; i < np; i++) {
       sf_cram_point3_close (cram_points[i]);
       sf_cram_rbranch3_close (cram_rbranches[i]);
    }
    free (cram_points);
    free (cram_rbranches);

    sf_cram_gather3_close (cram_gather);
    sf_cram_slowness3_close (cram_slowness);
    sf_cram_survey3_close (cram_survey);
    sf_cram_data2_close (cram_data);

    sf_fileclose (data);
    if (vz)
        sf_fileclose (vz);
    if (hits)
        sf_fileclose (hits);
    if (dimag)
        sf_fileclose (dimag);
    if (osmap)
        sf_fileclose (osmap);
    if (dsmap)
        sf_fileclose (dsmap);
    if (oimap)
        sf_fileclose (oimap);
    if (dimap)
        sf_fileclose (dimap);

    free (s);
    free (ext);

    free (esc[0][0]);
    free (esc[0]);
    free (esc);

    return 0;
}


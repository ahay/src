/* Prestack time migration (2-D/3-D) for irregular data. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "ktmig.h"

int main (int argc, char* argv[]) {
    /* Counters */
    int i = 0, j = 0, k = 0, l, n, m;

    /* Input data parameters */
    int nt, nx, ny = 1, nix = 1, nin = 1, osize, ntr, btr, dbtr;
    float ot, dt;

    /* Apperture parameters */
    int ix, iy, minix, miniy, maxix, maxiy;
    /* Apperture half-width in each direction */
    int apx, apy;

    /* Image(output) space parameters */
    int ont, onx, ony;
    float oot, oox, ooy;
    float odt, odx, ody;

    /* Antialias filter parameters */
    int maxtri;
    float trfact;

    /* Aperture corners */
    int el_cx1, el_cx2, el_cy1, el_cy2;
    float el_x, el_y;

    /* Input traces, output image, velocity */
    float *t, *img, *v;
    /* Aperture indices */
    int *ap;
    /* Coordinates: shot, receiver, and midpoint */
    sf_complex *sxy, *gxy, *cxy;

    /* Migration trace parameters */
    int oidx;
    float ox, oy, sx, sy, gx, gy;

    char *vrmsfile, *sxsyfile, *gxgyfile, *cxcyfile;
    sf_file data, image, vrms, sxsy, gxgy, cxcy;

    bool verb, time, aa, diff;

    sf_timer total_timer = NULL, aux_kernel_timer = NULL,
             main_kernel_timer = NULL;

    sf_init (argc, argv);

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* Verbosity flag */
    if (!sf_getbool ("time", &time)) time = false;
    /* Measure execution time */
    if (!sf_getbool ("aa", &aa)) aa = true;
    /* Antialiaing flag */
    if (!sf_getbool ("diff", &diff)) diff = true;
    /* Differentiation flag */

    data = sf_input ("in");
    image = sf_output ("out");

    if (SF_FLOAT != sf_gettype (data))
        sf_error ("Need float input");
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in input");
    if (!sf_histint (data, "n2", &nx)) sf_error ("No n2= in input");
    if (!sf_histint (data, "n3", &ny)) ny = 1;
    if (!sf_histint (data, "n4", &nin)) nin = 1;
    if (!sf_histint (data, "n5", &nix)) nix = 1;
    ntr = nx*ny*nin*nix;

    if (!sf_histfloat (data, "d1", &dt)) sf_error ("No d1= in input");
    if (!sf_histfloat (data, "o1", &ot)) ot = 0.;

    vrmsfile = sf_getstring ("vrms");
    /* File with RMS velocities */
    if (NULL == vrmsfile) sf_error ("Need vrms="); 
    vrms = sf_input ("vrms");
    if (SF_FLOAT != sf_gettype (vrms)) sf_error ("Need float vrms");

    if (!sf_histint (vrms, "n1", &ont)) sf_error ("No n1= in vrms");
    if (!sf_histint (vrms, "n2", &onx)) sf_error ("No n2= in vrms");
    if (!sf_histint (vrms, "n3", &ony)) ony = 1;
    osize = ont*onx*ony;

    if (!sf_histfloat (vrms, "d1", &odt)) sf_error ("No d1= in vrms");
    if (!sf_histfloat (vrms, "d2", &odx)) sf_error ("No d2= in vrms");
    if (!sf_histfloat (vrms, "d3", &ody)) ody = 1.0;

    if (!sf_histfloat (vrms, "o1", &oot)) oot = 0.;
    if (!sf_histfloat (vrms, "o2", &oox)) oox = 0.;
    if (!sf_histfloat (vrms, "o3", &ooy)) ooy = 0.;
    if (verb)
        sf_warning ("Image size: %d x %d x %d", ont, onx, ony);

    if (!sf_getint ("dbtr", &dbtr)) dbtr = 1000;
    /* Number of input traces to read at once */
    btr = dbtr;

    sxsyfile = sf_getstring ("sxsy");
    /* File with shot coordinates */
    if (NULL == sxsyfile) sf_error ("Need sxsy="); 
    sxsy = sf_input ("sxsy");
    if (SF_COMPLEX != sf_gettype (sxsy)) sf_error ("Need complex sxsy");
    if (!sf_histint (sxsy, "n2", &n)) sf_error ("No n2= in sxsy");
    if (n < ntr) sf_error ("Number of values in sxsy is less than number of input traces");

    gxgyfile = sf_getstring ("gxgy");
    /* File with receiver coordinates */
    if (NULL == gxgyfile) sf_error ("Need gxgy="); 
    gxgy = sf_input ("gxgy");
    if (SF_COMPLEX != sf_gettype (gxgy)) sf_error ("Need complex gxgy");
    if (!sf_histint (gxgy, "n2", &n)) sf_error ("No n2= in gxgy");
    if (n < ntr) sf_error ("Number of values in gxgy is less than number of input traces");

    cxcyfile = sf_getstring ("cxcy");
    /* File with midpoint coordinates */
    if (NULL == cxcyfile) sf_error ("Need cxcy="); 
    cxcy = sf_input ("cxcy");
    if (SF_COMPLEX != sf_gettype (cxcy)) sf_error ("Need complex cxcy");
    if (!sf_histint (cxcy, "n2", &n)) sf_error ("No n2= in cxcy");
    if (n < ntr) sf_error ("Number of values in cxcy is less than number of input traces");

    if (!sf_getint ("apx", &apx)) apx = onx/2;
    /* Apperture half-width in x direction */
    if (!sf_getint ("apy", &apy)) apy = ony/2;
    /* Apperture half-width in y direction */

    if (!sf_getint ("maxtri", &maxtri)) maxtri = 13;
    /* Maximum half-length of the antialias filter */
    if (!sf_getfloat ("trfact", &trfact)) trfact = 4.0*(0.5*(odx + ody)/dt);
    /* Trace factor for antialias filter length calculation */

    /* Initiate output */
    sf_putint (image, "n1", ont);
    sf_putint (image, "n2", onx);
    sf_putint (image, "n3", ony);
    sf_putint (image, "n4", 1);
    sf_putint (image, "n5", 1);
    sf_putfloat (image, "d1", odt);
    sf_putfloat (image, "d2", odx);
    sf_putfloat (image, "d3", ody);
    sf_putfloat (image, "d4", 0.0);
    sf_putfloat (image, "d5", 0.0);
    sf_putfloat (image, "o1", oot);
    sf_putfloat (image, "o2", oox);
    sf_putfloat (image, "o3", ooy);
    sf_putfloat (image, "o4", 0.0);
    sf_putfloat (image, "o5", 0.0);

    if (verb || time)
        total_timer = sf_timer_init ();
    if (verb) {
        main_kernel_timer = sf_timer_init ();
        aux_kernel_timer = sf_timer_init ();
    }

    if (verb || time)
        sf_timer_start (total_timer);

    v = sf_floatalloc (osize);
    img = sf_floatalloc (osize);
    memset (img, 0, osize*sizeof(float));

    sf_floatread (v, osize, vrms);

    t  = sf_floatalloc (btr*nt);
    sxy = sf_complexalloc (btr);
    gxy = sf_complexalloc (btr);
    cxy = sf_complexalloc (btr);

    /* Array of aperture indices */
    ap = sf_intalloc (onx*ony);

    if (verb)
        sf_warning ("Migrating traces in chunks of %d", btr);

    /* Loop over input traces */
    i = 0;
    while (i < ntr) {
        /* How many to read */
        k = ((i + btr) < ntr)
          ? btr
          : ntr - i;
        if (verb)
            sf_warning ("Processing traces %d-%d out of %d", i, i + k - 1, ntr);
        /* Read input data */
        sf_floatread (t, nt*k, data);
        sf_complexread (sxy, k, sxsy);
        sf_complexread (gxy, k, gxgy);
        sf_complexread (cxy, k, cxcy);

        /* Find CDP span */
        minix = onx - 1; maxix = 0;
        miniy = ony - 1; maxiy = 0;
        for (l = 0; l < k; l++) {
            ix = (int)((crealf (cxy[l]) - oox)/odx + 0.5f);
            iy = (int)((cimagf (cxy[l]) - ooy)/ody + 0.5f);
            if (ix < minix)
                minix = ix;
            if (ix > maxix)
                maxix = ix;
            if (iy < miniy)
                miniy = iy;
            if (iy > maxiy)
                maxiy = iy;
        }

        /* Aperture corners */
        el_cx1 = minix;
        el_cx2 = maxix;
        el_cy1 = miniy;
        el_cy2 = maxiy;
        /* Add apperture width */
        minix -= apx;
        if (minix < 0)
            minix = 0;
        miniy -= apy;
        if (miniy < 0)
            miniy = 0;
        maxix += apx;
        if (maxix >= onx)
            maxix = onx - 1;
        maxiy += apy;
        if (maxiy >= ony)
            maxiy = ony - 1;
        if (verb)
            sf_warning ("Rectangular aperture: %d-%d, %d-%d", minix, maxix, miniy, maxiy);
        /* Build aperture with rounded corners */
        l = 0;
        for (iy = miniy; iy <= maxiy; iy++) {
            for (ix = minix; ix <= maxix; ix++) {
                oidx = iy*onx + ix;
                if ((ix >= el_cx1 && ix <= el_cx2) ||
                    (iy >= el_cy1 && iy <= el_cy2)) {
                    ap[l] = oidx;
                    l++;
                    continue;
                }
                /* Distance to corners */
                if (ix < el_cx1)
                    el_x = ix - el_cx1;
                else
                    el_x = ix - el_cx2;
                if (iy < el_cy1)
                    el_y = iy - el_cy1;
                else
                    el_y = iy - el_cy2;
                /* Check if the point is within one of the ellipses */
                if ((el_x*el_x/(apx*apx) + el_y*el_y/(apy*apy)) < 1.0f) {
                    ap[l] = oidx;
                    l++;
                }
            }
        }

        /* Loop over input traces in the buffer */
        for (m = 0; m < k; m++) {
            sx = crealf (sxy[m]);
            sy = cimagf (sxy[m]);
            gx = crealf (gxy[m]);
            gy = cimagf (gxy[m]);
            /* Run antialiasing preparation */
            if (verb)
                sf_timer_start (aux_kernel_timer);
            if (diff)
                sf_ktmig_sbdiff (&t[m*nt], nt, dt);
            if (aa) {
                sf_ktmig_cint (&t[m*nt], nt);
                sf_ktmig_acint (&t[m*nt], nt);
            }
            if (verb) {
                sf_timer_stop (aux_kernel_timer);
                sf_timer_start (main_kernel_timer);
            }
            /* Loop over image traces within aperture */
            for (n = 0; n < l; n++) {
                oidx = ap[n];
                /* Skip positions outside of survey definition */
                if (v[oidx*ont] < odx)
                    continue;
                ox = oox + (oidx % onx)*odx;
                oy = ooy + (oidx / onx)*ody;
                sf_ktmig_kernel (&t[m*nt], &v[oidx*ont], &img[oidx*ont],
                                 ox, oy, sx, sy, gx, gy, nt, ont,
                                 ot, dt, oot, odt, maxtri, trfact, aa);
            }
            if (verb) {
                sf_timer_stop (main_kernel_timer);
                if ((i != 0 || m != 0) && (i + m) % 1000 == 0)
                    sf_warning ("Migrated %d traces", i + m);
            }
        }

        j++;
        i += k;
    } /* End of loop over input traces */

    sf_floatwrite (img, osize, image);

    if (verb || time)
        sf_timer_stop (total_timer);
    if (verb) {
        sf_warning ("*** Summary of wallclock time ***");
        sf_warning ("Auxillary kernels time: %f ms",
                    sf_timer_get_total_time (aux_kernel_timer));
        sf_warning ("Main kernel time: %f ms",
                    sf_timer_get_total_time (main_kernel_timer));
    }
    if (verb || time)
        sf_warning ("Total kernels + disk I/O time: %f ms",
                    sf_timer_get_total_time (total_timer));

    free (ap);
    free (v);
    free (img);
    free (t);
    free (sxy);
    free (gxy);
    free (cxy);

    if (verb) {
        free (main_kernel_timer);
        free (aux_kernel_timer);
    }
    if (verb || time)
        free (total_timer);

    exit (0);
}


/* Full angle-dependent slowness volume for 3-D reduced phase space. */
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

#include "esc_slow2.h"

int main (int argc, char* argv[]) {
    int nz, nx, na, nc;
    int iz, ix, ia;
    float dz, oz, dx, ox, da, oa, a;
    float savg = 0.0, smin = SF_HUGE, smax = -SF_HUGE;
    float corr_nterm = 0.0, corr = 0.0, new_savg = 0.0;
    float **s;
    bool verb;
    sf_file spdom, vspline = NULL, out;
    sf_esc_slowness2 esc_slow;

    sf_init (argc, argv);

    if (!sf_stdin ()) {
        spdom = NULL;
    } else {
        spdom = sf_input ("in");
        /* Spatial (z,x) domain */
    }

    out = sf_output ("out");
    /* Slowness values */

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
    if (!sf_getfloat ("dz", &dz) && !spdom) sf_error ("Need oz=");
    /* Sampling of z axis */
    if (!sf_getint ("nx", &nx) && !spdom) sf_error ("Need nx=");
    /* Number of samples in x axis */
    if (!sf_getfloat ("ox", &ox) && !spdom) sf_error ("Need ox=");
    /* Beginning of x axis */
    if (!sf_getfloat ("dx", &dx) && !spdom) sf_error ("Need ox=");
    /* Sampling of x axis */

    if (!sf_getint ("na", &na)) na = 360;
    /* Number of phase angles */

    da = 2.0*SF_PI/(float)na;
    oa = -SF_PI + 0.5*da;

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (spdom)
        sf_shiftdim (spdom, out, 1);

    /* Set up output */
    sf_putint (out, "n1", 4);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Slowness components");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", nz);
    sf_putfloat (out, "o2", oz);
    sf_putfloat (out, "d2", dz);
    if (!spdom) {
        sf_putstring (out, "label2", "Depth");
        sf_putstring (out, "unit2", "");
    }
    sf_putint (out, "n3", nx);
    sf_putfloat (out, "o3", ox);
    sf_putfloat (out, "d3", dx);
    if (!spdom) {
        sf_putstring (out, "label3", "Lateral");
        sf_putstring (out, "unit3", "");
    }
    sf_putint (out, "n4", na);
    sf_putfloat (out, "d4", da*180.0/SF_PI);
    sf_putfloat (out, "o4", oa*180.0/SF_PI);
    sf_putstring (out, "label4", "Angle");
    sf_putstring (out, "unit4", "Degrees");

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    if (!sf_histint (vspline, "Nc", &nc)) sf_error ("No Nc= in vspl");

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness2_init (vspline, verb);

    s = sf_floatalloc2 (4, nz);

    for (ia = 0; ia < na; ia++) {
        a = oa + ia*da;
        if (verb)
            sf_warning ("Computing angle plane %g;", a*180.0/SF_PI);
        for (ix = 0; ix < nx; ix++) {
            for (iz = 0; iz < nz; iz++) {
                 sf_esc_slowness2_get_components (esc_slow,
                                                  oz + iz*dz, ox + ix*dx, a,
                                                  &s[iz][0], &s[iz][3],
                                                  &s[iz][1], &s[iz][2]);
                if (s[iz][0] < smin)
                    smin = s[iz][0];
                if (s[iz][0] > smax)
                    smax = s[iz][0];
                /* Kahan summation algorithm to avoid roundoff errors
                   while accumulating the total sum */
                corr_nterm = s[iz][0] - corr;
                new_savg = savg + corr_nterm;
                corr = (new_savg - savg) - corr_nterm;
                savg = new_savg;
            }
            sf_floatwrite (&s[0][0], 4*nz, out);
        }
    }
    if (verb) {
        sf_warning (".");
        sf_warning ("Average velocity: %g [min=%g, max=%g]",
                    1.0/(savg/(float)((size_t)na*(size_t)nx*(size_t)nz)),
                    1.0/smax, 1.0/smin);
    }

    sf_esc_slowness2_close (esc_slow);

    free (s[0]);
    free (s);

    sf_fileclose (vspline);

    return 0;
}


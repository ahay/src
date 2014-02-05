/* Solution of escape equations by hybrid solver with narrow band for 2-D (an)isotropic media. */
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

#include "esc_nbgrid2.h"

int main (int argc, char* argv[]) {
    int nz, nx, na, morder;
    float dz, oz, dx, ox, mdist;
    sf_file spdom, vspline = NULL, out;

    bool verb, cmix, atraced, mtraced;
    sf_esc_slowness2 esc_slow;
    sf_esc_tracer2 esc_tracer;
    sf_esc_nbout2 esc_out;
    sf_esc_nbgrid2 esc_grid;

    sf_init (argc, argv);

    if (!sf_stdin ()) {
        spdom = NULL;
    } else {
        spdom = sf_input ("in");
        /* Spatial (z,x) domain */
    }

    out = sf_output ("out");
    /* Narrow band sweeps */

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

    if (!sf_getfloat ("mdist", &mdist)) mdist = SF_HUGE;
    /* Maximum distance between points in F-D stencil */
    
    if (!sf_getint ("morder", &morder)) morder = ESC2_MORDER;
    /* Maximum order in F-D stencil */

    if (!sf_getbool ("atraced", &atraced)) atraced = false;
    /* true - output map of all traced points */

    if (!sf_getbool ("mtraced", &mtraced)) mtraced = false;
    /* true - output map of points traced because of mdist criterion */

    if (!sf_getbool ("cmix", &cmix)) cmix = true;
    /* true - check for color mixing */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    /* Disk output module for the narrow band */
    esc_out = sf_esc_nbout2_init (nz, nx, na, 
                                  oz, ox, sf_esc_nbgrid2_get_oa (na),
                                  dz, dx, sf_esc_nbgrid2_get_da (na),
                                  atraced || mtraced, spdom, out);

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness2_init (vspline, verb);

    /* Ray tracer */
    esc_tracer = sf_esc_tracer2_init (esc_slow);

    /* Narrow band grid */
    esc_grid = sf_esc_nbgrid2_init (nz, nx, na, oz, ox, dz, dx, atraced, mtraced,
                                    esc_slow, esc_tracer, esc_out);
    sf_esc_nbgrid2_set_verb (esc_grid, verb);
    sf_esc_nbgrid2_set_cmix (esc_grid, cmix);
    sf_esc_nbgrid2_set_mdist (esc_grid, mdist);
    sf_esc_nbgrid2_set_morder (esc_grid, morder);

    /* Run hybrid solver computations in the narrow band */
    sf_esc_nbgrid2_compute (esc_grid);

    sf_esc_nbgrid2_close (esc_grid);
    sf_esc_tracer2_close (esc_tracer);
    sf_esc_slowness2_close (esc_slow);
    sf_esc_nbout2_close (esc_out);

    sf_fileclose (vspline);

    return 0;
}


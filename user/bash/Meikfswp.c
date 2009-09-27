/* Fast sweeping eikonal solver (2-D/3-D). */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "eikfswp.h"

int main (int argc,char* argv[]) {
    int n1, n2, n3, i, nshot, is, n123, ndim, niter;
    float o1, o2, o3, d1, d2, d3;
    float **s, *t, *v;
    char *sfile;
    bool isvel;
    sf_file vel, time, shots;

    sf_init (argc, argv);
    vel = sf_input ("in");
    time = sf_output ("out");

    if (SF_FLOAT != sf_gettype (vel))
        sf_error("Need float input");
    if (!sf_histint (vel, "n1", &n1)) sf_error ("No n1= in input");
    if (!sf_histint (vel, "n2", &n2)) sf_error ("No n2= in input");
    if (!sf_histint (vel, "n3", &n3)) n3 = 1;

    if (!sf_histfloat (vel, "d1", &d1)) sf_error ("No d1= in input");
    if (!sf_histfloat (vel, "d2", &d2)) sf_error ("No d2= in input");
    if (!sf_histfloat (vel, "d3", &d3)) d3 = d2;

    if (!sf_histfloat (vel, "o1", &o1)) o1 = 0.;
    if (!sf_histfloat (vel, "o2", &o2)) o2 = 0.;
    if (!sf_histfloat (vel, "o3", &o3)) o3 = 0.;

    if (!sf_getbool ("vel", &isvel)) isvel = true;
    /* if y, the input is velocity; n - slowness */

    if (!sf_getint ("niter", &niter)) niter = 2;
    /* number of sweeping iterations */

    sfile = sf_getstring ("shotfile");
    /* File with shot locations (n2=number of shots, n1=3) */

    if (NULL != sfile) {
        shots = sf_input ("shotfile");

        if (SF_FLOAT != sf_gettype (shots)) 
            sf_error ("Need float shotfile");
        if (!sf_histint (shots, "n2", &nshot))
            sf_error ("No n2= in shotfile");
        if (!sf_histint (shots, "n1", &ndim) || ndim != 3)
            sf_error ("Need n1=3 in shotfile");

        s = sf_floatalloc2 (ndim, nshot);
        sf_floatread (s[0], nshot * ndim, shots);

        sf_putint (time, 1 == n3 ? "n3" : "n4", nshot);
        free (sfile);
    } else {
        nshot = 1;
        ndim = 3;

        s = sf_floatalloc2 (ndim, nshot);

        if (!sf_getfloat ("zshot", &s[0][0])) s[0][0] = 0.; 
        /* Shot location (used if no shotfile) */
        if (!sf_getfloat ("yshot", &s[0][1])) s[0][1] = o2 + 0.5*(n2 - 1)*d2;
        if (!sf_getfloat ("xshot", &s[0][2])) s[0][2] = o3 + 0.5*(n3 - 1)*d3;

        sf_warning ("Shooting from zshot=%g yshot=%g xshot=%g",
                    s[0][0], s[0][1], s[0][2]);
    }

    n123 = n1 * n2 * n3;

    t  = sf_floatalloc (n123);
    v  = sf_floatalloc (n123);

    sf_floatread (v, n123, vel);
   /* transform velocity to slowness */
    if (isvel) {
        for (i = 0; i < n123; i++) {
            v[i] = 1. / v[i];
        }
    }

    sf_warning ("Performing %d-D sweeps", 2 + (n3 > 1));

    /* loop over shots */
    for (is = 0; is < nshot; is++) {
        if (nshot > 1)
            sf_warning ("Calculating shot %d of %d", is + 1, nshot);
        if (false == sf_init_fast_sweep (t,
                                         n3, n2, n1,
                                         o3, o2, o1,
                                         d3, d2, d1,
                                         s[is][0], s[is][1], s[is][2]))
            sf_error ("Incorrect shot location");

        sf_run_fast_sweep (t, v, niter,
                           n3, n2, n1,
                           o3, o2, o1,
                           d3, d2, d1);

        sf_floatwrite (t, n123, time);
    }
    sf_close();
    exit (0);
}

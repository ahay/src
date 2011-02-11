/* Adjoint state method for first-arrival traveltimes in 2-D. */
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

#include "adjfswp2.h"

int main (int argc,char* argv[]) {
    int n1, n2, nshot, is, niter, nh;
    float d1, d2, o1, o2, oh, dh;
    float **t1, **t2 = NULL, *dt, **l, *h = NULL;
    char *sfile = NULL;
    sf_eno horiz = NULL;
    sf_file time, lambda, horizon, deltat;

    sf_init (argc, argv);
    time = sf_input ("in");
    lambda = sf_output ("out");

    if (SF_FLOAT != sf_gettype (time))
        sf_error ("Need float input");
    if (!sf_histint (time, "n1", &n1)) sf_error ("No n1= in input");
    if (!sf_histint (time, "n2", &n2)) sf_error ("No n2= in input");
    if (!sf_histfloat (time, "d1", &d1)) sf_error ("No d1= in input");
    if (!sf_histfloat (time, "d2", &d2)) sf_error ("No d2= in input");
    if (!sf_histfloat (time, "o1", &o1)) o1 = 0.;
    if (!sf_histfloat (time, "o2", &o2)) o2 = 0.;

    if (!sf_getint ("niter", &niter)) niter = 2;
    /* number of sweeping iterations */

    sfile = sf_getstring ("deltat");
    /* File with traveltime differences */
    if (sfile) {
       deltat = sf_input (sfile);
       /* File with traveltime differences */
       free (sfile); sfile = NULL;
    } else
        sf_error ("Need deltat=");

    t1  = sf_floatalloc2 (n1, n2);
    t2  = sf_floatalloc2 (n1, n2);
    l  = sf_floatalloc2 (n1, n2);
    dt = sf_floatalloc (n2);

    sfile = sf_getstring ("horizon");
    /* File with a reflection interface */
    if (sfile) {
       horizon = sf_input (sfile);
       /* File with a reflection interface */
       if (!sf_histint (horizon, "n1", &nh)) nh = n2;
       if (!sf_histfloat (horizon, "o1", &oh)) oh = o2;
       if (!sf_histfloat (horizon, "d1", &dh)) dh = d2;

       if (!sf_histint (time, "n4", &nshot)) sf_error ("No n4= in input");

       h = sf_floatalloc (nh);
       sf_floatread (h, nh, horizon);

       horiz = sf_eno_init (5, nh);
       sf_eno_set (horiz, h);
       /* Collapse extra dimension */
       sf_putint (lambda, "n3", nshot);
       sf_putint (lambda, "n4", 1);

       free (sfile); sfile = NULL;
    } else
       if (!sf_histint (time, "n3", &nshot)) sf_error ("No n3= in input");

    /* loop over shots */
    for (is = 0; is < nshot; is++) {
        if (nshot > 1)
            sf_warning ("Calculating shot %d of %d", is + 1, nshot);

        sf_floatread (t1[0], n1*n2, time);
        sf_floatread (t2[0], n1*n2, time);
        sf_floatread (dt, n2, deltat);

        sf_init_lambda_2d_sweep (dt, l, t2, n2, n1, d2, d1);

        sf_run_lambda_2d_sweep2 (t2, l, niter, n2, n1, d2, d1);
/*
        if (horiz) {
        }
*/
        sf_floatwrite (l[0], n1*n2, lambda);
    }
    if (horiz)
        sf_eno_close (horiz);

    return 0;
}


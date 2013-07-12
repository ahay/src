/* Compute distance and traveltime difference between two escape tables. */
/*
  Copyright (C) 2013 University of Texas at Austin

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
#include "esc_point3.h"

int main (int argc, char* argv[]) {
    int i, ia, n, na;
    float maxd, maxt, mask;
    float **esc0, **esc1, **diff;
    bool is3d = false;
    sf_file esct0, esct1 = NULL, out;

    sf_init (argc, argv);

    esct0 = sf_input ("in");
    /* First set of escape tables */
    out = sf_output ("out");
    /* Difference */

    /* Phase space dimensions (also, imaging dimensions) */
    if (!sf_histint (esct0, "n1", &na)) sf_error ("No n1= in input");
    if (na != ESC2_NUM && na != ESC3_NUM) sf_error ("Need n1=%d or n1=%d in input", ESC2_NUM, ESC3_NUM);
    is3d = (ESC3_NUM == na);
    if (!sf_histint (esct0, "n2", &na)) sf_error ("No n2= in input");
    n = sf_leftsize (esct0, 2);

    if (sf_getstring ("esct")) {
        /* Second set of escape tables */
        esct1 = sf_input ("esct");
    } else {
        sf_error ("Need esct=");
    }

    if (!sf_getfloat ("maxd", &maxd)) maxd = SF_HUGE;
    /* Maximum allowed distance */
    if (!sf_getfloat ("maxt", &maxt)) maxt = SF_HUGE;
    /* Maximum allowed time */
    if (!sf_getfloat ("mask", &mask)) mask = SF_HUGE;
    /* Mask for values above maxd= and maxt= thresholds */

    sf_putint (out, "n1", 2);
    sf_putstring (out, "label1", "Diff");

    esc0 = sf_floatalloc2 (is3d ? ESC3_NUM : ESC2_NUM, na);
    esc1 = sf_floatalloc2 (is3d ? ESC3_NUM : ESC2_NUM, na);
    diff = sf_floatalloc2 (2, na);

    for (i = 0; i < n; i++) {
        sf_floatread (esc0[0], is3d ? ESC3_NUM*na : ESC2_NUM*na, esct0);
        sf_floatread (esc1[0], is3d ? ESC3_NUM*na : ESC2_NUM*na, esct1);
        for (ia = 0; ia < na; ia++) {
            if (is3d) {
                diff[ia][0] = sqrtf ((esc1[ia][ESC3_Z] - esc0[ia][ESC3_Z])*(esc1[ia][ESC3_Z] - esc0[ia][ESC3_Z]) +
                                     (esc1[ia][ESC3_X] - esc0[ia][ESC3_X])*(esc1[ia][ESC3_X] - esc0[ia][ESC3_X]) +
                                     (esc1[ia][ESC3_Y] - esc0[ia][ESC3_Y])*(esc1[ia][ESC3_Y] - esc0[ia][ESC3_Y]));
                diff[ia][1] = fabsf (esc1[ia][ESC3_T] - esc0[ia][ESC3_T]);
            } else {
                diff[ia][0] = sqrtf ((esc1[ia][ESC2_Z] - esc0[ia][ESC2_Z])*(esc1[ia][ESC2_Z] - esc0[ia][ESC2_Z]) +
                                     (esc1[ia][ESC2_X] - esc0[ia][ESC2_X])*(esc1[ia][ESC2_X] - esc0[ia][ESC2_X]));
                diff[ia][1] = fabsf (esc1[ia][ESC2_T] - esc0[ia][ESC2_T]);
            }
            if (diff[ia][0] > maxd)
                diff[ia][0] = mask;
            if (diff[ia][1] > maxt)
                diff[ia][1] = mask;
        }
        sf_floatwrite (diff[0], 2*na, out);
    }
    free (diff[0]);
    free (diff);
    free (esc1[0]);
    free (esc1);
    free (esc0[0]);
    free (esc0);

    sf_fileclose (esct1);

    return 0;
}


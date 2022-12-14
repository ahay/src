/* Calculate datum shift from elevation profile for 2-D shot gathers */
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

/* 2-D only for now, should add 3-D later */

#include <rsf.h>

int main (int argc, char* argv[]) {
    int n2, n3, i2, i3, nx, i;
    float o2, d2, o3, d3, ox, dx, v0;
    float sx, rx, sstat, rstat, f;
    float *profile, *recshift;
    sf_file in, out, elev;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype (in)) sf_error ("Need float input");

    /* Offset axis */
    if (!sf_histint (in, "n2", &n2)) sf_error ("Need n2= in input");
    if (!sf_histfloat (in,"d2", &d2)) sf_error ("Need d2= in input");
    if (!sf_histfloat (in,"o2", &o2)) sf_error ("Need o2= in input");
    /* Shot axis */
    if (!sf_histint (in, "n3", &n3)) sf_error ("Need n3= in input");
    if (!sf_histfloat (in,"d3", &d3)) sf_error ("Need d3= in input");
    if (!sf_histfloat (in,"o3", &o3)) sf_error ("Need o3= in input");

    if (NULL == sf_getstring ("elev")) sf_error ("Need elev=");
    elev = sf_input ("elev");
    /* Elevation profile */

    if (!sf_histint (elev, "n1", &nx)) sf_error ("Need n1= in elev=");
    if (!sf_histfloat (elev, "d1", &dx)) sf_error ("Need d1= in elev=");
    if (!sf_histfloat (elev, "o1", &ox)) sf_error ("Need o1= in elev=");

    if (!sf_getfloat ("v0", &v0)) sf_error ("Need v0=");

    recshift = sf_floatalloc (n2);
    profile = sf_floatalloc (nx);

    sf_floatread (profile, nx, elev);

    sf_unshiftdim (in, out, 1);

    for (i3 = 0; i3 < n3; i3++) { /* Shots */
        sx = o3 + i3*d3; /* Shot coordinate */
        f = (sx - ox)/dx; /* Elevation coordinate */
        i = (int)f; f -= (float)i;
        if (i < 0 || i > (nx - 2)) sf_error ("Shot #%d is out of profile", i3);
        /* Shot shift */
        sstat = (1.0 - f)*profile[i] + f*profile[i + 1];
        for (i2 = 0; i2 < n2; i2++) { /* Receivers */
            rx = sx + o2 + i2*d2; /* Receiver coordinate */
            f = (rx - ox)/dx; /* Elevation coordinate */
            i = (int)f; f -= (float)i;
            if (i < 0 || i > (nx - 2)) sf_error ("Receiver #%d for shot #%d is out of profile", i2, i3);
            /* Receiver shift */
            rstat = (1.0 - f)*profile[i] + f*profile[i + 1];
            recshift[i2] = (sstat + rstat)/v0;
        }
        sf_floatwrite (recshift, n2, out);
    }

    return 0;
}


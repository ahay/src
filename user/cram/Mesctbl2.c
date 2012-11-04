/* Esctape tables from solution of escape equations by the hybrid solver with narrow band. */
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
    int nz, nx, na, nc;
    int iz, ix;
    off_t off, off2;
    float dz, oz, dx, ox;
    float *buf1, *buf2;
    bool verb, otraced;
    sf_file in, out;

    sf_init (argc, argv);

    in = sf_input ("in");
    /* Narrow band sweeps */
    out = sf_output ("out");
    /* Escape values */

    if (!sf_histbool (in, "otraced", &otraced)) otraced = false;
    nc = otraced ? ESC2_NUM + 1 : ESC2_NUM;

    if (!sf_histint (in, "n1", &na)) sf_error ("No n1= in input");
    if (na != nc) sf_error ("Need n1=%d in input", nc);
    if (!sf_histint (in, "n5", &na)) sf_error ("No n5= in input");
    if (na != 2) sf_error ("Need n5=2 in input");

    /* Phase space dimensions */
    if (!sf_histint (in, "n2", &na)) sf_error ("No n2= in input");
    if (!sf_histint (in, "n3", &nx)) sf_error ("No n3= in input");
    if (!sf_histint (in, "n4", &nz)) sf_error ("No n4= in input");
    if (!sf_histfloat (in, "d3", &dx)) sf_error ("No d2= in input");
    if (!sf_histfloat (in, "o3", &ox)) ox = 0.0;
    if (!sf_histfloat (in, "d4", &dz)) sf_error ("No d1= in input");
    if (!sf_histfloat (in, "o4", &oz)) oz = 0.0;

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    /* Set up output */
    sf_putint (out, "n1", nc);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Escape variable");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", 2*na);
    sf_putfloat (out, "d2", sf_esc_nbgrid2_get_da (2*na)*180.0/SF_PI);
    sf_putfloat (out, "o2", sf_esc_nbgrid2_get_oa (2*na)*180.0/SF_PI);
    sf_oaxa (out, sf_iaxa (in, 4), 3);
    sf_oaxa (out, sf_iaxa (in, 3), 4);
    sf_putint (out, "n5", 1);
    sf_putfloat (out, "o5", 0.0);
    sf_putfloat (out, "d5", 1.0);
    sf_putstring (out, "label1", "");
    sf_putstring (out, "unit1", "");

    buf1 = sf_floatalloc (na*nc);
    buf2 = sf_floatalloc (na*nc);

    /* Starting offset in data */
    off = sf_tell (in);
    /* Offset to the upgoing part */ 
    off2 = (off_t)na*(off_t)nx*(off_t)nz;

    for (ix = 0; ix < nx; ix++) {
        if (verb)
            sf_warning ("Location %d of %d (%g);", ix + 1, nx, ox + ix*dx);
        for (iz = 0; iz < nz; iz++) {
            /* Read downgoing range - [-PI/2; PI/2) */
            sf_seek (in, off + ((off_t)na*(off_t)ix +
                                (off_t)na*(off_t)nx*(off_t)iz)*
                                (off_t)(nc*sizeof (float)), SEEK_SET);
            sf_floatread (buf1, (size_t)nc*(size_t)na, in);
            /* Read upgoing range - [PI/2; -PI/2) */
            sf_seek (in, off + (off2 + (off_t)na*(off_t)ix +
                                       (off_t)na*(off_t)nx*(off_t)(nz - iz - 1))*
                                (off_t)(nc*sizeof (float)), SEEK_SET);
            sf_floatread (buf2, (size_t)nc*(size_t)na, in);
            /* Combine output */
            sf_floatwrite (&buf2[nc*na/2], (size_t)nc*(size_t)na/2,
                           out); /* [-PI; -PI/2) */
            sf_floatwrite (buf1, (size_t)nc*(size_t)na,
                           out); /* [-PI/2; PI/2) */
            sf_floatwrite (&buf2[0], (size_t)nc*(size_t)na/2,
                           out); /* [PI/2; PI) */
        }
    }
    if (verb)
        sf_warning (".");

    free (buf2);
    free (buf1);

    return 0;
}


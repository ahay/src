/* 3-D phase-space escape variables output for the narrow band */
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

#ifndef _esc_nbout2_h

#include "esc_point2.h"
/*^*/

typedef struct EscNBOutput2 *sf_esc_nbout2;
/* abstract data type */
/*^*/

#endif

struct EscNBOutput2 {
    int         na;
    bool        otraced;
    float     **buf;
    sf_file     out;
};
/* concrete data type */

sf_esc_nbout2 sf_esc_nbout2_init (int nz, int nx, int na,
                                  float oz, float ox, float oa,
                                  float dz, float dx, float da,
                                  bool otraced, sf_file in, sf_file out)
/*< Initialize object >*/
{
    sf_esc_nbout2 esc_nbout = (sf_esc_nbout2)sf_alloc (1, sizeof (struct EscNBOutput2));

    esc_nbout->na = na/2;
    
    esc_nbout->buf = sf_floatalloc2 (otraced ? ESC2_NUM + 1 : ESC2_NUM, na/2);

    esc_nbout->otraced = otraced; /* flag: if true, output map of traced points */

    esc_nbout->out = out;

    sf_putint (out, "n1", otraced ? ESC2_NUM + 1 : ESC2_NUM);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Escape variable");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", na/2);
    sf_putfloat (out, "d2", da*180.0/SF_PI);
    sf_putfloat (out, "o2", oa*180.0/SF_PI);
    sf_putstring (out, "label2", "Angle");
    sf_putstring (out, "unit2", "Degrees");

    if (in) {
        sf_oaxa (out, sf_iaxa (in, 1), 4);
        sf_oaxa (out, sf_iaxa (in, 2), 3);
        sf_putint (out, "n3", nx);
        sf_putfloat (out, "o3", ox);
        sf_putint (out, "n4", nz);
        sf_putfloat (out, "o4", oz);
    } else {
        sf_putint (out, "n3", nx);
        sf_putfloat (out, "d3", dx);
        sf_putfloat (out, "o3", ox);
        sf_putstring (out, "label3", "Lateral");
        sf_putstring (out, "unit3", "");
        sf_putint (out, "n4", nz);
        sf_putfloat (out, "d4", dz);
        sf_putfloat (out, "o4", oz);
        sf_putstring (out, "label4", "Depth");
        sf_putstring (out, "unit4", "");
    }

    sf_putint (out, "n5", 2);
    sf_putfloat (out, "o5", (oa + da*na/4.0)*180/SF_PI);
    sf_putfloat (out, "d5", (da*na/2.0)*180/SF_PI);
    sf_putstring (out, "label1", "Narrow band sweeps");
    sf_putstring (out, "unit1", "");

    sf_putstring (out, "otraced", otraced ? "y" : "n");

    return esc_nbout;
}

void sf_esc_nbout2_close (sf_esc_nbout2 esc_nbout)
/*< Destroy object >*/
{
    free (esc_nbout->buf[0]);
    free (esc_nbout->buf);
    free (esc_nbout);
}

void sf_esc_nbout2_write_plane (sf_esc_nbout2 esc_nbout, sf_esc_point2 plane, unsigned long n)
/*< Write points in plane to disk >*/
{
    int i, j, sz;
    unsigned long nc = 0;
    unsigned char *points = (unsigned char*)plane;

    sz = sf_esc_point2_sizeof ();
    while (nc < n) {
        for (i = 0; i < esc_nbout->na; i++) {
            for (j = 0; j < ESC2_NUM; j++) {
                esc_nbout->buf[i][j] = sf_esc_point2_get_esc_var ((sf_esc_point2)&points[nc*sz], j);
            }
            if (esc_nbout->otraced)
                esc_nbout->buf[i][ESC2_NUM] = (float)sf_esc_point2_is_traced ((sf_esc_point2)&points[nc*sz]);
            nc++;
        }
        sf_floatwrite (esc_nbout->buf[0], esc_nbout->otraced ? 
                                          (size_t)(esc_nbout->na*(ESC2_NUM + 1)) :
                                          (size_t)(esc_nbout->na*ESC2_NUM),
                       esc_nbout->out);
    }
}


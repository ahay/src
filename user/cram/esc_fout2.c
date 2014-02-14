/* 3-D phase-space escape variables output */
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

#ifndef _esc_fout2_h

#include "esc_point2.h"
/*^*/

typedef struct EscFOutput2 *sf_esc_fout2;
/* abstract data type */
/*^*/

#endif

struct EscFOutput2 {
    int              nz, nx, na;
    float            oz, ox, oa;
    float            dz, dx, da;
    unsigned char ***t;
    bool             otraced;
    sf_file          out;
};
/* concrete data type */

sf_esc_fout2 sf_esc_fout2_init (int nz, int nx, int na,
                                float oz, float ox, float oa,
                                float dz, float dx, float da,
                                bool otraced, sf_file in, sf_file out)
/*< Initialize object >*/
{
    sf_esc_fout2 esc_fout = (sf_esc_fout2)sf_alloc (1, sizeof (struct EscFOutput2));

    esc_fout->nz = nz;
    esc_fout->nx = nx;
    esc_fout->na = na;
    esc_fout->oz = oz;
    esc_fout->ox = ox;
    esc_fout->oa = oa;
    esc_fout->dz = dz;
    esc_fout->dx = dx;
    esc_fout->da = da;

    esc_fout->t = sf_ucharalloc3 (nz*sf_esc_point2_sizeof (), nx, na);
    memset (esc_fout->t[0][0], 0, (size_t)nz*(size_t)nx*(size_t)na*
                                  (size_t)sf_esc_point2_sizeof ());

    esc_fout->out = out;

    esc_fout->otraced = otraced; /* flag: if true, output map of traced points */

    sf_putint (out, "n1", otraced ? ESC2_NUM + 1 : ESC2_NUM);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Escape variable");
    sf_putstring (out, "unit1", "");

    sf_putint (out, "n2", na);
    sf_putfloat (out, "d2", da*180.0/SF_PI);
    sf_putfloat (out, "o2", oa*180.0/SF_PI);
    sf_putstring (out, "label2", "Angle");
    sf_putstring (out, "unit2", "Degrees");

    if (in) {
        sf_oaxa (out, sf_iaxa (in, 1), 3);
        sf_oaxa (out, sf_iaxa (in, 2), 4);
        sf_putint (out, "n3", nz);
        sf_putfloat (out, "o3", oz);
        sf_putint (out, "n4", nx);
        sf_putfloat (out, "o4", ox);
    } else {
        sf_putint (out, "n3", nz);
        sf_putfloat (out, "d3", dz);
        sf_putfloat (out, "o3", oz);
        sf_putstring (out, "label3", "Depth");
        sf_putstring (out, "unit3", "");
        sf_putint (out, "n4", nx);
        sf_putfloat (out, "d4", dx);
        sf_putfloat (out, "o4", ox);
        sf_putstring (out, "label4", "Lateral");
        sf_putstring (out, "unit4", "");
    }

    sf_putstring (out, "otraced", otraced ? "y" : "n");

    return esc_fout;
}

void sf_esc_fout2_close (sf_esc_fout2 esc_fout)
/*< Destroy object >*/
{
    free (esc_fout->t[0][0]);
    free (esc_fout->t[0]);
    free (esc_fout->t);
    free (esc_fout);
}

sf_esc_point2 sf_esc_fout2_get_point (sf_esc_fout2 esc_fout, int iz, int ix, int ia)
/*< Return a point from the output array >*/
{
    if (iz < 0 || iz >= esc_fout->nz ||
        ix < 0 || ix >= esc_fout->nx)
        return NULL;

    while (ia < 0)
        ia += esc_fout->na;
    while (ia >= esc_fout->na)
        ia -= esc_fout->na;

    return (sf_esc_point2)&esc_fout->t[ia][ix][iz*sf_esc_point2_sizeof ()];
}

void sf_esc_fout2_write (sf_esc_fout2 esc_fout)
/*< Write the output array to disk >*/
{
    int sz, i, iz, ix, ia;
    float **buf;

    buf = sf_floatalloc2 (esc_fout->otraced ? ESC2_NUM + 1 : ESC2_NUM, esc_fout->na);
    sz = sf_esc_point2_sizeof ();

    for (ix = 0; ix < esc_fout->nx; ix++) {
        for (iz = 0; iz < esc_fout->nz; iz++) {
            for (ia = 0; ia < esc_fout->na; ia++) {
                for (i = 0; i < ESC2_NUM; i++) {
                    buf[ia][i] = sf_esc_point2_get_esc_var ((sf_esc_point2)&esc_fout->t[ia][ix][iz*sz], i);
                }
                if (esc_fout->otraced)
                    buf[ia][ESC2_NUM] = (float)sf_esc_point2_is_traced ((sf_esc_point2)&esc_fout->t[ia][ix][iz*sz]);
            }

            sf_floatwrite (buf[0], esc_fout->otraced ? 
                                   (size_t)(esc_fout->na*(ESC2_NUM + 1)) :
                                   (size_t)(esc_fout->na*ESC2_NUM),
                           esc_fout->out);
        }
    }
    free (buf[0]);
    free (buf);
}


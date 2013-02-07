/* 3-D slowness on request for angle migration */
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

#ifndef _cram_slowness3_h

typedef struct CRAMSlow3 *sf_cram_slowness3;
/* abstract data type */
/*^*/

#endif

struct CRAMSlow3 {
    int       nz, nx, ny;
    float     dz, dx, dy;
    float     z0, x0, y0;
    float     sconst;
    off_t     offs;
    float   **surf;
    sf_file   fvz;
};
/* concrete data type */

sf_cram_slowness3 sf_cram_slowness3_init (sf_file fvz, float vconst)
/*< Initialize object >*/
{
    int ix, iy;
    float vz;
    sf_cram_slowness3 cram_slowness = (sf_cram_slowness3)sf_alloc (1, sizeof (struct CRAMSlow3));

    cram_slowness->sconst = 1.0/vconst;

    cram_slowness->fvz = fvz;
    if (fvz) {
        if (!sf_histint (fvz, "n1", &cram_slowness->nz)) sf_error ("No n1= in slow");
        if (!sf_histint (fvz, "n2", &cram_slowness->nx)) sf_error ("No n2= in slow");
        if (!sf_histint (fvz, "n3", &cram_slowness->ny)) sf_error ("No n3= in slow");
        if (!sf_histfloat (fvz, "d1", &cram_slowness->dz)) sf_error ("No d1= in slow");
        if (!sf_histfloat (fvz, "o1", &cram_slowness->z0)) sf_error ("No o1= in slow");
        if (!sf_histfloat (fvz, "d2", &cram_slowness->dx)) sf_error ("No d2= in slow");
        if (!sf_histfloat (fvz, "o2", &cram_slowness->x0)) sf_error ("No o2= in slow");
        if (!sf_histfloat (fvz, "d3", &cram_slowness->dy)) sf_error ("No d3= in slow");
        if (!sf_histfloat (fvz, "o3", &cram_slowness->y0)) sf_error ("No o3= in slow");
    }
    cram_slowness->surf = cram_slowness->fvz ? sf_floatalloc2 (cram_slowness->nx,
                                                               cram_slowness->ny) : NULL;
    if (cram_slowness->fvz) { /* Read surface values */
        cram_slowness->offs = sf_tell (cram_slowness->fvz);
        for (iy = 0; iy < cram_slowness->ny; iy++) {
            for (ix = 0; ix < cram_slowness->nx; ix++) {
                sf_seek (cram_slowness->fvz, ((size_t)iy*(size_t)cram_slowness->nx*
                                                         (size_t)cram_slowness->nz +
                                              (size_t)ix*(size_t)cram_slowness->nz)*
                                             (size_t)sizeof(float) + cram_slowness->offs, SEEK_SET);
                sf_floatread (&vz, 1, cram_slowness->fvz);
                cram_slowness->surf[iy][ix] = 1.0/vz;
            }
        }
    }

    return cram_slowness;
}

void sf_cram_slowness3_close (sf_cram_slowness3 cram_slowness)
/*< Destroy object >*/
{
    if (cram_slowness->surf) {
        free (cram_slowness->surf[0]);
        free (cram_slowness->surf);
    }
    free (cram_slowness);
}

float sf_cram_slowness3_get_surf_value (sf_cram_slowness3 cram_slowness, float x, float y)
/*< Get slowness value for a surface location (x,y) >*/
{
    int ix, iy;
    float xf, yf;

    if (NULL == cram_slowness->surf)
        return cram_slowness->sconst;

    xf = (x - cram_slowness->x0)/cram_slowness->dx;
    ix = floorf (xf);
    xf -= (float)ix;
    if (ix < 0) {
        ix = 0;
        xf = 0.0;
    }
    if (ix >= (cram_slowness->nx - 1)) {
        ix = cram_slowness->nx - 2;
        xf = 1.0;
    }

    yf = (y - cram_slowness->y0)/cram_slowness->dy;
    iy = floorf (yf);
    yf -= (float)iy;
    if (iy < 0) {
        iy = 0;
        yf = 0.0;
    }
    if (iy >= (cram_slowness->ny - 1)) {
        iy = cram_slowness->ny - 2;
        yf = 1.0;
    }
    return (cram_slowness->surf[iy][ix]*(1.0 - yf)*(1.0 - xf)
            + cram_slowness->surf[iy][ix + 1]*(1.0 - yf)*xf
            + cram_slowness->surf[iy + 1][ix]*yf*(1.0 - xf)
            + cram_slowness->surf[iy + 1][ix + 1]*yf*xf);
}

float sf_cram_slowness3_get_value (sf_cram_slowness3 cram_slowness, float z, float x, float y)
/*< Get slowness value for a location (z,x,y) >*/
{
    float vz;
    int iz, ix, iy;

    if (cram_slowness->fvz) {
        /* Use nearest neighbor */
        iz = (int)((z - cram_slowness->z0)/cram_slowness->dz + 0.5);
        ix = (int)((x - cram_slowness->x0)/cram_slowness->dx + 0.5);
        iy = (int)((y - cram_slowness->y0)/cram_slowness->dy + 0.5);
        if (iz < 0)
            iz = 0;
        if (iz >= cram_slowness->nz)
            iz = cram_slowness->nz - 1;
        if (ix < 0)
            ix = 0;
        if (ix >= cram_slowness->nx)
            ix = cram_slowness->nx - 1;
        if (iy < 0)
            iy = 0;
        if (iy >= cram_slowness->ny)
            iy = cram_slowness->ny - 1;

        sf_seek (cram_slowness->fvz, ((size_t)iy*(size_t)cram_slowness->nx*
                                                 (size_t)cram_slowness->nz +
                                      (size_t)ix*(size_t)cram_slowness->nz + (size_t)iz)*
                                     (size_t)sizeof(float) + cram_slowness->offs, SEEK_SET);
        sf_floatread (&vz, 1, cram_slowness->fvz);
        return 1.0/vz;
    } else
        return cram_slowness->sconst;  
}


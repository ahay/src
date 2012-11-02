/* 2-D slowness on request for angle migration */
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

#include <rsf.h>

#ifndef _cram_slowness2_h

typedef struct CRAMSlow2 *sf_cram_slowness2;
/* abstract data type */
/*^*/

#endif

struct CRAMSlow2 {
    int      nz, nx;
    float    dz, dx;
    float    z0, x0;
    float    sconst;
    off_t    offs;
    float   *surf;
    sf_file  fvz;
};
/* concrete data type */

sf_cram_slowness2 sf_cram_slowness2_init (sf_file fvz, float vconst)
/*< Initialize object >*/
{
    int ix;
    float vz;
    sf_cram_slowness2 cram_slowness = (sf_cram_slowness2)sf_alloc (1, sizeof (struct CRAMSlow2));

    cram_slowness->sconst = 1.0/vconst;

    cram_slowness->fvz = fvz;
    if (fvz) {
        if (!sf_histint (fvz, "n1", &cram_slowness->nz)) sf_error ("No n1= in slow");
        if (!sf_histint (fvz, "n2", &cram_slowness->nx)) sf_error ("No n2= in slow");
        if (!sf_histfloat (fvz, "d1", &cram_slowness->dz)) sf_error ("No d1= in slow");
        if (!sf_histfloat (fvz, "o1", &cram_slowness->z0)) sf_error ("No o1= in slow");
        if (!sf_histfloat (fvz, "d2", &cram_slowness->dx)) sf_error ("No d2= in slow");
        if (!sf_histfloat (fvz, "o2", &cram_slowness->x0)) sf_error ("No o2= in slow");
    }
    cram_slowness->surf = cram_slowness->fvz ? sf_floatalloc (cram_slowness->nx) : NULL;
    if (cram_slowness->fvz) { /* Read surface values */
        cram_slowness->offs = sf_tell (cram_slowness->fvz);
        for (ix = 0; ix < cram_slowness->nx; ix++) {
            sf_seek (cram_slowness->fvz, (size_t)ix*(size_t)cram_slowness->nz*
                                         (size_t)sizeof(float) + cram_slowness->offs, SEEK_SET);
            sf_floatread (&vz, 1, cram_slowness->fvz);
            cram_slowness->surf[ix] = 1.0/vz;
        }
    }

    return cram_slowness;
}

void sf_cram_slowness2_close (sf_cram_slowness2 cram_slowness)
/*< Destroy object >*/
{
    if (cram_slowness->surf)
        free (cram_slowness->surf);
    free (cram_slowness);
}

float sf_cram_slowness2_get_surf_value (sf_cram_slowness2 cram_slowness, float x)
/*< Get slowness value for a surface location x >*/
{
    int ix;
    float xf;

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
    return (cram_slowness->surf[ix]*(1.0 - xf)
            + cram_slowness->surf[ix + 1]*xf);
}

float sf_cram_slowness2_get_value (sf_cram_slowness2 cram_slowness, float z, float x)
/*< Get slowness value for a location (z,x) >*/
{
    float vz;
    int iz, ix;

    if (cram_slowness->fvz) {
        /* Use nearest neighbor */
        iz = (int)((z - cram_slowness->z0)/cram_slowness->dz + 0.5);
        ix = (int)((x - cram_slowness->x0)/cram_slowness->dx + 0.5);
        if (iz < 0)
            iz = 0;
        if (iz >= cram_slowness->nz)
            iz = cram_slowness->nz - 1;
        if (ix < 0)
            ix = 0;
        if (ix >= cram_slowness->nx)
            ix = cram_slowness->nx - 1;

        sf_seek (cram_slowness->fvz, ((size_t)ix*(size_t)cram_slowness->nz +
                                      (size_t)iz)*(size_t)sizeof (float) + cram_slowness->offs, SEEK_SET);
        sf_floatread (&vz, 1, cram_slowness->fvz);
        return 1.0/vz;
    } else
        return cram_slowness->sconst;  
}


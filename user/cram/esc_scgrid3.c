/* 5-D phase-space escape values grid consisting of supercells */
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

#ifndef _esc_scgrid3_h

typedef struct EscSCgrid3 *sf_esc_scgrid3;
/* abstract data type */
/*^*/

#endif

#include "esc_scbl3.h"

struct EscSCgrid3 {
    size_t           ns, nsxz, nsz, nrt, nr;
    double           drt;
    int              nscz, nscx, nscy;
    float            oscz, oscx, oscy;
    float            dscz, dscx, dscy;
    sf_esc_scblock3 *cells;
};
/* concrete data type */

sf_esc_scgrid3 sf_esc_scgrid3_init (sf_file scgrid, sf_esc_tracer3 esc_tracer, bool verb)
/*< Initialize object >*/
{
    off_t sz = 0;
    size_t i, iscz, iscx, iscy;
    unsigned char sbuf[ESC3_BSC_BMAX];
    sf_file fcell;
    sf_simtab simtab;
    sf_esc_scgrid3 esc_scgrid = (sf_esc_scgrid3)sf_alloc (1, sizeof (struct EscSCgrid3));

    if (sf_gettype (scgrid) != SF_UCHAR)
        sf_error ("Wrong data type in supercell grid file");

    if (!sf_histint (scgrid, "n2", &esc_scgrid->nscz)) sf_error ("No n1= in supercell grid");
    /* Number of supercells in z */
    if (!sf_histint (scgrid, "n3", &esc_scgrid->nscx)) sf_error ("No n2= in supercell grid");
    /* Number of supercells in x */
    if (!sf_histint (scgrid, "n4", &esc_scgrid->nscy)) sf_error ("No n3= in supercell grid");
    /* Number of supercells in y */

    if (!sf_histfloat (scgrid, "Dscz", &esc_scgrid->dscz)) sf_error ("No Dscz= in supercell grid");
    /* Supercell size in z  */
    if (!sf_histfloat (scgrid, "Dscx", &esc_scgrid->dscx)) sf_error ("No Dscx= in supercell grid");
    /* Supercell size in x */
    if (!sf_histfloat (scgrid, "Dscy", &esc_scgrid->dscy)) sf_error ("No Dscy= in supercell grid");
    /* Supercell size in y */
    if (!sf_histfloat (scgrid, "o2", &esc_scgrid->oscz)) sf_error ("No Dscz= in supercell grid");
    /* Supercell start in z  */
    if (!sf_histfloat (scgrid, "o3", &esc_scgrid->oscx)) sf_error ("No Dscx= in supercell grid");
    /* Supercell start in x */
    if (!sf_histfloat (scgrid, "o4", &esc_scgrid->oscy)) sf_error ("No Dscy= in supercell grid");
    /* Supercell start in y */

    esc_scgrid->nsz = (size_t)esc_scgrid->nscz;
    esc_scgrid->nsxz = (size_t)esc_scgrid->nscx*(size_t)esc_scgrid->nscz;
    esc_scgrid->ns = (size_t)esc_scgrid->nscy*(size_t)esc_scgrid->nscx*(size_t)esc_scgrid->nscz;

    esc_scgrid->cells = (sf_esc_scblock3*)sf_alloc (esc_scgrid->ns,
                                                    sizeof(sf_esc_scblock3));
    simtab = sf_getpars ();

    /* Loop over supercells */
    for (iscy = 0; iscy < esc_scgrid->nscy; iscy++) {
        for (iscx = 0; iscx < esc_scgrid->nscx; iscx++) {
            for (iscz = 0; iscz < esc_scgrid->nscz; iscz++) {
                i = iscy*esc_scgrid->nsxz + iscx*esc_scgrid->nsz + iscz;
                /* Name of the file with this supercell */
                sf_ucharread (sbuf, ESC3_BSC_BMAX, scgrid);                    
                /* Add a key for the filename - the only way to
                   read a different RSF file without a command line parameter */
                sf_simtab_enter (simtab, "scell", (const char*)sbuf);
                fcell = sf_input ("scell");
                esc_scgrid->cells[i] = sf_esc_scblock3_init (fcell, esc_tracer);
                sz += sf_bytes (fcell);
                sf_fileclose (fcell);
            }
        }
    }

    if (verb) {
        sf_warning ("Grid - %dx%dx%d, sampling - %gx%gx%g", esc_scgrid->nscz,
                    esc_scgrid->nscx, esc_scgrid->nscy, esc_scgrid->dscz,
                    esc_scgrid->dscx, esc_scgrid->dscy);
        sf_warning ("Total number of supercells - %d", esc_scgrid->ns);
        sf_warning ("Memory used for all supercells - %g Mb", sz*1e-6);
    }

    esc_scgrid->nr = 0; /* Total number of trajectories traced */
    esc_scgrid->nrt = 0; /* Number of ray-traced segments */
    esc_scgrid->drt = 0.0; /* Total length of ray-traced segments */

    return esc_scgrid;
}

void sf_esc_scgrid3_close (sf_esc_scgrid3 esc_scgrid, bool verb)
/*< Destroy object >*/
{
    size_t i;

    if (verb) {
        sf_warning ("Number of ray-traced segments per trajectory - %g",
                    (double)esc_scgrid->nrt/(double)esc_scgrid->nr);
        sf_warning ("Average length of a ray-traced segment - %g",
                    esc_scgrid->nrt ? esc_scgrid->drt/(double)esc_scgrid->nrt : 0.0);
    }
    for (i = 0; i < esc_scgrid->ns; i++)
        sf_esc_scblock3_close (esc_scgrid->cells[i]);
    free (esc_scgrid->cells);

    free (esc_scgrid);
}

void sf_esc_scgrid3_compute (sf_esc_scgrid3 esc_scgrid,
                             float *z, float *x, float *y,
                             float *t, float *l, float *b, float *a)
/*< Compute escape values for a point with subsurface coordinates (z, x, y, b, a)
    by stitching local escape solutions in supercells of a phase-space grid >*/
{
    size_t i;
    int ic = 0, iscz, iscx, iscy;
    EscColor3 side = ESC3_SIDE_NUM;

    /* Initial supercell */
    iscz = floorf ((*z - esc_scgrid->oscz)/esc_scgrid->dscz);
    iscx = floorf ((*x - esc_scgrid->oscx)/esc_scgrid->dscx);
    iscy = floorf ((*y - esc_scgrid->oscy)/esc_scgrid->dscy);

    while (iscz >= 0 && iscz < esc_scgrid->nscz &&
           iscx >= 0 && iscx < esc_scgrid->nscx &&
           iscy >= 0 && iscy < esc_scgrid->nscy) {
        /* Supercell index in the grid */
        i = iscy*esc_scgrid->nsxz + iscx*esc_scgrid->nsz + iscz;
        /* Project trajectory across the supercell */
        side = sf_esc_scblock3_project_point (esc_scgrid->cells[i],
                                              z, x, y, t, l, b, a,
                                              ic != 0 ? &esc_scgrid->nrt : NULL,
                                              ic != 0 ? &esc_scgrid->drt : NULL);
        /* Move to the neighbor supercell */
        if (side & ESC3_TOP)
            iscz--;
        if (side & ESC3_BOTTOM)
            iscz++;
        if (side & ESC3_LEFT)
            iscx--;
        if (side & ESC3_RIGHT)
            iscx++;
        if (side & ESC3_NEAR)
            iscy--;
        if (side & ESC3_FAR)
            iscy++;
        if (ESC3_INSIDE == side)
            sf_error ("sf_esc_scgrid3_compute: point projection error");
        ic++;
    }
    esc_scgrid->nr++;
}


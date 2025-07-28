/* Angle gather at fixed lateral position in 3-D */
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

#ifndef _cram_gather3_h

#include "cram_point3.h"
/*^*/

typedef struct CRAMGather3 *sf_cram_gather3;
/* abstract data type */
/*^*/

#endif

struct CRAMGather3 {
    int    nb, na, nz, noa, nda, oda;
    float  b0, a0, db, da;
    float *stack;
    float *energy;
    float *illum;
    bool   outaz, inorm;
};
/* concrete data type */

sf_cram_gather3 sf_cram_gather3_init (int nb, int na, int nz,
                                      float b0, float db, float a0, float da,
                                      float oamax, float damax, bool outaz, bool inorm)
/*< Initialize object >*/
{
    sf_cram_gather3 cram_gather = (sf_cram_gather3)sf_alloc (1, sizeof (struct CRAMGather3));

    cram_gather->nb = nb; /* Inclination */
    cram_gather->b0 = b0;
    cram_gather->db = db;
    cram_gather->na = na; /* Azimuth */
    cram_gather->a0 = a0;
    cram_gather->da = da;
    cram_gather->nz = nz; /* Number of depth samples */

    cram_gather->outaz = outaz; /* Whether to output azimuth dimension */
    cram_gather->inorm = inorm; /* Whether to normalize gathers */

    /* Maximum number of scattering angles in the output */
    cram_gather->noa = (int)(oamax/db + 0.5) + 1;
    if (cram_gather->noa > nb)
        cram_gather->noa = nb;
    /* Maximum number of dip angles in the output */
    cram_gather->nda = (int)(2*damax/(0.5*db) + 0.5) + 1;
    if (cram_gather->nda > 4*nb)
        cram_gather->nda = 4*nb;
    /* Index of the first dip angle in the output */
    cram_gather->oda = (4*nb - cram_gather->nda)/2;
    if (cram_gather->oda < 0)
        cram_gather->oda = 0;

    /* Buffers for collapsing azimuth dimension */
    cram_gather->stack = sf_floatalloc (4*nb);
    cram_gather->energy = sf_floatalloc (4*nb);
    cram_gather->illum = sf_floatalloc (4*nb);

    return cram_gather;
}

void sf_cram_gather3_close (sf_cram_gather3 cram_gather)
/*< Destroy object >*/
{
    free (cram_gather->stack);
    free (cram_gather->energy);
    free (cram_gather->illum);
    free (cram_gather);
}

void sf_cram_gather3_image_output (sf_cram_gather3 cram_gather, sf_cram_point3 cram_point,
                                   sf_file imfile, sf_file ilfile)
/*< Write one spatial location from the image to disk >*/
{
    float smp, hits;

    smp = sf_cram_point3_get_image (cram_point, &hits);
    sf_floatwrite (&smp, 1, imfile);
    if (ilfile)
        sf_floatwrite (&hits, 1, ilfile);
}

void sf_cram_gather3_oangle_output (sf_cram_gather3 cram_gather, sf_cram_point3 cram_point,
                                    sf_file oafile, sf_file sqfile, sf_file ilfile)
/*< Add one spatial point to the scattering angle gather and write it to disk >*/
{
    int ib, ia;
    float **oimage, **osqimg, **ohits;

    oimage = sf_cram_point3_get_oimage (cram_point, &osqimg, &ohits);

    if (cram_gather->inorm) {
        for (ia = 0; ia < cram_gather->na/4; ia++) {
            for (ib = 0; ib < 2*cram_gather->nb; ib++) {
                if (ohits[ia][ib] >= 1.0) {
                    oimage[ia][ib] /= ohits[ia][ib];
                    osqimg[ia][ib] /= ohits[ia][ib];
                }
            }
        }
    }

    if (cram_gather->outaz) {
        for (ia = 0; ia < cram_gather->na/4; ia++) { /* [0; PI) */
            for (ib = 0; ib < cram_gather->noa; ib++) {
                /* Image */
                if (oafile)
                    sf_floatwrite (&oimage[ia][cram_gather->nb - ib], 1, oafile);
                /* Illumination map */
                if (ilfile)
                    sf_floatwrite (&ohits[ia][cram_gather->nb - ib], 1, ilfile);
                /* Energy image */ 
                if (sqfile) 
                    sf_floatwrite (&osqimg[ia][cram_gather->nb - ib], 1, sqfile);
            }
        }
        for (ia = 0; ia < cram_gather->na/4; ia++) { /* [PI; 2PI) */
            /* Image */
            if (oafile)
                sf_floatwrite (&oimage[ia][cram_gather->nb], cram_gather->noa, oafile);
            /* Illumination map */
            if (ilfile)
                sf_floatwrite (&ohits[ia][cram_gather->nb], cram_gather->noa, ilfile);
            /* Energy image */ 
            if (sqfile) 
                sf_floatwrite (&osqimg[ia][cram_gather->nb], cram_gather->noa, sqfile);
        }
    } else {
        memset (cram_gather->stack, 0, sizeof(float)*4*cram_gather->nb);
        memset (cram_gather->energy, 0, sizeof(float)*4*cram_gather->nb);
        memset (cram_gather->illum, 0, sizeof(float)*4*cram_gather->nb);
        for (ia = 0; ia < cram_gather->na/4; ia++) {
            for (ib = 0; ib < cram_gather->noa; ib++) {
                cram_gather->stack[ib] += oimage[ia][cram_gather->nb - ib];
                cram_gather->energy[ib] += osqimg[ia][cram_gather->nb - ib];
                cram_gather->illum[ib] += ohits[ia][cram_gather->nb - ib];
            }
        }
        for (ia = 0; ia < cram_gather->na/4; ia++) {
            for (ib = 1; ib < cram_gather->noa; ib++) {
                cram_gather->stack[ib] += oimage[ia][cram_gather->nb + ib];
                cram_gather->energy[ib] += osqimg[ia][cram_gather->nb + ib];
                cram_gather->illum[ib] += ohits[ia][cram_gather->nb + ib];
            }
        }
        /* Image */
        if (oafile)
            sf_floatwrite (cram_gather->stack, cram_gather->noa, oafile);
        /* Illumination map */
        if (ilfile)
            sf_floatwrite (cram_gather->illum, cram_gather->noa, ilfile);
        /* Energy image */ 
        if (sqfile) 
            sf_floatwrite (cram_gather->energy, cram_gather->noa, sqfile);
    }
}

void sf_cram_gather3_dangle_output (sf_cram_gather3 cram_gather, sf_cram_point3 cram_point,
                                    sf_file dafile, sf_file sqfile, sf_file ilfile)
/*< Add one spatial point to the dip angle gather and write it to disk >*/
{
    int ib, ia;
    float **dimage, **dsqimg, **dhits;

    dimage = sf_cram_point3_get_dimage (cram_point, &dsqimg, &dhits);

    if (cram_gather->inorm) {
        for (ia = 0; ia < cram_gather->na/4; ia++) {
            for (ib = 0; ib < 4*cram_gather->nb; ib++) {
                if (dhits[ia][ib] >= 1.0) {
                    dimage[ia][ib] /= dhits[ia][ib];
                    dsqimg[ia][ib] /= dhits[ia][ib];
                }
            }
        }
    }

    if (cram_gather->outaz) {
        for (ia = 0; ia < cram_gather->na/4; ia++) {
            /* Image */
            if (dafile)
                sf_floatwrite (&dimage[ia][cram_gather->oda], cram_gather->nda, dafile);
            /* Illumination map */
            if (ilfile)
                sf_floatwrite (&dhits[ia][cram_gather->oda], cram_gather->nda, ilfile);
            /* Energy image */ 
            if (sqfile) 
                sf_floatwrite (&dsqimg[ia][cram_gather->oda], cram_gather->nda, sqfile);
        }
    } else {
        memset (cram_gather->stack, 0, sizeof(float)*4*cram_gather->nb);
        memset (cram_gather->energy, 0, sizeof(float)*4*cram_gather->nb);
        memset (cram_gather->illum, 0, sizeof(float)*4*cram_gather->nb);
        for (ia = 0; ia < cram_gather->na/4; ia++) {
            for (ib = 0; ib < cram_gather->nda; ib++) {
                cram_gather->stack[ib] += dimage[ia][cram_gather->oda + ib];
                cram_gather->energy[ib] += dsqimg[ia][cram_gather->oda + ib];
                cram_gather->illum[ib] += dhits[ia][cram_gather->oda + ib];
            }
        }
        /* Image */
        if (dafile)
            sf_floatwrite (cram_gather->stack, cram_gather->nda, dafile);
        /* Illumination map */
        if (ilfile)
            sf_floatwrite (cram_gather->illum, cram_gather->nda, ilfile);
        /* Energy image */ 
        if (sqfile) 
            sf_floatwrite (cram_gather->energy, cram_gather->nda, sqfile);
    }
}

void sf_cram_gather3_setup_image_output (sf_cram_gather3 cram_gather, sf_file input,
                                          sf_file imfile)
/*< Setup dimensions in a file for the image >*/
{
    sf_oaxa (imfile, sf_iaxa (input, 4), 1);
    sf_oaxa (imfile, sf_iaxa (input, 5), 2);
    sf_oaxa (imfile, sf_iaxa (input, 6), 3);
    sf_putint (imfile, "n4", 1);
    sf_putint (imfile, "n5", 1);
    sf_putint (imfile, "n6", 1);
}

void sf_cram_gather3_setup_oangle_output (sf_cram_gather3 cram_gather, sf_file input,
                                          sf_file oafile)
/*< Setup dimensions in a file for scattering angle gathers >*/
{
    if (cram_gather->outaz)
        sf_unshiftdim (input, oafile, 0);
    else
        sf_unshiftdim2 (input, oafile, 0);
    sf_putint (oafile, "n1", cram_gather->noa);
    sf_putfloat (oafile, "d1", cram_gather->db);
    sf_putfloat (oafile, "o1", 0.0);
    sf_putstring (oafile, "label1", "Scattering angle");
    if (cram_gather->outaz) {
        sf_putint (oafile, "n2", cram_gather->na/2);
        sf_putfloat (oafile, "d2", 2.0*cram_gather->da);
        sf_putfloat (oafile, "o2", 0.0);
        sf_putstring (oafile, "label2", "Scattering angle azimuth");
    }
}

void sf_cram_gather3_setup_dangle_output (sf_cram_gather3 cram_gather, sf_file input,
                                          sf_file dafile)
/*< Setup dimensions in a file for dip angle gathers >*/
{
    if (cram_gather->outaz)
        sf_unshiftdim (input, dafile, 0);
    else
        sf_unshiftdim2 (input, dafile, 0);
    sf_putint (dafile, "n1", cram_gather->nda);
    sf_putfloat (dafile, "d1", 0.5*cram_gather->db);
    sf_putfloat (dafile, "o1", -180.0 + cram_gather->oda*0.5*cram_gather->db);
    sf_putstring (dafile, "label1", "Dip angle");
    if (cram_gather->outaz) {
        sf_putint (dafile, "n2", cram_gather->na/4);
        sf_putfloat (dafile, "d2", 2.0*cram_gather->da);
        sf_putfloat (dafile, "o2", 0.0);
        sf_putstring (dafile, "label2", "Dip angle azimuth");
    }
}


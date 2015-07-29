/* Angle gather at fixed lateral position in 2-D */
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

#ifndef _cram_gather2_h

#include "cram_point2.h"
/*^*/

typedef struct CRAMGather2 *sf_cram_gather2;
/* abstract data type */
/*^*/

#endif

struct CRAMGather2 {
    int     na, nz;
    int     noa, nda, oda;
    float   a0, da;
    bool    sqsmb;
    float **dstack;
    float **denergy;
    int   **dillum;
    float **rstack;
    float **renergy;
    int   **rillum;
    float  *smb;
};
/* concrete data type */

sf_cram_gather2 sf_cram_gather2_init (int na, int nz, float a0, float da,
                                      float oamax, float damax)
/*< Initialize object >*/
{
    sf_cram_gather2 cram_gather = (sf_cram_gather2)sf_alloc (1, sizeof (struct CRAMGather2));

    cram_gather->na = na; /* Number of angles */
    cram_gather->a0 = a0;
    cram_gather->da = da;
    cram_gather->nz = nz; /* Number of depth samples */

    cram_gather->sqsmb = false; /* Whether to output semblance or energy */

    cram_gather->dstack = sf_floatalloc2 (2*na, nz); /* Diff. stack image */
    cram_gather->denergy = sf_floatalloc2 (2*na, nz); /* Diff. energy image */
    cram_gather->dillum = sf_intalloc2 (2*na, nz); /* Diff. illumination map */
    cram_gather->rstack = sf_floatalloc2 (na/2, nz); /* Refl. stack image */
    cram_gather->renergy = sf_floatalloc2 (na/2, nz); /* Refl. energy image */
    cram_gather->rillum = sf_intalloc2 (na/2, nz); /* Refl. illumination map */
    cram_gather->smb = sf_floatalloc (2*na); /* Semblance */

    /* Maximum number of scattering angles in the output */
    cram_gather->noa = (int)(oamax/da + 0.5) + 1;
    if (cram_gather->noa > na/2)
        cram_gather->noa = na/2;
    /* Maximum number of dip angles in the output */
    cram_gather->nda = (int)(2*damax/(0.5*da) + 0.5) + 1;
    if (cram_gather->nda > 2*na)
        cram_gather->nda = 2*na;
    /* Index of the first dip angle in the output */
    cram_gather->oda = (2*na - cram_gather->nda)/2;
    if (cram_gather->oda < 0)
        cram_gather->oda = 0;

    return cram_gather;
}

void sf_cram_gather2_close (sf_cram_gather2 cram_gather)
/*< Destroy object >*/
{
    free (cram_gather->dstack[0]); free (cram_gather->dstack);
    free (cram_gather->denergy[0]); free (cram_gather->denergy);
    free (cram_gather->dillum[0]); free (cram_gather->dillum);
    free (cram_gather->rstack[0]); free (cram_gather->rstack);
    free (cram_gather->renergy[0]); free (cram_gather->renergy);
    free (cram_gather->rillum[0]); free (cram_gather->rillum);
    free (cram_gather->smb);
    free (cram_gather);
}

void sf_cram_gather2_set_sqsmb (sf_cram_gather2 cram_gather, bool sqsmb)
/*< Whether to output semblance or energy >*/
{
    cram_gather->sqsmb = sqsmb;
}

void sf_cram_gather2_set_point (sf_cram_gather2 cram_gather, int iz,
                                sf_cram_point2 cram_point)
/*< Add computed z location to the gathers >*/
{
    int ia, ja, h;
    float smp, amax;
    float **image, **sqimg;
    int **hits;

    if (0 == iz) {
        memset (cram_gather->dstack[0], 0, sizeof(float)*cram_gather->nz*2*cram_gather->na);
        memset (cram_gather->denergy[0], 0, sizeof(float)*cram_gather->nz*2*cram_gather->na);
        memset (cram_gather->dillum[0], 0, sizeof(int)*cram_gather->nz*2*cram_gather->na);
        memset (cram_gather->rstack[0], 0, sizeof(float)*cram_gather->nz*cram_gather->na/2);
        memset (cram_gather->renergy[0], 0, sizeof(float)*cram_gather->nz*cram_gather->na/2);
        memset (cram_gather->rillum[0], 0, sizeof(int)*cram_gather->nz*cram_gather->na/2);
    }

    /* Get pointers to the image */
    image = sf_cram_point2_get_image (cram_point, &sqimg, &hits, &amax);

    if (0.0 == amax)
        return; /* Empty image */

    for (ja = 0; ja < 2*cram_gather->na; ja++) {
        for (ia = 0; ia < cram_gather->na/2; ia++) {
            h = hits[ja][ia];
            if (0 == h)
                continue;
            smp = image[ja][ia];
            /* Normalize for the fold in the scattering angle direction */
            smp /= (float)(cram_gather->na/2);
            /* Store directional gather information */
            cram_gather->dstack[iz][ja] += smp;
            cram_gather->denergy[iz][ja] += smp*smp;
            cram_gather->dillum[iz][ja]++;
            /* There are 4 times as many dip angles as scattering angles,
               so renormalize the sample accordingly */
            smp *= 4.0;
            /* Store reflection gather information */
            cram_gather->rstack[iz][ia] += smp;
            cram_gather->renergy[iz][ia] += smp*smp;
            cram_gather->rillum[iz][ia]++;
        }
    }
}

void sf_cram_gather2_full_output (sf_cram_gather2 cram_gather, sf_cram_point2 cram_point,
                                  sf_file ofile, sf_file ifile)
/*< Compute scattering angle gather from the full gather and write it to disk >*/
{
    float **image, **sqimg, amax;
    int **hits;
    /* Get pointers to the image */
    image = sf_cram_point2_get_image (cram_point, &sqimg, &hits, &amax);

    /* Image */
    if (ofile)
        sf_floatwrite (image[0], cram_gather->na*cram_gather->na, ofile);
    /* Illumination map */
    if (ifile)
        sf_intwrite (hits[0], cram_gather->na*cram_gather->na, ifile);
}

void sf_cram_gather2_oangle_output (sf_cram_gather2 cram_gather, sf_file oafile,
                                    sf_file smfile, sf_file ilfile)
/*< Compute scattering angle gather from the full gather and write it to disk >*/
{
    const int zwidth = 5, awidth = 3;
    int ia, iz, iz0, iz1, h, ia0, ia1;
    float ave, rms;

    for (iz = 0; iz < cram_gather->nz; iz++) {
        /* Gather */
        if (oafile)
            sf_floatwrite (cram_gather->rstack[iz],
                           cram_gather->noa, oafile);
        /* Illumination map */
        if (ilfile)
            sf_intwrite (cram_gather->rillum[iz],
                         cram_gather->noa, ilfile);

        if (NULL == smfile)
            continue;

        if (cram_gather->sqsmb) /* Output squared traces instead of semblance */ 
            sf_floatwrite (cram_gather->renergy[iz],
                           cram_gather->noa, smfile);
    }
    if (NULL == smfile || cram_gather->sqsmb)
        return;

    /* Compute semblance gather */
    for (iz = 0; iz < cram_gather->nz; iz++) {
        for (ia = 0; ia < cram_gather->na/2; ia++) {
            iz0 = iz - zwidth;
            if (iz0 < 0)
                iz0 = 0;
            iz1 = iz + zwidth + 1;
            if (iz1 > cram_gather->nz)
                iz1 = cram_gather->nz;
            ave = rms = 0.0;
            for (; iz0 != iz1; iz0++) {
                ia0 = ia - awidth;
                if (ia0 < 0)
                    ia0 = 0;
                ia1 = ia + awidth + 1;
                if (ia1 > (cram_gather->na/2))
                    ia1 = cram_gather->na/2;
                for (; ia0 != ia1; ia0++) {
                    ave += cram_gather->rstack[iz0][ia0]/*cram_gather->rstack[iz0][ia0]*/;
                    rms += cram_gather->renergy[iz0][ia0];
                }
            }
            h = cram_gather->rillum[iz][ia];
            cram_gather->smb[ia] = (rms != 0.0 && h) ? ave/rms : 0.0;
        }
        /* Semblance map */
        sf_floatwrite (cram_gather->smb,
                       cram_gather->noa, smfile);
    }
}

void sf_cram_gather2_dangle_output (sf_cram_gather2 cram_gather, sf_file dafile,
                                    sf_file smfile, sf_file ilfile)
/*< Compute dip angle gather from the full gather and write it to disk >*/
{
    const int zwidth = 5, awidth = 1;
    int ia, iz, iz0, iz1, h, ia0, ia1;
    float ave, rms;

    for (iz = 0; iz < cram_gather->nz; iz++) {
        /* Gather */
        if (dafile)
            sf_floatwrite (&cram_gather->dstack[iz][cram_gather->oda],
                           cram_gather->nda, dafile);
        /* Illumination map */
        if (ilfile)
            sf_intwrite (&cram_gather->dillum[iz][cram_gather->oda],
                         cram_gather->nda, ilfile);

        if (NULL == smfile)
            continue;

        if (cram_gather->sqsmb) /* Output squared traces instead of semblance */ 
            sf_floatwrite (&cram_gather->denergy[iz][cram_gather->oda],
                           cram_gather->nda, smfile);
    }
    if (NULL == smfile || cram_gather->sqsmb)
        return;

    /* Compute semblance gather */
    for (iz = 0; iz < cram_gather->nz; iz++) {
        for (ia = 0; ia < 2*cram_gather->na; ia++) {
            iz0 = iz - zwidth;
            if (iz0 < 0)
                iz0 = 0;
            iz1 = iz + zwidth + 1;
            if (iz1 > cram_gather->nz)
                iz1 = cram_gather->nz;
            ave = rms = 0.0;
            for (; iz0 != iz1; iz0++) {
                ia0 = ia - awidth;
                if (ia0 < 0)
                    ia0 = 0;
                ia1 = ia + awidth + 1;
                if (ia1 > (2*cram_gather->na))
                    ia1 = 2*cram_gather->na;
                for (; ia0 != ia1; ia0++) {
                    ave += cram_gather->dstack[iz0][ia0]/*cram_gather->dstack[iz0][ia0]*/;
                    rms += cram_gather->denergy[iz0][ia0];
                }
            }
            h = cram_gather->dillum[iz][ia];
            cram_gather->smb[ia] = (rms != 0.0 && h) ? ave*ave/rms : 0.0;
        }
        /* Semblance map */
        sf_floatwrite (&cram_gather->smb[cram_gather->oda],
                       cram_gather->nda, smfile);
    }
}

void sf_cram_gather2_setup_full_output (sf_cram_gather2 cram_gather, sf_file ofile, bool imap)
/*< Setup dimensions in a file for the full image output >*/
{
    if (imap)
        sf_settype (ofile, SF_INT);
    sf_putint (ofile, "n1", cram_gather->na/2);
    sf_putfloat (ofile, "d1", cram_gather->da);
    sf_putfloat (ofile, "o1", 0.0);
    sf_putstring (ofile, "label1", "Scattering angle");
    sf_putint (ofile, "n2", 2*cram_gather->na);
    sf_putfloat (ofile, "d2", 0.5*cram_gather->da);
    sf_putfloat (ofile, "o2", cram_gather->a0 - 0.5*cram_gather->da);
    sf_putstring (ofile, "label2", "Dip angle");
}

void sf_cram_gather2_setup_oangle_output (sf_cram_gather2 cram_gather, sf_file input,
                                          sf_file oafile, bool imap)
/*< Setup dimensions in a file for scattering angle gathers >*/
{
    sf_unshiftdim (input, oafile, 0);
    if (imap)
        sf_settype (oafile, SF_INT);
    sf_putint (oafile, "n1", cram_gather->noa);
    sf_putfloat (oafile, "d1", cram_gather->da);
    sf_putfloat (oafile, "o1", 0.0);
    sf_putstring (oafile, "label1", "Scattering angle");
}

void sf_cram_gather2_setup_dangle_output (sf_cram_gather2 cram_gather, sf_file input,
                                          sf_file dafile, bool imap)
/*< Setup dimensions in a file for dip angle gathers >*/
{
    sf_unshiftdim (input, dafile, 0);
    if (imap)
        sf_settype (dafile, SF_INT);
    sf_putint (dafile, "n1", cram_gather->nda);
    sf_putfloat (dafile, "d1", 0.5*cram_gather->da);
    sf_putfloat (dafile, "o1", cram_gather->a0 - 0.5*cram_gather->da +
                               cram_gather->oda*0.5*cram_gather->da);
    sf_putstring (dafile, "label1", "Dip angle");
}


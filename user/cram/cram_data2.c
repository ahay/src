/* Common-shot data for angle migration (both 2-D and 3-D) */
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

#include <errno.h>
#include <sys/mman.h>

#include <rsf.h>

#ifndef _cram_data2_h

typedef struct CRAMData2 *sf_cram_data2;
/* abstract data type */
/*^*/

#endif

struct CRAMData2 {
    int            nt;
    size_t         n, offs;
    float          t0, dt;
    bool           kmah, filter, erefl;
    float         *data;
    unsigned char *mmaped;
};
/* concrete data type */

sf_cram_data2 sf_cram_data2_init (sf_file data)
/*< Initialize object >*/
{
    FILE *stream;
    size_t n;

    sf_cram_data2 cram_data = (sf_cram_data2)sf_alloc (1, sizeof (struct CRAMData2));

    /* Filter data according to slope */
    if (!sf_histbool (data, "filter", &cram_data->filter)) sf_error ("No filter= in data");
    /* Use KMAH phase shifts */
    if (!sf_histbool (data, "KMAH", &cram_data->kmah)) sf_error ("No KMAH= in data");
    /* Use exploding reflector assumption */
    if (!sf_histbool (data, "ExplRefl", &cram_data->erefl)) cram_data->erefl = false;

    if (!sf_histint (data, "n1", &cram_data->nt)) sf_error ("No n1= in data");
    if (!sf_histfloat (data, "d1", &cram_data->dt)) sf_error ("No d1= in data");
    if (!sf_histfloat (data, "o1", &cram_data->t0)) sf_error ("No o1= in data");

    cram_data->n = sf_leftsize (data, 1);

    /* Total size of the file memory mapping */
    n = (size_t)cram_data->nt*(size_t)cram_data->n*(size_t)sizeof(float);

    if (cram_data->kmah)
        sf_warning ("Using data prepared for KMAH shifts");
    if (cram_data->filter)
        sf_warning ("Using data prepared for anti-aliasing filter");
    if (cram_data->erefl)
        sf_warning ("Assuming data modeled by exploding reflector");
    sf_warning ("Total data size: %g Mb", 1e-6*(float)n);

    stream = sf_filestream (data);
    cram_data->offs = ftello (stream);
    cram_data->mmaped = (unsigned char*)mmap (NULL, n + cram_data->offs,
                                              PROT_READ, MAP_SHARED,
                                              fileno (stream), 0);
    if (cram_data->mmaped == MAP_FAILED)
        sf_error ("Data mmap failed: %s", strerror (errno));
    cram_data->data = (float*)(cram_data->mmaped + cram_data->offs);

    return cram_data;
}

void sf_cram_data2_close (sf_cram_data2 cram_data)
/*< Destroy object >*/
{
    size_t n = (size_t)cram_data->nt*(size_t)cram_data->n*(size_t)sizeof(float);

    munmap (cram_data->mmaped, n + cram_data->offs);
    free (cram_data);
}

bool sf_cram_data2_get_erefl (sf_cram_data2 cram_data)
/*< Return true if data is from exploding reflector modeling >*/
{
    return cram_data->erefl;
}

/* Return one sample from propely filtered trace according to requested slope */
static float sf_cram_data2_get_sample_from_trace (sf_cram_data2 cram_data, float *trace,
                                                  float trfact, float t, int *hits) {
    int it, k, ik1, ik2;
    float tf;

    tf = (t - cram_data->t0)/cram_data->dt;
    it = floorf (tf);
    tf -= (float)it;
    if (it < 0 || it >= (cram_data->nt - 1))
        return 0.0;

    if (cram_data->filter) {
        /* Filter length after Lumley-Claerbout-Bevc */
        k = 4.0*fabs (trfact) - 1.0;
        if (k < 1) k = 1;
        k = (k - 1)/2*2 + 1;
        ik1 = it - k/2 - 1;
        ik2 = it + k/2 + 1;
        if (ik1 < 0 || ik2 >= (cram_data->nt - 1))
            return 0.0;

        (*hits)++;
        return 1.0/((float)(k/2 + 1))*1.0/((float)(k/2 + 1))*(
                   2.0*(trace[it]*(1.0 - tf)
                        + trace[it + 1]*tf)
                  -1.0*(trace[ik1]*(1.0 - tf)
                        + trace[ik1 + 1]*tf)
                  -1.0*(trace[ik2]*(1.0 - tf)
                        + trace[ik2 + 1]*tf));
    } else {
        (*hits)++;
        return (trace[it]*(1.0 - tf)
                + trace[it + 1]*tf);
    }
}

float sf_cram_data2_get_sample (sf_cram_data2 cram_data, size_t i, float t,
                                float p, float dx, int kmah, int *hits)
/*< Get one sample for trace i, time t, slope p, and spatial sampling dx >*/
{
    size_t off;
    float trf, w = 1.0;
    float *trace = NULL;

    if (cram_data->erefl) {
        t *= 0.5;
        kmah /= 2;
    }

    if (cram_data->kmah)
        i *= (size_t)2;

    if (i >= cram_data->n)
        return 0.0;

    /* Choose largest trace factor for AA filter length */
    trf = fabs (p*dx/cram_data->dt);

    /* Choose trace according to KMAH index */
    if (cram_data->kmah) {
        kmah = kmah % 4; /* Enforce phase-shift periodicity */
    } else {
        kmah = 0;
    }
    off = (size_t)i*(size_t)cram_data->nt;

    switch (kmah) {
        case 0:
            trace = &cram_data->data[off];
            break;
        case 1: /* PI/2 */
            trace = &cram_data->data[off + (size_t)cram_data->nt];
            w = -1.0;
            break;
        case 2: /* PI */
            trace = &cram_data->data[off];
            w = -1.0;
            break;
        case 3: /* 3PI/4 */
            trace = &cram_data->data[off + (size_t)cram_data->nt];
            break;
    }

    return w*sf_cram_data2_get_sample_from_trace (cram_data, trace,
                                                  trf, t, hits);
}


/* Survey information (source and receiver locations) for 3-D angle migration */
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

#include <errno.h>
#include <sys/mman.h>

#include <rsf.h>

#ifndef _cram_survey3_h

typedef struct CRAMSurvey3 *sf_cram_survey3;
/* abstract data type */
/*^*/

#endif

struct CRAMSurvey3 {
    int            nst;
    size_t         nht;
    int           *sh;
    size_t        *ihtr;
    int           *hh;
    size_t         ns, nh;
    float          ds, dh;
    size_t        *ih;
    float         *s, *h, *mm;
    bool           tri;
    size_t         sz, offs;
    unsigned char *mmaped;
};
/* concrete data type */

sf_cram_survey3 sf_cram_survey3_init (sf_file survey, bool verb)
/*< Initialize object >*/
{
    FILE *stream;
    size_t offs;
    float dsx, dsy, dhx, dhy;
    sf_cram_survey3 cram_survey = (sf_cram_survey3)sf_alloc (1, sizeof (struct CRAMSurvey3));

    if (sf_gettype (survey) != SF_UCHAR)
        sf_error ("Wrong data type in the survey info file");

    if (!sf_histlargeint (survey, "n1", (off_t*)&cram_survey->sz)) sf_error ("No n1= in survey");
    if (!sf_histbool (survey, "Istri", &cram_survey->tri)) sf_error ("No Istri= in survey");
    if (!sf_histlargeint (survey, "Nshot", (off_t*)&cram_survey->ns)) sf_error ("No Nshot= in survey");
    if (!sf_histlargeint (survey, "Nrec", (off_t*)&cram_survey->nh)) sf_error ("No Nrec= in survey");
    if (cram_survey->tri) {
        if (!sf_histint (survey, "Nstri", &cram_survey->nst)) sf_error ("No Nstri= in survey");
        if (!sf_histlargeint (survey, "Nhtri", (off_t*)&cram_survey->nht)) sf_error ("No Nhtri= in survey");
    }

    stream = sf_filestream (survey);
    cram_survey->offs = ftello (stream);

    if (stream != stdin) {
        cram_survey->mmaped = (unsigned char*)mmap (NULL, cram_survey->offs + cram_survey->sz,
                                                    PROT_READ, MAP_SHARED,
                                                    fileno (stream), 0);
        if (cram_survey->mmaped == MAP_FAILED)
            sf_error ("Survey info file mmap failed: %s", strerror (errno));
        offs = cram_survey->offs;
        /* Array of number of receivers in each shot */
        cram_survey->ih = (size_t*)(cram_survey->mmaped + offs);
        offs += (size_t)cram_survey->ns*(size_t)sizeof(size_t);
        /* Shot coordinates */
        cram_survey->s = (float*)(cram_survey->mmaped + offs);
        offs += (size_t)cram_survey->ns*(size_t)(sizeof(float)*2);
        /* Min/max coordinates for every shot */
        cram_survey->mm = (float*)(cram_survey->mmaped + offs);
        offs += (size_t)cram_survey->ns*(size_t)(sizeof(float)*4);
        /* Receiver coordinates */
        cram_survey->h = (float*)(cram_survey->mmaped + offs);
        offs += (size_t)cram_survey->nh*(size_t)(sizeof(float)*2);
        if (cram_survey->tri) {
            /* Array of number of triangles in each receiver triangulation */
            cram_survey->ihtr = (size_t*)(cram_survey->mmaped + offs);
            offs += (size_t)cram_survey->ns*(size_t)sizeof(size_t);
            /* Shot triangulation */
            cram_survey->sh = (int*)(cram_survey->mmaped + offs);
            offs += (size_t)cram_survey->nst*(size_t)(sizeof(int)*3);
            /* Receiver triangulation */
            cram_survey->hh = (int*)(cram_survey->mmaped + offs);
        }
    } else {
        sf_warning ("Survey info file appears to be stdin");
        sf_warning ("mmap is not possible for the survey info");

        cram_survey->s = sf_floatalloc ((size_t)cram_survey->ns*(size_t)2);
        cram_survey->ih = (size_t*)sf_alloc ((size_t)cram_survey->ns, sizeof(size_t));
        cram_survey->mm = sf_floatalloc ((size_t)cram_survey->ns*(size_t)4);
        cram_survey->h = sf_floatalloc ((size_t)cram_survey->nh*(size_t)2);
        cram_survey->ihtr = (size_t*)sf_alloc ((size_t)cram_survey->ns, sizeof(size_t));
        cram_survey->sh = sf_intalloc ((size_t)cram_survey->nst*(size_t)3);
        cram_survey->hh = sf_intalloc ((size_t)cram_survey->nht*(size_t)3);
    
        sf_ucharread ((unsigned char*)cram_survey->ih,
                      cram_survey->ns*(size_t)sizeof(size_t), survey);
        sf_ucharread ((unsigned char*)cram_survey->s,
                      cram_survey->ns*(size_t)(sizeof(float)*2), survey);
        sf_ucharread ((unsigned char*)cram_survey->s,
                      cram_survey->ns*(size_t)(sizeof(float)*4), survey);
        sf_ucharread ((unsigned char*)cram_survey->h,
                      cram_survey->nh*(size_t)(sizeof(float)*2), survey);
        if (cram_survey->tri) {
            sf_ucharread ((unsigned char*)cram_survey->ihtr,
                          cram_survey->ns*(size_t)sizeof(size_t), survey);
            sf_ucharread ((unsigned char*)cram_survey->sh,
                          cram_survey->nst*(size_t)(sizeof(int)*3), survey);
            sf_ucharread ((unsigned char*)cram_survey->hh,
                          cram_survey->nht*(size_t)(sizeof(int)*3), survey);
        }
        cram_survey->mmaped = NULL;
    }

    if (verb) {
        sf_warning ("Using %d shots, %d receivers", cram_survey->ns, cram_survey->nh);
        if (cram_survey->tri)
            sf_warning ("Using %d shot triangles, %lu receiver triangles",
                        cram_survey->nst, cram_survey->nht);
    }

    dsx = dsy = dhx = dhy = 0.0125;
    /* Determine source and receiver samplings from loaded coordinates */
    if (cram_survey->ns > 1) {
        dsx = fabsf (cram_survey->s[2] - cram_survey->s[0]);
        dsy = fabsf (cram_survey->s[3] - cram_survey->s[1]);
    }
    if (cram_survey->nh > 1) {
        dhx = fabsf (cram_survey->h[2] - cram_survey->h[0]);
        dhy = fabsf (cram_survey->h[3] - cram_survey->h[1]);
    }
    if (1 == cram_survey->ns) {
        dsx = dhx;
        dsy = dhy;
    }
    cram_survey->ds = hypotf (dsx, dsy);
    cram_survey->dh = hypotf (dhx, dhy);

    if (verb)
        sf_warning ("Using %g shot sampling, %g receiver sampling",
                    cram_survey->ds, cram_survey->dh);

    return cram_survey;
}

void sf_cram_survey3_close (sf_cram_survey3 cram_survey)
/*< Destroy object >*/
{
    if (cram_survey->mmaped)
        munmap (cram_survey->mmaped, cram_survey->offs + cram_survey->sz);
    else {
        if (cram_survey->tri) {
            free (cram_survey->hh);
            free (cram_survey->sh);
            free (cram_survey->ihtr);
        }
        free (cram_survey->h);
        free (cram_survey->s);
        free (cram_survey->ih);
    }
    free (cram_survey);
}

float sf_cram_survey3_get_src_sampling (sf_cram_survey3 cram_survey)
/*< Return source spatial sampling >*/
{
    return cram_survey->ds;
}

float sf_cram_survey3_get_rcv_sampling (sf_cram_survey3 cram_survey)
/*< Return receiver spatial sampling >*/
{
    return cram_survey->dh;
}

int sf_cram_survey3_get_ns (sf_cram_survey3 cram_survey)
/*< Return number of sources >*/
{
    return cram_survey->ns;
}

int sf_cram_survey3_get_nstr (sf_cram_survey3 cram_survey)
/*< Return number of source triangles >*/
{
    return cram_survey->nst;
}

size_t sf_cram_survey3_get_first_source (sf_cram_survey3 cram_survey, float *sx, float *sy,
                                         size_t *nh, float *gxmin, float *gxmax,
                                         float *gymin, float *gymax)
/*< Initialize iterator through sources, return source index,
    number of receivers in this shot, number of sources >*/
{
    if (1 != cram_survey->ns)
        *nh = cram_survey->ih[1] - cram_survey->ih[0];
    else
        *nh = cram_survey->nh - cram_survey->ih[0];

    *sx = cram_survey->s[0];
    *sy = cram_survey->s[1];
    if (gxmin)
        *gxmin = cram_survey->mm[0];
    if (gxmax)
        *gxmax = cram_survey->mm[1];
    if (gymin)
        *gymin = cram_survey->mm[2];
    if (gymax)
        *gymax = cram_survey->mm[3];

    return 0;;
}

size_t sf_cram_survey3_get_next_source (sf_cram_survey3 cram_survey, size_t is,
                                        float *sx, float *sy, size_t *nh,
                                        float *gxmin, float *gxmax, float *gymin, float *gymax)
/*< Iterator through sources, return source index, source coordinates,
    number of receivers in this source >*/
{
    is += (size_t)1;

    if (is >= cram_survey->ns) {
        return (size_t)-1;
        *sx = SF_HUGE;
        *sy = SF_HUGE;
    }

    if (is != (cram_survey->ns - (size_t)1))
        *nh = cram_survey->ih[is + (size_t)1] - cram_survey->ih[is];
    else
        *nh = cram_survey->nh - cram_survey->ih[is];

    *sx = cram_survey->s[is*(size_t)2];
    *sy = cram_survey->s[is*(size_t)2 + (size_t)1];
    if (gxmin)
        *gxmin = cram_survey->mm[is*(size_t)4];
    if (gxmax)
        *gxmax = cram_survey->mm[is*(size_t)4 + (size_t)1];
    if (gymin)
        *gymin = cram_survey->mm[is*(size_t)4 + (size_t)2];
    if (gymax)
        *gymax = cram_survey->mm[is*(size_t)4 + (size_t)3];

    return is;
}

size_t sf_cram_survey3_get_first_receiver (sf_cram_survey3 cram_survey, size_t is, 
                                           size_t *i, float *gx, float *gy)
/*< Initialize iterator through receivers for a given source,
    return receiver index, number of receivers, trace index, and receiver coordinates >*/
{
    if (is >= cram_survey->ns) {
        *gx = SF_HUGE;
        *gy = SF_HUGE;
        return (size_t)-1;
    }

    *i = cram_survey->ih[is];

    *gx = cram_survey->h[(*i)*(size_t)2];
    *gy = cram_survey->h[(*i)*(size_t)2 + (size_t)1];

    return 0;
}

size_t sf_cram_survey3_get_next_receiver (sf_cram_survey3 cram_survey, size_t is,
                                          size_t ih, size_t *i, float *gx, float *gy)
/*< Iterator through receivers for a given source,
    return receiver index, trace index, and receiver coordinates >*/
{
    size_t nh;

    if (is >= cram_survey->ns) {
        *gx = SF_HUGE;
        *gy = SF_HUGE;
        return (size_t)-1;
    }

    if (is != (cram_survey->ns - (size_t)1))
        nh = cram_survey->ih[is + (size_t)1] - cram_survey->ih[is];
    else
        nh = cram_survey->nh - cram_survey->ih[is];

    ih += (size_t)1;

    if (ih >= nh) {
        *gx = SF_HUGE;
        *gy = SF_HUGE;
        return (size_t)-1;
    }

    *i = cram_survey->ih[is] + ih;

    *gx = cram_survey->h[(*i)*(size_t)2];
    *gy = cram_survey->h[(*i)*(size_t)2 + (size_t)1];

    return ih;
}

size_t sf_cram_survey3_get_first_trirec (sf_cram_survey3 cram_survey, size_t is, 
                                         size_t *nhtr, size_t *i, float *gx, float *gy)
/*< Initialize iterator through receiver triangles for a given source,
    return triangle index, number of triangles, receiver coordinates,
    and their corresponding trace indices >*/
{
    int j;
    size_t ii;
    if (is >= cram_survey->ns) {
        for (j = 0; j < 3; j++) {
            gx[j] = SF_HUGE;
            gy[j] = SF_HUGE;
        }
        return (size_t)-1;
    }

    /* Number of triangles */
    if (is != (cram_survey->ns - (size_t)1))
        *nhtr = cram_survey->ihtr[is + (size_t)1] - cram_survey->ihtr[is];
    else
        *nhtr = cram_survey->nht - cram_survey->ihtr[is];

    /* Extract three vertices and their coordinates */
    for (j = 0; j < 3; j++) {
        ii = cram_survey->ih[is] + (size_t)cram_survey->hh[cram_survey->ihtr[is]*(size_t)3 +
                                                           (size_t)j];  
        i[j] = ii;
        gx[j] = cram_survey->h[ii*(size_t)2];
        gy[j] = cram_survey->h[ii*(size_t)2 + (size_t)1];
    }

    return 0;
}

size_t sf_cram_survey3_get_next_trirec (sf_cram_survey3 cram_survey, size_t is,
                                        size_t ihtr, size_t *i, float *gx, float *gy)
/*< Iterator through receivers triangles for a given source,
    return triangle index, receiver coordinates, and their
    corresponding trace indices >*/
{
    int j;
    size_t ii, nhtr;

    if (is >= cram_survey->ns) {
        for (j = 0; j < 3; j++) {
            gx[j] = SF_HUGE;
            gy[j] = SF_HUGE;
        }
        return (size_t)-1;
    }

    /* Number of triangles */
    if (is != (cram_survey->ns - (size_t)1))
        nhtr = cram_survey->ihtr[is + (size_t)1] - cram_survey->ihtr[is];
    else
        nhtr = cram_survey->nht - cram_survey->ihtr[is];

    ihtr += (size_t)1;

    if (ihtr >= nhtr) {
        for (j = 0; j < 3; j++) {
            gx[j] = SF_HUGE;
            gy[j] = SF_HUGE;
        }
        return (size_t)-1;
    }

    /* Extract three vertices and their coordinates */
    for (j = 0; j < 3; j++) {
        ii = cram_survey->ih[is] + (size_t)cram_survey->hh[cram_survey->ihtr[is]*(size_t)3 +
                                                           ihtr*(size_t)3 + (size_t)j];  
        i[j] = ii;
        gx[j] = cram_survey->h[ii*(size_t)2];
        gy[j] = cram_survey->h[ii*(size_t)2 + (size_t)1];
    }

    return ihtr;
}


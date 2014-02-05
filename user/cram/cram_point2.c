/* Angle migration of one depth point in 2-D */
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

#ifndef _cram_point2_h

#include "cram_data2.h"
/*^*/
#include "cram_survey2.h"
/*^*/
#include "cram_slowness2.h"
/*^*/
#include "cram_rbranch2.h"
/*^*/

typedef struct CRAMPoint2 *sf_cram_point2;
/* abstract data type */
/*^*/

#endif

struct CRAMPoint2 {
    int                na;
    int                oaz; /* Maximum scattering angle */
    int                daz; /* Maximum dip angle (abs.value of) */
    int                ts, th; /* Tapering for source and receiver */
    float              a0, da, dx;
    bool               mute, amp;
    float              amax, smax, hmax;
    float            **image;
    float            **sqimg;
    int              **hits;
    sf_cram_data2      data;
    sf_cram_survey2    survey;
    sf_cram_slowness2  slowness;
    sf_cram_rbranch2   rbranch;
};
/* concrete data type */

sf_cram_point2 sf_cram_point2_init (int na, float a0, float da,
                                    sf_cram_data2 data, sf_cram_survey2 survey,
                                    sf_cram_slowness2 slowness, sf_cram_rbranch2 rbranch)
/*< Initialize object >*/
{
    float ds, dh;
    sf_cram_point2 cram_point = (sf_cram_point2)sf_alloc (1, sizeof (struct CRAMPoint2));

    cram_point->na = na; /* Number of angles */
    cram_point->a0 = a0;
    cram_point->da = da;

    cram_point->ts = 3;
    cram_point->th = 5;
   
    cram_point->mute = false;
    cram_point->oaz = 2*na;
    cram_point->daz = 2*na;

    cram_point->amax = 0.0;
    cram_point->smax = SF_HUGE;
    cram_point->hmax = SF_HUGE;

    cram_point->data = data;
    cram_point->survey = survey;
    cram_point->slowness = slowness;
    cram_point->rbranch = rbranch;

    cram_point->image = sf_floatalloc2 (na/2, 2*na); /* Regular image */
    cram_point->sqimg = sf_floatalloc2 (na/2, 2*na); /* RMS image */
    cram_point->hits = sf_intalloc2 (na/2, 2*na); /* Hit counts */

    ds = sf_cram_survey2_get_src_sampling (survey);
    dh = sf_cram_survey2_get_rcv_sampling (survey);
    cram_point->dx = fabs (ds) > fabs (dh) ? ds : dh; 

    return cram_point;
}

void sf_cram_point2_close (sf_cram_point2 cram_point)
/*< Destroy object >*/
{
    free (cram_point->image[0]); free (cram_point->image);
    free (cram_point->sqimg[0]); free (cram_point->sqimg);
    free (cram_point->hits[0]); free (cram_point->hits);
    free (cram_point);
}

void sf_cram_point2_set_mute (sf_cram_point2 cram_point, float oaz, float daz)
/*< Set mute limits in constant z plane >*/
{
    cram_point->mute = true;
    /* Maximum scattering angle */
    cram_point->oaz = (int)(oaz/cram_point->da + 0.5);
    /* Maximum dip angle */
    cram_point->daz = 2*(int)(daz/cram_point->da + 0.5);
}

void sf_cram_point2_set_taper (sf_cram_point2 cram_point, int ts, int th)
/*< Set tapering length >*/
{
    cram_point->ts = ts;
    cram_point->th = th;
}

void sf_cram_point2_set_shmax (sf_cram_point2 cram_point, float smax, float hmax)
/*< Set maximum allowed width of shot and receiver ray branches >*/
{
    cram_point->smax = smax;
    cram_point->hmax = hmax;
}

float** sf_cram_point2_get_image (sf_cram_point2 cram_point, float ***sqimg, int ***hits,
                                  float *amax)
/*< Return pointers to internal arrays of image samples >*/
{
    if (sqimg)
        *sqimg = cram_point->sqimg;
    if (hits)
        *hits = cram_point->hits;
    if (amax)
        *amax = cram_point->amax;
    return cram_point->image;
}

/* Check if angle is outside mute triangle */
static bool sf_cram_point2_idx_is_mute (float x, float y, float x1, float y1,
                                        float x2, float y2, float x3, float y3) {
    float AB = (y - y1)*(x2 - x1) - (x - x1)*(y2 - y1);
    float CA = (y - y3)*(x1 - x3) - (x - x3)*(y1 - y3);
    float BC = (y - y2)*(x3 - x2) - (x - x2)*(y3 - y2);

    return !(AB*BC >= 0.0 && BC*CA >= 0.0);
}

void sf_cram_point2_compute (sf_cram_point2 cram_point,
                             float z, float x, float **esc)
/*< Compute image for one point (z,x) and a correpsonding set of escape values >*/
{
    size_t i;
    int ies, ier, is, ih, ida, ioa, ht, ns, nh, ia, ja;
    float s, r, smp, w, ss = 1.0, sr, tp;
    sf_cram_surface_exit2 *srays, *rrays;

    /* Clean up */
    memset (cram_point->image[0], 0, sizeof(float)*cram_point->na*cram_point->na);
    memset (cram_point->sqimg[0], 0, sizeof(float)*cram_point->na*cram_point->na);
    memset (cram_point->hits[0], 0, sizeof(int)*cram_point->na*cram_point->na);

    /* Reinit escape branches */
    sf_cram_rbranch2_set_escapes (cram_point->rbranch, esc);

    /* Slowness for the source and receiver phases */
    if (cram_point->slowness)
        ss = sf_cram_slowness2_get_value (cram_point->slowness, z, x);
    sr = ss;

    /* Loop over known sources */
    s = sf_cram_survey2_get_first_source (cram_point->survey, &is, &ns);
    while (s != SF_HUGE) {
        srays = sf_cram_rbranch2_src_exits (cram_point->rbranch, s);
        ies = 0;
        /* Loop over found source ray branches */
        while (srays && srays[ies].ia >= 0) {
            if (srays[ies].d > cram_point->smax) {
                ies++;
                continue;
            }
            /* Loop over known receivers */
            r = sf_cram_survey2_get_first_receiver (cram_point->survey, is, &ih, &nh, &i);
            while (r != SF_HUGE) {
                rrays = sf_cram_rbranch2_rcv_exits (cram_point->rbranch, r);
                ier = 0;
                /* Loop over found receiver ray branches */
                while (rrays && rrays[ier].ia >= 0) {
                    if (rrays[ier].d > cram_point->smax) {
                        ier++;
                        continue;
                    }
                    ida = (srays[ies].ia + rrays[ier].ia); /* Dip angle index */
                    ioa = abs (srays[ies].ia - rrays[ier].ia); /* Scattering angle index */
                    if (cram_point->na/2 == ioa) {
                        ier++;
                        continue; /* Skip transmission */
                    }
                    if (ioa > cram_point->na/2) { /* Limit scattering angle to 180 */
                        ida += (ida < cram_point->na) /* Adjust dip */
                                ? cram_point->na
                                : -cram_point->na;
                        ioa = cram_point->na - ioa;
                    }
                    if (cram_point->mute &&
                        (abs (ida - cram_point->na) > cram_point->daz ||
                         ioa > cram_point->oaz ||
                         sf_cram_point2_idx_is_mute (ida - cram_point->na, ioa,
                                                     -cram_point->daz, 0,
                                                     cram_point->daz, 0,
                                                     0, cram_point->oaz))) {
                        ier++;
                        continue; /* Skip angles outside of limits */
                    }
                    ht = 0;
                    smp =  sf_cram_data2_get_sample (cram_point->data, i,
                                                     srays[ies].t + rrays[ier].t,
                                                     fabsf (srays[ies].p) + fabsf (rrays[ier].p)
                                                     /*srays[ies].p > rrays[ier].p ?
                                                     srays[ies].p : rrays[ier].p*/, cram_point->dx,
                                                     srays[ies].kmah + rrays[ier].kmah, &ht);
                    if (0 == ht) { /* No hit */
                        ier++;
                        continue;
                    }
                    /* Tapering */
                    tp = 1.0;
                    if (cram_point->ts) { /* Source direction */
                        if ((cram_point->ts - is) >= 0)
                            tp *= (float)(is + 1)/(float)cram_point->ts;
                        if ((ns - is) <= cram_point->ts)
                            tp *= (float)(ns - is)/(float)cram_point->ts;
                    }
                    if (cram_point->th) { /* Receiver direction */
                        if ((cram_point->th - ih) >= 0)
                            tp *= (float)(ih + 1)/(float)cram_point->th;
                        if ((nh - ih) <= cram_point->th)
                            tp *= (float)(nh - ih)/(float)cram_point->th;
                    }
                    /* Obliquity factor */
                    w = cos (SF_PI*ioa/(float)cram_point->na);
                    /* Multiply by geometric spreading */
                    smp *= tp*w*sqrt (ss*sr*srays[ies].j*rrays[ier].j);
                    cram_point->image[ida][ioa] += smp;
                    cram_point->sqimg[ida][ioa] += smp*smp;
                    cram_point->hits[ida][ioa]++;
                    ier++;
                }
                r = sf_cram_survey2_get_next_receiver (cram_point->survey, is, &ih, &i);
            } /* Loop over known receivers */
            ies++;
        }
        s = sf_cram_survey2_get_next_source (cram_point->survey, &is);
    } /* Loop over known sources */

    cram_point->amax = 0.0;
    /* Normalize for hits and find maximum absolute value */
    for (ja = 0; ja < 2*cram_point->na; ja++) {
        for (ia = 0; ia < cram_point->na/2; ia++) {
            ih = cram_point->hits[ja][ia];
            if (0 == ih)
                continue;
            cram_point->image[ja][ia] /= (float)ih;
            cram_point->sqimg[ja][ia] /= (float)ih;
            smp = fabs (cram_point->image[ja][ia]);
            if (smp > cram_point->amax)
                cram_point->amax = smp;
        }
    }
}


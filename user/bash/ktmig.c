/* Prestack time migration (2-D/3-D) kernel. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include <math.h>
#include <rsf.h>
/*^*/

#include "ktmig.h"

/****************************************************************
 *
 * Differentiation (backward second order)
 *
 ****************************************************************/

/*
  trace[nt] - input/output vector of differentiated trace
*/
void sf_ktmig_sbdiff (float *trace, const unsigned int nt, const float dt)
/*< Differentiation (backward second order) >*/
{
    int i;
    float val0, val1, val2;

    val2 = val1 = trace[0];
    for (i = 0; i < nt; i++) {
        val0 = trace[i];
        trace[i] = 0.5f*(3.0f*val0 - 4.0f*val1 + val2)/dt;
        val2 = val1;
        val1 = val0;
    }
}

/****************************************************************
 *
 * Causal integration
 *
 ****************************************************************/

/*
  trace[nt] - input/output vector of causally integrated trace
*/
void sf_ktmig_cint (float *trace, const unsigned int nt)
/*< Causal integration >*/
{
    int i;

    for (i = 1; i < nt; i++)
        trace[i] += trace[i - 1];
}

/****************************************************************
 *
 * Anti-causal integration
 *
 ****************************************************************/

/*
  trace[nt] - input/output vector of anti-causally integrated trace
*/
void sf_ktmig_acint (float *trace, const unsigned int nt)
/*< Anti-causal integration >*/
{
    int i;

    for (i = nt - 2; i >= 0; i--)
        trace[i] += trace[i + 1];
}

#define INTSMP(t, i) ((1.f - i + (float)((int)i))*t[(int)i] + (i - (float)((int)i))*t[(int)(i) + 1])

/****************************************************************
 *
 * PSTM kernel with anti-aliasing after Lumley-Claerbout
 *
 ****************************************************************/

/*
  trace[int] - input trace
  vrms[ont]  - RMS velocity vector at image location,
  image[ont] - output image vector at the same location,
  ox         - image coordinate in x,
  oy         - image coordinate in y,
  sx         - source coordinate in x,
  sy         - source coordinate in y,
  gx         - receiver coordinate in x,
  gy         - receiver coordinate in y,
  nt         - number of input samples,
  ont        - number of output samples,
  ot         - value of t0 in input,
  dt         - sample rate of input,
  oot        - value of t0 in output,
  odt        - sample rate of output,
  trm        - maximum AA filter length,
  trf        - trace factor for AA filter
  aa         - antialiasing (true/false)
*/
void sf_ktmig_kernel (float *trace, float *vrms, float *image,
                      const float ox, const float oy,
                      const float sx, const float sy,
                      const float gx, const float gy,
                      const unsigned int nt, const unsigned int ont,
                      const float ot, const float dt,
                      const float oot, const float odt,
                      const unsigned int trm, const float trf,
                      const bool aa)
/*< PSTM kernel with anti-aliasing after Lumley-Claerbout >*/
{
    float v, inv;
    float inv2trf, nf;
    int k;
    float j, scale, smp, so2, go2;
    float depth2, dx, dy, ts, tg;

    /* Loop over tau indices */
    for (k = 0; k < ont; k++) {
        /* RMS velocity at image location */
        v = vrms[k];
        /* Slowness at image location */
        inv = 1.0f/v;
        inv2trf = trf*inv*inv;
        depth2 = powf (0.5f*v*(oot + k*odt), 2.0f);
        /* squared distance to source from the image point on the surface */
        so2 = (sx - ox)*(sx - ox) + (sy - oy)*(sy - oy);
        /* squared distance to source from the image point on the surface */
        go2 = (gx - ox)*(gx - ox) + (gy - oy)*(gy - oy);
        /* Time from source to image point in pseudodepth */
        ts = sqrtf (so2 + depth2)*inv;
        /* Time from receiver to image point in pseudodepth */
        tg = sqrtf (go2 + depth2)*inv;
        /* double root square time = time to source + time to receiver */
        j = (ts + tg - ot)/dt; /* Input sample index */
        if (!aa) {
            if (j >= 0.f && j < nt)
                image[k] += INTSMP (trace, j);
            continue;
        }
        /* (distance to source.x)/(time to source) + (distance to receiver.x)/(time to receiver) */
        dx = (sx - ox)/ts + (gx - ox)/tg;
        /* (distance to source.y)/(time to source) + (distance to receiver.y)/(time to receiver) */
        dy = (sy - oy)/ts + (gy - oy)/tg;
        /* Filter length */
        nf = inv2trf*sqrtf (dx*dx + dy*dy);
        /* Truncate filter */
        if (nf > trm)
            nf = (float)trm;
        /* Check ranges */
        if ((j - nf - 1.0f) >= 0.0f && (j + nf + 1.0f) < nt) {
            /* Scaling factor */
            scale = 1.0f/(1.0f + nf);
            scale *= scale;
            /* Collect samples */
            smp = 2.0f*INTSMP (trace, j)
                  -INTSMP (trace, (j - nf - 1.0f))
                  -INTSMP (trace, (j + nf + 1.0f));
            /* Contribute to the image point */
            image[k] += scale*smp;
        }
    }
}


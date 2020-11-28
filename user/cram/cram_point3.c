/* Angle migration of one depth point in 3-D */
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

#ifndef _cram_point3_h

#include "cram_data2.h"
/*^*/
#include "cram_survey3.h"
/*^*/
#include "cram_slowness3.h"
/*^*/
#include "cram_rbranch3.h"
/*^*/

typedef struct CRAMPoint3 *sf_cram_point3;
/* abstract data type */
/*^*/

#endif

struct CRAMPoint3 {
    int                     nb, na;
    float                   oam; /* Maximum scattering angle */
    float                   dam; /* Maximum dip angle (abs.value of) */
    float                   img, hits;
    float                   b0, a0, db, da, ds, dh;
    float                   dym, dxm, atw;
    int                     mioz, midz;
    bool                    amp, mute, taper, extrap;
    bool                    agath, dipgath, erefl;
    float                 **oimage;
    float                 **osqimg;
    float                 **ohits;
    float                 **dimage;
    float                 **dsqimg;
    float                 **dhits;
    sf_cram_surface_exit3  *src_exits;
    sf_cram_surface_exit3  *rcv_exits;
    sf_cram_data2           data;
    sf_cram_survey3         survey;
    sf_cram_slowness3       slowness;
    sf_cram_rbranch3        rbranch;
};
/* concrete data type */

void sf_cram_point3_reset (sf_cram_point3 cram_point)
/*< Zero out internal image >*/
{
    if (cram_point->agath) {
        memset (cram_point->oimage[0], 0, sizeof(float)*cram_point->nb*cram_point->na/2);
        memset (cram_point->osqimg[0], 0, sizeof(float)*cram_point->nb*cram_point->na/2);
        memset (cram_point->ohits[0], 0, sizeof(float)*cram_point->nb*cram_point->na/2);
    }
    if (cram_point->dipgath) {
        memset (cram_point->dimage[0], 0, sizeof(float)*cram_point->nb*cram_point->na);
        memset (cram_point->dsqimg[0], 0, sizeof(float)*cram_point->nb*cram_point->na);
        memset (cram_point->dhits[0], 0, sizeof(float)*cram_point->nb*cram_point->na);
    }
    cram_point->hits = 0.0;
    cram_point->img = 0.0;
}

sf_cram_point3 sf_cram_point3_init (int nb, float b0, float db,
                                    int na, float a0, float da,
                                    bool agath, bool dipgath,
                                    sf_cram_data2 data, sf_cram_survey3 survey,
                                    sf_cram_slowness3 slowness, sf_cram_rbranch3 rbranch)
/*< Initialize object >*/
{
    sf_cram_point3 cram_point = (sf_cram_point3)sf_alloc (1, sizeof (struct CRAMPoint3));

    cram_point->na = na; /* Azimuth */
    cram_point->a0 = a0*SF_PI/180.0;
    cram_point->da = da*SF_PI/180.0;
    cram_point->nb = nb; /* Inclination */
    cram_point->b0 = b0*SF_PI/180.0;
    cram_point->db = db*SF_PI/180.0;

    cram_point->mute = true;

    cram_point->mute = false;
    cram_point->oam = SF_PI;
    cram_point->dam = SF_PI;
    cram_point->atw = 0.8;
    cram_point->mioz = nb;
    cram_point->midz = 2*nb;

    cram_point->taper = false;
    cram_point->dxm = 0.0;
    cram_point->dym = 0.0;

    cram_point->extrap = false;

    cram_point->data = data;
    cram_point->survey = survey;
    cram_point->slowness = slowness;
    cram_point->rbranch = rbranch;

    cram_point->erefl = sf_cram_data2_get_erefl (data);

    if (agath) {
        cram_point->oimage = sf_floatalloc2 (2*nb, na/4); /* Regular image, scattering angle */
        cram_point->osqimg = sf_floatalloc2 (2*nb, na/4); /* RMS image, scattering angle */
        cram_point->ohits = sf_floatalloc2 (2*nb, na/4); /* Hit counts, scattering angle */
    } else {
        cram_point->oimage = NULL;
        cram_point->osqimg = NULL;
        cram_point->ohits = NULL;
    }
    if (dipgath) {
        cram_point->dimage = sf_floatalloc2 (4*nb, na/4); /* Regular image, dip angle */
        cram_point->dsqimg = sf_floatalloc2 (4*nb, na/4); /* RMS image, dip angle */
        cram_point->dhits = sf_floatalloc2 (4*nb, na/4); /* Hit counts, dip angle */
    } else {
        cram_point->dimage = NULL;
        cram_point->dsqimg = NULL;
        cram_point->dhits = NULL;
    }

    cram_point->agath = agath; /* Compute scattering angle gather */
    cram_point->dipgath = dipgath; /* Compute dip angle gather */

    sf_cram_point3_reset (cram_point);

    cram_point->src_exits = (sf_cram_surface_exit3*)sf_alloc (na, sizeof (sf_cram_surface_exit3));
    cram_point->rcv_exits = (sf_cram_surface_exit3*)sf_alloc (na, sizeof (sf_cram_surface_exit3));

    cram_point->dh = sf_cram_survey3_get_rcv_sampling (survey);
    cram_point->ds = sf_cram_survey3_get_src_sampling (survey);

    return cram_point;
}

void sf_cram_point3_close (sf_cram_point3 cram_point)
/*< Destroy object >*/
{
    free (cram_point->src_exits);
    free (cram_point->rcv_exits);
    if (cram_point->agath) {
        free (cram_point->oimage[0]); free (cram_point->oimage);
        free (cram_point->osqimg[0]); free (cram_point->osqimg);
        free (cram_point->ohits[0]); free (cram_point->ohits);
    }
    if (cram_point->dipgath) {
        free (cram_point->dimage[0]); free (cram_point->dimage);
        free (cram_point->dsqimg[0]); free (cram_point->dsqimg);
        free (cram_point->dhits[0]); free (cram_point->dhits);
    }
    free (cram_point);
}

void sf_cram_point3_set_extrap (sf_cram_point3 cram_point, bool extrap)
/*< Set extrapolation flag for gathers >*/
{
    cram_point->extrap = extrap;
}

/* Tapering zone in degrees near the muting boundary */
#define OAM_TP 5.0
#define DAM_TP 10.0

void sf_cram_point3_set_mute (sf_cram_point3 cram_point, float oam, float dam)
/*< Set mute limits in constant z plane >*/
{
    cram_point->mute = true;
    /* Maximum scattering angle + tapering zone */
    cram_point->oam = (oam + OAM_TP)*SF_PI/180.0;
    if (cram_point->oam > SF_PI)
        cram_point->oam = SF_PI;
    /* Maximum dip angle + tapering zone */
    cram_point->dam = (dam + DAM_TP)*SF_PI/180.0;
    if (cram_point->dam > SF_PI)
        cram_point->dam = SF_PI;
    /* Zone without tapering (as a fraction of the full unmuted zone) */
    cram_point->atw = 0.5*((oam*SF_PI/180.0)/cram_point->oam + 
                           (dam*SF_PI/180.0)/cram_point->dam);
    /* Inidces of the maximum extent of the output samples
       in the scattering and dip angle domains */
    cram_point->mioz = (int)((oam*SF_PI/180.0)/cram_point->db + 0.5);
    if (cram_point->mioz > cram_point->nb)
        cram_point->mioz = cram_point->nb;
    cram_point->midz = (int)((dam*SF_PI/180.0)/(0.5*cram_point->db) + 0.5);
    if (cram_point->midz > 2*cram_point->nb)
        cram_point->midz = 2*cram_point->nb;
}

void sf_cram_point3_set_taper (sf_cram_point3 cram_point, float dxm, float dym)
/*< Set taper limits in receiver coordinate space >*/
{
    cram_point->taper = true;
    /* Taper length in x direction */
    cram_point->dxm = dxm;
    /* Taper length in y direction */
    cram_point->dym = dym;
}

void sf_cram_point3_set_amp (sf_cram_point3 cram_point, bool amp)
/*< Set amplitude correction flag >*/
{
    cram_point->amp = amp;
}

float sf_cram_point3_get_image (sf_cram_point3 cram_point, float *hits)
/*< Return current image sample and its hit counters >*/
{
    if (hits)
        *hits = cram_point->hits;
    return cram_point->img;
}

float** sf_cram_point3_get_oimage (sf_cram_point3 cram_point, float ***osqimg, float ***ohits)
/*< Return pointers to internal arrays of image samples >*/
{
    if (osqimg)
        *osqimg = cram_point->osqimg;
    if (ohits)
        *ohits = cram_point->ohits;
    return cram_point->oimage;
}

float** sf_cram_point3_get_dimage (sf_cram_point3 cram_point, float ***dsqimg, float ***dhits)
/*< Return pointers to internal arrays of image samples >*/
{
    if (dsqimg)
        *dsqimg = cram_point->dsqimg;
    if (dhits)
        *dhits = cram_point->dhits;
    return cram_point->dimage;
}

/* Compute scattering angle, scattering angle azimuth, dip angle, dip angle azimuth
   from source and receiver elevation/azimuth angles */
static void sf_cram_point3_angles (float sb, float sa, float hb, float ha,
                                   float *oa, float *oz, float *da, float *dz) {
    float sv[3], hv[3], ov[3], dv[3], pv[3], nv[3], mv[3];
    float a, l, dl, dll;

    /* Create shot and receiver vectors */
    sv[0] = cosf (sb);          /* z */
    sv[1] = sinf (sb)*cosf (sa); /* x */
    sv[2] = sinf (sb)*sinf (sa); /* y */
    hv[0] = cosf (hb);          /* z */
    hv[1] = sinf (hb)*cosf (ha); /* x */
    hv[2] = sinf (hb)*sinf (ha); /* y */

    /* Scattering angle vector */
    ov[0] = -sv[0] + hv[0];
    ov[1] = -sv[1] + hv[1];
    ov[2] = -sv[2] + hv[2];
    /* Dip vector */
    dv[0] = sv[0] + hv[0];
    dv[1] = sv[1] + hv[1];
    dv[2] = sv[2] + hv[2];
    dl = sqrtf (dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);

    if (fabsf (sb - hb) > 1e-4 ||
        fabsf (sa - ha) > 1e-4) {
        /* Scattering angle - dot product between vectors */
        a = sv[0]*hv[0] + sv[1]*hv[1] + sv[2]*hv[2];
        if (a > 1.0)
            a = 1.0;
        if (a < -1.0)
            a = -1.0;
        *oa = acosf (a);
        /* Cross-product between the two vectors - contained in the reflection plane */
        pv[0] = ov[1]*dv[2] - ov[2]*dv[1];
        pv[1] = ov[2]*dv[0] - ov[0]*dv[2];
        pv[2] = ov[0]*dv[1] - ov[1]*dv[0];
        /* Cross-product between the x direction (reference vector) and the dip -
           contained in the reflection plane */
        nv[0] = dv[2];  /* 1.0*dv[2] - 0.0*dv[1] */
        nv[1] = 0.0;    /* 0.0*dv[0] - 0.0*dv[2] */
        nv[2] = -dv[0]; /* 0.0*dv[1] - 1.0*dv[0] */
        /* Scattering angle azimuth - angle between the above two vectors in the
           reflection plane - see Fomel and Sava (2005) */
        dll = sqrtf (pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2])*
              sqrtf (nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
        a = (pv[0]*nv[0] + pv[1]*nv[1] + pv[2]*nv[2])/dll;
        if (a > 1.0)
            a = 1.0;
        if (a < -1.0)
            a = -1.0;
        *oz = acosf (a);
        /* Cross product between the above two vectors - should be parallel
           to the dip */
        mv[0] = nv[1]*pv[2] - nv[2]*pv[1];
        mv[1] = nv[2]*pv[0] - nv[0]*pv[2];
        mv[2] = nv[0]*pv[1] - nv[1]*pv[0];
        dll = dl*sqrtf (mv[0]*mv[0] + mv[1]*mv[1] + mv[2]*mv[2]);
        /* Check if it has the same direction with the dip */
        a = (dv[0]*mv[0] + dv[1]*mv[1] + dv[2]*mv[2])/dll;
        if (a < 0.0) /* Otherwise, flip the scattering azimuth */
            *oz = SF_PI - *oz;
        else
            *oa *= -1.0;
    } else {
        /* Scattering angle */
        *oa = 0.0;
        /* Scattering angle azimuth */
        *oz = 0.0;
    }

    /* Dip angle - angle between z axis and dv */
    a = dv[0]/dl;
    *da = acosf (a);
    /* Dip angle azimuth - angle between projection
       of dv onto x-y plane and x component of dv */
    l = hypotf (dv[1], dv[2]);
    if (l != 0.0) {
        a = dv[1]/l;
        if (a > 1.0)
            a = 1.0;
        if (a < -1.0)
            a = -1.0;
        *dz = acosf (a);
        if (dv[2] < 0.0)
            *dz = SF_PI - *dz;
        else
            *da *= -1.0;
    } else {
        *dz = 0.0;
    }
}

static int sf_cram_point3_angle_idx (float a, float a0, float da, float *fa) {
    int ia;
    float af;

    af = (a - a0)/da;
    ia = floorf (af);
    af -= (float)ia;
    if (fa)
        *fa = af;
    return ia;
}

/* Check if scattering and dip angles are in mute zone */
static bool sf_cram_point3_gangles_are_mute (sf_cram_point3 cram_point, float oa, float da) {
    /* Mute zone is comprised by an elliptic area in the
       scattering-dip anlge plane with oam and dam beign the
       two major semiaxes */
    return (oa*oa/(cram_point->oam*cram_point->oam) +
            da*da/(cram_point->dam*cram_point->dam)) > 1.0;
}

/* Return taper value for the current combination of dip and scattering angles;
   it varies from 1 to 0 as the point approaches edges of the mute zone */
static float sf_cram_point3_gataper (sf_cram_point3 cram_point, float oa, float da) {
    float d = oa*oa/(cram_point->oam*cram_point->oam) +
              da*da/(cram_point->dam*cram_point->dam);
    if (d < cram_point->atw)
        return 1.0;
    else if (d > 1.0)
        return 0.0;
    else
        return (1.0 - d)/(1.0 - cram_point->atw);
}

/* Compute how much azimuth and angle change from [ac, azc] to [a, az] pair,
   add difference in *dda, *ddz */
static void sf_cram_point3_aaz_spread (float ac, float azc, float a, float az,
                                       float *da, float *dz) {
    float rot[3]; 
    rot[0] = fabsf (azc - az);
    rot[1] = fabsf (azc - (az + (float) SF_PI));
    rot[2] = fabsf (azc - (az - (float) SF_PI));
    /* Find smallest rotation between the two azimuths */
    if (rot[0] <= rot[1] && rot[0] <= rot[2]) {
        *dz += rot[0];
        *da += fabsf (ac - a);
    } else if (rot[1] <= rot[0] && rot[1] <= rot[2]) {
        *dz += rot[1];
        *da += fabsf (ac + a);
    } else {
        *dz += rot[2];
        *da += fabsf (ac + a);
    }
}

static void sf_cram_point3_fill_abins (sf_cram_point3 cram_point, float smp, 
                                       float ss, float sr, int ies, int ier,
                                       float oac, float ozc, float dac, float dzc,
                                       bool zoffset) {
    float dtw = DAM_TP*SF_PI/180.0,
          otw = OAM_TP*SF_PI/180.0; /* Tapering width in degrees */
    int dnw = 2*(int)(dtw/cram_point->db + 0.5), onw = (int)(otw/cram_point->db + 0.5); /* Taper width in samples */
    int isxy, jsxy, ihxy, jhxy, ia, iz, ida, ioa, ioz, idz;
    int ioaf, ioal, iozf, iozl, idaf, idal, idzf, idzl;
    float doa = 0.0, doz = 0.0, dda = 0.0, ddz = 0.0;
    float oa, da, oz, dz;
    float ds, dh, sa, sb, ha, hb, dw, tw = 1.0;
    bool flip = false;

    /* Integral weight */
    if (cram_point->amp) {
        /* See Bleistein et al (2003), equation (21) for details */
        if (cram_point->erefl) {
            dw = sqrtf (ss);
            dw /= sqrtf (fabsf (cram_point->src_exits[ies].j));
            dw *= cram_point->src_exits[ies].cs;
            sb = sinf (cram_point->b0 + cram_point->src_exits[ies].ib*cram_point->db);
            dw *= sqrtf (fabsf (sb));
            smp *= dw;
        } else {
            dw = sqrtf (ss*sr);
            dw *= cosf (0.5*oac);
            dw /= sqrtf (fabsf (cram_point->src_exits[ies].j*cram_point->rcv_exits[ier].j));
            dw *= cram_point->src_exits[ies].cs*cram_point->rcv_exits[ier].cs;
            sb = sinf (cram_point->b0 + cram_point->src_exits[ies].ib*cram_point->db);
            hb = zoffset ? sb : sinf (cram_point->b0 + cram_point->rcv_exits[ier].ib*cram_point->db);
            dw *= sqrtf (fabsf (sb*hb));
            smp *= dw;
        }
    }

    /* Compute tapeting weights at the edges of the mute zone */
    dw = cram_point->oam - fabsf (oac);
    if (dw < 2.0*otw) {
        if (dw < otw)
            tw = 0.0;
        else
            tw *= (dw - otw)/otw;
    }
    dw = cram_point->dam - fabsf (dac);
    if (dw < 2.0*dtw) {
        if (dw < dtw)
            tw = 0.0;
        else
            tw *= (dw - dtw)/dtw;
    }

    /* Contribute sample to the stacked image */
    cram_point->img += tw*smp;
    cram_point->hits += 1.0;

    if (!cram_point->agath && !cram_point->dipgath)
        return;

    if (cram_point->extrap && !zoffset) {
        /* Compute change in inclination and azimuth angles for the source and receiver branches
           using db/dx,db/dy and da/dx,da/dy with respect to a displacement away from the
           source and receiver on the surface; then find dip and scattering angle deviation,
           which correspond to this displacement */
        for (isxy = 0; isxy < 2; isxy++) { /* d/dx, d/dy, source side */
            for (jsxy = 0; jsxy < 2; jsxy++) { /* +/- shift in x/y on the source side */
                ds = jsxy != 0 ? cram_point->ds : -cram_point->ds;
                sb = cram_point->b0 +
                     (cram_point->src_exits[ies].ib + cram_point->src_exits[ies].ibxy[isxy]*ds)*
                     cram_point->db;
                sa = cram_point->a0 +
                     (cram_point->src_exits[ies].ia + cram_point->src_exits[ies].iaxy[isxy]*ds)*
                     cram_point->da;
                for (ihxy = 0; ihxy < 2; ihxy++) { /* d/dx, d/dy, receiver side */
                    for (jhxy = 0; jhxy < 2; jhxy++) { /* +/- shift in x/y on the receiver side */
                        dh = jhxy != 0 ? cram_point->ds : -cram_point->ds;
                        hb = cram_point->b0 +
                             (cram_point->rcv_exits[ier].ib + cram_point->rcv_exits[ier].ibxy[ihxy]*dh)*
                             cram_point->db;
                        ha = cram_point->a0 +
                             (cram_point->rcv_exits[ier].ia + cram_point->rcv_exits[ier].iaxy[ihxy]*dh)*
                             cram_point->da;
                        /* New system of subsurface angles for the surface displacment */
                        sf_cram_point3_angles (sb, sa, hb, ha, &oa, &oz, &da, &dz);
                        /* Scattering angle spread with respect to the initial values */
                        sf_cram_point3_aaz_spread (oac, ozc, oa, oz, &doa, &doz);
                        /* Dip angle spread with respect to the initial values */
                        sf_cram_point3_aaz_spread (dac, dzc, da, dz, &dda, &ddz);
                    }
                }
            }
        }
        /* Average changes in dip and scattering angles and their azimuths */
        doa /= 16.0;
        doz /= 16.0;
        dda /= 16.0;
        ddz /= 16.0;
        /* Extrapolate azimuth further away */
        doa *= 0.75;
        /* Extrapolate dip only in the vicinity of a receiver */
        dw = cram_point->dh/cram_point->ds;
        dda *= 2.0*dw;
        ddz *= 2.0*dw;
    } else {
        doa = 0.5*cram_point->db;
        doz = cram_point->da;
        dda = 0.25*cram_point->db;
        ddz = cram_point->da;
    }

    /* Contribute sample to the scattering angle gather */
    if (cram_point->agath) {
        /* Scattering angle first index */
        ioaf = sf_cram_point3_angle_idx (oac - doa, -SF_PI, cram_point->db, NULL);
        if (ioaf < (cram_point->nb - cram_point->mioz))
            ioaf = cram_point->nb - cram_point->mioz;
        /* Scattering angle last index */
        ioal = sf_cram_point3_angle_idx (oac + doa, -SF_PI, cram_point->db, NULL);
        ioal += 1;
        if (ioal > (cram_point->nb + cram_point->mioz))
            ioal = cram_point->nb + cram_point->mioz;
        if (ioal >= 2*cram_point->nb)
            ioal = 2*cram_point->nb - 1;
        /* Scattering angle azimuth first index */
        iozf = sf_cram_point3_angle_idx (ozc - doz, 0.0, 2.0*cram_point->da, NULL);
        /* Scattering angle azimuth last index */
        iozl = sf_cram_point3_angle_idx (ozc + doz, 0.0, 2.0*cram_point->da, NULL);
        iozl += 1;
        dw = 1.0/*((ioal - ioaf + 1)*(iozl - iozf + 1))*/;
        for (ioz = iozf; ioz <= iozl; ioz++) { /* Scattering angle azimuths */
            iz = ioz;
            if (iz < 0) {
                while (iz < 0)
                    iz += cram_point->na/4;
                flip = true;
            } else if (iz >= cram_point->na/4) {
                while (iz >= cram_point->na/4)
                    iz -= cram_point->na/4;
                flip = true;
            }
            for (ioa = ioaf; ioa <= ioal; ioa++) { /* Scattering angles */
                tw = dw;
                ia = ioa;
                if (flip && ia != 0) {
                    ia = 2*cram_point->nb - ia;
                }
                /* Tapering at the edges */
                if ((cram_point->mioz - abs (ia - cram_point->nb)) < onw) {
                    tw *= (cram_point->mioz - abs (ia - cram_point->nb))/(float)onw;
                }
                cram_point->oimage[iz][ia] += tw*smp;
                cram_point->osqimg[iz][ia] += (tw*smp)*(tw*smp);
                cram_point->ohits[iz][ia] += 1.0;
            }
            flip = false;
        }
    }

    /* Contribute sample to the dip angle gather */
    if (cram_point->dipgath) {
        /* Dip angle first index */
        idaf = sf_cram_point3_angle_idx (dac - dda, -SF_PI, 0.5*cram_point->db, NULL);
        if (idaf < (2*cram_point->nb - cram_point->midz))
            idaf = 2*cram_point->nb - cram_point->midz;
        /* Dip angle last index */
        idal = sf_cram_point3_angle_idx (dac + dda, -SF_PI, 0.5*cram_point->db, NULL);
        idal += 1;
        if (idal > (2*cram_point->nb + cram_point->midz))
            idal = 2*cram_point->nb + cram_point->midz;
        if (idal >= 4*cram_point->nb)
            idal = 4*cram_point->nb - 1;
        /* Dip angle azimuth first index */
        idzf = sf_cram_point3_angle_idx (dzc - ddz, 0.0, 2.0*cram_point->da, NULL);
        /* Dip angle azimuth last index */
        idzl = sf_cram_point3_angle_idx (dzc + ddz, 0.0, 2.0*cram_point->da, NULL);
        idzl += 1;
        if (cram_point->amp)
            dw = 1.0/(sinf (0.5*cram_point->oam)*sinf (0.5*cram_point->oam));
        else
            dw = 1.0;
        for (idz = idzf; idz <= idzl; idz++) { /* Dip angle azimuths */
            iz = idz;
            if (iz < 0) {
                while (iz < 0)
                    iz += cram_point->na/4;
                flip = true;
            } else if (iz >= cram_point->na/4) {
                while (iz >= cram_point->na/4)
                    iz -= cram_point->na/4;
                flip = true;
            }
            for (ida = idaf; ida <= idal; ida++) { /* Dip angles */
                tw = dw;
                ia = ida;
                if (flip && ia != 0) {
                    ia = 4*cram_point->nb - ia;
                }
                /* Tapering at the edges */
                if ((cram_point->midz - abs (ia - 2*cram_point->nb)) < dnw) {
                    tw *= (cram_point->midz - abs (ia - 2*cram_point->nb))/(float)dnw;
                }
                cram_point->dimage[iz][ia] += tw*smp;
                cram_point->dsqimg[iz][ia] += (tw*smp)*(tw*smp);
                cram_point->dhits[iz][ia] += 1.0;
            }
            flip = false;
        }
    }
}

/* Process all source and receiver branches according to the exit locations:
   get data samples and fill dip and scattering angle bins in the image */
static void sf_cram_point3_process_branches (sf_cram_point3 cram_point, float gx, float gy,
                                             size_t i, int nes, int ner,
                                             float ss, float sr, float tw,
                                             bool zoffset) {
    int ies, ier, ht;
    float oa, da, oz, dz, aw;
    float smp, sa, sb, ha, hb;

    ies = 0;
    /* Loop over found source ray branches */
    while (ies < nes) {
        /* Shot ray angles */
        sb = cram_point->b0 + cram_point->src_exits[ies].ib*cram_point->db;
        sa = cram_point->a0 + cram_point->src_exits[ies].ia*cram_point->da;
        if (zoffset) {
            /* Prevent double summation in the zero offset case */
            ier = ies;
            ner = ies + 1;
        } else
            ier = 0;
        /* Loop over found receiver ray branches */
        while (ier < ner) {
            /* Receiver ray angles */
            if (zoffset) {
                hb = sb;
                ha = sa;
            } else {
                hb = cram_point->b0 + cram_point->rcv_exits[ier].ib*cram_point->db;
                ha = cram_point->a0 + cram_point->rcv_exits[ier].ia*cram_point->da;
            }
            /* Convert to dip and scattering angles */
            sf_cram_point3_angles (sb, sa, hb, ha, &oa, &oz, &da, &dz);
            if (sf_cram_point3_gangles_are_mute (cram_point, oa, da)) {
                ier++;
                continue;
            }
            ht = 0;
            smp = sf_cram_data2_get_sample (cram_point->data, i,
                                            zoffset ? 2.0*cram_point->src_exits[ies].t
                                                    : cram_point->src_exits[ies].t +
                                                      cram_point->rcv_exits[ier].t,
                                            cram_point->rcv_exits[ier].p[0],
                                            cram_point->dh,
                                            zoffset ? 2*cram_point->src_exits[ies].kmah
                                                    : cram_point->src_exits[ies].kmah +
                                                      cram_point->rcv_exits[ier].kmah,
                                            &ht);
            if (0 == ht || 0.0 == smp) { /* No data samples have been extracted */
                ier++;
                continue;
            }
            /* Taper associated with the edges of the mute zone */
            aw = sf_cram_point3_gataper (cram_point, oa, da);
            /* Add sample to the scattering and dip angle gathers */
            sf_cram_point3_fill_abins (cram_point, aw*tw*smp, ss, sr, ies, ier,
                                       oa, oz, da, dz, zoffset);
            ier++;
        } /* Loop over known receiver branches */
        ies++;
    } /* Loop over known source branches */
}

/* Compute taper function(weight) depending on the distance from the
   receiver to the edges of the current shot point */
static float sf_cram_point3_staper (sf_cram_point3 cram_point, float gx, float gy,
                                    float gxmin, float gxmax, float gymin, float gymax) {
    float w = 1.0;

    if (gx < gxmin || gx > gxmax ||
        gy < gymin || gy > gymax)
        return 0.0;

    if ((gx - gxmin) < cram_point->dxm)
        w *= 1.0 - (gx - gxmin)/cram_point->dxm;
    if ((gxmax - gx) < cram_point->dxm)
        w *= (gxmax - gx)/cram_point->dxm;
    if ((gy - gymin) < cram_point->dym)
        w *= 1.0 - (gy - gymin)/cram_point->dym;
    if ((gymax - gy) < cram_point->dym)
        w *= (gymax - gy)/cram_point->dym;

    return w;
}

void sf_cram_point3_compute (sf_cram_point3 cram_point,
                             float z, float x, float y, float ***esc)
/*< Compute image for one point (z,x,y) and a correpsonding set of escape values >*/
{
    int nes, ner;
    size_t i, is, ih, nh;
    float sx, sy, gx, gy, ss, sr, w;
    float gxmin, gxmax, gymin, gymax;

    sf_cram_point3_reset (cram_point);

    /* Reinit escape branches */
    sf_cram_rbranch3_set_escapes (cram_point->rbranch, esc);

    /* Slowness for the source and receiver phases */
    ss = sf_cram_slowness3_get_value (cram_point->slowness, z, x, y);
    sr = ss;

    /* Loop over known sources */
    is = sf_cram_survey3_get_first_source (cram_point->survey, &sx, &sy, &nh,
                                           &gxmin, &gxmax, &gymin, &gymax);
    while (is != (size_t)-1) {
        nes = sf_cram_rbranch3_find_exits (cram_point->rbranch, sx, sy,
                                           cram_point->na, cram_point->src_exits);
        if (nes) { /* Proceed to receivers, if shot branches are found */
            /* Loop over known receivers */
            ih = sf_cram_survey3_get_first_receiver (cram_point->survey, is,
                                                     &i, &gx, &gy);
            while (ih != (size_t)-1) {
                ner = sf_cram_rbranch3_find_exits (cram_point->rbranch, gx, gy,
                                                   cram_point->na, cram_point->rcv_exits);
                if (ner) {
                    if (nh > 1 && cram_point->taper)
                        w = sf_cram_point3_staper (cram_point, gx, gy, gxmin, gxmax, gymin, gymax);
                    else
                        w = 1.0;
                    /* Extract samples from data and contribute to the image */
                    sf_cram_point3_process_branches (cram_point, gx, gy, i, nes, ner, ss, sr, w,
                                                     1 == nh && nes == ner && /* Zero offset */
                                                     hypotf (sx - gx, sy - gy) < 0.1*cram_point->dh);
                }
                ih = sf_cram_survey3_get_next_receiver (cram_point->survey, is, ih,
                                                        &i, &gx, &gy);
            } /* Loop over known receivers */
        } /* If there are shot branches */
        is = sf_cram_survey3_get_next_source (cram_point->survey, is, &sx, &sy, &nh,
                                              &gxmin, &gxmax, &gymin, &gymax);
    } /* Loop over known sources */
}

bool sf_cram_point3_compute_one_trace (sf_cram_point3 cram_point, float sx, float sy,
                                       float gx, float gy, float gxmin, float gxmax,
                                       float gymin, float gymax, float s,
                                       size_t i, size_t nh)
/*< Compute image for one trace with source position (sx,sy) and receiver position(gx,gy);
    s - phase slowness at the image point, nh - number of traces in the current gather,
    i - current trace index according to cram_survey3 >*/
{
    int nes, ner;
    float w;

    nes = sf_cram_rbranch3_find_exits (cram_point->rbranch, sx, sy,
                                       cram_point->na, cram_point->src_exits);
    if (nes) { /* Proceed to receivers, if shot branches are found */
        ner = sf_cram_rbranch3_find_exits (cram_point->rbranch, gx, gy,
                                           cram_point->na, cram_point->rcv_exits);
        if (ner) {
            if (nh > 1 && cram_point->taper)
                w = sf_cram_point3_staper (cram_point, gx, gy, gxmin, gxmax, gymin, gymax);
            else
                w = 1.0;
            /* Extract samples from data and contribute to the image */
            sf_cram_point3_process_branches (cram_point, gx, gy, i, nes, ner, s, s, w,
                                             1 == nh && nes == ner && /* Zero offset */
                                             hypotf (sx - gx, sy - gy) < 0.1*cram_point->dh);
            return true;
       }
    }
    return false;
}


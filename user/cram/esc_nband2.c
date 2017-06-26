/* Narrow-band for escape equations in 3-D reduced phase space */
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

#ifndef _esc_nband2_h

#include "esc_point2.h"
/*^*/

#define ESC2_NQUADS 4
/* Number of distinct phase space direction */
#define ESC2_MORDER 2
/* Highest F-D stencil order */
/*^*/

typedef struct EscNBand2 *sf_esc_nband2;
/* abstract data type */
/*^*/

typedef struct EscNBChildIter2 *sf_esc_nbchild2_iter;
/* abstract data iterator type */
/*^*/

#endif

/* Iterator direction */
typedef enum { NB2_BACK = 0, NB2_FORW = 1 } NBIter2;


struct EscNBand2 {
    int              nz, nx, na; /* Dimensions of phase space */
    unsigned long    nc; /* Number of points in the band */
    int              iz; /* Current band position */
    EscDirection2  zdir; /* Direction of the band movement in z */
    unsigned char ***nb; /* Values in the narrow band */
    unsigned char  **nb2; /* Pointer for the middle NB level */ 
};
/* concrete data type */

struct EscNBChildIter2 {
    int           ix, ia; /* Current angle and x position */ 
    int           nz, nx, na;
    NBIter2       xdir, adir; /* Direction of the iterations */
    sf_esc_nband2 nb; /* Link to the narrow band */
};
/* concrete data iterator type */

sf_esc_nband2 sf_esc_nband2_init (int nz, int nx, int na)
/*< Initialize object >*/
{
    sf_esc_nband2 esc_nband = (sf_esc_nband2)sf_alloc (1, sizeof (struct EscNBand2));

    esc_nband->nz = nz;
    esc_nband->nx = nx;
    esc_nband->na = na;

    esc_nband->iz = 0;
    esc_nband->zdir = ESC2_FORW;
    esc_nband->nc = na/2*nx*(ESC2_MORDER + 1);

    esc_nband->nb = sf_ucharalloc3 (na/2*sf_esc_point2_sizeof (), nx, ESC2_MORDER + 1);
    memset (esc_nband->nb[0][0], 0, (size_t)esc_nband->nc*(size_t)sf_esc_point2_sizeof ());
    /* Contents of nb[...] will be moved later as the narrow band progress,
       nb2 preserves the original pointer for a proper free () in the end */
    esc_nband->nb2 = esc_nband->nb[0];

    return esc_nband;
}

void sf_esc_nband2_close (sf_esc_nband2 esc_nband)
/*< Destroy object >*/
{
    free (esc_nband->nb2[0]);
    free (esc_nband->nb2);
    free (esc_nband->nb);
    free (esc_nband);
}

/* Convert narrow band index iia into global phase angle index */
static int sf_esc_nband2_iia_to_ia (sf_esc_nband2 esc_nband, int iia) {
    if (ESC2_BACK == esc_nband->zdir) {
        if (iia < esc_nband->na/4)
            return iia + 3*esc_nband->na/4;
        else
            return iia - esc_nband->na/4;
    } else
        return iia + esc_nband->na/4;
}

sf_esc_point2 sf_esc_nband2_get_child (sf_esc_nband2 esc_nband, unsigned long i,
                                       int *iz, int *ix, int *ia)
/*< Returns child point for index i and its 3-D coordinate indices >*/
{
    int iia = i % (unsigned long)(esc_nband->na/2);

    *ia = sf_esc_nband2_iia_to_ia (esc_nband, iia);
    *ix = i / (unsigned long)(esc_nband->na/2);
    *iz = esc_nband->iz;

    if ((*ix) >= esc_nband->nx)
        return NULL;

    return (sf_esc_point2)&esc_nband->nb[0][(*ix)][iia*sf_esc_point2_sizeof ()];
}

sf_esc_point2 sf_esc_nband2_get_used_point (sf_esc_nband2 esc_nband, unsigned long i,
                                            int *iz, int *ix, int *ia)
/*< Returns parent point (outermost layer) for index i and its 3-D coordinate indices >*/
{
    unsigned long nc;
    int ic, iia;

    /* Top/bottom at the beginnig of z sweeps, no used points get discarded */
    if ((esc_nband->iz < ESC2_MORDER && ESC2_FORW == esc_nband->zdir) ||
        ((esc_nband->nz - ESC2_MORDER) <= esc_nband->iz && ESC2_BACK == esc_nband->zdir))
        return NULL;

    /* Top/bottom at the end of z sweeps, all points are marked as used */
    if (((esc_nband->nz - 1) == esc_nband->iz && ESC2_FORW == esc_nband->zdir) ||
        (0 == esc_nband->iz && ESC2_BACK == esc_nband->zdir)) {
        ic = i/((unsigned long)(esc_nband->nx*esc_nband->na/2));

        if (ic > ESC2_MORDER)
            return false;
        nc = ic*((unsigned long)(esc_nband->nx*esc_nband->na/2));

        iia = (i - nc) % (unsigned long)(esc_nband->na/2);
        *ia = sf_esc_nband2_iia_to_ia (esc_nband, iia);
        *ix = (i - nc) / (unsigned long)(esc_nband->na/2);
        *iz = (ESC2_FORW == esc_nband->zdir) ? (esc_nband->iz - ic)
                                             : (esc_nband->iz + ic);

        return (sf_esc_point2)&esc_nband->nb[ic][(*ix)][iia*sf_esc_point2_sizeof ()];
    } else { /* All other cases - only outermost layer of parents gets discarded */
        iia = i % (unsigned long)(esc_nband->na/2);
        *ia = sf_esc_nband2_iia_to_ia (esc_nband, iia);
        *ix = i / (unsigned long)(esc_nband->na/2);
        *iz = (ESC2_FORW == esc_nband->zdir) ? (esc_nband->iz - ESC2_MORDER)
                                             : (esc_nband->iz + ESC2_MORDER);

        if ((*ix) >= esc_nband->nx)
            return NULL;

        return (sf_esc_point2)&esc_nband->nb[ESC2_MORDER][(*ix)][iia*sf_esc_point2_sizeof ()];
    }
}

sf_esc_point2 sf_esc_nband2_get_used_plane (sf_esc_nband2 esc_nband, unsigned long i,
                                            unsigned long *size)
/*< Returns pointer to a plane of used points; size is a number of points in the plane >*/
{
    /* Top/bottom at the beginnig of z sweeps, there are no used planes to be discarded */
    if ((esc_nband->iz < ESC2_MORDER && ESC2_FORW == esc_nband->zdir) ||
        ((esc_nband->nz - ESC2_MORDER) <= esc_nband->iz && ESC2_BACK == esc_nband->zdir))
        return NULL;

    /* Top/bottom at the end of z sweeps, return all used planes one by one */
    if (((esc_nband->nz - 1) == esc_nband->iz && ESC2_FORW == esc_nband->zdir) ||
        (0 == esc_nband->iz && ESC2_BACK == esc_nband->zdir)) {
        if (i > ESC2_MORDER)
            return NULL;
        *size = esc_nband->nx*esc_nband->na/2;
        return (sf_esc_point2)&esc_nband->nb[ESC2_MORDER - i][0][0];
    } else { /* All other cases - only outermost plane of points gets returned */
        if (i != 0)
            return NULL;
        *size = esc_nband->nx*esc_nband->na/2;
        return (sf_esc_point2)&esc_nband->nb[ESC2_MORDER][0][0];
    }
}

bool sf_esc_nband2_next_generation (sf_esc_nband2 esc_nband)
/*< Discards previous generation of parents, make room for a new
    generation of child points >*/
{
    int i;
    unsigned char **tmp;

    if (0 == esc_nband->iz && ESC2_BACK == esc_nband->zdir)
        return false;

    /* Check if the narrow band direction has to be reversed */
    if ((esc_nband->nz - 1) == esc_nband->iz && ESC2_FORW == esc_nband->zdir)
        esc_nband->zdir = ESC2_BACK;
    else
        esc_nband->iz += (ESC2_FORW == esc_nband->zdir) ? 1 : -1;

    /* Rotate layers in the narrow band */
    tmp = esc_nband->nb[ESC2_MORDER];
    for (i = (ESC2_MORDER - 1); i >= 0; i--) {
        esc_nband->nb[i + 1] = esc_nband->nb[i];
    }
    esc_nband->nb[0] = tmp;
    /* Reset points in the child layer */
    memset (&esc_nband->nb[0][0][0], 0, (size_t)esc_nband->nc/(size_t)(ESC2_MORDER + 1)*
                                        (size_t)sf_esc_point2_sizeof ());

    return true;
}

unsigned long sf_esc_nband2_child_count (sf_esc_nband2 esc_nband)
/*< Returns number of child points in the band >*/
{
    return (unsigned long)(esc_nband->nx*esc_nband->na/2);
}

unsigned long sf_esc_nband2_point_count (sf_esc_nband2 esc_nband)
/*< Returns number of points in the band >*/
{
    return (unsigned long)(esc_nband->nc);
}

/*
 *
 *  Iterator over child points
 *
 */

void sf_esc_nbchild2_iter_reset (sf_esc_nbchild2_iter esc_citer)
/*< Destroy iterator object >*/
{
    esc_citer->ia = -1;
    esc_citer->ix = 0;
    esc_citer->adir = NB2_FORW;
    esc_citer->xdir = NB2_FORW;
}

sf_esc_nbchild2_iter sf_esc_nbchild2_iter_init (sf_esc_nband2 esc_nband)
/*< Initialize iterator object >*/
{
    sf_esc_nbchild2_iter esc_citer = (sf_esc_nbchild2_iter)sf_alloc (1, sizeof (struct EscNBChildIter2));

    esc_citer->nb = esc_nband;
    esc_citer->na = esc_nband->na/4;
    esc_citer->nx = esc_nband->nx;
    esc_citer->nz = esc_nband->nz;

    sf_esc_nbchild2_iter_reset (esc_citer);

    return esc_citer;
}

void sf_esc_nbchild2_iter_close (sf_esc_nbchild2_iter esc_citer)
/*< Destroy iterator object >*/
{
    free (esc_citer);
}

sf_esc_point2 sf_esc_nbchild2_iter_next (sf_esc_nbchild2_iter esc_citer,
                                         int *iz, int *ix, int *ia,
                                         sf_esc_point2 *zp, sf_esc_point2 *xp,
                                         sf_esc_point2 *apf, sf_esc_point2 *apb)
/*< Find next child in the narrow band >*/
{
    int i, iia, iix, iiz, sz = sf_esc_point2_sizeof ();

    for (i = 0; i < ESC2_MORDER; i++) {
        zp[i] = NULL;
        xp[i] = NULL;
        apf[i] = NULL;
        apb[i] = NULL;
    }

    esc_citer->ia += (NB2_FORW == esc_citer->adir) ? 1 : -1;

    /* Check for limits, switch sweeping directions */
    if (NB2_FORW == esc_citer->adir && esc_citer->na == esc_citer->ia) {
        esc_citer->ia = esc_citer->na - 2;
        esc_citer->adir = NB2_BACK;
    } else if (NB2_BACK == esc_citer->adir && -1 == esc_citer->ia) {
        esc_citer->ix += (NB2_FORW == esc_citer->xdir) ? 1 : -1;
        if (NB2_FORW == esc_citer->xdir && esc_citer->nx == esc_citer->ix) {
            esc_citer->ix = esc_citer->nx - 1;
            esc_citer->xdir = NB2_BACK;
        } else if (NB2_BACK == esc_citer->xdir && -1 == esc_citer->ix) {
            return NULL;
        }
        esc_citer->ia = 0;
        esc_citer->adir = NB2_FORW;
    }
    iix =  esc_citer->ix;

    /* Translate from sweeping indices ia to array indices iia
       and then to global phase-space indices depending on the sweeping direction */
    if (ESC2_BACK == esc_citer->nb->zdir) { /* Moving upward */
        if (NB2_BACK == esc_citer->xdir) { /* Right->left */
            iia = esc_citer->na + esc_citer->ia;
            *ia = esc_citer->ia;
        } else { /* Left->right */
            iia = esc_citer->na - esc_citer->ia - 1;
            *ia = iia + 3*esc_citer->na;
        }
    } else { /* Moving downward */
        if (NB2_BACK == esc_citer->xdir) { /* Right->left */
            iia = esc_citer->na - esc_citer->ia - 1;
        } else { /* Left->right */
            iia = esc_citer->na + esc_citer->ia;
        }
        *ia = iia + esc_citer->na;
    }
    *ix = iix;
    *iz = esc_citer->nb->iz;

    /* Fill in neighbor arrays */
    iiz = *iz;
    if (ESC2_FORW == esc_citer->nb->zdir) { /* Moving upward */
        for (i = 0; i < ESC2_MORDER; i++) {
            if ((iiz - i - 1) >= 0)
               zp[i] = (sf_esc_point2)&esc_citer->nb->nb[i + 1][iix][iia*sz];
        }
    } else { /* Moving downward */
        for (i = 0; i < ESC2_MORDER; i++) {
            if ((iiz + i + 1) < esc_citer->nz)
               zp[i] = (sf_esc_point2)&esc_citer->nb->nb[i + 1][iix][iia*sz];
        }
    }
    /* Now, neighbors in x */
    if ((*ia) >= 2*esc_citer->na) {
        for (i = 0; i < ESC2_MORDER; i++) {
            if ((iix - i - 1) >= 0)
               xp[i] = (sf_esc_point2)&esc_citer->nb->nb[0][iix - i - 1][iia*sz];
        }
    } else {
        for (i = 0; i < ESC2_MORDER; i++) {
            if ((iix + i + 1) < esc_citer->nx)
               xp[i] = (sf_esc_point2)&esc_citer->nb->nb[0][iix + i + 1][iia*sz];
        }
    }
    /* Finally, angle neighbors */
    for (i = 0; i < ESC2_MORDER; i++) {
        if ((iia - i - 1) >= 0) /* Backward */
            apb[i] = (sf_esc_point2)&esc_citer->nb->nb[0][iix][(iia - i - 1)*sz];
        if ((iia + i + 1) < 2*esc_citer->na) /* Forward */
            apf[i] = (sf_esc_point2)&esc_citer->nb->nb[0][iix][(iia + i + 1)*sz];
    }

    return (sf_esc_point2)&esc_citer->nb->nb[0][iix][iia*sz];
}


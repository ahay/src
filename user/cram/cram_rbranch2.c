/* Determine surface-exiting ray branches from escape values in 2-D */
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

#include "esc_point2.h"

#include "cram_slowness2.h"

#ifndef _cram_rbranch2_h

typedef struct CRAMRBranch2 *sf_cram_rbranch2;
/* abstract data type */
/*^*/

typedef struct {
    int   ia, kmah;
    float t, j, p, d;
} sf_cram_surface_exit2;
/*^*/

#endif

typedef struct {
    int   kmah, state;
    float x0, x1;
    float t0, t1;
    float j, p;
} sf_cram_surface_branch2;

struct CRAMRBranch2 {
    int                      na, ne;
    float                    z0, tmax, mwidth;
    float                    xmin, xmax;
    int                      ib0, ib1;
    sf_cram_slowness2        slowness;
    sf_cram_surface_branch2 *surface_branches;
    sf_cram_surface_exit2   *src_exits;
    sf_cram_surface_exit2   *rcv_exits;
};
/* concrete data type */

typedef enum { STATE_NORMAL = 0, STATE_REMOVED = 1 } branch_state;

sf_cram_rbranch2 sf_cram_rbranch2_init (int na, float z0, float tmax, sf_cram_slowness2 slowness)
/*< Initialize object >*/
{
    sf_cram_rbranch2 cram_rbranch = (sf_cram_rbranch2)sf_alloc (1, sizeof (struct CRAMRBranch2));

    cram_rbranch->na = na; /* Number of escape angles */
    cram_rbranch->z0 = z0; /* Escape surface (where data are) */
    cram_rbranch->tmax = tmax; /* Maximum time along a ray */
    cram_rbranch->slowness = slowness; /* For velocity at the surface */

    cram_rbranch->xmin = SF_HUGE;
    cram_rbranch->xmax = -SF_HUGE;
    cram_rbranch->mwidth = SF_HUGE;

    cram_rbranch->ib0 = na;
    cram_rbranch->ib1 = -1;
    cram_rbranch->surface_branches = (sf_cram_surface_branch2*)sf_alloc (na, sizeof (sf_cram_surface_branch2));
    cram_rbranch->src_exits = (sf_cram_surface_exit2*)sf_alloc (na + 1, sizeof (sf_cram_surface_exit2));
    cram_rbranch->rcv_exits = (sf_cram_surface_exit2*)sf_alloc (na + 1, sizeof (sf_cram_surface_exit2));

    return cram_rbranch;
}

void sf_cram_rbranch2_close (sf_cram_rbranch2 cram_rbranch)
/*< Destroy object >*/
{
    free (cram_rbranch->src_exits);
    free (cram_rbranch->rcv_exits);
    free (cram_rbranch->surface_branches);
    free (cram_rbranch);
}

void sf_cram_rbranch2_set_maxwidth (sf_cram_rbranch2 cram_rbranch, float mwidth)
/*< Set maximum allowed width for an exit ray on the surface >*/
{
    cram_rbranch->mwidth = mwidth;
}

void sf_cram_rbranch2_set_escapes (sf_cram_rbranch2 cram_rbranch, float **esc)
/*< Set new escape variables >*/
{
    int ia0, ia1, na = cram_rbranch->na;
    float x0, x1, t0, t1, vsurf, p, sn;
    float xmin = SF_HUGE, xmax = -SF_HUGE;

    cram_rbranch->ib0 = na;
    cram_rbranch->ib1 = -1;

    /* Scan ray branches and store those exiting on the surface */
    for (ia1 = 0; ia1 < na; ia1++) {
        ia0 = (ia1 != 0) ? ia1 - 1 : cram_rbranch->na - 1;
        if (esc[ia0][ESC2_Z] > cram_rbranch->z0 ||
            esc[ia1][ESC2_Z] > cram_rbranch->z0) {
            cram_rbranch->surface_branches[ia0].state = STATE_REMOVED;
            continue; /* Rays do not exit on the surface */
        }
        t0 = esc[ia0][ESC2_T];
        t1 = esc[ia1][ESC2_T];
        if (t0 < 0.0 || t1 < 0.0) {
            cram_rbranch->surface_branches[ia0].state = STATE_REMOVED;
            continue; /* Ray branch is outside of time constraints */
        }
        if (t0 > cram_rbranch->tmax && t1 > cram_rbranch->tmax) {
            cram_rbranch->surface_branches[ia0].state = STATE_REMOVED;
            continue; /* Ray branch is outside of time constraints */
        }
        x0 = esc[ia0][ESC2_X];
        x1 = esc[ia1][ESC2_X];
        if (x0 == x1 || fabsf (x0 - x1) > cram_rbranch->mwidth) {
            cram_rbranch->surface_branches[ia0].state = STATE_REMOVED;
            continue; /* Ray branch is too narrow or too wide */
        }
        /* Slope */
        p = (t0 - t1)/(x1 - x0);
        vsurf = 1.0/sf_cram_slowness2_get_surf_value (cram_rbranch->slowness, 0.5*(x1 + x0));
        /* Sine of exit angle */
        sn = fabsf (vsurf*p);
        if (sn >= 1.0) {
            cram_rbranch->surface_branches[ia0].state = STATE_REMOVED;
            continue; /* Too wide ray tube */
        }
        /* Store branch */
        cram_rbranch->surface_branches[ia0].t0 = t0;
        cram_rbranch->surface_branches[ia0].t1 = t1;
        cram_rbranch->surface_branches[ia0].x0 = x0;
        cram_rbranch->surface_branches[ia0].x1 = x1;
        cram_rbranch->surface_branches[ia0].p = p;
        /* Geometric spreading */
        cram_rbranch->surface_branches[ia0].j = fabsf ((x1 - x0)*sqrtf (1.0f - sn*sn));
        cram_rbranch->surface_branches[ia0].kmah = (x0 - x1) < 0.0 ? 1 : 0;
        cram_rbranch->surface_branches[ia0].state = STATE_NORMAL;
        if (na == cram_rbranch->ib0)
            cram_rbranch->ib0 = ia1;
        cram_rbranch->ib1 = ia1;
        if (x0 < xmin)
            xmin = x0;
        if (x0 > xmax)
            xmax = x0;
        if (x1 < xmin)
            xmin = x1;
        if (x1 > xmax)
            xmax = x1;
    }
    cram_rbranch->xmin = xmin;
    cram_rbranch->xmax = xmax;
}

/* Find all exit positions and return array with exit information */
static sf_cram_surface_exit2* sf_cram_rbranch2_find_exits (sf_cram_rbranch2 cram_rbranch, float x,
                                                           sf_cram_surface_exit2* exits)
{
    int ib, ie = 0;
    float dx0, dx1, f, t;
    sf_cram_surface_branch2 branch;

    if (cram_rbranch->na == cram_rbranch->ib0 ||
        x > cram_rbranch->xmax || x < cram_rbranch->xmin)
        return NULL;

    ib = cram_rbranch->ib0;
    /* Search for normal branches and check for hits */
    while (ib <= cram_rbranch->ib1) {
        if (STATE_NORMAL == cram_rbranch->surface_branches[ib].state) {
            branch = cram_rbranch->surface_branches[ib];
            dx0 = branch.x0 - x;
            dx1 = branch.x1 - x;
            if ((dx0 <= 0.0 && dx1 >= 0.0) ||
                (dx0 >= 0.0 && dx1 <= 0.0)) {
                f = fabs ((x - branch.x0)/(branch.x1 - branch.x0));
                /* Exit time */
                t = branch.t0 + f*(branch.t1 - branch.t0);
                if (t <= cram_rbranch->tmax) {
                    exits[ie].t = t;
                    exits[ie].ia = ib;
                    exits[ie].kmah = branch.kmah;
                    exits[ie].p = branch.p;
                    exits[ie].j = branch.j;
                    exits[ie].d = fabsf (branch.x0 - branch.x1);
                    ie++;
                }
            }
        }
        ib++;
    }
    exits[ie].t = SF_HUGE;
    exits[ie].ia = -1;

    return (ie != 0) ? exits : NULL;
}

sf_cram_surface_exit2* sf_cram_rbranch2_src_exits (sf_cram_rbranch2 cram_rbranch, float x)
/*< Find all source exit positions and return array with exit information >*/
{
    return sf_cram_rbranch2_find_exits (cram_rbranch, x, cram_rbranch->src_exits);
}

sf_cram_surface_exit2* sf_cram_rbranch2_rcv_exits (sf_cram_rbranch2 cram_rbranch, float x)
/*< Find all receiver exit positions and return array with exit information >*/
{
    return sf_cram_rbranch2_find_exits (cram_rbranch, x, cram_rbranch->rcv_exits);
}


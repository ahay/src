/* 3-D phase-space grid for 2-D medium with narrow-band support */
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

#ifndef _esc_nbgrid2_h

#include "esc_slow2.h"
/*^*/

#include "esc_nband2.h"
/*^*/

#include "esc_nbout2.h"
/*^*/

#include "esc_tracer2.h"
/*^*/

typedef struct EscNBGrid2 *sf_esc_nbgrid2;
/* abstract data type */
/*^*/

#endif

struct EscNBGrid2 {
    int                  nz, nx, na, jz, jx, morder;
    float                oz, ox, oa;
    float                dz, dx, da;
    float                xmax, zmax;
    float                d[ESC2_DIMS], mdist;
    size_t               nct, ndt, fdt, nlt, nbt, pgt, bct;
    bool                 verb, cmix, atraced, mtraced;
    sf_timer             ttime, rtime;
    sf_esc_nband2        nband;
    sf_esc_nbchild2_iter nbciter;
    sf_esc_slowness2     esc_slow;
    sf_esc_tracer2       esc_tracer;
    sf_esc_nbout2        esc_out;
};
/* concrete data type */

float sf_esc_nbgrid2_get_oa (int na)
/*< Returns origin of the angular axis >*/
{
    return -SF_PI + SF_PI/(float)na;
}

float sf_esc_nbgrid2_get_da (int na)
/*< Returns sampling of the angular axis >*/
{
    return 2.0*SF_PI/(float)na;
}

sf_esc_nbgrid2 sf_esc_nbgrid2_init (int nz, int nx, int na,
                                    float oz, float ox, float dz, float dx,
                                    bool atraced, bool mtraced,
                                    sf_esc_slowness2 esc_slow, sf_esc_tracer2 esc_tracer,
                                    sf_esc_nbout2 esc_out)
/*< Initialize object >*/
{
    sf_esc_nbgrid2 esc_nbgrid = (sf_esc_nbgrid2)sf_alloc (1, sizeof (struct EscNBGrid2));

    esc_nbgrid->nz = nz;
    esc_nbgrid->nx = nx;

    if (na % ESC2_NQUADS)
        na += ESC2_NQUADS - (na % ESC2_NQUADS);
    esc_nbgrid->na = na;

    esc_nbgrid->dz = dz;
    esc_nbgrid->d[ESC2_AXIS_Z] = dz;
    esc_nbgrid->dx = dx;
    esc_nbgrid->d[ESC2_AXIS_X] = dx;
    esc_nbgrid->da = sf_esc_nbgrid2_get_da (na);
    esc_nbgrid->d[ESC2_AXIS_A] = esc_nbgrid->da;

    esc_nbgrid->oz = oz;
    esc_nbgrid->ox = ox;
    esc_nbgrid->oa = sf_esc_nbgrid2_get_oa (na);

    esc_nbgrid->zmax = oz + (nz - 1)*dz;
    esc_nbgrid->xmax = ox + (nx - 1)*dx;

    esc_nbgrid->mdist = SF_HUGE;
    esc_nbgrid->morder = ESC2_MORDER;
    esc_nbgrid->cmix = true;

    esc_nbgrid->nct = 0; /* Number of ray traced because of color mixing */
    esc_nbgrid->ndt = 0; /* Number of ray traced because of maximum allowed exit distance spread */
    esc_nbgrid->fdt = 0; /* Number of ray traced because of F-D overshooting */
    esc_nbgrid->nlt = 0; /* Number of ray traced because of deadlocks */
    esc_nbgrid->pgt = 0; /* Number of ray traced because of phase/group velocity deviation */
    esc_nbgrid->nbt = 0; /* Number of ray traced for the narrow band edges */
    esc_nbgrid->bct = 0; /* Number of ray traced for the boundary conditions */

    esc_nbgrid->verb = false;

    esc_nbgrid->atraced = atraced;
    esc_nbgrid->mtraced = mtraced;

    esc_nbgrid->ttime = sf_timer_init (); /* Total time */
    esc_nbgrid->rtime = sf_timer_init (); /* Ray tracing time */

    esc_nbgrid->esc_slow = esc_slow;
    esc_nbgrid->esc_tracer = esc_tracer;
    esc_nbgrid->esc_out = esc_out;

    esc_nbgrid->nband = sf_esc_nband2_init (nz, nx, na);
    esc_nbgrid->nbciter = sf_esc_nbchild2_iter_init (esc_nbgrid->nband);

    return esc_nbgrid;
}

void sf_esc_nbgrid2_close (sf_esc_nbgrid2 esc_nbgrid)
/*< Destroy object >*/
{
    sf_esc_nbchild2_iter_close (esc_nbgrid->nbciter);
    sf_esc_nband2_close (esc_nbgrid->nband);
    sf_timer_close (esc_nbgrid->ttime);
    sf_timer_close (esc_nbgrid->rtime);
    free (esc_nbgrid);
}

void sf_esc_nbgrid2_set_verb (sf_esc_nbgrid2 esc_nbgrid, bool verb)
/*< Set verbatim flag >*/
{
    esc_nbgrid->verb = verb;
}

void sf_esc_nbgrid2_set_mdist (sf_esc_nbgrid2 esc_nbgrid, float mdist)
/*< Set maximum distance between escape values in F-D stencil >*/
{
    esc_nbgrid->mdist = mdist;
}

void sf_esc_nbgrid2_set_morder (sf_esc_nbgrid2 esc_nbgrid, int morder)
/*< Set maximum order in F-D stencil >*/
{
    if (morder < 1)
        morder = 1;
    if (morder > ESC2_MORDER)
        morder = ESC2_MORDER;
    esc_nbgrid->morder = morder;
}

void sf_esc_nbgrid2_set_cmix (sf_esc_nbgrid2 esc_nbgrid, bool cmix)
/*< Set color mixing check flag >*/
{
    esc_nbgrid->cmix = cmix;
}

/* Do simple Gauss-Seidel update for one point and one type of variable */
static float sf_esc_nbgrid2_gs_update (sf_esc_nbgrid2 esc_nbgrid,
                                       float *n, float *e,
                                       float f[ESC2_DIMS][ESC2_MORDER], float rhs) {
    int j;
    float v = rhs, dd = 0.0;

    for (j = 0; j < ESC2_DIMS; j++) {
        if (n[j] > 1) { /* Second order */
            v += e[j]*(2.0*f[j][0] - 0.5*f[j][1]);
            dd += 1.5*e[j];
        } else if (1 == n[j]) { /* First order */
            v += e[j]*f[j][0];
            dd += e[j];
        }
    }

    return v/dd;
}

/* Apply finite-difference stencil for three neighbors */
static void sf_esc_nbgrid2_compute_point (sf_esc_nbgrid2 esc_nbgrid, int iz, int ix, int ia,
                                          sf_esc_point2 child, sf_esc_point2 parents[ESC2_DIMS][ESC2_MORDER]) {
    int j, k, l;
    float a, rhs, e[ESC2_DIMS], n[ESC2_DIMS], f[ESC2_DIMS][ESC2_MORDER];
    float s, sa = 0.0, sz = 0.0, sx = 0.0;
    EscColor2 col = 0;

    /* Get coefficients of the escape equations */
    a = esc_nbgrid->oa + ia*esc_nbgrid->da;
    sf_esc_slowness2_get_components (esc_nbgrid->esc_slow,
                                     esc_nbgrid->oz + iz*esc_nbgrid->dz,
                                     esc_nbgrid->ox + ix*esc_nbgrid->dx,
                                     a, &s,
                                     parents[ESC2_AXIS_Z] ||
                                     parents[ESC2_AXIS_X] ? &sa : NULL,
                                     parents[ESC2_AXIS_A] ? &sz : NULL,
                                     parents[ESC2_AXIS_A] ? &sx : NULL);
    sf_esc_slowness2_get_coefs (esc_nbgrid->esc_slow, a, s, sa, sz, sx,
                                &e[ESC2_AXIS_Z], &e[ESC2_AXIS_X], &e[ESC2_AXIS_A]); 

    for (j = 0; j < ESC2_DIMS; j++) {
        e[j] = (parents[j][0] != NULL) ? fabsf (e[j]/esc_nbgrid->d[j])
                                       : 0.0;
    }

    /* Apply regular finite-difference for all escape variables */
    for (k = 0; k < ESC2_NUM; k++) {
        for (j = 0; j < ESC2_DIMS; j++) {
            n[j] = 0;
            for (l = 0; (l < ESC2_MORDER && l < esc_nbgrid->morder); l++) {
                if (parents[j][l]) {
                    f[j][l] = sf_esc_point2_get_esc_var (parents[j][l], k);
                    n[j]++;
                } else
                    break;
            }
        }
        if (ESC2_T == k)
            rhs = s*s;
#ifdef ESC_EQ_WITH_L
        else if (ESC2_L == k)
            rhs = s;
#endif
        else
            rhs = 0.0;

        sf_esc_point2_set_esc_var (child, k,
                                   sf_esc_nbgrid2_gs_update (esc_nbgrid, n, e, f, rhs));
    }

    for (j = 0; j < ESC2_DIMS; j++) {
        for (k = 0; k < ESC2_MORDER; k++) {
            if (parents[j][k])
                col |= sf_esc_point2_get_col (parents[j][k]);
        }
    }

    sf_esc_point2_set_col (child, col);
}

/* Returns true, if point belongs to boundary conditions */
static bool sf_esc_nbgrid2_is_bc (sf_esc_nbgrid2 esc_nbgrid, int iz, int ix, int ia) {
    if ((0 == iz && ia >= esc_nbgrid->na/4 && ia < 3*esc_nbgrid->na/4) ||
        ((esc_nbgrid->nz - 1) == iz && (ia < esc_nbgrid->na/4 || ia >= 3*esc_nbgrid->na/4)) ||
        (0 == ix && ia >= esc_nbgrid->na/2) ||
        ((esc_nbgrid->nx - 1) == ix && ia < esc_nbgrid->na/2))
        return true;
    return false;
}

/* Returns true, if assumed phase direction (ia) coincides with
   direction of parent points */
static bool sf_esc_nbgrid2_same_group_phase (sf_esc_nbgrid2 esc_nbgrid, bool *hasp,
                                             EscDirection2 *dirs, int ia) {
    if (hasp[ESC2_AXIS_Z]) {
        if (((ia >= 3*esc_nbgrid->na/4 || ia < esc_nbgrid->na/4) &&
             ESC2_BACK == dirs[ESC2_AXIS_Z]) ||
            (ia < 3*esc_nbgrid->na/4 && ia >= esc_nbgrid->na/4 &&
             ESC2_FORW == dirs[ESC2_AXIS_Z]))
            return false;
    }
    if (hasp[ESC2_AXIS_X]) {
        if ((ia < esc_nbgrid->na/2 && ESC2_BACK == dirs[ESC2_AXIS_X]) ||
            (ia >= esc_nbgrid->na/2 && ESC2_FORW == dirs[ESC2_AXIS_X]))
            return false;
    }
    return true;
}

/* Compute a child point, if all of its parents are present */
static bool sf_esc_nbgrid2_evaluate_child (sf_esc_nbgrid2 esc_nbgrid, int iz, int ix, int ia,
                                           sf_esc_point2 child, sf_esc_point2 parents[ESC2_DIMS][ESC2_MORDER]) {
    int j, k, l, m, nb = esc_nbgrid->na/4;
    bool trace = false, trmdist = false;
    bool hasp[ESC2_DIMS] = { false, false, false };
    EscDirection2 pdir, dirs[ESC2_DIMS];
    EscColor2 col1, col2;
    float x1, x2, x, z;

    /* Boundary condition? */
    if (sf_esc_nbgrid2_is_bc (esc_nbgrid, iz, ix, ia)) {
        trace = true;
        esc_nbgrid->bct++;
    } else {
        /* See in which directions parent points are expected */
        for (j = 0; j < ESC2_DIMS; j++) {
            hasp[j] = sf_esc_point2_has_parent_link (child, j, &dirs[j]);
        }
    }

    /* Check if phase and group direction are in the same quadrants */
    if (!trace && !sf_esc_nbgrid2_same_group_phase (esc_nbgrid, hasp, dirs, ia)) {
        trace = true;
        esc_nbgrid->pgt++;
    }

    if (!trace && hasp[ESC2_AXIS_A]) {
        /* Check if child point and its parent(s) point in the same angular direction */
        for (j = 0; j < ESC2_MORDER; j++) {
            if (parents[ESC2_AXIS_A][j] && 
                sf_esc_point2_has_parent_link (parents[ESC2_AXIS_A][j], ESC2_AXIS_A, &pdir) &&
                pdir != dirs[ESC2_AXIS_A]) {
                trace = true;
                esc_nbgrid->nlt++;
                break;
            }
        }
    }

    if (!trace && hasp[ESC2_AXIS_A]) {
        /* Check if current point is on the edge of the narrow band (angle direction) */
        if (NULL == parents[ESC2_AXIS_A][0] ||
            (!sf_esc_point2_is_parent (parents[ESC2_AXIS_A][0]) &&
             ((0 == ia % nb && ESC2_BACK == dirs[ESC2_AXIS_A]) ||
             ((nb - 1) == ia % nb && ESC2_FORW == dirs[ESC2_AXIS_A])))) {
            trace = true;
            esc_nbgrid->nbt++;
        } else if (parents[ESC2_AXIS_A][0] && !sf_esc_point2_is_parent (parents[ESC2_AXIS_A][0]))
            return false;
        /* If one of the parents is not ready - fall back to first order in angle direction */
        if (parents[ESC2_AXIS_A][1] && !sf_esc_point2_is_parent (parents[ESC2_AXIS_A][1]))
            parents[ESC2_AXIS_A][1] = NULL;
    }

    if (!trace) {
        /* Check if all parents in z and x are present */
        for (j = 0; j < (ESC2_DIMS - 1); j++) {
            if (hasp[j]) {
                if (NULL == parents[j][0] || !sf_esc_point2_is_parent (parents[j][0]))
                    sf_error ("Missing spatial parent");
            }
        }
    }

    if (!trace && (esc_nbgrid->cmix || esc_nbgrid->mdist != SF_HUGE)) {
        /* Check for color mixing */
        for (j = 0; (j < (ESC2_DIMS - 1) && !trace); j++) {
            if (hasp[j]) {
                for (l = 0; (l < ESC2_MORDER && !trace && parents[j][l]); l++) {
                    x1 = sf_esc_point2_get_esc_var (parents[j][l], ESC2_X);
                    col1 = sf_esc_point2_get_col (parents[j][l]);
                    for (k = j + 1; (k < ESC2_DIMS && !trace); k++) {
                        if (hasp[k]) {
                            for (m = 0; (m < ESC2_MORDER && !trace && parents[k][m]); m++) {
                                x2 = sf_esc_point2_get_esc_var (parents[k][m], ESC2_X);
                                col2 = sf_esc_point2_get_col (parents[k][m]);
                                /* Hitting corner */
                                if (col1 != col2 && esc_nbgrid->cmix) {
                                    trace = true;
                                    esc_nbgrid->nct++;
                                    break;
                                }
                                /* Exceeding maximum allowed distance between parents */
                                if (ESC2_TOP == col1 && fabs (x1 - x2) > esc_nbgrid->mdist) {
                                    trace = true;
                                    esc_nbgrid->ndt++;
                                    trmdist = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (!trace) { /* F-D */
        sf_esc_nbgrid2_compute_point (esc_nbgrid, iz, ix, ia, child,
                                      parents);
        /* Check for F-D overshooting */
        z =  sf_esc_point2_get_esc_var (child, ESC2_Z);
        x =  sf_esc_point2_get_esc_var (child, ESC2_X);
        if (z < esc_nbgrid->oz || z > esc_nbgrid->zmax ||
            x < esc_nbgrid->ox || x > esc_nbgrid->xmax) {
            esc_nbgrid->fdt++;
            trace = true;
        }
    }

    if (trace) { /* Do explicit ray tracing */
        if (esc_nbgrid->verb)
            sf_timer_start (esc_nbgrid->rtime);
        sf_esc_tracer2_compute (esc_nbgrid->esc_tracer,
                                esc_nbgrid->oz + iz*esc_nbgrid->dz,
                                esc_nbgrid->ox + ix*esc_nbgrid->dx,
                                esc_nbgrid->oa + ia*esc_nbgrid->da,
                                0.0, 0.0, child, NULL);
        if ((esc_nbgrid->mtraced && trmdist) || esc_nbgrid->atraced)
            sf_esc_point2_set_traced (child, true);
        if (esc_nbgrid->verb)
            sf_timer_stop (esc_nbgrid->rtime);
    }

    sf_esc_point2_become_parent (child);

    return true;
}

/* Sweep across points in the narrow band and compute those children,
   which have all of their parents present; return number of computed points */
static unsigned long sf_esc_nbgrid2_evaluate_children (sf_esc_nbgrid2 esc_nbgrid) {
    unsigned long nc = 0;
    int j, k, iz, ix, ia;
    sf_esc_point2 child = NULL;
    sf_esc_point2 apf[ESC2_MORDER] = { NULL, NULL }, apb[ESC2_MORDER] = { NULL, NULL };
    sf_esc_point2 parents[ESC2_DIMS][ESC2_MORDER] = { { NULL, NULL },
                                                      { NULL, NULL },
                                                      { NULL, NULL } };
    EscDirection2 cdir;

    while ((child = (sf_esc_point2)sf_esc_nbchild2_iter_next (esc_nbgrid->nbciter,
                                                              &iz, &ix, &ia,
                                                              parents[ESC2_AXIS_Z],
                                                              parents[ESC2_AXIS_X],
                                                              apf, apb))) {
        if (!sf_esc_point2_is_child (child))
            continue;
        /* Choose between a + 1 and a - 1 parent point according to what
           child point actually expects */
        if (sf_esc_point2_has_parent_link (child, ESC2_AXIS_A, &cdir)) {
            for (j = 0; j < ESC2_MORDER; j++) {
                 parents[ESC2_AXIS_A][j] = (ESC2_FORW == cdir) ? apf[j] : apb[j];
            }
        }
        nc += sf_esc_nbgrid2_evaluate_child (esc_nbgrid, iz, ix, ia, child, parents);
        for (j = 0; j < ESC2_DIMS; j++) {
            for (k = 0; k < ESC2_MORDER; k++) {
                parents[j][k] = NULL;
            }
        }
        for (j = 0; j < ESC2_MORDER; j++) {
            apf[j] = NULL; apb[j] = NULL;
        }
    }

    if (nc != sf_esc_nband2_child_count (esc_nbgrid->nband))
        sf_error ("sf_esc_nbgrid2: not all points have been evaluated [%lu, %lu]",
                  nc, sf_esc_nband2_child_count (esc_nbgrid->nband));

    sf_esc_nbchild2_iter_reset (esc_nbgrid->nbciter);

    return nc;
}

/* Find parents in the narrow band, save their escape values,
   return number of found points */
static unsigned long sf_esc_nbgrid2_save_used_points (sf_esc_nbgrid2 esc_nbgrid) {
    unsigned long i = 0, size, nc = 0;
    sf_esc_point2 plane = NULL;

    while ((plane = sf_esc_nband2_get_used_plane (esc_nbgrid->nband, i, &size))) {
        sf_esc_nbout2_write_plane (esc_nbgrid->esc_out, plane, size);
        nc += size;
        i++;
    }

    return nc;
}

/* Remove previous generation of parents, turn present
   generation of child points into parents */
static bool sf_esc_nbgrid2_create_new_generation (sf_esc_nbgrid2 esc_nbgrid) {
    int iz, ix, ia, i = 0;
    float a, s, sa, sz, sx, fz, fx, fa;
    sf_esc_point2 point = NULL;

    if (sf_esc_nband2_next_generation (esc_nbgrid->nband)) {
        while ((point = (sf_esc_point2)sf_esc_nband2_get_child (esc_nbgrid->nband, i,
                                                                &iz, &ix, &ia))) {
            i++;
            a = esc_nbgrid->oa + ia*esc_nbgrid->da;
            sf_esc_slowness2_get_components (esc_nbgrid->esc_slow,
                                             esc_nbgrid->oz + iz*esc_nbgrid->dz,
                                             esc_nbgrid->ox + ix*esc_nbgrid->dx,
                                             a, &s, &sa, &sz, &sx);
            sf_esc_slowness2_get_coefs (esc_nbgrid->esc_slow, a, s, sa, sz, sx,
                                        &fz, &fx, &fa); 
            if (fz != 0.0)
                sf_esc_point2_add_parent_link (point, ESC2_AXIS_Z,
                                               (fz < 0.0) ? ESC2_FORW : ESC2_BACK); 
            if (fx != 0.0)
                sf_esc_point2_add_parent_link (point, ESC2_AXIS_X,
                                               (fx < 0.0) ? ESC2_FORW : ESC2_BACK); 
            if (fa != 0.0)
                sf_esc_point2_add_parent_link (point, ESC2_AXIS_A,
                                               (fa < 0.0) ? ESC2_FORW : ESC2_BACK);
        }
        return true;
    }

    return false;
}

void sf_esc_nbgrid2_compute (sf_esc_nbgrid2 esc_nbgrid)
/*< Run escape values computations  >*/
{
    int i = 0, sz = sf_esc_point2_sizeof ();
    unsigned long nc, nt = 0, nr = 0, nr0, nr1;

    if (esc_nbgrid->verb) {
        sf_warning ("na = %d, oa = %g", esc_nbgrid->na, esc_nbgrid->oa*180.0/SF_PI);
        sf_warning ("Size of one entry: %d", sz);
    }

    nc = sf_esc_nband2_point_count (esc_nbgrid->nband);

    if (esc_nbgrid->verb) {
        sf_warning ("%lu points in the band", nc);
        sf_warning ("%g Mb occupied", 1e-6*nc*sf_esc_point2_sizeof ());
        sf_timer_start (esc_nbgrid->ttime);
    }

    do {
        nr0 = esc_nbgrid->nlt + esc_nbgrid->pgt + esc_nbgrid->nbt;
        nc = sf_esc_nbgrid2_evaluate_children (esc_nbgrid);
        nr1 = esc_nbgrid->nlt + esc_nbgrid->pgt + esc_nbgrid->nbt;
        if (esc_nbgrid->verb)
            sf_warning ("Iteration %d: %lu points computed (%g %% traced)",
                        i + 1, nc, (nr1 - nr0)/(float)nc*100.0);
        nt += nc;
        nc = sf_esc_nbgrid2_save_used_points (esc_nbgrid);
        nr += nc;
        if (esc_nbgrid->verb)
            sf_warning ("Iteration %d: %g %% completed", i + 1,
                        nr/(float)((size_t)esc_nbgrid->na*(size_t)esc_nbgrid->nx*
                                   (size_t)esc_nbgrid->nz)*100.0);
        i++;
    } while (sf_esc_nbgrid2_create_new_generation (esc_nbgrid));
    if (esc_nbgrid->verb) {
        sf_warning ("Total: %lu points computed, %lu points removed", nt, nr);
        sf_warning ("Total: %lu points had colors mixed (%g %%)",
                    esc_nbgrid->nct, esc_nbgrid->nct/(float)nr*100.0);
        sf_warning ("Total: %lu points exceeded maximum exit spread (%g %%)",
                    esc_nbgrid->ndt, esc_nbgrid->ndt/(float)nr*100.0);
        sf_warning ("Total: %lu points had F-D overshooting (%g %%)",
                    esc_nbgrid->fdt, esc_nbgrid->fdt/(float)nr*100.0);
        sf_warning ("Total: %lu points had deadlocks (%g %%)",
                    esc_nbgrid->nlt, esc_nbgrid->nlt/(float)nr*100.0);
        sf_warning ("Total: %lu points used for boundary conditions (%g %%)",
                    esc_nbgrid->bct, esc_nbgrid->bct/(float)nr*100.0);
        sf_warning ("Total: %lu points ray-traced for narrow band edges (%g %%)",
                    esc_nbgrid->nbt, esc_nbgrid->nbt/(float)nr*100.0);
        sf_warning ("Total: %lu points ray-traced for phase/group velocity deviation (%g %%)",
                    esc_nbgrid->pgt, esc_nbgrid->pgt/(float)nr*100.0);
        nc = esc_nbgrid->pgt + esc_nbgrid->nbt;
        sf_warning ("Total: %lu points ray-traced (%g %%)", nc, nc/(float)nr*100.0);
        sf_timer_stop (esc_nbgrid->ttime);
        sf_warning ("Total runtime: %g mins",
                    sf_timer_get_total_time (esc_nbgrid->ttime)/60000.0);
        sf_warning ("Ray-tracing runtime: %g mins (%g %%)",
                    sf_timer_get_total_time (esc_nbgrid->rtime)/60000.0,
                    100.0*sf_timer_get_total_time (esc_nbgrid->rtime)/
                          sf_timer_get_total_time (esc_nbgrid->ttime));
    }
}


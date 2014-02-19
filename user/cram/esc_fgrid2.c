/* 3-D full phase-space grid for 2-D medium, uses finite differences with Gauss-Seidel updates */
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

#ifndef _esc_fgrid2_h

#include "esc_slow2.h"
/*^*/

#include "esc_fout2.h"
/*^*/

#include "esc_tracer2.h"
/*^*/

#include "esc_nbgrid2.h"

typedef struct EscFGrid2 *sf_esc_fgrid2;
/* abstract data type */
/*^*/

#endif

struct EscFGrid2 {
    int                  nz, nx, na, morder, niter;
    float                oz, ox, oa;
    float                dz, dx, da;
    float                gzmin, gxmin, gzmax, gxmax;
    float                d[ESC2_DIMS], mdist, thresh;
    unsigned long        nct, ndt, bct;
    bool                 verb, cmix, tracebc, atraced, mtraced;
    sf_timer             ttime, rtime;
    sf_esc_slowness2     esc_slow;
    sf_esc_tracer2       esc_tracer;
    sf_esc_fout2         esc_out;
};
/* concrete data type */

typedef enum { ESC2_FGRID_FORW = 0, ESC2_FGRID_BACK = 1 } sf_esc_fgrid2_dir;

void sf_esc_fgrid2_set_tracebc (sf_esc_fgrid2 esc_fgrid, bool tracebc)
/*< Set trace boundary points flag >*/
{
    esc_fgrid->tracebc = tracebc;
    if (!tracebc) {
        esc_fgrid->gzmin = esc_fgrid->oz;
        esc_fgrid->gzmax = esc_fgrid->oz + (esc_fgrid->nz - 1)*esc_fgrid->dz;
        esc_fgrid->gxmin = esc_fgrid->ox;
        esc_fgrid->gxmax = esc_fgrid->ox + (esc_fgrid->nx - 1)*esc_fgrid->dx;
    } else {
        esc_fgrid->gzmin = sf_esc_slowness2_oz (esc_fgrid->esc_slow);
        esc_fgrid->gzmax = esc_fgrid->gzmin +
                          (sf_esc_slowness2_nz (esc_fgrid->esc_slow) - 1)*sf_esc_slowness2_dz (esc_fgrid->esc_slow);
        esc_fgrid->gxmin = sf_esc_slowness2_ox (esc_fgrid->esc_slow);
        esc_fgrid->gxmax = esc_fgrid->gxmin +
                          (sf_esc_slowness2_nx (esc_fgrid->esc_slow) - 1)*sf_esc_slowness2_dx (esc_fgrid->esc_slow);
    }
    sf_esc_tracer2_set_zmin (esc_fgrid->esc_tracer, esc_fgrid->gzmin);
    sf_esc_tracer2_set_zmax (esc_fgrid->esc_tracer, esc_fgrid->gzmax);
    sf_esc_tracer2_set_xmin (esc_fgrid->esc_tracer, esc_fgrid->gxmin);
    sf_esc_tracer2_set_xmax (esc_fgrid->esc_tracer, esc_fgrid->gxmax);
}

sf_esc_fgrid2 sf_esc_fgrid2_init (int nz, int nx, int na,
                                  float oz, float ox, float dz, float dx,
                                  bool atraced, bool mtraced,
                                  sf_esc_slowness2 esc_slow, sf_esc_tracer2 esc_tracer,
                                  sf_esc_fout2 esc_out)
/*< Initialize object >*/
{
    sf_esc_fgrid2 esc_fgrid = (sf_esc_fgrid2)sf_alloc (1, sizeof (struct EscFGrid2));

    esc_fgrid->nz = nz;
    esc_fgrid->nx = nx;
    esc_fgrid->na = na;

    esc_fgrid->dz = dz;
    esc_fgrid->d[ESC2_AXIS_Z] = dz;
    esc_fgrid->dx = dx;
    esc_fgrid->d[ESC2_AXIS_X] = dx;
    esc_fgrid->da = sf_esc_nbgrid2_get_da (na);
    esc_fgrid->d[ESC2_AXIS_A] = esc_fgrid->da;

    esc_fgrid->oz = oz;
    esc_fgrid->ox = ox;
    esc_fgrid->oa = sf_esc_nbgrid2_get_oa (na);

    esc_fgrid->thresh = 1e-6;
    esc_fgrid->morder = ESC2_MORDER;
    esc_fgrid->niter = 100;
    esc_fgrid->cmix = false;
    esc_fgrid->atraced = atraced;
    esc_fgrid->mtraced = mtraced;

    esc_fgrid->nct = 0; /* Number of ray traced because of color mixing */
    esc_fgrid->ndt = 0; /* Number of ray traced because of maximum allowed exit distance spread */
    esc_fgrid->bct = 0; /* Number of ray traced for the boundary conditions */

    esc_fgrid->verb = false;

    esc_fgrid->ttime = sf_timer_init (); /* Total time */
    esc_fgrid->rtime = sf_timer_init (); /* Ray tracing time */

    esc_fgrid->esc_slow = esc_slow;
    esc_fgrid->esc_tracer = esc_tracer;
    esc_fgrid->esc_out = esc_out;

    sf_esc_fgrid2_set_tracebc (esc_fgrid, true);

    return esc_fgrid;
}

void sf_esc_fgrid2_close (sf_esc_fgrid2 esc_fgrid)
/*< Destroy object >*/
{
    sf_timer_close (esc_fgrid->ttime);
    sf_timer_close (esc_fgrid->rtime);
    free (esc_fgrid);
}

void sf_esc_fgrid2_set_verb (sf_esc_fgrid2 esc_fgrid, bool verb)
/*< Set verbatim flag >*/
{
    esc_fgrid->verb = verb;
}

void sf_esc_fgrid2_set_threshold (sf_esc_fgrid2 esc_fgrid, float thresh)
/*< Set convergence threshold for Gauss-Seidel iterations >*/
{
    esc_fgrid->thresh = thresh;
}

void sf_esc_fgrid2_set_mdist (sf_esc_fgrid2 esc_fgrid, float mdist)
/*< Set maximum distance between escape values in F-D stencil >*/
{
    esc_fgrid->mdist = mdist;
}

void sf_esc_fgrid2_set_morder (sf_esc_fgrid2 esc_fgrid, int morder)
/*< Set maximum order in F-D stencil >*/
{
    if (morder < 1)
        morder = 1;
    if (morder > ESC2_MORDER)
        morder = ESC2_MORDER;
    esc_fgrid->morder = morder;
}

void sf_esc_fgrid2_set_niter (sf_esc_fgrid2 esc_fgrid, float niter)
/*< Set maximum number of Gauss-Seidel iterations >*/
{
    esc_fgrid->niter = niter;
}

void sf_esc_fgrid2_set_cmix (sf_esc_fgrid2 esc_fgrid, bool cmix)
/*< Set color mixing check flag >*/
{
    esc_fgrid->cmix = cmix;
}

/* Do simple Gauss-Seidel update for one point and one type of variable */
static float sf_esc_fgrid2_gs_update (sf_esc_fgrid2 esc_fgrid,
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
static void sf_esc_fgrid2_compute_point (sf_esc_fgrid2 esc_fgrid, int iz, int ix, int ia,
                                         float fz, float fx, float fa, float s, sf_esc_point2 point,
                                         sf_esc_point2 neighbors[ESC2_DIMS][ESC2_MORDER]) {
    int j, k, l;
    float val, rhs, e[ESC2_DIMS], n[ESC2_DIMS], f[ESC2_DIMS][ESC2_MORDER];
    EscColor2 col = 0;

    /* Coefficients of the escape equations */
    e[ESC2_AXIS_Z] = fz;
    e[ESC2_AXIS_X] = fx;
    e[ESC2_AXIS_A] = fa;

    for (j = 0; j < ESC2_DIMS; j++) {
        e[j] = (neighbors[j][0] != NULL) ? fabsf (e[j]/esc_fgrid->d[j])
                                         : 0.0;
    }

    /* Apply regular finite-difference for all escape variables */
    for (k = 0; k < ESC2_NUM; k++) {
        for (j = 0; j < ESC2_DIMS; j++) {
            n[j] = 0;
            for (l = 0; (l < ESC2_MORDER && l < esc_fgrid->morder); l++) {
                if (neighbors[j][l]) {
                    f[j][l] = sf_esc_point2_get_esc_var (neighbors[j][l], k);
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

        val = sf_esc_fgrid2_gs_update (esc_fgrid, n, e, f, rhs);
        /* Check limits */
        if (ESC2_Z == k) {
            if (val < esc_fgrid->gzmin)
                val = esc_fgrid->gzmin;
            if (val > esc_fgrid->gzmax)
                val = esc_fgrid->gzmax;
        } else if (ESC2_X == k) {
            if (val < esc_fgrid->gxmin)
                val = esc_fgrid->gxmin;
            if (val > esc_fgrid->gxmax)
                val = esc_fgrid->gxmax;
        }
        sf_esc_point2_set_esc_var (point, k, val);
    }

    for (j = 0; j < ESC2_DIMS; j++) {
        for (k = 0; k < ESC2_MORDER; k++) {
            if (neighbors[j][k])
                col |= sf_esc_point2_get_col (neighbors[j][k]);
        }
    }

    sf_esc_point2_set_col (point, col);
}

/* Returns true, if point belongs to boundary conditions */
static bool sf_esc_fgrid2_is_bc (sf_esc_fgrid2 esc_fgrid, int iz, int ix,
                                 float fz, float fx, EscColor2 *col) {
    EscColor2 color = 0;

    if (0 == iz)
        color |= ESC2_TOP;
    if ((esc_fgrid->nz - 1) == iz)
        color |= ESC2_BOTTOM;
    if (0 == ix)
        color |= ESC2_LEFT;
    if ((esc_fgrid->nx - 1) == ix)
        color |= ESC2_RIGHT;

    *col = color;

    return ((color & ESC2_TOP) && fz > 0.0) ||
           ((color & ESC2_BOTTOM) && fz < 0.0) ||
           ((color & ESC2_LEFT) && fx > 0.0) ||
           ((color & ESC2_RIGHT) && fx < 0.0);
}

/* Evaluate whether a point has to be computed (and how - ray tracing or F-D */
static void sf_esc_fgrid2_evaluate_point (sf_esc_fgrid2 esc_fgrid, int iz, int ix, int ia,
                                          double *cvt) {
    int i, j, k, l, m;
    float z, x, a, s, sa, sz, sx, fz, fx, fa;
    bool trace = false, hasn[ESC2_DIMS] = { false, false, false }, trmdist = false;
    sf_esc_point2 neighbors[ESC2_DIMS][ESC2_MORDER] = { { NULL, NULL },
                                                        { NULL, NULL },
                                                        { NULL, NULL } };
    sf_esc_point2 point;
    EscColor2 col1, col2;
    float x1, x2;

    point = sf_esc_fout2_get_point (esc_fgrid->esc_out, iz, ix, ia);

    if (sf_esc_point2_is_child (point)) {
        z = esc_fgrid->oz + iz*esc_fgrid->dz;
        x = esc_fgrid->ox + ix*esc_fgrid->dx;
        a = esc_fgrid->oa + ia*esc_fgrid->da;
        sf_esc_slowness2_get_components (esc_fgrid->esc_slow, z, x, a,
                                         &s, &sa, &sz, &sx);
        sf_esc_slowness2_get_coefs (esc_fgrid->esc_slow, a, s, sa, sz, sx,
                                    &fz, &fx, &fa);
        hasn[ESC2_AXIS_Z] = fz != 0.0;
        hasn[ESC2_AXIS_X] = fx != 0.0;
        hasn[ESC2_AXIS_A] = fa != 0.0;
        /* Collect neighbors for F-D stencil */
        if (hasn[ESC2_AXIS_Z]) {
            neighbors[ESC2_AXIS_Z][0] = sf_esc_fout2_get_point (esc_fgrid->esc_out,
                                                                fz < 0.0 ? iz + 1 : iz - 1, ix, ia);
            if (esc_fgrid->morder > 1)
                neighbors[ESC2_AXIS_Z][1] = sf_esc_fout2_get_point (esc_fgrid->esc_out,
                                                                    fz < 0.0 ? iz + 2 : iz - 2, ix, ia);
        }
        if (hasn[ESC2_AXIS_X]) {
            neighbors[ESC2_AXIS_X][0] = sf_esc_fout2_get_point (esc_fgrid->esc_out,
                                                                iz, fx < 0.0 ? ix + 1 : ix - 1, ia);
            if (esc_fgrid->morder > 1)
                neighbors[ESC2_AXIS_X][1] = sf_esc_fout2_get_point (esc_fgrid->esc_out,
                                                                    iz, fx < 0.0 ? ix + 2 : ix - 2, ia);
        }
        if (hasn[ESC2_AXIS_A]) {
            neighbors[ESC2_AXIS_A][0] = sf_esc_fout2_get_point (esc_fgrid->esc_out,
                                                                iz, ix, fa < 0.0 ? ia + 1 : ia - 1);
            if (esc_fgrid->morder > 1)
                neighbors[ESC2_AXIS_A][1] = sf_esc_fout2_get_point (esc_fgrid->esc_out,
                                                                    iz, ix, fa < 0.0 ? ia + 2 : ia - 2);
        }

        /* Boundary condition? */
        if (sf_esc_fgrid2_is_bc (esc_fgrid, iz, ix, fz, fx, &col1)) {
            if (esc_fgrid->tracebc) {
                trace = true;
                esc_fgrid->bct++;
            } else {
                sf_esc_point2_set_esc_var (point, ESC2_Z, z);                
                sf_esc_point2_set_esc_var (point, ESC2_X, x);                
#ifdef ESC_EQ_WITH_L
                sf_esc_point2_set_esc_var (point, ESC2_L, 0.0);                
#endif
                sf_esc_point2_set_esc_var (point, ESC2_T, 0.0);                
                sf_esc_point2_set_col (point, col1);
                sf_esc_point2_become_parent (point);               
            }
        }

        if (!trace && (esc_fgrid->cmix || esc_fgrid->mdist != SF_HUGE)) {
            /* Check for color mixing */
            for (j = 0; (j < (ESC2_DIMS - 1) && !trace); j++) {
                if (hasn[j]) {
                    for (l = 0; (l < ESC2_MORDER && !trace && neighbors[j][l]); l++) {
                        x1 = sf_esc_point2_get_esc_var (neighbors[j][l], ESC2_X);
                        col1 = sf_esc_point2_get_col (neighbors[j][l]);
                        for (k = j + 1; (k < ESC2_DIMS && !trace); k++) {
                            if (hasn[k]) {
                                for (m = 0; (m < ESC2_MORDER && !trace && neighbors[k][m]); m++) {
                                    x2 = sf_esc_point2_get_esc_var (neighbors[k][m], ESC2_X);
                                    col2 = sf_esc_point2_get_col (neighbors[k][m]);
                                    /* Hitting corner */
                                    if (col1 != col2 && esc_fgrid->cmix) {
                                        trace = true;
                                        esc_fgrid->nct++;
                                        break;
                                    }
                                    /* Exceding maximum allowed distance between parents */
                                    if (ESC2_TOP == col1 && fabs (x1 - x2) > esc_fgrid->mdist) {
                                        trace = true;
                                        esc_fgrid->ndt++;
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
        
        if (trace) { /* Do explicit ray tracing */
            if (esc_fgrid->verb)
                sf_timer_start (esc_fgrid->rtime);
            sf_esc_tracer2_compute (esc_fgrid->esc_tracer,
                                    esc_fgrid->oz + iz*esc_fgrid->dz,
                                    esc_fgrid->ox + ix*esc_fgrid->dx,
                                    esc_fgrid->oa + ia*esc_fgrid->da,
                                    0.0, 0.0, point, NULL);
            if (esc_fgrid->verb)
                sf_timer_stop (esc_fgrid->rtime);
            sf_esc_point2_become_parent (point);
            if ((esc_fgrid->mtraced && trmdist) || esc_fgrid->atraced)
                sf_esc_point2_set_traced (point, true);
        } else if (sf_esc_point2_is_child (point)) { /* F-D */
            sf_esc_fgrid2_compute_point (esc_fgrid, iz, ix, ia, fz, fx, fa, s,
                                         point, neighbors);
        }
    }

    for (i = 0; i < ESC2_NUM; i++)
        cvt[i] += fabsf (sf_esc_point2_get_esc_var (point, i));
}

/* Compute one Gauss-Seidel iteration */
static void sf_esc_fgrid2_compute_one_iter (sf_esc_fgrid2 esc_fgrid, sf_esc_fgrid2_dir dir,
                                            double *cvt) {
    int iz, ix, ia, ix0, ix1, ixs;
    int iz0, iz1, izs, ia0, ia1, ias;
    float a, cs, sn;

    ia0 = (ESC2_FGRID_FORW == dir) ? 0 : esc_fgrid->na - 1;
    ias = (ESC2_FGRID_FORW == dir) ? 1 : -1;
    ia1 = (ESC2_FGRID_FORW == dir) ? esc_fgrid->na : -1;

    for (ia = ia0; ia != ia1; ia += ias) {
        a = esc_fgrid->oa + ia*esc_fgrid->da;
        cs = cosf(a); sn = sinf(a);
        if (fabsf (cs) == 1.0) sn = 0.0;
        if (fabsf (sn) == 1.0) cs = 0.0;
        /* Determine sweep direction from angle */
        if (sn > 0.0) { /* Left->Right */
            ix0 = 0; ix1 = esc_fgrid->nx; ixs = 1;
        } else { /* Right->Left */
            ix0 = esc_fgrid->nx - 1; ix1 = -1; ixs = -1;
        }
        if (cs > 0.0) { /* Up->Down */
            iz0 = 0; iz1 = esc_fgrid->nz; izs = 1;
        } else { /* Down->Up */
            iz0 = esc_fgrid->nz - 1; iz1 = -1; izs = -1;
        }
        for (ix = ix0; ix != ix1; ix += ixs) {
            for (iz = iz0; iz != iz1; iz += izs) {
                sf_esc_fgrid2_evaluate_point (esc_fgrid, iz, ix, ia, cvt);
            }
        }
    }
}

void sf_esc_fgrid2_compute (sf_esc_fgrid2 esc_fgrid)
/*< Run escape values computations  >*/
{
    int i, iter;
    size_t np = (size_t)esc_fgrid->na*(size_t)esc_fgrid->nx*(size_t)esc_fgrid->nz;
    double old_gcvt[ESC2_NUM], gcvt[ESC2_NUM];
    bool stop = true, cmix = esc_fgrid->cmix;
    float mdist = esc_fgrid->mdist;

    if (esc_fgrid->verb)
        sf_timer_start (esc_fgrid->ttime);

    /* Gauss-Seidel iterations until convergence or
       maximum number of iterations */
    esc_fgrid->mdist = SF_HUGE;
    esc_fgrid->cmix = false;
    for (i = 0; i < ESC2_NUM; i++) {
        old_gcvt[i] = SF_HUGE;
        gcvt[i] = 0.0;
    }
    /* Iterate without injected tracing first to obtain a starting solution */
    for (iter = 0; iter < esc_fgrid->niter; iter++) {
        sf_esc_fgrid2_compute_one_iter (esc_fgrid, ESC2_FGRID_FORW, gcvt);
        sf_esc_fgrid2_compute_one_iter (esc_fgrid, ESC2_FGRID_BACK, gcvt);
        for (i = 0; i < ESC2_NUM; i++) {
            gcvt[i] *= 0.5; /* Average between two angle iterations */
            gcvt[i] /= (double)np;
            old_gcvt[i] -= gcvt[i];
            if (esc_fgrid->verb)
                sf_warning ("Iteration %d: L1(%s)=%g, change=%g", iter + 1,
                            sf_esc_point2_str[i], gcvt[i], old_gcvt[i]);
            if (fabsf (old_gcvt[i]) > esc_fgrid->thresh)
                stop = false;
            old_gcvt[i] = gcvt[i];
            gcvt[i] = 0.0;
        }
        if (stop)
            break;
        stop = true;
    }
    if (cmix || mdist != SF_HUGE) {
        /* Iterate with injected tracing first to obtain final solution */
        for (i = 0; i < ESC2_NUM; i++) {
            old_gcvt[i] = SF_HUGE;
            gcvt[i] = 0.0;
        }
        esc_fgrid->cmix = cmix;
        esc_fgrid->mdist = mdist;
        for (iter = 0; iter < esc_fgrid->niter; iter++) {
            sf_esc_fgrid2_compute_one_iter (esc_fgrid, ESC2_FGRID_FORW, gcvt);
            sf_esc_fgrid2_compute_one_iter (esc_fgrid, ESC2_FGRID_BACK, gcvt);
            for (i = 0; i < ESC2_NUM; i++) {
                gcvt[i] *= 0.5; /* Average between two angle iterations */
                gcvt[i] /= (double)np;
                old_gcvt[i] -= gcvt[i];
                if (esc_fgrid->verb)
                    sf_warning ("Iteration %d: L1(%s)=%g, change=%g", iter + 1,
                                sf_esc_point2_str[i], gcvt[i], old_gcvt[i]);
                if (fabsf (old_gcvt[i]) > esc_fgrid->thresh)
                    stop = false;
                old_gcvt[i] = gcvt[i];
                gcvt[i] = 0.0;
            }
            if (stop)
                break;
            stop = true;
        }
    }

    if (esc_fgrid->verb) {
        if (esc_fgrid->cmix)
            sf_warning ("Total: %lu points had colors mixed (%g %%)",
                        esc_fgrid->nct, esc_fgrid->nct/(float)np*100.0);
        sf_warning ("Total: %lu points exceeded maximum exit spread (%g %%)",
                    esc_fgrid->ndt, esc_fgrid->ndt/(float)np*100.0);
        sf_warning ("Total: %lu points used for boundary conditions (%g %%)",
                    esc_fgrid->bct, esc_fgrid->bct/(float)np*100.0);
        sf_timer_stop (esc_fgrid->ttime);
        sf_warning ("Total runtime: %g mins",
                    sf_timer_get_total_time (esc_fgrid->ttime)/60000.0);
        sf_warning ("Ray-tracing runtime: %g mins (%g %%)",
                    sf_timer_get_total_time (esc_fgrid->rtime)/60000.0,
                    100.0*sf_timer_get_total_time (esc_fgrid->rtime)/
                          sf_timer_get_total_time (esc_fgrid->ttime));
    }
    if (esc_fgrid->verb)
        sf_warning ("Writing escape variables");
    sf_esc_fout2_write (esc_fgrid->esc_out);
}


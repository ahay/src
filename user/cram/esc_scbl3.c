/* 5-D phase-space supercell block */
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

#include "esc_scbl3.h"
#include "einspline.h"

#ifndef _esc_scbl3_h

#include "esc_bgrid3.h"
/*^*/

#include "esc_tracer3.h"
/*^*/

typedef struct EscSCBlock3 *sf_esc_scblock3;
/* abstract data type */
/*^*/

#endif

struct EscSCBlock3 {
    int                  nz, nx, ny, na, nb;
    float                oz, ox, oy, oa, ob;
    float                dz, dx, dy, da, db;
    float                zmin, zmax;
    float                xmin, xmax;
    float                ymin, ymax;
    size_t               offs, n;
    unsigned char       *mmaped;
    multi_UBspline_3d_s *scsplines[ESC3_SIDE_NUM];
    sf_esc_point3        esc_point;
    sf_esc_tracer3       esc_tracer;
};
/* concrete data type */

sf_esc_scblock3 sf_esc_scblock3_init (sf_file scblock, sf_esc_tracer3 esc_tracer)
/*< Initialize object >*/
{
    int i, ia;
    size_t nc = 0, offs;
    FILE *stream;
    sf_esc_scblock3 esc_scbl = (sf_esc_scblock3)sf_alloc (1, sizeof (struct EscSCBlock3));

    if (sf_gettype (scblock) != SF_UCHAR)
        sf_error ("Wrong data type in supercell file");

    if (!sf_histlargeint (scblock, "n1", (off_t*)&esc_scbl->n)) sf_error ("No n1= in supercell file");

    if (!sf_histint (scblock, "Nz", &esc_scbl->nz)) sf_error ("No Nz= in supercell file");
    if (!sf_histfloat (scblock, "Dz", &esc_scbl->dz)) sf_error ("No Dz= in supercell file");
    if (!sf_histfloat (scblock, "Oz", &esc_scbl->oz)) sf_error ("No Oz= in supercell file");
    if (!sf_histint (scblock, "Nx", &esc_scbl->nx)) sf_error ("No Nx= in supercell file");
    if (!sf_histfloat (scblock, "Dx", &esc_scbl->dx)) sf_error ("No Dx= in supercell file");
    if (!sf_histfloat (scblock, "Ox", &esc_scbl->ox)) sf_error ("No Ox= in supercell file");
    if (!sf_histint (scblock, "Ny", &esc_scbl->ny)) sf_error ("No Ny= in supercell file");
    if (!sf_histfloat (scblock, "Dy", &esc_scbl->dy)) sf_error ("No Dy= in supercell file");
    if (!sf_histfloat (scblock, "Oy", &esc_scbl->oy)) sf_error ("No Oy= in supercell file");
    if (!sf_histint (scblock, "Na", &esc_scbl->na)) sf_error ("No Na= in supercell file");
    if (!sf_histint (scblock, "Nb", &esc_scbl->nb)) sf_error ("No Nb= in supercell file");

    esc_scbl->zmin = esc_scbl->oz;
    esc_scbl->zmax = esc_scbl->oz + (esc_scbl->nz - 1)*esc_scbl->dz;
    esc_scbl->xmin = esc_scbl->ox;
    esc_scbl->xmax = esc_scbl->ox + (esc_scbl->nx - 1)*esc_scbl->dx;
    esc_scbl->ymin = esc_scbl->oy;
    esc_scbl->ymax = esc_scbl->oy + (esc_scbl->ny - 1)*esc_scbl->dy;

    stream = sf_filestream (scblock);
    if (stdin == stream)
        sf_error ("Can not mmap supercell file from stdin");

    esc_scbl->offs = ftello (stream);

    /* Read spline structures */
    for (i = 0; i < ESC3_SIDE_NUM; i++) {
        esc_scbl->scsplines[i] = (multi_UBspline_3d_s*)sf_alloc ((size_t)esc_scbl->na,
                                                                 sizeof(multi_UBspline_3d_s));
        for (ia = 0; ia < esc_scbl->na; ia++) {
            sf_ucharread ((unsigned char*)&(esc_scbl->scsplines[i][ia]),
                          sizeof(multi_UBspline_3d_s), scblock);
        }
    }

    offs = ftello (stream);
    esc_scbl->mmaped = (unsigned char*)mmap (NULL, (size_t)esc_scbl->offs +
                                                   (size_t)esc_scbl->n,
                                             PROT_READ, MAP_SHARED,
                                             fileno (stream), 0);
    if (esc_scbl->mmaped == MAP_FAILED)
        sf_error ("Supercell spline coefficients mmap failed: %s", strerror (errno));
 
    /* Mmap spline coefficients, put correct pointers into the structures */
    for (i = 0; i < ESC3_SIDE_NUM; i++) {
        for (ia = 0; ia < esc_scbl->na; ia++) {
            esc_scbl->scsplines[i][ia].coefs = (float*)(esc_scbl->mmaped + offs);
            offs += (size_t)esc_scbl->scsplines[i][ia].nc;
        }
    }

    esc_scbl->oa = sf_esc_bgrid3_get_oa (esc_scbl->na);
    esc_scbl->da = sf_esc_bgrid3_get_da (esc_scbl->na);

    esc_scbl->esc_tracer = esc_tracer;
    esc_scbl->esc_point = sf_esc_point3_init ();

    return esc_scbl;
}

void sf_esc_scblock3_close (sf_esc_scblock3 esc_scbl)
/*< Destroy object >*/
{
    int i;
    munmap (esc_scbl->mmaped, (size_t)esc_scbl->offs +
                              (size_t)esc_scbl->n);
    for (i = 0; i < ESC3_SIDE_NUM; i++)
        free (esc_scbl->scsplines[i]);
    sf_esc_point3_close (esc_scbl->esc_point);
    free (esc_scbl);
}

#define SPEPS 1e-3

/* Return boundary (side of the supercell) index correpsonding to z,x,y coordinate */
static EscSide3 sf_esc_scblock3_xyz_side (sf_esc_scblock3 esc_scbl,
                                          float z, float x, float y) {
    if (z <= esc_scbl->zmin + SPEPS*esc_scbl->dz)
        return ESC3_SIDE_TOP;
    if (z >= esc_scbl->zmax - SPEPS*esc_scbl->dz)
        return ESC3_SIDE_BOTTOM;
    if (x <= esc_scbl->xmin + SPEPS*esc_scbl->dx)
        return ESC3_SIDE_LEFT;
    if (x >= esc_scbl->xmax - SPEPS*esc_scbl->dx)
        return ESC3_SIDE_RIGHT;
    if (y <= esc_scbl->ymin + SPEPS*esc_scbl->dy)
        return ESC3_SIDE_NEAR;
    if (y >= esc_scbl->ymax - SPEPS*esc_scbl->dy)
        return ESC3_SIDE_FAR;
    /* Not on any of the six sides */
    return ESC3_SIDE_NUM;
}

/* Convert normalized phase vector components into polar angles */
static void sf_esc_scblock3_p_to_ab (sf_esc_scblock3 esc_scbl,
                                     float pz, float px, float py,
                                     float *a, float *b) {
    if (pz < -1.0) pz = -1.0;
    if (pz > 1.0) pz = 1.0;
    if (px < -1.0) px = -1.0;
    if (px > 1.0) px = 1.0;
    if (py < -1.0) py = -1.0;
    if (py > 1.0) py = 1.0;

    *b = acosf (pz);
    *a = acosf (px);
    if (py < 0.0)
        *a = 2.0*SF_PI - *a;

    /* Make sure that a is in [0;2PI] and b is in [0;PI] */
    if (*a < 0.0)
        *a += 2.0*SF_PI;
    else if (*a > 2.0*SF_PI)
        *a -= 2.0*SF_PI;
    if (*b < 0.0) {
        *a += SF_PI;
        *b = fabsf (*b);
    } else if (*b > SF_PI) {
        *a -= SF_PI;
        *b = 2.0*SF_PI - *b;
    }  
    if (*a < 0.0)
        *a += 2.0*SF_PI;
    else if (*a > 2.0*SF_PI)
        *a -= 2.0*SF_PI;
}

static void sf_esc_scblock3_interpolate (sf_esc_scblock3 esc_scbl,
                                         float *z, float *x, float *y,
                                         float *t, float *l, float *a, float *b) {
    int i, ia, iap;
    float af;
    float vals[2][ESC3_NUM + 3], p[3];

    /* Make small shifts to account for the lack of boundary checks
       in the spline interpolation routine */
    if (*z <= esc_scbl->zmin)
        *z = esc_scbl->zmin + SPEPS*esc_scbl->dz;
    if (*z >= esc_scbl->zmax)
        *z = esc_scbl->zmax - SPEPS*esc_scbl->dz;
    if (*x <= esc_scbl->xmin)
        *x = esc_scbl->xmin + SPEPS*esc_scbl->dx;
    if (*x >= esc_scbl->xmax)
        *x = esc_scbl->xmax - SPEPS*esc_scbl->dx;
    if (*y <= esc_scbl->ymin)
        *y = esc_scbl->ymin + SPEPS*esc_scbl->dy;
    if (*y >= esc_scbl->ymax)
        *y = esc_scbl->ymax - SPEPS*esc_scbl->dy;

    /* Azimuth plane */
    af = (*a - esc_scbl->oa)/esc_scbl->da;
    ia = floorf (af);
    af -= (float)ia;

    /* Next azimuth plane */
    iap = ia + 1;
    if (ia == esc_scbl->na)
        iap = 0;

    /* Get values for the both planes */
    eval_multi_UBspline_3d_s (esc_scbl->scsplines[ia], *y, *x, *z, &vals[0][0]);
    eval_multi_UBspline_3d_s (esc_scbl->scsplines[iap], *y, *x, *z, &vals[1][0]);

    /* Interpolate phase components */
    for (i = 0; i < 3; i++)
        p[i] = (1.0 - af)*vals[0][ESC3_NUM + i] +
                       af*vals[1][ESC3_NUM + i];
    /* Convert to polar angles */
    sf_esc_scblock3_p_to_ab (esc_scbl, p[0], p[1], p[2], a, b);

    /* Do linear interpolation for the spatial components */
    *z = (1.0 - af)*vals[0][ESC3_Z] +
                 af*vals[1][ESC3_Z];
    *x = (1.0 - af)*vals[0][ESC3_X] +
                 af*vals[1][ESC3_X];
    *y = (1.0 - af)*vals[0][ESC3_Y] +
                 af*vals[1][ESC3_Y];
    *t += (1.0 - af)*vals[0][ESC3_T] +
                  af*vals[1][ESC3_T];
#ifdef ESC_EQ_WITH_L
    *l += (1.0 - af)*vals[0][ESC3_T] +
                  af*vals[1][ESC3_T];
#endif
}

EscSide3 sf_esc_scblock3_project_point (sf_esc_scblock3 esc_scbl,
                                        float *ze, float *xe, float *ye,
                                        float *te, float *le, float *be, float *ae,
                                        size_t *rt, double *rd)
/*< Project point (z, x, y, b, a) to the boundary of the supercell,
    return new interpolated exit values >*/
{
    float z = *ze, x = *xe, y = *ye, b = *be, a = *ae;
    float zp = z, xp = x, yp = y, t = *te, l = *le;
    EscSide3 side = sf_esc_scblock3_xyz_side (esc_scbl, z, x, y);

    if (z < (esc_scbl->zmin - SPEPS*esc_scbl->dz) ||
        z > (esc_scbl->zmax + SPEPS*esc_scbl->dz) ||
        x < (esc_scbl->xmin - SPEPS*esc_scbl->dx) ||
        x > (esc_scbl->xmax + SPEPS*esc_scbl->dx) ||
        y < (esc_scbl->ymin - SPEPS*esc_scbl->dy) ||
        y > (esc_scbl->ymax + SPEPS*esc_scbl->dy))
        sf_error ("sf_esc_scblock3_project_point: point is outside of this supercell");

    /* If point is on any of the side, do interpolation
       and get a new exit location */ 
    if (ESC3_SIDE_NUM != side)
        sf_esc_scblock3_interpolate (esc_scbl, &z, &x, &y, &t, &l, &a, &b);

    /* Check what side we are landed on */
    side = sf_esc_scblock3_xyz_side (esc_scbl, z, x, y);

    if (ESC3_SIDE_NUM == side) {
        /* Point is not on the boundary, project it;
           this might happen due to inaccuracies in the interpolation,
           or if the original point is not on the boundary */
        sf_esc_tracer3_set_zmin (esc_scbl->esc_tracer, esc_scbl->zmin);
        sf_esc_tracer3_set_zmax (esc_scbl->esc_tracer, esc_scbl->zmax);
        sf_esc_tracer3_set_xmin (esc_scbl->esc_tracer, esc_scbl->xmin);
        sf_esc_tracer3_set_xmax (esc_scbl->esc_tracer, esc_scbl->xmax);
        sf_esc_tracer3_set_ymin (esc_scbl->esc_tracer, esc_scbl->ymin);
        sf_esc_tracer3_set_ymax (esc_scbl->esc_tracer, esc_scbl->ymax);
        sf_esc_tracer3_compute (esc_scbl->esc_tracer, z, x, y, b, a,
                                t, l, esc_scbl->esc_point, ae, be);
        b = *be;
        a = *ae;
        z = sf_esc_point3_get_esc_var (esc_scbl->esc_point, ESC3_Z);
        x = sf_esc_point3_get_esc_var (esc_scbl->esc_point, ESC3_X);
        y = sf_esc_point3_get_esc_var (esc_scbl->esc_point, ESC3_Y);
        t = sf_esc_point3_get_esc_var (esc_scbl->esc_point, ESC3_T);
#ifdef ESC_EQ_WITH_L
        l = sf_esc_point3_get_esc_var (esc_scbl->esc_point, ESC3_L);
#endif
        side = sf_esc_scblock3_xyz_side (esc_scbl, z, x, y);
        if (rt)
            *rt += 1;
        /* Crude approximation to the distance covered by
           the ray tracing */
        if (rd)
            *rd += sqrtf ((z - zp)*(z - zp) +
                          (x - xp)*(x - xp) +
                          (y - yp)*(y - yp));
    }

    *ze = z;
    *xe = x;
    *ye = y;
    *te = t;
    *le = l;
    *ae = a;
    *be = b;

    return side;
}


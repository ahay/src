/* 5-D boundary phase-space grid for 3-D medium */
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

#ifndef _esc_bgrid3_h

#include "einspline.h"
/*^*/

#include "esc_point3.h"
/*^*/

#include "esc_tracer3.h"
/*^*/

typedef struct EscBGrid3 *sf_esc_bgrid3;
/* abstract data type */
/*^*/

#define ESC3_BSC_BMAX 128
#define ESC3_BGRID3_APRON 2
/*^*/

/* Colors for different boundaries */
typedef enum { ESC3_SIDE_TOP = 0, ESC3_SIDE_BOTTOM = 1,
               ESC3_SIDE_LEFT = 2, ESC3_SIDE_RIGHT = 3,
               ESC3_SIDE_NEAR = 4, ESC3_SIDE_FAR = 5,
               ESC3_SIDE_NUM = 6 } EscSide3;
/*^*/

#endif

struct EscBGrid3 {
    int                   nz, nx, ny, na, nb;
    float                 oz, ox, oy, oa, ob;
    float                 dz, dx, dy, da, db;
    float                 zmax, xmax, ymax;
    bool                  verb;
    sf_timer              ttime;
    sf_esc_tracer3        esc_tracer;
};
/* concrete data type */

float sf_esc_bgrid3_get_db (int nb)
/*< Returns sampling of the angular axis >*/
{
    return SF_PI/(float)nb;
}

float sf_esc_bgrid3_get_ob (int nb)
/*< Returns origin of the angular axis >*/
{
    return 0.5*sf_esc_bgrid3_get_db (nb);
}

float sf_esc_bgrid3_get_da (int na)
/*< Returns sampling of the angular axis >*/
{
    return 2.0*SF_PI/(float)na;
}

float sf_esc_bgrid3_get_oa (int na)
/*< Returns origin of the angular axis >*/
{
    return 0.5*sf_esc_bgrid3_get_da (na);
}

sf_esc_bgrid3 sf_esc_bgrid3_init (int nz, int nx, int ny, int na, int nb,
                                  float oz, float ox, float oy,
                                  float dz, float dx, float dy, 
                                  sf_esc_tracer3 esc_tracer)
/*< Initialize object >*/
{
    sf_esc_bgrid3 esc_bgrid = (sf_esc_bgrid3)sf_alloc (1, sizeof (struct EscBGrid3));

    esc_bgrid->nz = nz;
    esc_bgrid->nx = nx;
    esc_bgrid->ny = ny;
    esc_bgrid->na = na;
    esc_bgrid->nb = nb;

    esc_bgrid->dz = dz;
    esc_bgrid->dx = dx;
    esc_bgrid->dy = dy;
    esc_bgrid->da = sf_esc_bgrid3_get_da (na);
    esc_bgrid->db = sf_esc_bgrid3_get_db (nb);

    esc_bgrid->zmax = oz + (nz - 1)*dz;
    esc_bgrid->xmax = ox + (nx - 1)*dx;
    esc_bgrid->ymax = oy + (ny - 1)*dy;

    esc_bgrid->oz = oz;
    esc_bgrid->ox = ox;
    esc_bgrid->oy = oy;
    esc_bgrid->oa = sf_esc_bgrid3_get_oa (na);
    esc_bgrid->ob = sf_esc_bgrid3_get_ob (nb);

    esc_bgrid->verb = false;

    esc_bgrid->ttime = sf_timer_init (); /* Total time */

    esc_bgrid->esc_tracer = esc_tracer;

    return esc_bgrid;
}

void sf_esc_bgrid3_close (sf_esc_bgrid3 esc_bgrid)
/*< Destroy object >*/
{
    sf_timer_close (esc_bgrid->ttime);
    free (esc_bgrid);
}

void sf_esc_bgrid3_set_verb (sf_esc_bgrid3 esc_bgrid, bool verb)
/*< Set verbatim flag >*/
{
    esc_bgrid->verb = verb;
}

static void sf_esc_bgrid3_set_tracer_limits (sf_esc_bgrid3 esc_bgrid)
{
    sf_esc_tracer3_set_zmin (esc_bgrid->esc_tracer, esc_bgrid->oz);
    sf_esc_tracer3_set_zmax (esc_bgrid->esc_tracer, esc_bgrid->zmax);
    sf_esc_tracer3_set_xmin (esc_bgrid->esc_tracer, esc_bgrid->ox);
    sf_esc_tracer3_set_xmax (esc_bgrid->esc_tracer, esc_bgrid->xmax);
    sf_esc_tracer3_set_ymin (esc_bgrid->esc_tracer, esc_bgrid->oy);
    sf_esc_tracer3_set_ymax (esc_bgrid->esc_tracer, esc_bgrid->ymax);
}

static void sf_esc_bgrid3_free5 (float *****v)
{
    free (v[0][0][0][0]);
    free (v[0][0][0]);
    free (v[0][0]);
    free (v[0]);
    free (v);
}

multi_UBspline_3d_s** sf_esc_bgrid3_compute_topbottom (sf_esc_bgrid3 esc_bgrid, bool istop)
/*< Run escape values computations for the top or bottom boundary >*/
{
    int ix, iy, ib, ia, iap, i, j;
    float z, x, y, b, a, be, ae;
    float *****top;
    sf_esc_point3 esc_point;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s **splines;

    sf_esc_bgrid3_set_tracer_limits (esc_bgrid);
    top = sf_floatalloc5 (esc_bgrid->nx, esc_bgrid->ny, esc_bgrid->nb + 2*ESC3_BGRID3_APRON,
                          ESC3_NUM + 3, esc_bgrid->na);

    if (esc_bgrid->verb) {
         sf_timer_reset (esc_bgrid->ttime);
         sf_timer_start (esc_bgrid->ttime);
    }

    esc_point = sf_esc_point3_init ();

    z = istop ? esc_bgrid->oz : esc_bgrid->zmax;
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        a = esc_bgrid->oa + ia*esc_bgrid->da;
        for (ib = 0; ib < esc_bgrid->nb; ib++) {
            b = esc_bgrid->ob + ib*esc_bgrid->db;
            for (iy = 0; iy < esc_bgrid->ny; iy++) {
                y = esc_bgrid->oy + iy*esc_bgrid->dy;
                for (ix = 0; ix < esc_bgrid->nx; ix++) {
                    x = esc_bgrid->ox + ix*esc_bgrid->dx;
                    sf_esc_tracer3_compute (esc_bgrid->esc_tracer, z, x, y, b, a,
                                            0.0, 0.0, esc_point, &ae, &be);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3_NUM; i++)
                        top[ia][i][ib + ESC3_BGRID3_APRON][iy][ix] =
                        sf_esc_point3_get_esc_var (esc_point, i);
                    top[ia][ESC3_NUM][ib + ESC3_BGRID3_APRON][iy][ix] = cosf (be);
                    top[ia][ESC3_NUM + 1][ib + ESC3_BGRID3_APRON][iy][ix] = /*sinf (be)*/cosf (ae);
                    top[ia][ESC3_NUM + 2][ib + ESC3_BGRID3_APRON][iy][ix] = /*sinf (be)*/sinf (ae);
                }
            }
        }
    }
    /* Copy boundary values in b direction to the apron area -
       this will later help with the proper boundary conditions
       during interpolation */
    for (i = 0; i < (ESC3_NUM + 3); i++) {
        for (ia = 0; ia < esc_bgrid->na; ia++) {
            /* Boundary condition on b - as we go past 0 or PI,
               we add +/-PI to a */
            if (ia < esc_bgrid->na/2)
                iap = esc_bgrid->na/2 + ia;
            else
                iap = ia - esc_bgrid->na/2;
            for (iy = 0; iy < esc_bgrid->ny; iy++) {
                for (ix = 0; ix < esc_bgrid->nx; ix++) {
                    for (j = 0; j < ESC3_BGRID3_APRON; j++) {
                        top[ia][i][ESC3_BGRID3_APRON - 1 - j][iy][ix] =
                        top[iap][i][ESC3_BGRID3_APRON + j][iy][ix];
                        top[ia][i][esc_bgrid->nb + ESC3_BGRID3_APRON + j][iy][ix] =
                        top[iap][i][esc_bgrid->nb + ESC3_BGRID3_APRON - 1 - j][iy][ix];
                    }
                }
            }
        }
    }

    splines = (multi_UBspline_3d_s**)sf_alloc ((size_t)esc_bgrid->na, sizeof(multi_UBspline_3d_s*));
    z_grid.start = esc_bgrid->ox;
    z_grid.end = esc_bgrid->ox + (esc_bgrid->nx - 1)*esc_bgrid->dx;
    z_grid.num = esc_bgrid->nx;
    x_grid.start = esc_bgrid->oy;
    x_grid.end = esc_bgrid->oy + (esc_bgrid->ny - 1)*esc_bgrid->dy;
    x_grid.num = esc_bgrid->ny;
    y_grid.start = esc_bgrid->ob - ESC3_BGRID3_APRON*esc_bgrid->db;
    y_grid.end = esc_bgrid->ob + (ESC3_BGRID3_APRON + esc_bgrid->nb - 1)*esc_bgrid->db;
    y_grid.num = esc_bgrid->nb + 2*ESC3_BGRID3_APRON;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    yBC.lCode = yBC.rCode = NATURAL;
    /* Create splines for each constant azimuth hyperplane */
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        splines[ia] = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3_NUM + 3);
        for (i = 0; i < (ESC3_NUM + 3); i++) {
            set_multi_UBspline_3d_s (splines[ia], i, &top[ia][i][0][0][0]);
        }
    }
    sf_esc_bgrid3_free5 (top);

    if (esc_bgrid->verb) {
        sf_timer_stop (esc_bgrid->ttime);
        sf_warning ("Runtime: %g mins",
                   sf_timer_get_total_time (esc_bgrid->ttime)/60000.0);
    }

    sf_esc_point3_close (esc_point);

    return splines;
}

multi_UBspline_3d_s** sf_esc_bgrid3_compute_leftright (sf_esc_bgrid3 esc_bgrid, bool isleft)
/*< Run escape values computations for the left of right boundary >*/
{
    int iz, iy, ib, ia, iap, i, j;
    float z, x, y, b, a, be, ae;
    float *****left;
    sf_esc_point3 esc_point;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s **splines;

    sf_esc_bgrid3_set_tracer_limits (esc_bgrid);
    left = sf_floatalloc5 (esc_bgrid->nz, esc_bgrid->ny, esc_bgrid->nb + 2*ESC3_BGRID3_APRON,
                           ESC3_NUM + 3, esc_bgrid->na);

    if (esc_bgrid->verb) {
         sf_timer_reset (esc_bgrid->ttime);
         sf_timer_start (esc_bgrid->ttime);
    }

    esc_point = sf_esc_point3_init ();

    x = isleft ? esc_bgrid->ox : esc_bgrid->xmax;
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        a = esc_bgrid->oa + ia*esc_bgrid->da;
        for (ib = 0; ib < esc_bgrid->nb; ib++) {
            b = esc_bgrid->ob + ib*esc_bgrid->db;
            for (iy = 0; iy < esc_bgrid->ny; iy++) {
                y = esc_bgrid->oy + iy*esc_bgrid->dy;
                for (iz = 0; iz < esc_bgrid->nz; iz++) {
                    z = esc_bgrid->oz + iz*esc_bgrid->dz;
                    sf_esc_tracer3_compute (esc_bgrid->esc_tracer, z, x, y, b, a,
                                            0.0, 0.0, esc_point, &ae, &be);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3_NUM; i++)
                        left[ia][i][ib + ESC3_BGRID3_APRON][iy][iz] =
                        sf_esc_point3_get_esc_var (esc_point, i);
                    left[ia][ESC3_NUM][ib + ESC3_BGRID3_APRON][iy][iz] = cosf (be);
                    left[ia][ESC3_NUM + 1][ib + ESC3_BGRID3_APRON][iy][iz] = /*sinf (be)*/cosf (ae);
                    left[ia][ESC3_NUM + 2][ib + ESC3_BGRID3_APRON][iy][iz] = /*sinf (be)*/sinf (ae);
                }
            }
        }
    }
    /* Copy boundary values in b direction to the apron area -
       this will later help with the proper boundary conditions
       during interpolation */
    for (i = 0; i < (ESC3_NUM + 3); i++) {
        for (ia = 0; ia < esc_bgrid->na; ia++) {
            /* Boundary condition on b - as we go past 0 or PI,
               we add +/-PI to a */
            if (ia < esc_bgrid->na/2)
                iap = esc_bgrid->na/2 + ia;
            else
                iap = ia - esc_bgrid->na/2;
            for (iy = 0; iy < esc_bgrid->ny; iy++) {
                for (iz = 0; iz < esc_bgrid->nz; iz++) {
                    for (j = 0; j < ESC3_BGRID3_APRON; j++) {
                        left[ia][i][ESC3_BGRID3_APRON - 1 - j][iy][iz] =
                        left[iap][i][ESC3_BGRID3_APRON + j][iy][iz];
                        left[ia][i][esc_bgrid->nb + ESC3_BGRID3_APRON + j][iy][iz] =
                        left[iap][i][esc_bgrid->nb + ESC3_BGRID3_APRON - 1 - j][iy][iz];
                    }
                }
            }
        }
    }

    splines = (multi_UBspline_3d_s**)sf_alloc ((size_t)esc_bgrid->na, sizeof(multi_UBspline_3d_s*));
    z_grid.start = esc_bgrid->oz;
    z_grid.end = esc_bgrid->oz + (esc_bgrid->nz - 1)*esc_bgrid->dz;
    z_grid.num = esc_bgrid->nz;
    x_grid.start = esc_bgrid->oy;
    x_grid.end = esc_bgrid->oy + (esc_bgrid->ny - 1)*esc_bgrid->dy;
    x_grid.num = esc_bgrid->ny;
    y_grid.start = esc_bgrid->ob - ESC3_BGRID3_APRON*esc_bgrid->db;
    y_grid.end = esc_bgrid->ob + (ESC3_BGRID3_APRON + esc_bgrid->nb - 1)*esc_bgrid->db;
    y_grid.num = esc_bgrid->nb + 2*ESC3_BGRID3_APRON;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    yBC.lCode = yBC.rCode = NATURAL;
    /* Create splines for each constant azimuth hyperplane */
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        splines[ia] = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3_NUM + 3);
        for (i = 0; i < (ESC3_NUM + 3); i++) {
            set_multi_UBspline_3d_s (splines[ia], i, &left[ia][i][0][0][0]);
        }
    }
    sf_esc_bgrid3_free5 (left);
    if (esc_bgrid->verb) {
        sf_timer_stop (esc_bgrid->ttime);
        sf_warning ("Runtime: %g mins",
                   sf_timer_get_total_time (esc_bgrid->ttime)/60000.0);
    }

    sf_esc_point3_close (esc_point);

    return splines;
}

multi_UBspline_3d_s** sf_esc_bgrid3_compute_nearfar (sf_esc_bgrid3 esc_bgrid, bool isnear)
/*< Run escape values computations for the near of far boundary >*/
{
    int iz, ix, ib, ia, iap, i, j;
    float z, x, y, b, a, be, ae;
    float *****near;
    sf_esc_point3 esc_point;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s **splines;

    sf_esc_bgrid3_set_tracer_limits (esc_bgrid);
    near = sf_floatalloc5 (esc_bgrid->nz, esc_bgrid->nx, esc_bgrid->nb + 2*ESC3_BGRID3_APRON,
                           ESC3_NUM + 3, esc_bgrid->na);

    if (esc_bgrid->verb) {
         sf_timer_reset (esc_bgrid->ttime);
         sf_timer_start (esc_bgrid->ttime);
    }

    esc_point = sf_esc_point3_init ();

    y = isnear ? esc_bgrid->oy : esc_bgrid->ymax;
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        a = esc_bgrid->oa + ia*esc_bgrid->da;
        for (ib = 0; ib < esc_bgrid->nb; ib++) {
            b = esc_bgrid->ob + ib*esc_bgrid->db;
            for (ix = 0; ix < esc_bgrid->nx; ix++) {
                x = esc_bgrid->ox + ix*esc_bgrid->dx;
                for (iz = 0; iz < esc_bgrid->nz; iz++) {
                    z = esc_bgrid->oz + iz*esc_bgrid->dz;
                    sf_esc_tracer3_compute (esc_bgrid->esc_tracer, z, x, y, b, a,
                                            0.0, 0.0, esc_point, &ae, &be);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3_NUM; i++)
                        near[ia][i][ib + ESC3_BGRID3_APRON][ix][iz] =
                        sf_esc_point3_get_esc_var (esc_point, i);
                    near[ia][ESC3_NUM][ib + ESC3_BGRID3_APRON][ix][iz] = cosf (be);
                    near[ia][ESC3_NUM + 1][ib + ESC3_BGRID3_APRON][ix][iz] = /*sinf (be)*/cosf (ae);
                    near[ia][ESC3_NUM + 2][ib + ESC3_BGRID3_APRON][ix][iz] = /*sinf (be)*/sinf (ae);
                }
            }
        }
    }
    /* Copy boundary values in b direction to the apron area -
       this will later help with the proper boundary conditions
       during interpolation */
    for (i = 0; i < (ESC3_NUM + 3); i++) {
        for (ia = 0; ia < esc_bgrid->na; ia++) {
            /* Boundary condition on b - as we go past 0 or PI,
               we add +/-PI to a */
            if (ia < esc_bgrid->na/2)
                iap = esc_bgrid->na/2 + ia;
            else
                iap = ia - esc_bgrid->na/2;
            for (ix = 0; ix < esc_bgrid->nx; ix++) {
                for (iz = 0; iz < esc_bgrid->nz; iz++) {
                    for (j = 0; j < ESC3_BGRID3_APRON; j++) {
                        near[ia][i][ESC3_BGRID3_APRON - 1 - j][ix][iz] =
                        near[iap][i][ESC3_BGRID3_APRON + j][ix][iz];
                        near[ia][i][esc_bgrid->nb + ESC3_BGRID3_APRON + j][ix][iz] =
                        near[iap][i][esc_bgrid->nb + ESC3_BGRID3_APRON - 1 - j][ix][iz];
                    }
                }
            }
        }
    }

    splines = (multi_UBspline_3d_s**)sf_alloc ((size_t)esc_bgrid->na, sizeof(multi_UBspline_3d_s*));
    z_grid.start = esc_bgrid->oz;
    z_grid.end = esc_bgrid->oz + (esc_bgrid->nz - 1)*esc_bgrid->dz;
    z_grid.num = esc_bgrid->nz;
    x_grid.start = esc_bgrid->ox;
    x_grid.end = esc_bgrid->ox + (esc_bgrid->nx - 1)*esc_bgrid->dx;
    x_grid.num = esc_bgrid->nx;
    y_grid.start = esc_bgrid->ob - ESC3_BGRID3_APRON*esc_bgrid->db;
    y_grid.end = esc_bgrid->ob + (ESC3_BGRID3_APRON + esc_bgrid->nb - 1)*esc_bgrid->db;
    y_grid.num = esc_bgrid->nb + 2*ESC3_BGRID3_APRON;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    yBC.lCode = yBC.rCode = NATURAL;
    /* Create splines for each constant azimuth hyperplane */
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        splines[ia] = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3_NUM + 3);
        for (i = 0; i < (ESC3_NUM + 3); i++) {
            set_multi_UBspline_3d_s (splines[ia], i, &near[ia][i][0][0][0]);
        }
    }
    sf_esc_bgrid3_free5 (near);

    if (esc_bgrid->verb) {
        sf_timer_stop (esc_bgrid->ttime);
        sf_warning ("Runtime: %g mins",
                   sf_timer_get_total_time (esc_bgrid->ttime)/60000.0);
    }

    sf_esc_point3_close (esc_point);

    return splines;
}


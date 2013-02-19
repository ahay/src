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

#include "esc_slow3.h"
/*^*/

#include "esc_point3.h"
/*^*/

#include "esc_tracer3.h"
/*^*/

typedef struct EscBGrid3 *sf_esc_bgrid3;
/* abstract data type */
/*^*/

#define ESC3_BGRID3_APRON 2
/*^*/

#endif

struct EscBGrid3 {
    int                   nz, nx, ny, na, nb;
    float                 oz, ox, oy, oa, ob;
    float                 dz, dx, dy, da, db;
    float                 zmax, xmax, ymax;
    bool                  verb;
    sf_timer              ttime;
    sf_esc_slowness3      esc_slow;
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
                                  sf_esc_slowness3 esc_slow, sf_esc_tracer3 esc_tracer)
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

    esc_bgrid->esc_slow = esc_slow;
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
    float z, x, y, b, a;
    float *****top;
    sf_esc_point3f esc_pointf;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s **splines;

    sf_esc_bgrid3_set_tracer_limits (esc_bgrid);
    top = sf_floatalloc5 (esc_bgrid->nx, esc_bgrid->ny, esc_bgrid->nb + 2*ESC3_BGRID3_APRON,
                          esc_bgrid->na, ESC3F_NUM);

    if (esc_bgrid->verb) {
         sf_timer_reset (esc_bgrid->ttime);
         sf_timer_start (esc_bgrid->ttime);
    }

    esc_pointf = sf_esc_point3f_init ();

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
                                            0.0, 0.0, NULL, esc_pointf);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3F_NUM; i++)
                        top[i][ia][ib + ESC3_BGRID3_APRON][iy][ix] =
                        sf_esc_point3f_get_esc_var (esc_pointf, i);
                }
            }
        }
    }
    /* Copy boundary values in b direction to the apron area -
       this will later help with the proper boundary conditions
       during interpolation */
    for (i = 0; i < ESC3F_NUM; i++) {
        for (ia = 0; ia < esc_bgrid->na; ia++) {
            /* Boundary condition on b - as we go past 0 or PI,
               we add +/-PI to a */
            if (ia < esc_bgrid->na/2)
                iap = esc_bgrid->na/2 + ia;
            else
                iap = ia - esc_bgrid->na/2;
            for (iy = 0; iy < esc_bgrid->ny; iy++) {
                for (ix = 0; ix < esc_bgrid->ny; ix++) {
                    for (j = 0; j < ESC3_BGRID3_APRON; j++) {
                        top[i][ia][ESC3_BGRID3_APRON - 1 - j][iy][ix] =
                        top[i][iap][ESC3_BGRID3_APRON + j][iy][ix];
                        top[i][ia][esc_bgrid->nb + ESC3_BGRID3_APRON + j][iy][ix] =
                        top[i][iap][esc_bgrid->nb + ESC3_BGRID3_APRON - 1 - j][iy][ix];
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
        splines[ia] = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3F_NUM);
        for (i = 0; i < ESC3F_NUM; i++) {
            set_multi_UBspline_3d_s (splines[ia], i, &top[i][ia][0][0][0]);
        }
    }
    sf_esc_bgrid3_free5 (top);

    if (esc_bgrid->verb) {
        sf_timer_stop (esc_bgrid->ttime);
        sf_warning ("Runtime: %g mins",
                   sf_timer_get_total_time (esc_bgrid->ttime)/60000.0);
    }

    sf_esc_point3f_close (esc_pointf);

    return splines;
}

multi_UBspline_3d_s** sf_esc_bgrid3_compute_leftright (sf_esc_bgrid3 esc_bgrid, bool isleft)
/*< Run escape values computations for the left of right boundary >*/
{
    int iz, iy, ib, ia, iap, i, j;
    float z, x, y, b, a;
    float *****left;
    sf_esc_point3f esc_pointf;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s **splines;

    sf_esc_bgrid3_set_tracer_limits (esc_bgrid);
    left = sf_floatalloc5 (esc_bgrid->nz, esc_bgrid->ny, esc_bgrid->nb + 2*ESC3_BGRID3_APRON,
                           esc_bgrid->na, ESC3F_NUM);

    if (esc_bgrid->verb) {
         sf_timer_reset (esc_bgrid->ttime);
         sf_timer_start (esc_bgrid->ttime);
    }

    esc_pointf = sf_esc_point3f_init ();

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
                                            0.0, 0.0, NULL, esc_pointf);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3F_NUM; i++)
                        left[i][ia][ib + ESC3_BGRID3_APRON][iy][iz] =
                        sf_esc_point3f_get_esc_var (esc_pointf, i);
                }
            }
        }
    }
    /* Copy boundary values in b direction to the apron area -
       this will later help with the proper boundary conditions
       during interpolation */
    for (i = 0; i < ESC3F_NUM; i++) {
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
                        left[i][ia][ESC3_BGRID3_APRON - 1 - j][iy][iz] =
                        left[i][iap][ESC3_BGRID3_APRON + j][iy][iz];
                        left[i][ia][esc_bgrid->nb + ESC3_BGRID3_APRON + j][iy][iz] =
                        left[i][iap][esc_bgrid->nb + ESC3_BGRID3_APRON - 1 - j][iy][iz];
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
        splines[ia] = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3F_NUM);
        for (i = 0; i < ESC3F_NUM; i++) {
            set_multi_UBspline_3d_s (splines[ia], i, &left[i][ia][0][0][0]);
        }
    }
    sf_esc_bgrid3_free5 (left);
    if (esc_bgrid->verb) {
        sf_timer_stop (esc_bgrid->ttime);
        sf_warning ("Runtime: %g mins",
                   sf_timer_get_total_time (esc_bgrid->ttime)/60000.0);
    }

    sf_esc_point3f_close (esc_pointf);

    return splines;
}

multi_UBspline_3d_s** sf_esc_bgrid3_compute_nearfar (sf_esc_bgrid3 esc_bgrid, bool isnear)
/*< Run escape values computations for the near of far boundary >*/
{
    int iz, ix, ib, ia, iap, i, j;
    float z, x, y, b, a;
    float *****near;
    sf_esc_point3f esc_pointf;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s **splines;

    sf_esc_bgrid3_set_tracer_limits (esc_bgrid);
    near = sf_floatalloc5 (esc_bgrid->nz, esc_bgrid->nx, esc_bgrid->nb + 2*ESC3_BGRID3_APRON,
                           esc_bgrid->na, ESC3F_NUM);

    if (esc_bgrid->verb) {
         sf_timer_reset (esc_bgrid->ttime);
         sf_timer_start (esc_bgrid->ttime);
    }

    esc_pointf = sf_esc_point3f_init ();

    y = isnear ? esc_bgrid->oy : esc_bgrid->ymax;
    for (ia = 0; ia < esc_bgrid->na; ia++) {
        a = esc_bgrid->oa + ia*esc_bgrid->da;
        for (ib = 0; ib < esc_bgrid->nb; ib++) {
            b = esc_bgrid->ob + ib*esc_bgrid->db;
            for (ix = 0; ix < esc_bgrid->ny; ix++) {
                x = esc_bgrid->ox + ix*esc_bgrid->dx;
                for (iz = 0; iz < esc_bgrid->nz; iz++) {
                    z = esc_bgrid->oz + iz*esc_bgrid->dz;
                    sf_esc_tracer3_compute (esc_bgrid->esc_tracer, z, x, y, b, a,
                                            0.0, 0.0, NULL, esc_pointf);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3F_NUM; i++)
                        near[i][ia][ib + ESC3_BGRID3_APRON][ix][iz] =
                        sf_esc_point3f_get_esc_var (esc_pointf, i);
                }
            }
        }
    }
    /* Copy boundary values in b direction to the apron area -
       this will later help with the proper boundary conditions
       during interpolation */
    for (i = 0; i < ESC3F_NUM; i++) {
        for (ia = 0; ia < esc_bgrid->na; ia++) {
            /* Boundary condition on b - as we go past 0 or PI,
               we add +/-PI to a */
            if (ia < esc_bgrid->na/2)
                iap = esc_bgrid->na/2 + ia;
            else
                iap = ia - esc_bgrid->na/2;
            for (ix = 0; ix < esc_bgrid->ny; ix++) {
                for (iz = 0; iz < esc_bgrid->nz; iz++) {
                    for (j = 0; j < ESC3_BGRID3_APRON; j++) {
                        near[i][ia][ESC3_BGRID3_APRON - 1 - j][ix][iz] =
                        near[i][iap][ESC3_BGRID3_APRON + j][ix][iz];
                        near[i][ia][esc_bgrid->nb + ESC3_BGRID3_APRON + j][ix][iz] =
                        near[i][iap][esc_bgrid->nb + ESC3_BGRID3_APRON - 1 - j][ix][iz];
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
        splines[ia] = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3F_NUM);
        for (i = 0; i < ESC3F_NUM; i++) {
            set_multi_UBspline_3d_s (splines[ia], i, &near[i][ia][0][0][0]);
        }
    }
    sf_esc_bgrid3_free5 (near);

    if (esc_bgrid->verb) {
        sf_timer_stop (esc_bgrid->ttime);
        sf_warning ("Runtime: %g mins",
                   sf_timer_get_total_time (esc_bgrid->ttime)/60000.0);
    }

    sf_esc_point3f_close (esc_pointf);

    return splines;
}


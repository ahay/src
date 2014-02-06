/* Ray tracer in 5-D phase-space for 3-D medium with escape equations */
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

#ifndef _esc_tracer3_h

#include "esc_slow3.h"
#include "esc_point3.h"
/*^*/

typedef struct EscTracer3 *sf_esc_tracer3;
/* abstract data type */
/*^*/

#endif

typedef void (*sf_esc_tracer3_traj)(float z, float x, float y, float b,
                                    float a, int it, void *ud)
/*< Callback to output points along a ray trajectory >*/
;

struct EscTracer3 {
    int                  nz, nx, ny, na, nb;
    float                oz, ox, oy, oa, ob;
    float                dz, dx, dy, da, db;
    float                dt, df, md;
    float                zmin, zmax;
    float                xmin, xmax;
    float                ymin, ymax;
    float                zgmin, zgmax;
    float                xgmin, xgmax;
    float                ygmin, ygmax;
    bool                 parab;
    sf_esc_slowness3     esc_slow;
    sf_esc_tracer3_traj  traj;
    void                *ud;
};
/* concrete data type */

void sf_esc_tracer3_reset_bounds (sf_esc_tracer3 esc_tracer)
/*< Reset spatial bounds >*/
{
    esc_tracer->zmin = esc_tracer->zgmin;
    esc_tracer->zmax = esc_tracer->zgmax;
    esc_tracer->xmin = esc_tracer->xgmin;
    esc_tracer->xmax = esc_tracer->xgmax;
    esc_tracer->ymin = esc_tracer->ygmin;
    esc_tracer->ymax = esc_tracer->ygmax;
}

sf_esc_tracer3 sf_esc_tracer3_init (sf_esc_slowness3 esc_slow)
/*< Initialize object >*/
{
    sf_esc_tracer3 esc_tracer = (sf_esc_tracer3)sf_alloc (1, sizeof (struct EscTracer3));

    esc_tracer->nz = sf_esc_slowness3_nz (esc_slow);
    esc_tracer->nx = sf_esc_slowness3_nx (esc_slow);
    esc_tracer->ny = sf_esc_slowness3_ny (esc_slow);
    esc_tracer->oz = sf_esc_slowness3_oz (esc_slow);
    esc_tracer->dz = sf_esc_slowness3_dz (esc_slow);
    esc_tracer->ox = sf_esc_slowness3_ox (esc_slow);
    esc_tracer->dx = sf_esc_slowness3_dx (esc_slow);
    esc_tracer->oy = sf_esc_slowness3_oy (esc_slow);
    esc_tracer->dy = sf_esc_slowness3_dy (esc_slow);
    /* Use 1/4 of a degree as acceptable angle displacement
       in the straght ray approximation below */
    esc_tracer->da = 0.25*SF_PI/180.0;
    esc_tracer->db = 0.25*SF_PI/180.0;

    esc_tracer->md = SF_HUGE; /* Maximum allowed distance along a ray */

    /* Use smaller spatial steps than those in the velocity model;
       analytical solution (either parabolic or straight line) tends
       to be inaccurate with large steps; typical velocity volumes
       have 20-50 meter sampling with velocity gradients changing
       noticeably across one cell; this should probably be adaptive */
    esc_tracer->df = 0.25;

    /* Global limits */
    esc_tracer->zgmin = esc_tracer->oz;
    esc_tracer->zgmax = esc_tracer->oz + (esc_tracer->nz - 1)*esc_tracer->dz;
    esc_tracer->xgmin = esc_tracer->ox;
    esc_tracer->xgmax = esc_tracer->ox + (esc_tracer->nx - 1)*esc_tracer->dx;
    esc_tracer->ygmin = esc_tracer->oy;
    esc_tracer->ygmax = esc_tracer->oy + (esc_tracer->ny - 1)*esc_tracer->dy;

    sf_esc_tracer3_reset_bounds (esc_tracer);

    esc_tracer->parab = true;

    esc_tracer->esc_slow = esc_slow;

    esc_tracer->traj = NULL;
    esc_tracer->dt = 0.001;
    esc_tracer->ud = NULL;

    return esc_tracer;
}

void sf_esc_tracer3_close (sf_esc_tracer3 esc_tracer)
/*< Destroy object >*/
{
    free (esc_tracer);
}

void sf_esc_tracer3_set_trajcb (sf_esc_tracer3 esc_tracer,
                                sf_esc_tracer3_traj traj, float dt, void *ud)
/*< Set trajectory callback  >*/
{
    /* Callback to output points along a ray trajectory */
    esc_tracer->traj = traj;
    esc_tracer->dt = dt; /* Output points every dt time intervals */
    esc_tracer->ud = ud;
}

void sf_esc_tracer3_set_zmin (sf_esc_tracer3 esc_tracer, float zmin)
/*< Set spatial bound >*/
{
    if (zmin > esc_tracer->zgmin)
        esc_tracer->zmin = zmin;
    else
        esc_tracer->zmin = esc_tracer->zgmin;
}

void sf_esc_tracer3_set_zmax (sf_esc_tracer3 esc_tracer, float zmax)
/*< Set spatial bound >*/
{
    if (zmax < esc_tracer->zgmax)
        esc_tracer->zmax = zmax;
    else
        esc_tracer->zmax = esc_tracer->zgmax;
}

void sf_esc_tracer3_set_xmin (sf_esc_tracer3 esc_tracer, float xmin)
/*< Set spatial bound >*/
{
    if (xmin > esc_tracer->xgmin)
        esc_tracer->xmin = xmin;
    else
        esc_tracer->xmin = esc_tracer->xgmin;
}

void sf_esc_tracer3_set_xmax (sf_esc_tracer3 esc_tracer, float xmax)
/*< Set spatial bound >*/
{
    if (xmax < esc_tracer->xgmax)
        esc_tracer->xmax = xmax;
    else
        esc_tracer->xmax = esc_tracer->xgmax;
}

void sf_esc_tracer3_set_ymin (sf_esc_tracer3 esc_tracer, float ymin)
/*< Set spatial bound >*/
{
    if (ymin > esc_tracer->ygmin)
        esc_tracer->ymin = ymin;
    else
        esc_tracer->ymin = esc_tracer->ygmin;
}

void sf_esc_tracer3_set_ymax (sf_esc_tracer3 esc_tracer, float ymax)
/*< Set spatial bound >*/
{
    if (ymax < esc_tracer->ygmax)
        esc_tracer->ymax = ymax;
    else
        esc_tracer->ymax = esc_tracer->ygmax;
}

void sf_esc_tracer3_set_parab (sf_esc_tracer3 esc_tracer, bool parab)
/*< Set parabolic/straight ray flag >*/
{
    esc_tracer->parab = parab;
}

void sf_esc_tracer3_set_mdist (sf_esc_tracer3 esc_tracer, float md)
/*< Set maximum allowed distance to travel for a ray >*/
{
    esc_tracer->md = md;
}

void sf_esc_tracer3_set_df (sf_esc_tracer3 esc_tracer, float df)
/*< Set maximum allowed distance to travel per step (fraction of the cell size) >*/
{
    esc_tracer->df = df;
}

float sf_esc_tracer3_sintersect (sf_esc_tracer3 esc_tracer, float *z, float *x, float *y,
                                 float *b, float *a, float dz, float dx, float dy,
                                 float db, float da, float fz, float fx, float fy,
                                 float fb, float fa)
/*< Compute intersection of a straight trajectory from (z, x, y, b, a) with 
    the nearest wall defined by (dz, dx, dy, db, da), return pseudotime along the trajectory >*/
{
    float tz, tx, ty, tb, ta, sigma, rp = 1e-6;
    /* Loop until curvature approximation by a straight line is adequate */
    do {
        /* Time to phase-space cell walls */
        tz = fabsf (fz) > rp ? fabsf (dz/fz) : SF_HUGE;
        tx = fabsf (fx) > rp ? fabsf (dx/fx) : SF_HUGE;
        ty = fabsf (fy) > rp ? fabsf (dy/fy) : SF_HUGE;
        tb = fabsf (fb) > rp ? fabsf (db/fb) : SF_HUGE;
        ta = fabsf (fa) > rp ? fabsf (da/fa) : SF_HUGE;
        /* Hitting the angle wall first - too much curvature on the ray,
           reduce distance by half */
        if ((ta < tz && ta < tx && ta < ty) ||
            (tb < tz && tb < tx && tb < ty)) {
            dz *= 0.5;
            dx *= 0.5;
            dy *= 0.5;
        }
    } while ((ta < tz && ta < tx && ta < ty) ||
             (tb < tz && tb < tx && tb < ty));
    sigma = SF_MIN (ta, tb);
    sigma = SF_MIN (ty, sigma);
    sigma = SF_MIN (tx, sigma);
    sigma = SF_MIN (tz, sigma);
    *z -= fz*sigma;
    *x -= fx*sigma;
    *y -= fy*sigma;
    *b -= fb*sigma;
    *a -= fa*sigma;
    return sigma;
}

float sf_esc_tracer3_pintersect (sf_esc_tracer3 esc_tracer, float *z, float *x, float *y,
                                 float *b, float *a, float *t, float *dd, float dz, float dx, float dy,
                                 float fz, float fx, float fy, float s, float sz, float sx, float sy,
                                 float md)
/*< Compute intersection of a parabolic trajectory from (z, x, y, b, a) with 
    the nearest wall defined by (dz, dx, dy), return pseudotime along the trajectory >*/
{
    float A, B, C, D, s1, s2, sigma, az, ax, ay,
          pz, px, py, pz0, px0, py0, l, ddz, ddx, ddy;
    *dd = 0.0;
    /* Assume locally constant slowness and slowness gradients */
    /* Parabola - dz = -v_z*sigma + 0.5*a_z*sigma^2 */

    fz = -fz;
    fx = -fx;
    fy = -fy;
    az = s*sz;
    ax = s*sx;
    ay = s*sy;
    pz0 = -s*cosf (*b);
    px0 = s*sinf (*b)*cosf (*a);
    py0 = s*sinf (*b)*sinf (*a);

    do {
        if (*dd > md) {
            dz *= md/(*dd);
            dx *= md/(*dd);
            dy *= md/(*dd);
        }
        if (*b <= SF_PI/4.0 || *b >= 3.0*SF_PI/4.0) {
            /* Intersection with z */
            A = 0.5*az;
            B = pz0;
            C = -dz;
        } else if ((*a <= SF_PI/4.0 || *a >= 7.0*SF_PI/4.0) ||
                   (*a >= 3.0*SF_PI/4.0 && *a <= 5.0*SF_PI/4.0)) {
            /* Intersection with x */
            A = 0.5*ax;
            B = px0;
            C = -dx;
        } else {
            /* Intersection with y */
            A = 0.5*ay;
            B = py0;
            C = -dy;
        }
        if (4.0*A*C > B*B)
/*          sf_error ("Parabola miss");*/
            return SF_HUGE;
        /* Solve the parabolic equation */
        D = sqrt (B*B - 4.0*A*C);
        if (fabsf (A) > 1e-7 && (-B + D) > 1e-7 && (-B - D) > 1e-7) {
            s1 = (-B + D)/(2.0*A);
            s2 = (-B - D)/(2.0*A);
            sigma = fabsf (s1) < fabsf (s2) ? s1 : s2;
        } else
            sigma = -C/B;
        ddz = fz*sigma + 0.5*az*sigma*sigma;
        ddx = fx*sigma + 0.5*ax*sigma*sigma;
        ddy = fy*sigma + 0.5*ay*sigma*sigma;
        *dd = sqrtf (ddz*ddz + ddx*ddx + ddy*ddy);
    } while (*dd > md);

    *z += ddz;
    *x += ddx;
    *y += ddy;
    pz = pz0 + az*sigma;
    px = px0 + ax*sigma;
    py = py0 + ay*sigma;
    /* Find new phase angle */
    l = sqrt (px*px + py*py);
    if (l > 1e-6) {
        if (py >= 0.0)
            *a = acosf (px/l);
        else
            *a = 2.0*SF_PI - acosf (px/l);
    }
    l = sqrt (pz*pz + px*px + py*py);
    *b = acosf (-pz/l);
    *t += (pz0*pz0 + px0*px0 + py0*py0)*sigma + (pz0*az + px0*ax + py0*ay)*sigma*sigma +
          (az*az + ax*ax + ay*ay)*sigma*sigma*sigma/3.0;
    return sigma;
}

bool sf_esc_tracer3_inside (sf_esc_tracer3 esc_tracer, float *z, float *x, float *y,
                            bool snap)
/*< Return true, if point (z,x,y) is inside the current limits;
    snap point to the boundary otherwise, if snap=true >*/
{
    float eps = 1e-2;
    float ezmin, ezmax, exmin, exmax, eymin, eymax;

    /* Bounding box + epsilon */
    ezmin = esc_tracer->zmin + eps*esc_tracer->dz;
    ezmax = esc_tracer->zmax - eps*esc_tracer->dz;
    exmin = esc_tracer->xmin + eps*esc_tracer->dx;
    exmax = esc_tracer->xmax - eps*esc_tracer->dx;
    eymin = esc_tracer->ymin + eps*esc_tracer->dy;
    eymax = esc_tracer->ymax - eps*esc_tracer->dy;

    if (*y > eymin && *y < eymax &&
        *x > exmin && *x < exmax &&
        *z > ezmin && *z < ezmax)
        return true;

    if (snap) {
        if (*z <= ezmin) {
            *z = esc_tracer->zmin;
        } else if (*z >= ezmax) {
            *z = esc_tracer->zmax;
        }
        if (*x <= exmin) {
            *x = esc_tracer->xmin;
        } else if (*x >= exmax) {
            *x = esc_tracer->xmax;
        }
        if (*y <= eymin) {
            *y = esc_tracer->ymin;
        } else if (*y >= eymax) {
            *y = esc_tracer->ymax;
        }
    }

    return false;
}

void sf_esc_tracer3_compute (sf_esc_tracer3 esc_tracer, float z, float x, float y,
                             float b, float a, float t, float l,
                             sf_esc_point3 point, float *ae, float *be)
/*< Compute escape values for a point with subsurface coordinates (z, x, y, b, a) >*/
{
    int pit = -1, it = 0;
    float eps = 1e-2;
    float s, sp, sb, sa, sz, sx, sy, dd;
    float dz, dx, dy, db, da, fz, fx, fy, fb, fa, ll = 0.0, sigma;
    float ezmin, ezmax, exmin, exmax, eymin, eymax;
    EscColor3 col = 0;

    /* Bounding box + epsilon */
    ezmin = esc_tracer->zmin + eps*esc_tracer->dz;
    ezmax = esc_tracer->zmax - eps*esc_tracer->dz;
    exmin = esc_tracer->xmin + eps*esc_tracer->dx;
    exmax = esc_tracer->xmax - eps*esc_tracer->dx;
    eymin = esc_tracer->ymin + eps*esc_tracer->dy;
    eymax = esc_tracer->ymax - eps*esc_tracer->dy;

    /* Get slowness and derivatives */
    sf_esc_slowness3_get_components (esc_tracer->esc_slow, z, x, y, b, a,
                                     &s, &sb, &sa, &sz, &sx, &sy);
    do {
        /* Call trajectory point callback */
        if (esc_tracer->traj) {
            it = t/esc_tracer->dt;
            if (pit != it) {
                pit = it;
                esc_tracer->traj (z, x, y, b*180.0/SF_PI, a*180.0/SF_PI,
                                  it, esc_tracer->ud);
            }
        }
        /* Advection coefficients */
        sf_esc_slowness3_get_coefs (esc_tracer->esc_slow,
                                    b, a, s, sz, sx, sy, sa, sb,
                                    &fz, &fx, &fy, &fa, &fb);
        /* Displacements */
        dz = fz < 0.0 ? esc_tracer->dz : -esc_tracer->dz;
        dx = fx < 0.0 ? esc_tracer->dx : -esc_tracer->dx;
        dy = fy < 0.0 ? esc_tracer->dy : -esc_tracer->dy;
        da = fa < 0.0 ? esc_tracer->da : -esc_tracer->da;
        db = fb < 0.0 ? esc_tracer->db : -esc_tracer->db;
        dz *= esc_tracer->df;
        dx *= esc_tracer->df;
        dy *= esc_tracer->df;
        /* Adjust if near boundaries */
        if ((z + dz) < esc_tracer->zmin)
            dz = esc_tracer->zmin - z;
        if ((z + dz) > esc_tracer->zmax)
            dz = esc_tracer->zmax - z;
        if ((x + dx) < esc_tracer->xmin)
            dx = esc_tracer->xmin - x;
        if ((x + dx) > esc_tracer->xmax)
            dx = esc_tracer->xmax - x;
        if ((y + dy) < esc_tracer->ymin)
            dy = esc_tracer->ymin - y;
        if ((y + dy) > esc_tracer->ymax)
            dy = esc_tracer->ymax - y;
        sp = s;
        if (esc_tracer->parab) { /* Intersection with a parabolic trajectory */
            sigma = SF_HUGE;
            while (SF_HUGE == sigma) {
                sigma = sf_esc_tracer3_pintersect (esc_tracer, &z, &x, &y, &b, &a, &t, &dd,
                                                   dz, dx, dy, fz, fx, fy, s, sz, sx, sy,
                                                   esc_tracer->md != SF_HUGE ?
                                                   esc_tracer->md - ll : esc_tracer->md);
                if (SF_HUGE == sigma) {
                    /* Adjust cell sizes, if the analytic solution failed;
                       this happens in case of overturning rays */
                    dz *= 0.5;
                    dx *= 0.5;
                    dy *= 0.5;
                }
            }
        } else /* Intersection with a straight trajectory */
            sigma = sf_esc_tracer3_sintersect (esc_tracer, &z, &x, &y, &b, &a,
                                               dz, dx, dy, db, da, fz, fx, fy, fb, fa);
        /* Keep b in [0; pi] and a in [0; 2*pi] range */
        if (a < 0.0)
            a += 2.0*SF_PI;
        else if (a > 2.0*SF_PI)
            a -= 2.0*SF_PI;
        if (b < 0.0) {
            a += SF_PI;
            b = fabsf (b);
        } else if (b > SF_PI) {
            a -= SF_PI;
            b = 2.0*SF_PI - b;
        }  
        if (a < 0.0)
            a += 2.0*SF_PI;
        else if (a > 2.0*SF_PI)
            a -= 2.0*SF_PI;
        /* Update slowness and it components at the new location */
        sf_esc_slowness3_get_components (esc_tracer->esc_slow, z, x, y, b, a,
                                         &s, &sb, &sa, &sz, &sx, &sy);
        if (!esc_tracer->parab) {
            /* Length of this segment of the characteristic */
            dd = fabsf (sigma*sqrt (fz*fz + fx*fx + fy*fy));
            t += dd*(s + sp)*0.5;
        }
        ll += dd;
#ifdef ESC_EQ_WITH_L
        l += dd;
#endif
    } while (y > eymin && y < eymax &&
             x > exmin && x < exmax &&
             z > ezmin && z < ezmax &&
             ll < esc_tracer->md);
    /* Snap to boundary */
    if (z <= ezmin) {
        z = esc_tracer->zmin;
        col |= ESC3_TOP;
    } else if (z >= ezmax) {
        z = esc_tracer->zmax;
        col |= ESC3_BOTTOM;
    }
    if (x <= exmin) {
        x = esc_tracer->xmin;
        col |= ESC3_LEFT;
    } else if (x >= exmax) {
        x = esc_tracer->xmax;
        col |= ESC3_RIGHT;
    }
    if (y <= eymin) {
        y = esc_tracer->ymin;
        col |= ESC3_LEFT;
    } else if (y >= eymax) {
        y = esc_tracer->ymax;
        col |= ESC3_RIGHT;
    }

    /* Last trajectory point */
    if (esc_tracer->traj) {
        it = t/esc_tracer->dt;
        if (pit != it) {
            pit = it;
            esc_tracer->traj (z, x, y, b*180.0/SF_PI, a*180.0/SF_PI,
                              it, esc_tracer->ud);
        }
    }

    if (point) {
        sf_esc_point3_set_esc_var (point, ESC3_Z, z);
        sf_esc_point3_set_esc_var (point, ESC3_X, x);
        sf_esc_point3_set_esc_var (point, ESC3_Y, y);
        sf_esc_point3_set_esc_var (point, ESC3_T, t);
#ifdef ESC_EQ_WITH_L
        sf_esc_point3_set_esc_var (point, ESC3_L, l);
#endif
        sf_esc_point3_set_col (point, col);
    }
    if (ae)
        *ae = a;
    if (be)
        *be = b;
}


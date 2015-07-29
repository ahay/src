/* Ray tracer in 3-D phase-space for 2-D medium with escape equations */
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

#ifndef _esc_tracer2_h

#include "esc_slow2.h"
#include "esc_point2.h"
/*^*/

typedef struct EscTracer2 *sf_esc_tracer2;
/* abstract data type */
/*^*/

#endif

typedef void (*sf_esc_tracer2_traj)(float z, float x, float a, int it, void *ud)
/*< Callback to output points along a ray trajectory >*/
;

struct EscTracer2 {
    int                  nz, nx, na;
    float                oz, ox, oa;
    float                dz, dx, da;
    float                dt, df, md;
    float                zmin, zmax;
    float                xmin, xmax;
    float                zgmin, zgmax;
    float                xgmin, xgmax;
    bool                 parab;
    sf_esc_slowness2     esc_slow;
    sf_esc_tracer2_traj  traj;
    void                *ud;
};
/* concrete data type */

void sf_esc_tracer2_reset_bounds (sf_esc_tracer2 esc_tracer)
/*< Reset spatial bounds >*/
{
    esc_tracer->zmin = esc_tracer->zgmin;
    esc_tracer->zmax = esc_tracer->zgmax;
    esc_tracer->xmin = esc_tracer->xgmin;
    esc_tracer->xmax = esc_tracer->xgmax;
}

sf_esc_tracer2 sf_esc_tracer2_init (sf_esc_slowness2 esc_slow)
/*< Initialize object >*/
{
    sf_esc_tracer2 esc_tracer = (sf_esc_tracer2)sf_alloc (1, sizeof (struct EscTracer2));

    esc_tracer->nz = sf_esc_slowness2_nz (esc_slow);
    esc_tracer->nx = sf_esc_slowness2_nx (esc_slow);
    esc_tracer->oz = sf_esc_slowness2_oz (esc_slow);
    esc_tracer->ox = sf_esc_slowness2_ox (esc_slow);
    esc_tracer->dz = sf_esc_slowness2_dz (esc_slow);
    esc_tracer->dx = sf_esc_slowness2_dx (esc_slow);
    /* Use 1/4 of a degree as acceptable angle displacement
       in the straght ray approximation below */
    esc_tracer->da = 0.25*SF_PI/180.0;

    esc_tracer->md = SF_HUGE; /* Maximum allowed distance along a ray */

    /* Use smaller spatial steps than those in the velocity model;
       analytical solution (either parabolic or straight line) tends
       to be inaccurate with large steps; typical velocity volumes
       have 20-50 meter sampling with velocity gradients changing
       noticeably across one cell; this should probably be adaptive */
    esc_tracer->df = 0.5;

    /* Global limits */
    esc_tracer->zgmin = esc_tracer->oz;
    esc_tracer->zgmax = esc_tracer->oz + (esc_tracer->nz - 1)*esc_tracer->dz;
    esc_tracer->xgmin = esc_tracer->ox;
    esc_tracer->xgmax = esc_tracer->ox + (esc_tracer->nx - 1)*esc_tracer->dx;

    sf_esc_tracer2_reset_bounds (esc_tracer);

    esc_tracer->parab = true;

    esc_tracer->esc_slow = esc_slow;

    esc_tracer->traj = NULL;
    esc_tracer->dt = 0.001;
    esc_tracer->ud = NULL;

    return esc_tracer;
}

void sf_esc_tracer2_close (sf_esc_tracer2 esc_tracer)
/*< Destroy object >*/
{
    free (esc_tracer);
}

void sf_esc_tracer2_set_trajcb (sf_esc_tracer2 esc_tracer,
                                sf_esc_tracer2_traj traj, float dt, void *ud)
/*< Set trajectory callback  >*/
{
    /* Callback to output points along a ray trajectory */
    esc_tracer->traj = traj;
    esc_tracer->dt = dt; /* Output points every dt time intervals */
    esc_tracer->ud = ud;
}

void sf_esc_tracer2_set_zmin (sf_esc_tracer2 esc_tracer, float zmin)
/*< Set spatial bound >*/
{
    if (zmin > esc_tracer->zgmin)
        esc_tracer->zmin = zmin;
    else
        esc_tracer->zmin = esc_tracer->zgmin;
}

void sf_esc_tracer2_set_zmax (sf_esc_tracer2 esc_tracer, float zmax)
/*< Set spatial bound >*/
{
    if (zmax < esc_tracer->zgmax)
        esc_tracer->zmax = zmax;
    else
        esc_tracer->zmax = esc_tracer->zgmax;
}

void sf_esc_tracer2_set_xmin (sf_esc_tracer2 esc_tracer, float xmin)
/*< Set spatial bound >*/
{
    if (xmin > esc_tracer->xgmin)
        esc_tracer->xmin = xmin;
    else
        esc_tracer->xmin = esc_tracer->xgmin;
}

void sf_esc_tracer2_set_xmax (sf_esc_tracer2 esc_tracer, float xmax)
/*< Set spatial bound >*/
{
    if (xmax < esc_tracer->xgmax)
        esc_tracer->xmax = xmax;
    else
        esc_tracer->xmax = esc_tracer->xgmax;
}

void sf_esc_tracer2_set_parab (sf_esc_tracer2 esc_tracer, bool parab)
/*< Set parabolic/straight ray flag >*/
{
    esc_tracer->parab = parab;
}

void sf_esc_tracer2_set_mdist (sf_esc_tracer2 esc_tracer, float md)
/*< Set maximum allowed distance to travel for a ray >*/
{
    esc_tracer->md = md;
}

void sf_esc_tracer2_set_df (sf_esc_tracer2 esc_tracer, float df)
/*< Set maximum allowed distance to travel per step (fraction of the cell size) >*/
{
    esc_tracer->df = df;
}

float sf_esc_tracer2_sintersect (sf_esc_tracer2 esc_tracer, float *z, float *x, float *a,
                                 float dz, float dx, float da, float fz, float fx, float fa)
/*< Compute intersection of a straight trajectory from (z, x, a) with 
    the nearest wall defined by (dz, dx, da), return pseudotime along the trajectory >*/
{
    float tz, tx, ta, sigma, rp = 1e-6;
    /* Loop until curvature approximation by a straight line is adequate */
    do {
        /* Time to phase-space cell walls */
        tz = fabsf (fz) > rp ? fabsf (dz/fz) : SF_HUGE;
        tx = fabsf (fx) > rp ? fabsf (dx/fx) : SF_HUGE;
        ta = fabsf (fa) > rp ? fabsf (da/fa) : SF_HUGE;
        /* Hitting the angle wall first - too much curvature on the ray,
           reduce distance by half */
        if (ta < tz && ta < tx) {
            dz *= 0.5;
            dx *= 0.5;
        }
    } while (ta < tz && ta < tx);
    sigma = SF_MIN (tz, SF_MIN (tx, ta));
    *z -= fz*sigma;
    *x -= fx*sigma;
    *a -= fa*sigma;
    return sigma;
}

float sf_esc_tracer2_pintersect (sf_esc_tracer2 esc_tracer, float *z, float *x, float *a, float *t, float *dd,
                                 float dz, float dx, float fz, float fx, float s, float sz, float sx, float md)
/*< Compute intersection of a parabolic trajectory from (z, x, a) with 
    the nearest wall defined by (dz, dx), return pseudotime along the trajectory >*/
{
    float A, B, C, D, s1, s2, sigma, az, ax, pz, px, pz0, px0, l, ddz, ddx;
    *dd = 0.0;
    /* Assume locally constant slowness and slowness gradients */
    /* Parabola - dz = -v_z*sigma + 0.5*a_z*sigma^2 */

    fz = -fz;
    fx = -fx;
    az = s*sz;
    ax = s*sx;
    pz0 = -s*cos (*a);
    px0 = -s*sin (*a);

    do {
        if (*dd > md) {
            dz *= md/(*dd);
            dx *= md/(*dd);
        }
        if ((*a >= -SF_PI/4.0 && *a <= SF_PI/4.0) ||
            (*a <= -3.0*SF_PI/4.0 || *a >= 3.0*SF_PI/4.0)) {
            /* Intersection with z */
            A = 0.5*az;
            B = pz0;
            C = -dz;
        } else {
            /* Intersection with x */
            A = 0.5*ax;
            B = px0;
            C = -dx;
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
        *dd = sqrtf (ddz*ddz + ddx*ddx);
    } while (*dd > md);

    *z += ddz;
    *x += ddx;
    pz = pz0 + az*sigma;
    px = px0 + ax*sigma;

    /* Find new phase angle */
    l = hypotf (pz, px);
    if (px < 0.0)
        *a = acos (-pz/l);
    else
        *a = -acos (-pz/l);
    *t += (pz0*pz0 + px0*px0)*sigma + (pz0*az + px0*ax)*sigma*sigma +
          (az*az + ax*ax)*sigma*sigma*sigma/3.0;
    return sigma;
}

bool sf_esc_tracer2_inside (sf_esc_tracer2 esc_tracer, float *z, float *x,
                            bool snap)
/*< Return true, if point (z,x) is inside the current limits;
    snap point to the boundary otherwise, if snap=true >*/
{
    float eps = 1e-2;
    float ezmin, ezmax, exmin, exmax;

    /* Bounding box + epsilon */
    ezmin = esc_tracer->zmin + eps*esc_tracer->dz;
    ezmax = esc_tracer->zmax - eps*esc_tracer->dz;
    exmin = esc_tracer->xmin + eps*esc_tracer->dx;
    exmax = esc_tracer->xmax - eps*esc_tracer->dx;

    if (*x > exmin && *x < exmax &&
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
    }

    return false;
}

void sf_esc_tracer2_compute (sf_esc_tracer2 esc_tracer, float z, float x, float a,
                             float t, float l, sf_esc_point2 point, float *ae)
/*< Compute escape values for a point with subsurface coordinates (z, x, a) >*/
{
    int pit = -1, it = 0;
    float eps = 1e-2;
    float s, sp, sa, sz, sx, dd, ll = 0;
    float dz, dx, da, fz, fx, fa, sigma;
    float ezmin, ezmax, exmin, exmax;
    EscColor2 col = 0;

    /* Bounding box + epsilon */
    ezmin = esc_tracer->zmin + eps*esc_tracer->dz;
    ezmax = esc_tracer->zmax - eps*esc_tracer->dz;
    exmin = esc_tracer->xmin + eps*esc_tracer->dx;
    exmax = esc_tracer->xmax - eps*esc_tracer->dx;

    /* Get slowness and derivatives */
    sf_esc_slowness2_get_components (esc_tracer->esc_slow, z, x, a,
                                     &s, &sa, &sz, &sx);
    do {
        /* Call trajectory point callback */
        if (esc_tracer->traj) {
            it = t/esc_tracer->dt;
            if (pit != it) {
                pit = it;
                esc_tracer->traj (z, x, a*180.0/SF_PI, it, esc_tracer->ud);
            }
        }
        /* Advection coefficients */
        sf_esc_slowness2_get_coefs (esc_tracer->esc_slow, a,
                                    s, sa, sz, sx, &fz, &fx, &fa); 
        /* Displacements */
        dz = fz < 0.0 ? esc_tracer->dz : -esc_tracer->dz;
        dx = fx < 0.0 ? esc_tracer->dx : -esc_tracer->dx;
        da = fa < 0.0 ? esc_tracer->da : -esc_tracer->da;
        /* Use half steps from velocity model */
        dz *= esc_tracer->df;
        dx *= esc_tracer->df;
        /* Adjust if near boundaries */
        if ((z + dz) < esc_tracer->zmin)
            dz = esc_tracer->zmin - z;
        if ((z + dz) > esc_tracer->zmax)
            dz = esc_tracer->zmax - z;
        if ((x + dx) < esc_tracer->xmin)
            dx = esc_tracer->xmin - x;
        if ((x + dx) > esc_tracer->xmax)
            dx = esc_tracer->xmax - x;
        sp = s;
        if (esc_tracer->parab)  { /* Intersection with a parabolic trajectory */
            sigma = SF_HUGE;
            while (SF_HUGE == sigma) {
                sigma = sf_esc_tracer2_pintersect (esc_tracer, &z, &x, &a, &t, &dd,
                                                   dz, dx, fz, fx, s, sz, sx,
                                                   esc_tracer->md != SF_HUGE ?
                                                   esc_tracer->md - ll : esc_tracer->md);
                if (SF_HUGE == sigma) {
                    dz *= 0.5;
                    dx *= 0.5;
                }
            }
        } else /* Intersection with a straight trajectory */
            sigma = sf_esc_tracer2_sintersect (esc_tracer, &z, &x, &a, dz, dx, da,
                                               fz, fx, fa);
        if (a < -SF_PI)
            a += 2.0*SF_PI;
        else if (a > SF_PI)
            a -= 2.0*SF_PI;
        /* Update slowness and it components at the new location */
        sf_esc_slowness2_get_components (esc_tracer->esc_slow, z, x, a,
                                         &s, &sa, &sz, &sx);
        /* Length of this segment of the characteristic */
        dd = fabs (sigma*sqrt (fz*fz + fx*fx));
        if (!esc_tracer->parab) {
            /* Length of this segment of the characteristic */
            dd = fabsf (sigma*sqrt (fz*fz + fx*fx));
            t += dd*(s + sp)*0.5;
        }
        ll += dd;
#ifdef ESC_EQ_WITH_L
        l += dd;
#endif
    } while (x > exmin && x < exmax &&
             z > ezmin && z < ezmax);
    /* Snap to boundary */
    if (z <= ezmin) {
        z = esc_tracer->zmin;
        col |= ESC2_TOP;
    } else if (z >= ezmax) {
        z = esc_tracer->zmax;
        col |= ESC2_BOTTOM;
    }
    if (x <= exmin) {
        x = esc_tracer->xmin;
        col |= ESC2_LEFT;
    } else if (x >= exmax) {
        x = esc_tracer->xmax;
        col |= ESC2_RIGHT;
    }

    /* Last trajectory point */
    if (esc_tracer->traj) {
        it = t/esc_tracer->dt;
        if (pit != it) {
            pit = it;
            esc_tracer->traj (z, x, a*180.0/SF_PI, it, esc_tracer->ud);
        }
    }

    if (point) {
        sf_esc_point2_set_esc_var (point, ESC2_Z, z);
        sf_esc_point2_set_esc_var (point, ESC2_X, x);
        sf_esc_point2_set_esc_var (point, ESC2_T, t);
#ifdef ESC_EQ_WITH_L
        sf_esc_point2_set_esc_var (point, ESC2_L, l);
#endif
        sf_esc_point2_set_col (point, col);
    }
    if (ae)
        *ae = a;
}


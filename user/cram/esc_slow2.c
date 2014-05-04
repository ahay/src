/* 3-D phase-space slowness container for (an)isotropic 2-D medium */
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

#include "einspline.h"

#ifndef _esc_slow2_h

typedef struct EscSlowness2 *sf_esc_slowness2;
/* abstract data type */
/*^*/

#endif

struct EscSlowness2 {
    int                  nz, nx;
    float                oz, ox;
    float                dz, dx;
    float                zmin, zmax;
    float                xmin, xmax;
    bool                 aniso;
    size_t               offs;
    unsigned char       *mmaped;
    multi_UBspline_2d_s  velspline;
};
/* concrete data type */

sf_esc_slowness2 sf_esc_slowness2_init (sf_file vspline, bool verb)
/*< Initialize object >*/
{
    FILE *stream;
    bool spl;
    sf_esc_slowness2 esc_slow = (sf_esc_slowness2)sf_alloc (1, sizeof (struct EscSlowness2));

    if (sf_gettype (vspline) != SF_UCHAR)
        sf_error ("Wrong data type in velocity spline file");
    if (!sf_histbool (vspline, "splines", &spl) || !spl)
        sf_error ("No spline data in velocity spline file");

    if (!sf_histint (vspline, "Nz", &esc_slow->nz)) sf_error ("No Nz= in velocity spline file");
    if (!sf_histfloat (vspline, "Dz", &esc_slow->dz)) sf_error ("No Dz= in velocity spline file");
    if (!sf_histfloat (vspline, "Oz", &esc_slow->oz)) sf_error ("No Oz= in velocity spline file");
    if (!sf_histint (vspline, "Nx", &esc_slow->nx)) sf_error ("No Nx= in velocity spline file");
    if (!sf_histfloat (vspline, "Dx", &esc_slow->dx)) sf_error ("No Dx= in velocity spline file");
    if (!sf_histfloat (vspline, "Ox", &esc_slow->ox)) sf_error ("No Ox= in velocity spline file");

    esc_slow->zmin = esc_slow->oz;
    esc_slow->zmax = esc_slow->oz + (esc_slow->nz - 1)*esc_slow->dz;
    esc_slow->xmin = esc_slow->ox;
    esc_slow->xmax = esc_slow->ox + (esc_slow->nx - 1)*esc_slow->dx;

    if (verb) {
        sf_warning ("Splined velocity model dimensions: nz=%d, z=[%g, %g]", esc_slow->nz,
                    esc_slow->zmin, esc_slow->zmax);
        sf_warning ("Splined velocity model dimensions: nx=%d, x=[%g, %g]", esc_slow->nx,
                    esc_slow->xmin, esc_slow->xmax);
    }

    sf_ucharread ((unsigned char*)&(esc_slow->velspline),  sizeof(multi_UBspline_2d_s),
                  vspline);

    esc_slow->aniso = (esc_slow->velspline.num_splines != 1);

    if (verb) {
        sf_warning ("Number of velocity spline coefficients: %lu",
                    esc_slow->velspline.nc/(size_t)sizeof(float));
        if (esc_slow->aniso) {
            if (1 == esc_slow->velspline.num_splines)
                sf_warning ("Assuming isotropic velocity model");
            else if (3 == esc_slow->velspline.num_splines)
                sf_warning ("Assuming anisotropic VTI velocity model");
            else if (4 == esc_slow->velspline.num_splines)
                sf_warning ("Assuming anisotropic TTI velocity model");
            else
                sf_error ("Incorrect splined velocity model");
        }
    }

    stream = sf_filestream (vspline);
    esc_slow->offs = ftello (stream);
#ifdef HAVE_SSE
    if (sizeof(multi_UBspline_2d_s) % 64)
        esc_slow->offs += 64 - (sizeof(multi_UBspline_2d_s) % 64);
#endif

    if (stream != stdin
#ifdef HAVE_SSE
        && 0 == (esc_slow->offs % 64)
#endif
        ) {
        esc_slow->mmaped = (unsigned char*)mmap (NULL, (size_t)esc_slow->offs +
                                                       (size_t)esc_slow->velspline.nc,
                                                 PROT_READ, MAP_SHARED,
                                                 fileno (stream), 0);
        if (esc_slow->mmaped == MAP_FAILED)
            sf_error ("Velocity spline coefficients mmap failed: %s", strerror (errno));
        esc_slow->velspline.coefs = (float*)(esc_slow->mmaped + esc_slow->offs);
    } else {
        if (stream == stdin)
            sf_warning ("Velocity coefficients file appears to be stdin");
        else
            sf_warning ("Velocity coefficients file does not have the right padding");
        sf_warning ("mmap is not possible for velocity coefficients");
#ifndef HAVE_SSE
        esc_slow->velspline.coefs = (float*)sf_ucharalloc (esc_slow->velspline.nc);
#else
        posix_memalign ((void**)&esc_slow->velspline.coefs, 64, esc_slow->velspline.nc);
#endif
        sf_ucharread ((unsigned char*)esc_slow->velspline.coefs,
                      esc_slow->velspline.nc, vspline);
        esc_slow->mmaped = NULL;
    }
    init_einspline ();

    if (verb)
        sf_warning ("Size of splined velocity model: %g Mb", 1e-6*(float)esc_slow->velspline.nc);

    return esc_slow;
}

void sf_esc_slowness2_close (sf_esc_slowness2 esc_slow)
/*< Destroy object >*/
{
    if (esc_slow->mmaped)
        munmap (esc_slow->mmaped, (size_t)esc_slow->offs +
                                  (size_t)esc_slow->velspline.nc);
    else
        free (esc_slow->velspline.coefs);
    free (esc_slow);
}

int sf_esc_slowness2_nz (sf_esc_slowness2 esc_slow)
/*< Number of samples in depth >*/
{
    return esc_slow->nz;
}

float sf_esc_slowness2_oz (sf_esc_slowness2 esc_slow)
/*< Beginning of depth axis >*/
{
    return esc_slow->oz;
}

float sf_esc_slowness2_dz (sf_esc_slowness2 esc_slow)
/*< Depth axis sampling >*/
{
    return esc_slow->dz;
}

int sf_esc_slowness2_nx (sf_esc_slowness2 esc_slow)
/*< Number of samples in lateral direction >*/
{
    return esc_slow->nx;
}

float sf_esc_slowness2_ox (sf_esc_slowness2 esc_slow)
/*< Beginning of lateral axis >*/
{
    return esc_slow->ox;
}

float sf_esc_slowness2_dx (sf_esc_slowness2 esc_slow)
/*< Lateral axis sampling >*/
{
    return esc_slow->dx;
}

/* Return phase slowness and derivatives for phase angle a */
static void sf_esc_slowness2_compute_components (sf_esc_slowness2 esc_slow, float a,
                                                 /* Anistoropic components: V_z^2, V_x^2, et, tilt */
                                                 float vz2, float vx2, float et, float ti,
                                                 /* Derivatives of V_z^2 and V_x^2 */
                                                 float vz2z, float vz2x, float vx2z, float vx2x,
                                                 /* Derivatives of et and tilt */
                                                 float etz, float etx, float tiz, float tix,
                                                 float *s, float *sa, float *sz, float *sx)
{
    float cs, sn2, sqcs, sqsn, prodsq, m2cs, v2;
    float A, B, C, AB, dAdz, dAdx, dAda, dBdz, dBdx, dBda,
          dCdz, dCdx, dCda, dDdz, dDdx, dDda;

    /* For anisotropic case, use the following relations for phase slowness:
       V = 1/S - phase velocity
       V^2 = 1/2*([V_z*cos(c)]^2 + [V_x*sin(c)]^2) +
             1/2*sqrt([[V_z*cos(c)]^2 + [V_x*sin(c)]^2]^2 -
                      8.0*eta/(1.0 + 2.0*eta)*[V_x*V_z*cos(c)*sin(c)]^2) */
    /* Let
       et = -8.0*eta/(1.0 + 2.0*eta)*[V_x*V_z]^2
       c = (a - ti) - angle between tilt axis and the phase vector */
    /* Also
       cos(c) = cos(a)*cos(ti) + sin(a)*sin(ti) = cos(a - ti),
       sin(c)^2 = 1 - cos(c)^2 */
    /* If D = cos(c)^2, (1 - D) = sin(c)^2, E = V_z^2, F = V_x^2, then
       V^2 = 1/2*(E*D + F*(1 - D) +
             1/2*sqrt([E*D + F*(1 - D)]^2 + et*D*(1 - D)) */
    /* For 2-D,
       dD/dx = sin(2*(a - ti))*d(ti)/dx,
       d(1-D)/dx = -sin(2*(a - ti))*d(ti)/dx,
       dD/da = -sin(2*(a - ti)),
       d(1-D)/da = sin(2*(a - ti)) */
    /* Then, V^2 can be written as
       0.5*(A + B + C) */

    if (esc_slow->aniso) {
        /* Find derivatives of A, B, C; restore V^2 with its derivatives, and then
           find values for the phase slowness S and its derivatives */
        cs = cosf (a - ti);
        /* sn = sinf (a - ti); */
        sn2 = sinf (2.0*(a - ti));
        sqcs = cs*cs; /* D */
        sqsn = 1.0 - sqcs;
        prodsq = sqcs*sqsn; /* D*(1 - D) */
        dDdz = sn2*tiz;
        dDdx = sn2*tix;
        dDda = -sn2;
        /* d(D*(1 - D)]/dx = (1 - 2*D)*dD/dx = m2cs*dDdx */
        m2cs = 1.0 - 2.0*sqcs;
        /* Derivatives of A = [V_z*cos(c)]^2 = E*D */
        A = vz2*sqcs;
        dAdz = vz2z*sqcs + vz2*dDdz;
        dAdx = vz2x*sqcs + vz2*dDdx;
        dAda = vz2*dDda;
        /* Derivatives of B = [V_x*sin(c)]^2 = F*(1 - D)*/
        B = vx2*sqsn;
        dBdz = vx2z*sqsn - vx2*dDdz;
        dBdx = vx2x*sqsn - vx2*dDdx;
        dBda = -vx2*dDda;
        AB = A + B; /* E*D + F*(1 - D) */
        /* Derivatives of C = sqrt([E*D + F*(1 - D)]^2 + et*D*(1 - D)) */
        C = sqrt (AB*AB + et*prodsq);
        dCdz = (2.0*AB*(dAdz + dBdz) + etz*prodsq + et*dDdz*m2cs)/(2.0*C); 
        dCdx = (2.0*AB*(dAdx + dBdx) + etx*prodsq + et*dDdx*m2cs)/(2.0*C); 
        dCda = (2.0*AB*(dAda + dBda) + et*dDda*m2cs)/(2.0*C);
        /* V^2 - velocity squared */
        v2 = 0.5*(A + B + C);
        if (s)
            *s = 1.0/sqrt (v2);
        /* Slowness derivatives; use the following equality - 
                    S' = (1/V^2)' = -2/V^3*(V^2)', and
                    V^3 = [V^2]^(1.5) */
        /* (V^2)' = 1/2*(A' + B' + C')
           S' = 1/2*(A' + B' + C')/(-2/[V^2]^(1.5)) =
              = -1/4*[V^2]^(-1.5)*(A' + B' + C') */
        v2 = -0.25*pow (v2, -1.5);
        if (sa)
            *sa = v2*(dAda + dBda + dCda);
        if (sz)
            *sz = v2*(dAdz + dBdz + dCdz);
        if (sx)
            *sx = v2*(dAdx + dBdx + dCdx);
    } else {
        /* For isotropic case, V_z has been converted to
           slowness beforehand */
        if (s)
            *s = vz2;
        if (sa)
            *sa = 0.0;
        if (sz)
            *sz = vz2z;
        if (sx)
            *sx = vz2x;
    }
}

#define SPEPS 1e-3

void sf_esc_slowness2_get_svg (sf_esc_slowness2 esc_slow,
                               float z, float x,
                               float *vel, float *grad)
/*< Set slowness and its derivatives for position z, x, and phase angle a >*/
{
    if (z <= esc_slow->zmin)
        z = esc_slow->zmin + SPEPS*esc_slow->dz;
    if (z >= esc_slow->zmax)
        z = esc_slow->zmax - SPEPS*esc_slow->dz;
    if (x <= esc_slow->xmin)
        x = esc_slow->xmin + SPEPS*esc_slow->dx;
    if (x >= esc_slow->xmax)
        x = esc_slow->xmax - SPEPS*esc_slow->dx;

    eval_multi_UBspline_2d_s_vg (&esc_slow->velspline, x, z, vel, grad);

    /* Zero out angle components for VTI */
    if (esc_slow->velspline.num_splines > 1 &&
        esc_slow->velspline.num_splines < 4) {
        vel[3] = 0.0;
        grad[7] = 0.0;
        grad[6] = 0.0;
    }
}

void sf_esc_slowness2_get_components (sf_esc_slowness2 esc_slow,
                                      float z, float x, float a,
                                      float *s, float *sa, float *sz, float *sx)
/*< Set slowness and its derivatives for position z, x, and phase angle a >*/
{
    float vel[4] = { 0.0, 0.0, 0.0, 0.0 },
          grad[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    sf_esc_slowness2_get_svg (esc_slow, z, x, vel, grad);

    sf_esc_slowness2_compute_components (esc_slow, a,
                                         vel[0], vel[1], vel[2], vel[3], /* vz, vx, et, ti */
                                         grad[1], grad[0], grad[3], grad[2], /* vz/dz, vz/dx, vx/dz, vx/dx */
                                         grad[5], grad[4], grad[7], grad[6], /* et/dz, et/dx, ti/dz, ti/dx */
                                         s, sa, sz, sx);
}

void sf_esc_slowness2_get_coefs (sf_esc_slowness2 esc_slow,
                                 float a, float s, float sa, float sz, float sx,
                                 float *fz, float *fx, float *fa)
/*< Return coefficients in front of dT/dz, dT/dx, and dT/da terms of escape equations >*/
{
    float cs, sn;

    cs = cosf (a);
    sn = sinf (a);
    if (fabsf (cs) == 1.0) sn = 0.0;
    if (fabsf (sn) == 1.0) cs = 0.0;

    if (fz)
        *fz = s*cs + sa*sn;
    if (fx)
        *fx = s*sn - sa*cs;
    if (fa)
        *fa = sx*cs - sz*sn;
}


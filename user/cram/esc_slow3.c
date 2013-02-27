/* 5-D phase-space slowness container for (an)isotropic 3-D medium */
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

#ifndef _esc_slow3_h

typedef struct EscSlowness3 *sf_esc_slowness3;
/* abstract data type */
/*^*/

#endif

struct EscSlowness3 {
    int                  nz, nx, ny;
    float                oz, ox, oy;
    float                dz, dx, dy;
    float                zmin, zmax;
    float                xmin, xmax;
    float                ymin, ymax;
    bool                 aniso;
    size_t               offs;
    unsigned char       *mmaped;
    multi_UBspline_3d_s  velspline;
};
/* concrete data type */

sf_esc_slowness3 sf_esc_slowness3_init (sf_file vspline, bool verb)
/*< Initialize object >*/
{
    FILE *stream;
    bool spl;
    sf_esc_slowness3 esc_slow = (sf_esc_slowness3)sf_alloc (1, sizeof (struct EscSlowness3));

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
    if (!sf_histint (vspline, "Ny", &esc_slow->ny)) sf_error ("No Ny= in velocity spline file");
    if (!sf_histfloat (vspline, "Dy", &esc_slow->dy)) sf_error ("No Dy= in velocity spline file");
    if (!sf_histfloat (vspline, "Oy", &esc_slow->oy)) sf_error ("No Oy= in velocity spline file");

    esc_slow->zmin = esc_slow->oz;
    esc_slow->zmax = esc_slow->oz + (esc_slow->nz - 1)*esc_slow->dz;
    esc_slow->xmin = esc_slow->ox;
    esc_slow->xmax = esc_slow->ox + (esc_slow->nx - 1)*esc_slow->dx;
    esc_slow->ymin = esc_slow->oy;
    esc_slow->ymax = esc_slow->oy + (esc_slow->ny - 1)*esc_slow->dy;

    if (verb) {
        sf_warning ("Splined velocity model dimensions: nz=%d, z=[%g, %g]", esc_slow->nz,
                    esc_slow->zmin, esc_slow->zmax);
        sf_warning ("Splined velocity model dimensions: nx=%d, x=[%g, %g]", esc_slow->nx,
                    esc_slow->xmin, esc_slow->xmax);
        sf_warning ("Splined velocity model dimensions: ny=%d, y=[%g, %g]", esc_slow->ny,
                    esc_slow->ymin, esc_slow->ymax);
    }

    sf_ucharread ((unsigned char*)&(esc_slow->velspline),  sizeof(multi_UBspline_3d_s),
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
            else if (5 == esc_slow->velspline.num_splines)
                sf_warning ("Assuming anisotropic TTI velocity model");
            else
                sf_error ("Incorrect splined velocity model");
        }
    }

    stream = sf_filestream (vspline);
    esc_slow->offs = ftello (stream);
#ifdef HAVE_SSE
    if (sizeof(multi_UBspline_3d_s) % 64)
        esc_slow->offs += 64 - (sizeof(multi_UBspline_3d_s) % 64);
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

void sf_esc_slowness3_close (sf_esc_slowness3 esc_slow)
/*< Destroy object >*/
{
    if (esc_slow->mmaped)
        munmap (esc_slow->mmaped, (size_t)esc_slow->offs +
                                  (size_t)esc_slow->velspline.nc);
    else
        free (esc_slow->velspline.coefs);
    free (esc_slow);
}

int sf_esc_slowness3_nz (sf_esc_slowness3 esc_slow)
/*< Number of samples in depth >*/
{
    return esc_slow->nz;
}

float sf_esc_slowness3_oz (sf_esc_slowness3 esc_slow)
/*< Beginning of depth axis >*/
{
    return esc_slow->oz;
}

float sf_esc_slowness3_dz (sf_esc_slowness3 esc_slow)
/*< Depth axis sampling >*/
{
    return esc_slow->dz;
}

int sf_esc_slowness3_nx (sf_esc_slowness3 esc_slow)
/*< Number of samples in lateral x direction >*/
{
    return esc_slow->nx;
}

float sf_esc_slowness3_ox (sf_esc_slowness3 esc_slow)
/*< Beginning of lateral x axis >*/
{
    return esc_slow->ox;
}

float sf_esc_slowness3_dx (sf_esc_slowness3 esc_slow)
/*< Lateral x axis sampling >*/
{
    return esc_slow->dx;
}

int sf_esc_slowness3_ny (sf_esc_slowness3 esc_slow)
/*< Number of samples in lateral y direction >*/
{
    return esc_slow->ny;
}

float sf_esc_slowness3_oy (sf_esc_slowness3 esc_slow)
/*< Beginning of lateral y axis >*/
{
    return esc_slow->oy;
}

float sf_esc_slowness3_dy (sf_esc_slowness3 esc_slow)
/*< Lateral y axis sampling >*/
{
    return esc_slow->dy;
}

/* Return phase slowness and derivatives for phase angle a */
static void sf_esc_slowness3_compute_components (sf_esc_slowness3 esc_slow, float b, float a,
                                                 /* Anistoropic components: V_z^2, V_x^2, et, tilt */
                                                 float vz2, float vx2, float et, float bti, float ati,
                                                 /* Derivatives of V_z^2 and V_x^2 */
                                                 float vz2z, float vz2x, float vz2y,
                                                 float vx2z, float vx2x, float vx2y,
                                                 /* Derivatives of et and tilt angles */
                                                 float etz, float etx, float ety,
                                                 float btiz, float btix, float btiy,
                                                 float atiz, float atix, float atiy,
                                                 float *s, float *sb, float *sa,
                                                 float *sz, float *sx, float *sy)
{
    float v2=0.0;
/*  float cs, sn, sn2, sqcs, sqsn, prodsq, m2cs, v2;
    float A, B, C, AB, dAdz, dAdx, dAda, dBdz, dBdx, dBda,
          dCdz, dCdx, dCda, dDdz, dDdx, dDda;*/

    /* For anisotropic case, use the following relations for phase slowness:
       V = 1/S - phase velocity
       V^2 = 1/2*([V_z*cos(c)]^2 + [V_x*sin(c)]^2) +
             1/2*sqrt([[V_z*cos(c)]^2 + [V_x*sin(c)]^2]^2 -
                      8.0*eta/(1.0 + 2.0*eta)*[V_x*V_z*cos(c)*sin(c)]^2) */
    /* Let
       et = -8.0*eta/(1.0 + 2.0*eta)*[V_x*V_z]^2
       c - angle between tilt axis and the phase vector */
    /* Also, from the dot product between the phase vector and the tilt
       cos(c) = cos(b)*cos(bt) + sin(b)*cos(a)*sin(bt)*cos(at) +
                               + sin(b)*sin(a)*sin(bt)*sin(at) =
              = cos(b)*cos(bt) + sin(b)*sin(bt)*cos(a - at)
       sin(c)^2 = 1 - cos(c)^2 */
    /* For 3-D,
       d(cos(c))/dx = sin(b)*cos(bt)*cos(a - at)*d(bt)/dx +
                      sin(b)*sin(bt)*sin(a - at)*d(at)/dx - cos(b)*sin(bt)*d(bt)/dx
       d(cos(c))/db = cos(b)*sin(bt)*cos(a - at) - sin(b)*cos(bt)
       d(cos(c))/da = -sin(b)*sin(bt)*sin(a - at) */
    /* If D = cos(c)^2, (1 - D) = sin(c)^2, E = V_z^2, F = V_x^2, then
       V^2 = 1/2*(E*D + F*(1 - D) +
             1/2*sqrt([E*D + F*(1 - D)]^2 + et*D*(1 - D)) */
    /* Then, V^2 can be written as
       0.5*(A + B + C) */
    /* A = E*D, B = F*(1 - D),
       D = sqrt([E*D + F*(1 - D)]^2 + et*D*(1 - D)) */

    if (esc_slow->aniso) {
        /* Find derivatives of A, B, C; restore V^2 with its derivatives, and then
           find values for the phase slowness S and its derivatives */
        /* Code up 3-D anisotropic case based on 2-D code and the relations
           above */
        if (s)
            *s = 1.0/sqrt (v2);
/*      if (sb)
            *sb = v2*(dAdb + dBdb + dCdb);
        if (sa)
            *sa = v2*(dAda + dBda + dCda);
        if (sz)
            *sz = v2*(dAdz + dBdz + dCdz);
        if (sx)
            *sx = v2*(dAdx + dBdx + dCdx);
        if (sy)
            *sy = v2*(dAdy + dBdy + dCdy);*/
    } else {
        /* For isotropic case, V_z has been converted to
           slowness beforehand */
        if (s)
            *s = vz2;
        if (sb)
            *sb = 0.0;
        if (sa)
            *sa = 0.0;
        if (sz)
            *sz = vz2z;
        if (sx)
            *sx = vz2x;
        if (sy)
            *sy = vz2y;
    }
}

#define SPEPS 1e-3

void sf_esc_slowness3_get_svg (sf_esc_slowness3 esc_slow,
                               float z, float x, float y,
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
    if (y <= esc_slow->ymin)
        y = esc_slow->ymin + SPEPS*esc_slow->dy;
    if (y >= esc_slow->ymax)
        y = esc_slow->ymax - SPEPS*esc_slow->dy;

    eval_multi_UBspline_3d_s_vg (&esc_slow->velspline, y, x, z, vel, grad);

    /* Zero out angle components for VTI */
    if (esc_slow->velspline.num_splines > 1 &&
        esc_slow->velspline.num_splines < 4) {
        vel[4] = 0.0;
        vel[5] = 0.0;
        grad[9] = 0.0; grad[10] = 0.0; grad[11] = 0.0;
        grad[12] = 0.0; grad[13] = 0.0; grad[14] = 0.0;
    }

}

void sf_esc_slowness3_get_components (sf_esc_slowness3 esc_slow,
                                      float z, float x, float y, float b, float a,
                                      float *s, float *sb, float *sa,
                                      float *sz, float *sx, float *sy)
/*< Set slowness and its derivatives for position z, x, and phase angle a >*/
{
    float vel[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 },
          grad[15] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    sf_esc_slowness3_get_svg (esc_slow, z, x, y, vel, grad);

    sf_esc_slowness3_compute_components (esc_slow, b, a,
                                         vel[0], vel[1], vel[2], vel[3], vel[4], /* vz, vx, et, bti, ati */
                                         grad[2], grad[1], grad[0], /* vz/dz, vz/dx, vz/dy */
                                         grad[5], grad[4], grad[3], /* vx/dz, vx/dz, vx/dy */
                                         grad[8], grad[7], grad[6], /* et/dz, et/dz, et/dy */
                                         grad[11], grad[10], grad[9], /* bti/dz, bti/dz, bti/dy */
                                         grad[14], grad[13], grad[12], /* ati/dz, ati/dz, ati/dy */
                                         s, sb, sa, sz, sx, sy);
}

void sf_esc_slowness3_get_coefs (sf_esc_slowness3 esc_slow,
                                 float b, float a, float s,
                                 float sz, float sx, float sy, float sa, float sb,
                                 float *fz, float *fx, float *fy, float *fa, float *fb)
/*< Return coefficients in front of dT/dz, dT/dx, dT/dy, and dT/da, dT/db terms of escape equations 
>*/
{
    float csa, csb, sna, snb;

    csa = cosf (a);
    sna = sinf (a);
    csb = cosf (b);
    snb = sinf (b);
    if (fabsf (csa) == 1.0) sna = 0.0;
    if (fabsf (sna) == 1.0) csa = 0.0;
    if (fabsf (csb) == 1.0) snb = 0.0;
    if (fabsf (snb) == 1.0) csb = 0.0;
    if (snb < 1e-5 && snb >= 0.0)
        snb += 1e-5;
    if (snb > -1e-5 && snb <= 0.0)
        snb -= 1e-5;

    if (fz)
        *fz = s*csb + sb*snb;
    if (fx)
        *fx = -s*csa*snb + sb*csa*csb - sa*sna/snb;
    if (fy)
        *fy = -s*sna*snb + sb*sna*csb + sa*csa/snb;
    if (fb)
        *fb = -csb*(sx*csa + sy*sna) - sz*snb;
    if (fa)
        *fa = -(sy*csa - sx*sna)/snb;
}


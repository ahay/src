/* Adjoint state method for 2-D eikonal equation. */
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

#include <math.h>
#include <rsf.h>

void sf_init_lambda_2d_sweep (float *dt, float **l, float **grtz, float **grtx,
                              int n1, int n2, float d1, float d2)
/*< Initialize lambda on the top of the model >*/
{
    int i, j;
    float grz;

    for (i = 0; i < n2; i++) {
        grz = grtz[i][0];
        if (grz != 0.0)
            l[i][0] = -dt[i] / grz;
        else
            l[i][0] = 0.0;
    }

    for (j = 0; j < n2; j++) {
        for (i = 1; i < n1; i++) {
            l[j][i] = 0.0;
        }
    }
}

static void sf_fast_lsweep_2d_stencil (float **l, float **grtz, float **grtx,
                                       int i, int j, int n1, int n2,
                                       float d1, float d2) {
    float fz, fx, dd, grz, grx;
    float fzmax, fzmin, fxmax, fxmin;
    int im, jm, ip, jp;

    grz = grtz[j][i];
    grx = grtx[j][i];

    fz = grz/d1;
    fx = grx/d2;
    fzmin = SF_MAX(-fz,0.);
    fzmax = SF_MAX(fz,0.);
    fxmin = SF_MAX(-fx,0.);
    fxmax = SF_MAX(fx,0.);

    /* Diagonal term */
    dd = (fxmax + fxmin + fzmax + fzmin);

    ip = i + 1; im = i - 1;
    jp = j + 1; jm = j - 1;
    if (i == 0) {
        if (fzmin != 0.0) return; /* Inflow */
        im = 0;
    } else if ((n1 - 1) == i) {
        if (fzmax != 0.0) return; /* Inflow */
        ip = n1 - 1;
    }
    if (j == 0) {
        if (fxmin != 0.0) return; /* Inflow */
        jm = 0;
    } else if ((n2 - 1) == j) {
        if (fxmax != 0.0) return; /* Inflow */
        jp = n2 - 1;
    }

    /* Choose between left-sided and right-sided second
       derivative depending on the coefficients */
    if (fzmin != 0.0 && i != 0) {
        dd += (grz - grtz[j][im])/d1;
    } else if (fzmax != 0.0 && i != (n1 - 1)) {
        dd += (grtz[j][ip] - grz)/d1;
    }
    if (fxmin != 0.0 && j != 0) {
        dd += (grx - grtx[jm][i])/d2;
    } else if (fxmax != 0.0 && j != (n1 - 1)) {
        dd += (grtx[jp][i] - grx)/d2;
    }

    l[j][i] = (fzmax*l[j][ip] + fzmin*l[j][im] +
               fxmax*l[jp][i] + fxmin*l[jm][i])/dd;
}

static int sf_get_horiz2d_iz (sf_eno horiz, float oh, float dh, int nh, int ix,
                              int n1, int n2, float o1, float o2, float d1, float d2) {
    int iz;
    float x;
    x = o2 + ix*d2;
    if (x < oh) x = oh;
    if (x > (oh + (nh - 1)*dh)) x = oh + (nh - 1)*dh;
    x = (x - oh)/dh;
    /* Find horizon location in z for this x */
    sf_eno_apply (horiz, (int)x, x - (int)x, &x, NULL, FUNC);
    iz = (x - o1)/d1;
    iz -= 2;
    return iz < 0 ? 0 : iz;
}

void sf_run_lambda_2d_sweep (float **l, float **grtz, float **grtx, int niter,
                             int n1, int n2, float o1, float o2, float d1, float d2,
                             sf_eno horiz, float oh, float dh, int nh)
/*< Run lambda calculation over whole domain >*/
{
    int i, j, k;
    int imin, imax=0, jmin, jmax;
    float *hor = NULL;

    if (horiz) {
        hor = sf_floatalloc (n2);
        for (j = 0; j < n2; j++) {
            hor[j] = sf_get_horiz2d_iz (horiz, oh, dh, nh, j, n1, n2,
                                        o1, o2, d1, d2);
            imax = hor[j];
            for (i = 0; i < imax; i++)
                l[j][i] = 0.0;
        }
    }

    /* Sweeps */
    for (k = 0; k < niter; k++) {
        jmin = 0; jmax = n2 - 1;
        imin = 0; imax = n1 - 1;
        for (j = jmin; j <= jmax; j++) {
            if (horiz)
                imax = hor[j];
            for (i = imin; i <= imax; i++) {
                sf_fast_lsweep_2d_stencil (l, grtz, grtx, i, j,
                                           n1, n2, d1, d2);
            }
        }
        for (j = jmax; j >= jmin; j--) {
            if (horiz)
                imax = hor[j];
            for (i = imin; i <= imax; i++) {
                sf_fast_lsweep_2d_stencil (l, grtz, grtx, i, j,
                                           n1, n2, d1, d2);
            }
        }
        for (j = jmax; j >= jmin; j--) {
            if (horiz)
                imax = hor[j];
            for (i = imax; i >= imin; i--) {
                sf_fast_lsweep_2d_stencil (l, grtz, grtx, i, j,
                                           n1, n2, d1, d2);
            }
        }
        for (j = jmin; j <= jmax; j++) {
            if (horiz)
                imax = hor[j];
            for (i = imax; i >= imin; i--) {
                sf_fast_lsweep_2d_stencil (l, grtz, grtx, i, j,
                                           n1, n2, d1, d2);
            }
        }
    }

    if (horiz)
        free (hor);
}


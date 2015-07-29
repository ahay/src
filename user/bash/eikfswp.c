/* Fast sweeping eikonal solver (2-D/3-D). */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
/*^*/

#include "eikfswp.h"

bool sf_init_fast_sweep (float *t,
                         int n3, int n2, int n1,
                         float o3, float o2, float o1,
                         float d3, float d2, float d1,
                         float shotz, float shoty, float shotx)
/*< initialize >*/
{
    int i, n123;
    int x = (int)((shotx - o3) / d3 + 0.5f);
    int y = (int)((shoty - o2) / d2 + 0.5f);
    int z = (int)((shotz - o1) / d1 + 0.5f);

    if (NULL == t)
        return false;

    if (x < 0 || x >= n3 ||
        y < 0 || y >= n2 ||
        z < 0 || z >= n1)
        return false;

    n123 = n1 * n2 * n3;

    for (i = 0; i < n123; i++)
        t[i] = SF_HUGE;

    t[x * n1 * n2 + y * n1 + z] = 0.0;

    return true;
}

static void sf_fast_sweep_2d_stencil (float *t, float *s,
                                      int i, int j, int n1, int n2,
                                      float d1, float d2) {

    float Ujmin = 0 == j ? t[n1 + i] : ((n2 - 1) == j ? t[(n2 - 2) * n1 + i] :
                   (t[(j - 1) * n1 + i] < t[(j + 1) * n1 + i] ?
                    t[(j - 1) * n1 + i] : t[(j + 1) * n1 + i]));
    float Uimin = 0 == i ? t[j * n1 + 1] : ((n1 - 1) == i ? t[j * n1 + n1 - 2] :
                   (t[j * n1 + i - 1] < t[j * n1 + i + 1] ?
                    t[j * n1 + i - 1] : t[j * n1 + i + 1]));
    float Uj, Ui, d12, d22;
    float sl = s[j * n1 + i];
    float tij = t[j * n1 + i], ttemp = tij;

    if (SF_HUGE == Ujmin && SF_HUGE == Uimin)
        return;

    Uj = Ujmin + d2 * sl;
    Ui = Uimin + d1 * sl;

    if (Ui <= Ujmin) {
        ttemp = Ui;
    } else if (Uj <= Uimin) {
        ttemp = Uj;
    } else {
       d12 = d1 * d1;
       d22 = d2 * d2;
       ttemp = (d1 * d2 * sqrt ((d12 + d22) * sl * sl - (Ujmin - Uimin) * (Ujmin - Uimin)) +
                Uimin * d22 + Ujmin * d12) / (d12 + d22);
    }

    if (ttemp < tij)
        t[j * n1 + i] = ttemp;
}

#define SWAPF(a,b,t) t = b; b = a; a = t;

static void sf_fast_sweep_3d_stencil (float *t, float *s,
                                      int i, int j, int k,
                                      int n1, int n2, int n3,
                                      float d1, float d2, float d3) {
    int n12 = n1 * n2;
    float Ukmin = 0 == k ? t[n12 + j * n1 + i] : ((n3 - 1) == k ? t[(n3 - 2) * n12 + j * n1 + i] :
                   (t[(k - 1) * n12 + n1 * j + i] < t[(k + 1) * n12 + n1 * j + i] ?
                    t[(k - 1) * n12 + n1 * j + i] : t[(k + 1) * n12 + n1 * j + i]));
    float Ujmin = 0 == j ? t[k * n12 + n1 + i] : ((n2 - 1) == j ? t[k * n12 + (n2 - 2) * n1 + i] :
                   (t[k * n12 + (j - 1) * n1 + i] < t[k * n12 + (j + 1) * n1 + i] ?
                    t[k * n12 + (j - 1) * n1 + i] : t[k * n12 + (j + 1) * n1 + i]));
    float Uimin = 0 == i ? t[k * n12 + j * n1 + 1] : ((n1 - 1) == i ? t[k * n12 + j * n1 + n1 - 2] :
                   (t[k * n12 + j * n1 + i - 1] < t[k * n12 + j * n1 + i + 1] ?
                    t[k * n12 + j * n1 + i - 1] : t[k * n12 + j * n1 + i + 1]));
    float d12, d22, d32, u[3], d[3];
    float sl = s[k * n12 + j * n1 + i];
    float tij = t[k * n12 + j * n1 + i], ttemp = tij, swap;

    if (SF_HUGE == Ukmin && SF_HUGE == Ujmin && SF_HUGE == Uimin)
        return;

    u[0] = Uimin;
    u[1] = Ujmin;
    u[2] = Ukmin;
    d[0] = d1;
    d[1] = d2;
    d[2] = d3;

    if (u[0] > u[2]) {
        SWAPF (u[0], u[2], swap)
        SWAPF (d[0], d[2], swap)
    }
    if (u[0] > u[1]) {
        SWAPF (u[0], u[1], swap)
        SWAPF (d[0], d[1], swap)
    }
    if (u[1] > u[2]) {
        SWAPF (u[1], u[2], swap)
        SWAPF (d[1], d[2], swap)
    }
    ttemp = u[0] + d[0] * sl;
    if (ttemp > u[1]) {
        d12 = d[0] * d[0];
        d22 = d[1] * d[1];
        ttemp = (d[0] * d[1] * sqrt ((d12 + d22) * sl * sl - (u[0] - u[1]) * (u[0] - u[1])) +
                 u[0] * d22 + u[1] * d12) / (d12 + d22);
        if (ttemp > u[2]) {
            d32 = d[2] * d[2];
            ttemp = (d[0] * d[1] * d[2] * sqrt ((d22 * d32 + d12 * d32 + d12 * d22) * sl * sl
                                                 + 2.0 * d32 * u[0] * u[1]
                                                 + 2.0 * d22 * u[0] * u[2]
                                                 + 2.0 * d12 * u[1] * u[2] 
                                                 - (d32 + d22) * u[0] * u[0]
                                                 - (d32 + d12) * u[1] * u[1]
                                                 - (d22 + d12) * u[2] * u[2])
                     + d12 * d22 * u[2] + d12 * d32 * u[1] + d22 * d32 * u[0])
                    / (d22 * d32 + d12 * d32 + d12 * d22);
        }
    }

    if (ttemp < tij)
        t[k * n12 + j * n1 + i] = ttemp;
}

void sf_run_fast_sweep (float *t, float *s, int niter,
                        int n3, int n2, int n1,
                        float o3, float o2, float o1,
                        float d3, float d2, float d1)
/*< run sweeps >*/
{
    int i, j ,k = 0, l = 0;

    if (n3 <= 1) { /* 2-D */
        for (l = 0; l < niter; l++) {
            for (j = 0; j < n2; j++) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil (t, s, i, j, n1, n2,
                                              d1, d2);
                }
            }
            for (j = n2 - 1; j >= 0; j--) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil (t, s, i, j, n1, n2,
                                              d1, d2);
                }
            }
            for (j = n2 - 1; j >= 0; j--) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil (t, s, i, j, n1, n2,
                                              d1, d2);
                }
            }
            for (j = 0; j < n2; j++) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil (t, s, i, j, n1, n2,
                                              d1, d2);
                }
            }
        }
    } else { /* 3-D */
        for (l = 0; l < niter; l++) {
            if (niter > 1)
                sf_warning ("Iteration %d of %d", l + 1, niter);
            for (k = 0; k < n3; k++) {
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 1 */
            for (k = n3 - 1; k >= 0; k--) {
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 2 */
            for (k = 0; k < n3; k++) {
                for (j = n2 - 1; j >= 0; j--) {
                    for (i = 0; i < n1; i++) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 3 */
            for (k = n3 - 1; k >= 0; k--) {
                for (j = n2 - 1; j >= 0; j--) {
                    for (i = 0; i < n1; i++) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 4 */
            for (k = 0; k < n3; k++) {
                for (j = 0; j < n2; j++) {
                    for (i = n1 - 1; i >= 0; i--) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 5 */
            for (k = n3 - 1; k >= 0; k--) {
                for (j = 0; j < n2; j++) {
                    for (i = n1 - 1; i >= 0; i--) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 6 */
            for (k = 0; k < n3; k++) {
                for (j = n2 - 1; j >= 0; j--) {
                    for (i = n1 - 1; i >= 0; i--) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 7 */
            for (k = n3 - 1; k >= 0; k--) {
                for (j = n2 - 1; j >= 0; j--) {
                    for (i = n1 - 1; i >= 0; i--) {
                        sf_fast_sweep_3d_stencil (t, s, i, j, k,
                                                  n1, n2, n3,
                                                  d1, d2, d3);
                    }
                }
            } /* 8 */
        }
    }
}

void sf_insert_horiz_2d (float *v, int n2, int n1,
                         float o2, float o1, float d2, float d1,
                         sf_eno horiz, float oh, float dh, int nh)
/*< modify velocity model acording to horizon interface >*/
{
    int ix, iz;
    float x;
    for (ix = 0; ix < n2; ix++) {
        x = o2 + ix*d2;
        if (x < oh) x = oh;
        if (x > (oh + (nh - 1)*dh)) x = oh + (nh - 1)*dh;
        x = (x - oh)/dh;
        /* Find horizon location in z for this x */
        sf_eno_apply (horiz, (int)x, x - (int)x, &x, NULL, FUNC);
        iz = (x - o1)/d1;
        if (iz < 0 || iz >= n1) continue;
        x = v[ix*n1 + iz];
        iz++;
        /* Spread velocity value fronm above the horizon below it */
        while (iz < n1) {
            v[ix*n1 + iz] = x;
            iz++;
        }
    }
}

void sf_extract_horiz_2d (float *t, float *s, int n2, int n1,
                          float o2, float o1, float d2, float d1,
                          sf_eno horiz, float oh, float dh, int nh)
/*< modify velocity model acording to horizon interface >*/
{
    int ix, iz, iiz;
    float x, th, izf;

    for (ix = 0; ix < n2; ix++) { /* Horizon points */
        x = o2 + ix*d2;
        if (x < oh) x = oh;
        if (x > (oh + (nh - 1)*dh)) x = oh + (nh - 1)*dh;
        x = (x - oh)/dh;
        /* Find horizon location in z for this x */
        sf_eno_apply (horiz, (int)x, x - (int)x, &x, NULL, FUNC);
        iz = (x - o1)/d1;
        if (iz < 0 || iz >= n1) continue;
        izf = (x - o1)/d1 - (float)iz;
        /* Traveltime on the horizon at x */
        if (iz < (n1 - 1))
            th = (1.0 - izf)*t[ix*n1 + iz] + izf*t[ix*n1 + iz + 1];
        else
            th = t[ix*n1 + iz];
        /* Advance to the nearest grid location */
        t[ix*n1 + iz] = th + s[ix*n1 + iz]*d1*izf;
        iiz = iz - 1;
        /* Unset all value above and below */
        while (iiz >= 0) {
            t[ix*n1 + iiz] = SF_HUGE;
            iiz--;
        }
        iiz = iz + 1;
        while (iiz < n1) {
            t[ix*n1 + iiz] = SF_HUGE;
            iiz++;
        }
    }
}


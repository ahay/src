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

#include<math.h>
#include<rsf.h>

bool sf_init_fast_sweep (float *t,
                         int n3, int n2, int n1,
                         float o3, float o2, float o1,
                         float d3, float d2, float d1,
                         float shotz, float shoty, float shotx) {
    if (NULL == t)
        return false;

    int i, n123;
    int x = (int)((shotx - o3) / d3 + 0.5f);
    int y = (int)((shoty - o2) / d2 + 0.5f);
    int z = (int)((shotz - o1) / d1 + 0.5f);

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
                        float d3, float d2, float d1) {
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

int main (int argc,char* argv[]) {
    int n1, n2, n3, i, nshot, is, n123, ndim, niter;
    float o1, o2, o3, d1, d2, d3;
    float **s, *t, *v;
    char *sfile;
    bool isvel;
    sf_file vel, time, shots;

    sf_init (argc, argv);
    vel = sf_input ("in");
    time = sf_output ("out");

    if (SF_FLOAT != sf_gettype (vel))
        sf_error("Need float input");
    if (!sf_histint (vel, "n1", &n1)) sf_error ("No n1= in input");
    if (!sf_histint (vel, "n2", &n2)) sf_error ("No n2= in input");
    if (!sf_histint (vel, "n3", &n3)) n3 = 1;

    if (!sf_histfloat (vel, "d1", &d1)) sf_error ("No d1= in input");
    if (!sf_histfloat (vel, "d2", &d2)) sf_error ("No d2= in input");
    if (!sf_histfloat (vel, "d3", &d3)) d3 = d2;

    if (!sf_histfloat (vel, "o1", &o1)) o1 = 0.;
    if (!sf_histfloat (vel, "o2", &o2)) o2 = 0.;
    if (!sf_histfloat (vel, "o3", &o3)) o3 = 0.;

    if (!sf_getbool ("vel", &isvel)) isvel = true;
    /* if y, the input is velocity; n - slowness */

    if (!sf_getint ("niter", &niter)) niter = 2;
    /* number of sweeping iterations */

    sfile = sf_getstring ("shotfile");
    /* File with shot locations (n2=number of shots, n1=3) */

    if (NULL != sfile) {
        shots = sf_input ("shotfile");

        if (SF_FLOAT != sf_gettype (shots)) 
            sf_error ("Need float shotfile");
        if (!sf_histint (shots, "n2", &nshot))
            sf_error ("No n2= in shotfile");
        if (!sf_histint (shots, "n1", &ndim) || ndim != 3)
            sf_error ("Need n1=3 in shotfile");

        s = sf_floatalloc2 (ndim, nshot);
        sf_floatread (s[0], nshot * ndim, shots);

        sf_putint (time, "n4", nshot);
        free (sfile);
    } else {
        nshot = 1;
        ndim = 3;

        s = sf_floatalloc2 (ndim, nshot);

        if (!sf_getfloat ("zshot", &s[0][0])) s[0][0] = 0.; 
        /* Shot location (used if no shotfile) */
        if (!sf_getfloat ("yshot", &s[0][1])) s[0][1] = o2 + 0.5*(n2 - 1)*d2;
        if (!sf_getfloat ("xshot", &s[0][2])) s[0][2] = o3 + 0.5*(n3 - 1)*d3;

        sf_warning ("Shooting from zshot=%g yshot=%g xshot=%g",
                    s[0][0], s[0][1], s[0][2]);
    }

    n123 = n1 * n2 * n3;

    t  = sf_floatalloc (n123);
    v  = sf_floatalloc (n123);

    sf_floatread (v, n123, vel);
   /* transform velocity to slowness squared */
    if (isvel) {
        for (i = 0; i < n123; i++) {
            v[i] = 1. / v[i];
        }
    }

    sf_warning ("Performing %d-D sweeps", 2 + (n3 > 1));

    /* loop over shots */
    for (is = 0; is < nshot; is++) {
        if (nshot > 1)
            sf_warning ("Calculating shot %d of %d", is + 1, nshot);
        if (false == sf_init_fast_sweep (t,
                                         n3, n2, n1,
                                         o3, o2, o1,
                                         d3, d2, d1,
                                         s[is][0], s[is][1], s[is][2]))
            sf_error ("Incorrect shot location");

        sf_run_fast_sweep (t, v, niter,
                           n3, n2, n1,
                           o3, o2, o1,
                           d3, d2, d1);

        sf_floatwrite (t, n123, time);
    }
    
    exit (0);
}

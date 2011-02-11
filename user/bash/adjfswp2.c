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

void sf_init_lambda_2d_sweep (float *dt, float **l, float **t,
                              int n2, int n1, float d2, float d1)
/*< Initialize lambda on the top of the model >*/
{
    int i, j;

    for (i = 0; i < n2; i++) {
        l[i][0] = dt[i] /
                  -(-t[i][2] + 4.0*t[i][1] - 3.0*t[i][0])*d1*2.0;
        if (isinf (l[i][0]) || isnan (l[i][0]))
            l[i][0] = 0.0;
    }

    for (j = 0; j < n2; j++) {
        for (i = 1; i < n1; i++) {
            l[j][i] = 0.0;
        }
    }
}

static void sf_fast_lsweep_2d_stencil2 (float **l, sf_eno2 dt,
                                        int i, int j, int n1, int n2,
                                        float d1, float d2) {
    float gr[2], gr2[2];
    float fz, fx, dd;
    float fzmax, fzmin, fxmax, fxmin;
    int im, jm, ip, jp;

    sf_eno2_apply (dt, i, j, 0., 0., NULL, gr, DER);

    fz = gr[0]/d1/d1;
    fx = gr[1]/d2/d2;

    fzmin = SF_MAX(-fz,0.);
    fzmax = SF_MAX(fz,0.);

    fxmin = SF_MAX(-fx,0.);
    fxmax = SF_MAX(fx,0.);

    /* Diagonal term */
    dd = (fxmax + fxmin + fzmax + fzmin);

    ip = i + 1; im = i - 1;
    jp = j + 1; jm = j - 1;
    if (i == 0) {
        im = 0;
    } else if ((n1 - 1) == i) {
        ip = n1 - 1;
    }
    if (j == 0) {
        jm = 0;
    } else if ((n2 - 1) == j) {
        jp = n2 - 1;
    }

    /* Choose between left-sided and right-sided second
       derivative depending on the coefficients */
    if (fzmin != 0.0 && i != 0) {
        sf_eno2_apply (dt, i - 1, j, 0., 0., NULL, gr2, DER);
        dd += (gr[0] - gr2[0])/d1/d1;
    } else if (fzmax != 0.0 && i != (n1 - 1)) {
        sf_eno2_apply (dt, i + 1, j, 0., 0., NULL, gr2, DER);
        dd += (gr2[0] - gr[0])/d1/d1;
    }
    if (fxmin != 0.0 && j != 0) {
        sf_eno2_apply (dt, i, j - 1, 0., 0., NULL, gr2, DER);
        dd += (gr[1] - gr2[1])/d2/d2;
    } else if (fxmax != 0.0 && j != (n1 - 1)) {
        sf_eno2_apply (dt, i, j + 1, 0., 0., NULL, gr2, DER);
        dd += (gr2[1] - gr[1])/d2/d2;
    }

    l[j][i] = (fzmax*l[j][ip] + fzmin*l[j][im] +
               fxmax*l[jp][i] + fxmin*l[jm][i])/dd;
}

/* Apply lambda calculation to one point */
static void sf_fast_lsweep_2d_stencil (float **t, float **l,
                                       int i, int j, int n1, int n2,
                                       float d1, float d2) {
    float ap, am, bp, bm, app, apm, amp, amm, bpp, bpm, bmp, bmm;
    float tij, tim, tip, tjm, tjp, lij, lim, lip, ljm, ljp;

    tij = t[j][i];
    tim = t[j][i - 1];
    tip = t[j][i + 1];
    tjm = t[j - 1][i];
    tjp = t[j + 1][i];

    lij = l[j][i];
    lim = l[j][i - 1];
    lip = l[j][i + 1];
    ljm = l[j - 1][i];
    ljp = l[j + 1][i];

    ap = -(tip - tij) / d1;
    am = -(tij - tim) / d1;
    bp = -(tjp - tij) / d2;
    bm = -(tij - tjm) / d2;

    app = 0.5*(ap + (ap < 0 ? -ap : ap));
    apm = 0.5*(am + (am < 0 ? -am : am));
    amp = 0.5*(ap - (ap < 0 ? -ap : ap));
    amm = 0.5*(am - (am < 0 ? -am : am));
    bpp = 0.5*(bp + (bp < 0 ? -bp : bp));
    bpm = 0.5*(bm + (bm < 0 ? -bm : bm));
    bmp = 0.5*(bp - (bp < 0 ? -bp : bp));
    bmm = 0.5*(bm - (bm < 0 ? -bm : bm));

    lij = ((apm*lim - amp*lip)/d1 + (bpm*ljm - bmp*ljp)/d2)  /
          ((app - amm)/d1 + (bpp - bmm)/d2);
    l[j][i] = lij;
}

void sf_run_lambda_2d_sweep (float **t, float **l,
                             int niter, int n2, int n1,
                             float d2, float d1)
/*< Run lambda calculation over whole domain >*/
{
    int i, j, k = 0, is = -1, js = -1;
    float sl = 0.0;

    /* Find shot location */
    for (j = 0; j < n2; j++) {
        for (i = 0; i < n1; i++) {
            if (t[j][i] == 0.0) {
                js = j;
                is = i;
                j = n2 - 1;
                break;
            }
        }
    }

    /* Sweeps */
    for (k = 0; k < niter; k++) {
        for (j = 1; j < (n2 - 1); j++) {
            for (i = 1; i < (n1 - 1); i++) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil (t, l, i, j, n1, n2,
                                               d1, d2);
            }
        }
        for (j = n2 - 2; j > 0; j--) {
            for (i = 1; i < (n1 - 1); i++) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil (t, l, i, j, n1, n2,
                                               d1, d2);
            }
        }
        for (j = n2 - 2; j > 0; j--) {
            for (i = n1 - 2; i > 0; i--) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil (t, l, i, j, n1, n2,
                                               d1, d2);
            }
        }
        for (j = 1; j < (n2 - 1); j++) {
            for (i = n1 - 2; i > 0; i--) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil (t, l, i, j, n1, n2,
                                               d1, d2);
            }
        }
    }

    /* Fill the shot location */
    if (is != -1 && js != -1) {
        i = 0;
        if (is > 0) {
            sl = l[js][is - 1];
            i++;
        }
        if (is < (n1 - 1)) {
            sl = l[js][is + 1];
            i++;
        }
        if (js > 0) {
            sl = l[js - 1][is];
            i++;
        }
        if (js < (n2 - 1)) {
            sl = l[js + 1][is];
            i++;
        }
        sl /= (float)i;
        l[js][is] = sl;
    }
}

void sf_run_lambda_2d_sweep2 (float **t, float **l,
                              int niter, int n2, int n1,
                              float d2, float d1)
/*< Run lambda calculation over whole domain >*/
{
    int i, j, k = 0, is = -1, js = -1;
    float sl = 0.0;
    sf_eno2 dt;

    dt = sf_eno2_init (3, n1, n2);
    sf_eno2_set (dt, t);

    /* Find shot location */
    for (j = 0; j < n2; j++) {
        for (i = 0; i < n1; i++) {
            if (t[j][i] == 0.0) {
                js = j;
                is = i;
                j = n2 - 1;
                break;
            }
        }
    }

    /* Sweeps */
    for (k = 0; k < niter; k++) {
        for (j = 1; j < n2; j++) {
            for (i = 1; i < n1; i++) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil2 (l, dt, i, j, n1, n2,
                                               d1, d2);
            }
        }
        for (j = n2 - 1; j >= 0; j--) {
            for (i = 1; i < n1; i++) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil2 (l, dt, i, j, n1, n2,
                                               d1, d2);
            }
        }
        for (j = n2 - 1; j >= 0; j--) {
            for (i = n1 - 1; i > 0; i--) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil2 (l, dt, i, j, n1, n2,
                                               d1, d2);
            }
        }
        for (j = 0; j < n2; j++) {
            for (i = n1 - 1; i > 0; i--) {
                if (is != i || js != j)
                    sf_fast_lsweep_2d_stencil2 (l, dt, i, j, n1, n2,
                                               d1, d2);
            }
        }
    }

    sf_eno2_close (dt);

    /* Fill the shot location */
    if (is != -1 && js != -1) {
        i = 0;
        if (is > 0) {
            sl = l[js][is - 1];
            i++;
        }
        if (is < (n1 - 1)) {
            sl = l[js][is + 1];
            i++;
        }
        if (js > 0) {
            sl = l[js - 1][is];
            i++;
        }
        if (js < (n2 - 1)) {
            sl = l[js + 1][is];
            i++;
        }
        sl /= (float)i;
        l[js][is] = sl;
    }
}


/* Angle migration of one depth point */
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

#include <rsf.h>

#define DHORDER 3


/* F-D coefficients (FDo9p) from Bogey and Bailly, 2004, JCP */
#define C1 0.841570125482
#define C2 -0.244678631765
#define C3 0.059463584768
#define C4 -0.007650904064

static float sf_cram_bbfd_stencil (float am4, float am3, float am2, float am1,
                                   float ap1, float ap2, float ap3, float ap4) {
    return -C4*am4 -C3*am3 -C2*am2 -C1*am1
           +C1*ap1 +C2*ap2 +C3*ap3 +C4*ap4;
}

void sf_cram_trace_deriv (float *tracein, float *traceout, int n, float d)
/*< First order derivative of a trace[n] >*/
{
    int i;
    float invd = 1.0/d;
    float b = tracein[0], e = tracein[n - 1];
    float bd = tracein[0] - tracein[1];
    float ed = tracein[n - 1] - tracein[n - 2];

    for (i = 0; i < n; i++) {
        if (i < 4)
            traceout[i] = sf_cram_bbfd_stencil (b + (4.0 - (float)i)*bd,
                                                i > 2 ? tracein[i - 3] : b + (3.0 - (float)i)*bd,
                                                i > 1 ? tracein[i - 2] : b + (2.0 - (float)i)*bd,
                                                i > 0 ? tracein[i - 1] : b + (1.0 - (float)i)*bd,
                                                tracein[i + 1], tracein[i + 2],
                                                tracein[i + 3], tracein[i + 4]);
        else if (i > (n - 5))
            traceout[i] = sf_cram_bbfd_stencil (tracein[i - 4], tracein[i - 3],
                                                tracein[i - 2], tracein[i - 1],
                                                i < (n - 1) ? tracein[i + 1] :
                                                e + (1.0 - (float)(n - i - 1))*ed,
                                                i < (n - 2) ? tracein[i + 2] :
                                                e + (2.0 - (float)(n - i - 1))*ed,
                                                i < (n - 3) ? tracein[i + 3] :
                                                e + (3.0 - (float)(n - i - 1))*ed,
                                                e + (4.0 - (float)(n - i - 1))*ed);
        else
            traceout[i] = sf_cram_bbfd_stencil (tracein[i - 4], tracein[i - 3],
                                                tracein[i - 2], tracein[i - 1],
                                                tracein[i + 1], tracein[i + 2],
                                                tracein[i + 3], tracein[i + 4]);
        traceout[i] *= invd;
    }
}

void sf_cram_trace_cint (float *trace, int n)
/*< Causal integration of a trace[n] >*/
{
    int i;

    for (i = 1; i < n; i++)
        trace[i] += trace[i - 1];
}

void sf_cram_trace_acint (float *trace, int n)
/*< Anticausal integrations of a trace[n] >*/
{
    int i;

    for (i = n - 2; i >= 0; i--)
        trace[i] += trace[i + 1];
}

/* Convolution operator - borrowed from SU */
static void sf_cram_trace_conv (int lx, int ifx, float *x,
                                int ly, int ify, float *y,
                                int lz, int ifz, float *z) {
    int ilx = ifx + lx - 1,
        ily = ify + ly - 1,
        ilz = ifz + lz - 1,
        i, j, jlow, jhigh;
    float sum;

    x -= ifx;  y -= ify;  z -= ifz;
    for (i=ifz; i<=ilz; ++i) {
        jlow = i - ily;
        if (jlow < ifx)
            jlow = ifx;
        jhigh = i-ify;
        if (jhigh > ilx)
            jhigh = ilx;
        for (j = jlow, sum = 0.0; j <= jhigh; ++j)
            sum += x[j]*y[i-j];
        z[i] = sum;
    }
}

/* Hilbert transform - borrowed from SU */
#define LHHALF 30       /* half-length of Hilbert transform filter*/
#define LH 2*LHHALF+1   /* filter length must be odd */
void sf_cram_trace_hilbert (int n, float *x, float *y)

/*< Hilbert tracnform of a trace x[n] -> y[n] >*/
{
    static int madeh=0;
    static float h[LH];
    int i;
    float taper;

    /* if not made, make Hilbert transform filter; use Hamming window */
    if (!madeh) {
        h[LHHALF] = 0.0;
        for (i=1; i<=LHHALF; i++) {
            taper = 0.54 + 0.46*cos (M_PI*(float)i/(float)(LHHALF));
            h[LHHALF+i] = taper*(-(float)(i % 2)*2.0/
                         (M_PI*(float)(i)));
            h[LHHALF-i] = -h[LHHALF + i];
        }
        madeh = 1;
    }

    /* convolve Hilbert transform with input array */
    sf_cram_trace_conv (LH, -LHHALF, h, n, 0, x, n, 0, y);
}


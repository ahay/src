/*
  Copyright (C) 2012 KAUST

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

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/* fast fft size */
int fft_size(int n)
/*< >*/
{
    int m;

    while (1) {
        m = n;
        while (!(n % 2)) n /= 2;
        while (!(n % 3)) n /= 3;
        while (!(n % 5)) n /= 5;
        while (!(n % 7)) n /= 7;
        if (1 == n) break;
        else n = m + 1;
    }
    return m;
}

/* update boundary to fit fast fft size */
void fft_expand(int n,
                int *l, /* [1] */
                int *h  /* [1] */)
/*< >*/
{
    int N,d;

    N = n + *l + *h;
    d = fft_size(N) - N;

    *l += d / 2;
    *h += d - d / 2;
}

/* expand 2D model */
void expand(float *q, /* [N2][N1]    */
            const float *p,    /* [n2][n1]    */
            const int *n,      /* {n1,n2} */
            const int *l,      /* {l1,l2} */
            const int *h       /* {h1,h2} */)
/*< >*/
{
    int i1,i2,n1,n2,l1,l2,h1,h2,N1,N2,j;
    float aa,bb;

    N1 = (n1 = n[0]) + (l1 = l[0]) + (h1 = h[0]);
    N2 = (n2 = n[1]) + (l2 = l[1]) + (h2 = h[1]);

    for (i2=0; i2 < n2; i2++) { /* extend axis 1 */
        aa = p[i2*n1];
        bb = p[i2*n1 + n1-1];
        for (i1=0; i1 < N1; i1++) {
            j = (l2+i2) * N1 + i1;
            if      (i1 < l1) 
                q[j] = aa;
            else if (i1 < l1 + n1)
                q[j] = p[i2*n1 + i1-l1];
            else
                q[j] = bb;
        }
    }

    for (i1=0; i1 < N1; i1++) { /* extend axis 2 */
        for (aa = q[l2*N1 + i1], i2=0; i2 < l2; i2++)
            q[i2*N1 + i1] = aa;
        for (bb = q[(l2+n2-1)*N1 + i1], i2=l2+n2; i2 < N2; i2++)
            q[i2*N1 + i1] = bb;
    }
}
/* todo : set corners to zero */

/* compute wavenumbers
   0     to N/2 positive
   N/2+1 to N-1 negative */
void compute_k(float *k, /* [nk] */
               int nk)
/*< >*/
{
    int i;
    float dk;

    for (k[0]=0, dk=1./nk, i=1; i < nk; i++) {
        k[i] = k[i-1] + dk;
        if (nk/2+1 == i) k[i] -= 1.; 
    }
}


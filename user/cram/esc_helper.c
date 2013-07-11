/* Auxiliary functions for escape solvers */
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

/* ================================= Quaternions ================================================ */

void sf_quat_rotmat (float *q, float *M)
/*< Build rotation matrix[3][3] from quaternion q[4] >*/
{
    M[0] = 1.0 - 2.0*q[2]*q[2] - 2.0*q[3]*q[3];
    M[1] = 2.0*q[1]*q[2] - 2.0*q[3]*q[0];
    M[2] = 2.0*q[1]*q[3] + 2.0*q[2]*q[0];
    M[3] = 2.0*q[1]*q[2] + 2.0*q[3]*q[0];
    M[4] = 1.0 - 2.0*q[1]*q[1] - 2.0*q[3]*q[3];
    M[5] = 2.0*q[3]*q[2] - 2.0*q[1]*q[0];
    M[6] = 2.0*q[1]*q[3] - 2.0*q[2]*q[0];
    M[7] = 2.0*q[3]*q[2] + 2.0*q[1]*q[0];
    M[8] = 1.0 - 2.0*q[1]*q[1] - 2.0*q[2]*q[2];
}

void sf_quat_norm (float *q)
/*< Normalize quaternion q[4] >*/
{
    float f = sqrtf (q[0]*q[0] + q[1]*q[1] + 
                     q[2]*q[2] + q[3]*q[3]);
    q[0] /= f;
    q[1] /= f;
    q[2] /= f;
    q[3] /= f;
    return;
}

void sf_quat_vecrot (float *vf, float *vt, float *q)
/*< Find quaternion q[4] from rotation of vector vf[3] to vt[3] >*/
{
    float s, invs, d, c[3];

    /* Dot product */
    d = vf[0]*vt[0] + vf[1]*vt[1] + vf[2]*vt[2];

    if (d >= 1.0) { /* Same vectors */
        q[0] = 1.0;
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 0.0;
        return;
    } else if (d <= -1.0) { /* Opposite vectors */
        /* Cross product with the vertical direction */
        c[0] = 0.0;
        c[1] = vf[2];
        c[2] = -vf[1];
        q[0] = 0.0;
    } else {
        /* If phi is rotation angle, then d = cos(phi),
           d = 2*cos(0.5*phi)^2 - 1, 2*cos(0.5*phi)^2 = d + 1,
           4*cos(0.5*phi)^2 = 2*(d + 1),
           s = 2*cos(0.5*phi) = sqrt(2*(d + 1)), and
           sin(phi) = 2*sin(0.5*phi)*cos(0.5*phi),
           sin(phi) = sin(0.5*phi)*s,
           s = sin(phi)/sin(0.5*phi), and
           vfxvt cross product is sin(phi)*n [n - unit vector],
           then (vfxvt)/s = n*sin(phi)/(sin(phi)/sin(0.5*phi)) = n*sin(0.5*phi)
         */
        s = sqrtf ((1.0 + d)*2.0);
        invs = 1.0/s;
        c[0] = (vf[1]*vt[2] - vf[2]*vt[1])*invs;
        c[1] = (vf[2]*vt[0] - vf[0]*vt[2])*invs;
        c[2] = (vf[0]*vt[1] - vf[1]*vt[0])*invs;
        q[0] = 0.5*s;
    }
    q[1] = c[0];
    q[2] = c[1];
    q[3] = c[2];
    /* Normalize */
    sf_quat_norm (q);
    return;
}

/* ================================= LU matrix decomposition ================================================ */

/* LU decomposition with Doolittle method, original code can
   be found at http://mymathlib.com */

int sf_ludlt_decomposition (float *A, int *pivot, int n)
/*< Find LU decomposition of matrix A[nxn], return
    -1 for singular matrix, otherwise - 0 >*/
{
    int i, j, k;
    float *p_k, *p_row, *p_col = NULL;
    float max;

    /* For each row and column, k = 0, ..., n-1, */
    for (k = 0, p_k = A; k < n; p_k += n, k++) {
        /* Find the pivot row */
        pivot[k] = k;
        max = fabsf (*(p_k + k));
        for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
            if (max < fabsf (*(p_row + k))) {
                max = fabsf (*(p_row + k));
                pivot[k] = j;
                p_col = p_row;
            }
        }
        /* and if the pivot row differs from the current row, then
           interchange the two rows.*/
        if (pivot[k] != k) {
            for (j = 0; j < n; j++) {
                max = *(p_k + j);
                *(p_k + j) = *(p_col + j);
                *(p_col + j) = max;
            }
        }
        /* and if the matrix is singular, return error */
        if (*(p_k + k) == 0.0)
            return -1;
        /* otherwise find the lower triangular matrix
           elements for column k.*/
        for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
            *(p_row + k) /= *(p_k + k);
        }
        /* update remaining matrix */
        for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
            for (j = k+1; j < n; j++) {
                *(p_row + j) -= *(p_row + k) * *(p_k + j);
            }
        }
    }
    return 0;
}

/* Matrix vector multiplication for a LU-decomposed system with Doolittle method,
   original code can be found at http://mymathlib.com */

int sf_ludtl_solve (float *A, float *b, int *pivot, float *x, int n)
/*< Find vecotr x[n] for RHS b[n] and prveiously
    LU-processed matrix A[nxn], return -1 for improper
    matrix, otherwise - 0 >*/
{
    int i, k;
    float *p_k;
    float dum;
    
    /* Solve the linear equation Lx = B for x, where L is a lower
       triangular matrix with an implied 1 along the diagonal.*/
    for (k = 0, p_k = A; k < n; p_k += n, k++) {
        if (pivot[k] != k) {
            dum = b[k];
            b[k] = b[pivot[k]];
            b[pivot[k]] = dum;
        }
        x[k] = b[k];
        for (i = 0; i < k; i++) {
            x[k] -= x[i] * *(p_k + i);
        }
    }
    /* Solve the linear equation Ux = y, where y is the solution
       obtained above of Lx = B and U is an upper triangular matrix. */
    for (k = n-1, p_k = A + n*(n-1); k >= 0; k--, p_k -= n) {
        if (pivot[k] != k) {
            dum = b[k];
            b[k] = b[pivot[k]];
            b[pivot[k]] = dum;
        }
        for (i = k + 1; i < n; i++) {
            x[k] -= x[i] * *(p_k + i);
        }
        if (*(p_k + k) == 0.0)
            return -1;
        x[k] /= *(p_k + k);
    }
    return 0;
}

/* Thin plate spline function */
static double sf_tps_base_func (double r) {
    if (r == 0.0)
        return 0.0;
    else
        return r*r*log (r);
}

/* ================================= TPS interpolation ================================================ */

/* Thin plate spline interpolation, original code by Jarno Elonen */

void sf_tps_build_matrix (float **L, float *x, float *y, int nd,
                          float eps)
/*< Build thin plate spline matrix L[nd+3][nd+3] for coordinate
    locations x[nd] and y[nd]; eps is bending energy >*/
{
    int i, j;
    double a = 0.0, elen;

    if (nd < 3)
        return;

    /* Fill K (p x p, upper left of L) and calculate
       mean edge length from control points.
       K is symmetrical, thus only half of the 
       coefficients has to be calculated. */
    for (i = 0; i < nd; i++) {
        for (j = i + 1; j < nd; j++) {
            elen = sqrt ((x[i] - x[j])*(x[i] - x[j]) + 
                         (y[i] - y[j])*(y[i] - y[j]));
            L[i][j] = sf_tps_base_func (elen);
            L[j][i] = L[i][j];
            a += elen*2;
        }
    }
    a /= (double)(nd*nd);

    /* Fill the rest of L */
    for (i = 0; i < nd; i++) {
        /* Diagonal: reqularization parameters (lambda * a^2) */
        L[i][i] = eps*(a*a);
        /* P (p x 3, upper right) */
        L[i][nd] = 1.0;
        L[i][nd + 1] = x[i];
        L[i][nd + 2] = y[i];
        /* P transposed (3 x p, bottom left) */
        L[nd][i] = 1.0;
        L[nd + 1][i] = x[i];
        L[nd + 2][i] = y[i];
    }
    /* O (3 x 3, lower right) */
    for (i = nd; i < nd + 3; i++) {
        for (j = nd; j < nd + 3; j++) {
            L[i][j] = 0.0;
        }
    }
}

float sf_tps_compute_point (float *w, float *x, float *y, int nd,
                            float xc, float yc)
/*< Compute point at xc,yc using thin-plate spline coefficients w[nd+3]
    and control nodes at x[nd], y[nd] >*/
{
    int i;
    float h;

    /* Do interpolation */
    h = w[nd] + w[nd + 1]*xc + w[nd + 2]*yc;
    for (i = 0; i < nd; i++) {
        h += w[i]*sf_tps_base_func (sqrtf ((x[i] - xc)*(x[i] - xc) + 
                                           (y[i] - yc)*(y[i] - yc)));
    }
    return h;
}


/* ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// Modified by Jiubing Cheng under Madagascar enveriment (2012.12)
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//
// Reference:
// Kopp J., 2008, Efficient numerical diagonalization of hermitian 3 × 3 matrices.
//
// ----------------------------------------------------------------------------*/
#include <rsf.h>

/* Constants */
#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151   /* sqrt(3) */
#endif

/* Macros */
#define SQR(x)      ((x)*(x))                        /* x^2 */ 
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */


void dsyev2(double A, double B, double C, double *rt1, double *rt2,
                   double *cs, double *sn)
/*< ----------------------------------------------------------------------------
// Calculates the eigensystem of a real symmetric 2x2 matrix
//    [ A  B ]
//    [ B  C ]
// in the form
//    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
//    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
// where rt1 >= rt2. Note that this convention is different from the one used
// in the LAPACK routine DLAEV2, where |rt1| >= |rt2|.
// --------------------------------------------------------------------------->*/
{
  double sm = A + C;
  double df = A - C;
  double rt = sqrt(SQR(df) + 4.0*B*B);
  double t;

  if (sm > 0.0)
  {
    *rt1 = 0.5 * (sm + rt);
    t = 1.0/(*rt1);
    *rt2 = (A*t)*C - (B*t)*B;
  }
  else if (sm < 0.0)
  {
    *rt2 = 0.5 * (sm - rt);
    t = 1.0/(*rt2);
    *rt1 = (A*t)*C - (B*t)*B;
  }
  else       /* This case needs to be treated separately to avoid div by 0 */
  {
    *rt1 = 0.5 * rt;
    *rt2 = -0.5 * rt;
  }

  /* Calculate eigenvectors */
  if (df > 0.0)
    *cs = df + rt;
  else
    *cs = df - rt;

  if (fabs(*cs) > 2.0*fabs(B))
  {
    t   = -2.0 * B / *cs;
    *sn = 1.0 / sqrt(1.0 + SQR(t));
    *cs = t * (*sn);
  }
  else if (fabs(B) == 0.0)
  {
    *cs = 1.0;
    *sn = 0.0;
  }
  else
  {
    t   = -0.5 * (*cs) / B;
    *cs = 1.0 / sqrt(1.0 + SQR(t));
    *sn = t * (*cs);
  }

  if (df > 0.0)
  {
    t   = *cs;
    *cs = -(*sn);
    *sn = t;
  }
}

void slvsec3(double d[3], double z[3], double w[3],
                    double R[3][3], int i0, int i1, int i2)
/*< ----------------------------------------------------------------------------
// Finds the three roots w_j of the secular equation
//   f(w_j) = 1 + Sum[ z_i / (d_i - w_j) ]  ==  0.
// It is assumed that d_0 <= d_1 <= d_2, and that all z_i have the same sign.
// The arrays P_i will contain the information required for the calculation
// of the eigenvectors:
//   P_ij = d_i - w_j.
// These differences can be obtained with better accuracy from intermediate
// results.
// ---------------------------------------------------------------------------->*/
{
  int i, j, nIter;
  double a[4];            /* Bounds of the intervals bracketing the roots */
  double delta;           /* Shift of the d_i which ensures better accuracy */
  double dd[3];           /* Shifted coefficients dd_i = d_i - delta */
  double xl, xh;          /* Interval which straddles the current root. f(xl) < 0, f(xh) > 0 */
  double x;               /* Current estimates for the root */
  double x0[3];           /* Analytically calculated roots, used as starting values */
  double F, dF;           /* Function value f(x) and derivative f'(x) */
  double dx, dxold;       /* Current and last stepsizes */
  double error;           /* Numerical error estimate, used for termination condition */
  double t[3];            /* Temporary storage used for evaluating f */
  double alpha, beta, gamma;       /* Coefficients of polynomial f(x) * Product [ d_i - x ] */
  double p, sqrt_p, q, c, s, phi;  /* Intermediate results of analytical calculation */
  
  /* Determine intervals which must contain the roots */
  if (z[0] > 0)
  {
    a[0] = d[i0];
    a[1] = d[i1];
    a[2] = d[i2];
    a[3] = fabs(d[0] + 3.0*z[0]) + fabs(d[1] + 3.0*z[1]) + fabs(d[2] + 3.0*z[2]);
  }
  else
  {    
    a[0] = -fabs(d[0] + 3.0*z[0]) - fabs(d[1] + 3.0*z[1]) - fabs(d[2] + 3.0*z[2]);
    a[1] = d[i0];
    a[2] = d[i1];
    a[3] = d[i2];
  }

  /* Calculate roots of f(x) = 0 analytically (analogous to ZHEEVC3) */
  t[0]  = d[1]*d[2];
  t[1]  = d[0]*d[2];
  t[2]  = d[0]*d[1];
  gamma = t[0]*d[0] + (z[0]*t[0] + z[1]*t[1] + z[2]*t[2]);    /* Coefficients */
  beta  = (z[0]*(d[1]+d[2]) + z[1]*(d[0]+d[2]) + z[2]*(d[0]+d[1]))
           + (t[0] + t[1] + t[2]);
  alpha = (z[0] + z[1] + z[2]) + (d[0] + d[1] + d[2]);
  
  p = SQR(alpha) - 3.0*beta;    /* Transformation that removes the x^2 term */
  q = alpha*(p - (3.0/2.0)*beta) + (27.0/2.0)*gamma;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(beta)*(p - beta) - gamma*(q - 27.0/4.0*gamma));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*fabs(sin(phi));

  x0[0] = x0[1] = x0[2] = (1.0/3.0)*(alpha - c);
  if (c > s)             /* Make sure the roots are in ascending order */
  {
    x0[0] -= s;
    x0[1] += s;
    x0[2] += c;
  }
  else if (c < -s)
  {
    x0[0] += c;
    x0[1] -= s;
    x0[2] += s;
  }
  else
  {
    x0[0] -= s;
    x0[1] += c;
    x0[2] += s;
  }

  /* Refine roots with a combined Bisection/Newton-Raphson method */
  for (i=0; i < 3; i++)
  {
      xl = a[i];               /* Lower bound of bracketing interval */
      xh = a[i+1];             /* Upper bound of bracketing interval */
    dx = dxold = 0.5 * (xh - xl);

    /* Make sure that xl != xh */
    if (dx == 0.0)
    {
      w[i] = xl;
      for (j=0; j < 3; j++)
        R[j][i] = d[j] - xl;
      continue;
    }
    
    /* Shift the root close to zero to achieve better accuracy */
    if (x0[i] >= xh)
    {
      delta = xh;
      x     = -dx;
      for (j=0; j < 3; j++)
      {
        dd[j]   = d[j] - delta;
        R[j][i] = dd[j] - x;
      }
    }
    else if (x0[i] <= xl)
    {
      delta = xl;
      x     = dx;
      for (j=0; j < 3; j++)
      {
        dd[j]   = d[j] - delta;
        R[j][i] = dd[j] - x;
      }
    }
    else
    {
      delta = x0[i];
      x     = 0.0;
      for (j=0; j < 3; j++)
        R[j][i] = dd[j] = d[j] - delta;
    }
    xl -= delta;
    xh -= delta;
   
    /* Make sure that f(xl) < 0 and f(xh) > 0 */
    if (z[0] < 0.0)
    {
      double t = xh;
      xh = xl;
      xl = t;
    }

    /* Main iteration loop */
    for (nIter=0; nIter < 500; nIter++)
    {
	/* Evaluate f and f', and calculate an error estimate */
      F     = 1.0;
      dF    = 0.0;
      error = 1.0;
      for (j=0; j < 3; j++)
      {
        t[0]   = 1.0 / R[j][i];
        t[1]   = z[j] * t[0];
        t[2]   = t[1] * t[0];
        F     += t[1];
        error += fabs(t[1]);
        dF    += t[2];
      }

      /* Check for convergence */ 
      if (fabs(F) <= DBL_EPSILON * (8.0 * error + fabs(x*dF)))
        break;

      /* Adjust interval boundaries */
      if (F < 0.0)
        xl   = x;
      else
        xh   = x;

      /* Check, whether Newton-Raphson would converge fast enough. If so,
	 give it a try. If not, or if it would run out of bounds, use bisection */
      if (fabs(2.0 * F) < fabs(dxold * dF))
      {
        dxold = dx;
        dx    = F / dF;
        x     = x - dx;
        if ((x - xh) * (x - xl) >= 0.0)
        {
          dx = 0.5 * (xh - xl);
          x  = xl + dx;
        }
      }
      else
      {
        dx = 0.5 * (xh - xl);
        x  = xl + dx;
      }

      /* Prepare next iteration */
      for (j=0; j < 3; j++)
        R[j][i] = dd[j] - x;
    }
     
    /* Un-shift result */
    w[i] = x + delta;
  }
}

void dsytrd3(double A[3][3], double Q[3][3], double d[3], double e[2])
/*<  ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// --------------------------------------------------------------------------->*/
{
  int i, j;
#define n 3
  double u[n], q[n];
  double omega, f;
  double K, h, g;
  
  /* Initialize Q to the identitity matrix */
#ifndef EVALS_ONLY
  for (i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  /* Bring first row and column to the desired form */
  h = SQR(A[0][1]) + SQR(A[0][2]);
  if (A[0][1] > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[0][1];
  u[1] = A[0][1] - g;
  u[2] = A[0][2];
  
  omega = h - f;
  if (omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    for (i=1; i < n; i++)
    {
      f    = A[1][i] * u[1] + A[i][2] * u[2];
      q[i] = omega * f;                  /* p */
      K   += u[i] * f;                   /* u* A u */
    }
    K *= 0.5 * SQR(omega);

    for (i=1; i < n; i++)
      q[i] = q[i] - K * u[i];
    
    d[0] = A[0][0];
    d[1] = A[1][1] - 2.0*q[1]*u[1];
    d[2] = A[2][2] - 2.0*q[2]*u[2];
    
    /* Store inverse Householder transformation in Q */
#ifndef EVALS_ONLY
    for (j=1; j < n; j++)
    {
      f = omega * u[j];
      for (i=1; i < n; i++)
        Q[i][j] = Q[i][j] - f*u[i];
    }
#endif

    /* Calculate updated A[1][2] and store it in e[1] */
    e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
  }
  else
  {
    for (i=0; i < n; i++)
      d[i] = A[i][i];
    e[1] = A[1][2];
  }
}

int dsyevc3(double A[3][3], double w[3])
/*< ----------------------------------------------------------------------------
// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ---------------------------------------------------------------------------->*/
{
  double m, c1, c0;
  
  /* Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  | */
  double de = A[0][1] * A[1][2];                                    /* d * e */
  double dd = SQR(A[0][1]);                                         /* d^2 */
  double ee = SQR(A[1][2]);                                         /* e^2 */
  double ff = SQR(A[0][2]);                                         /* f^2 */
  double p, sqrt_p, q, c, s, phi;

  m  = A[0][0] + A[1][1] + A[2][2];
  c1 = (A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2])        /* a*b + a*c + b*c - d^2 - e^2 - f^2 */
          - (dd + ee + ff);
  c0 = A[2][2]*dd + A[0][0]*ee + A[1][1]*ff - A[0][0]*A[1][1]*A[2][2]
      - 2.0 * A[0][2]*de;                                     /* c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e) */

  p = SQR(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
  
  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);

  w[1]  = (1.0/3.0)*(m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;

  return 0;
}

int dsyevv3(double A[3][3], double Q[3][3], double w[3])
/*< ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using Cardano's method for the eigenvalues and an analytical
// method based on vector cross products for the eigenvectors.
// Only the diagonal and upper triangular parts of A need to contain meaningful
// values. However, all of A may be used as temporary storage and may hence be
// destroyed.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   dsyevc3()
// ----------------------------------------------------------------------------
// Version history:
//   v1.1 (12 Mar 2012): Removed access to lower triangualr part of A
//     (according to the documentation, only the upper triangular part needs
//     to be filled)
//   v1.0: First released version
// ---------------------------------------------------------------------------->*/
{
#ifndef EVALS_ONLY
    double norm;          /* Squared norm or inverse norm of current eigenvector */
  double n0, n1;        /* Norm of first and second columns of A */
  double n0tmp, n1tmp;  /* "Templates" for the calculation of n0/n1 - saves a few FLOPS */
  double thresh;        /* Small number used as threshold for floating point comparisons */
  double error;         /* Estimated maximum roundoff error in some steps */
  double wmax;          /* The eigenvalue of maximum modulus */
  double f, t;          /* Intermediate storage */
  int i, j;             /* Loop counters */
#endif

  /* Calculate eigenvalues */
  dsyevc3(A, w);

#ifndef EVALS_ONLY
  wmax = fabs(w[0]);
  if ((t=fabs(w[1])) > wmax)
    wmax = t;
  if ((t=fabs(w[2])) > wmax)
    wmax = t;
  thresh = SQR(8.0 * DBL_EPSILON * wmax);

  /* Prepare calculation of eigenvectors */
  n0tmp   = SQR(A[0][1]) + SQR(A[0][2]);
  n1tmp   = SQR(A[0][1]) + SQR(A[1][2]);
  Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
  Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
  Q[2][1] = SQR(A[0][1]);

  /* Calculate first eigenvector by the formula
  //   v[0] = (A - w[0]).e1 x (A - w[0]).e2 */
  A[0][0] -= w[0];
  A[1][1] -= w[0];
  Q[0][0] = Q[0][1] + A[0][2]*w[0];
  Q[1][0] = Q[1][1] + A[1][2]*w[0];
  Q[2][0] = A[0][0]*A[1][1] - Q[2][1];
  norm    = SQR(Q[0][0]) + SQR(Q[1][0]) + SQR(Q[2][0]);
  n0      = n0tmp + SQR(A[0][0]);
  n1      = n1tmp + SQR(A[1][1]);
  error   = n0 * n1;
  
  if (n0 <= thresh)         /* If the first column is zero, then (1,0,0) is an eigenvector */
  {
    Q[0][0] = 1.0;
    Q[1][0] = 0.0;
    Q[2][0] = 0.0;
  }
  else if (n1 <= thresh)    /* If the second column is zero, then (0,1,0) is an eigenvector */
  {
    Q[0][0] = 0.0;
    Q[1][0] = 1.0;
    Q[2][0] = 0.0;
  }
  else if (norm < SQR(64.0 * DBL_EPSILON) * error)
  {                         /* If angle between A[0] and A[1] is too small, don't use */
      t = SQR(A[0][1]);       /* cross product, but calculate v ~ (1, -A0/A1, 0) */
    f = -A[0][0] / A[0][1];
    if (SQR(A[1][1]) > t)
    {
      t = SQR(A[1][1]);
      f = -A[0][1] / A[1][1];
    }
    if (SQR(A[1][2]) > t)
      f = -A[0][2] / A[1][2];
    norm    = 1.0/sqrt(1 + SQR(f));
    Q[0][0] = norm;
    Q[1][0] = f * norm;
    Q[2][0] = 0.0;
  }
  else                      /* This is the standard branch */
  {
    norm = sqrt(1.0 / norm);
    for (j=0; j < 3; j++)
      Q[j][0] = Q[j][0] * norm;
  }

  
  /* Prepare calculation of second eigenvector */
  t = w[0] - w[1];
  if (fabs(t) > 8.0 * DBL_EPSILON * wmax)
  {
      /* For non-degenerate eigenvalue, calculate second eigenvector by the formula
      //   v[1] = (A - w[1]).e1 x (A - w[1]).e2 */
    A[0][0] += t;
    A[1][1] += t;
    Q[0][1]  = Q[0][1] + A[0][2]*w[1];
    Q[1][1]  = Q[1][1] + A[1][2]*w[1];
    Q[2][1]  = A[0][0]*A[1][1] - Q[2][1];
    norm     = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1]);
    n0       = n0tmp + SQR(A[0][0]);
    n1       = n1tmp + SQR(A[1][1]);
    error    = n0 * n1;
 
    if (n0 <= thresh)       /* If the first column is zero, then (1,0,0) is an eigenvector */
    {
      Q[0][1] = 1.0;
      Q[1][1] = 0.0;
      Q[2][1] = 0.0;
    }
    else if (n1 <= thresh)  /* If the second column is zero, then (0,1,0) is an eigenvector */
    {
      Q[0][1] = 0.0;
      Q[1][1] = 1.0;
      Q[2][1] = 0.0;
    }
    else if (norm < SQR(64.0 * DBL_EPSILON) * error)
    {                       /* If angle between A[0] and A[1] is too small, don't use */
	t = SQR(A[0][1]);     /* cross product, but calculate v ~ (1, -A0/A1, 0) */
      f = -A[0][0] / A[0][1];
      if (SQR(A[1][1]) > t)
      {
        t = SQR(A[1][1]);
        f = -A[0][1] / A[1][1];
      }
      if (SQR(A[1][2]) > t)
        f = -A[0][2] / A[1][2];
      norm    = 1.0/sqrt(1 + SQR(f));
      Q[0][1] = norm;
      Q[1][1] = f * norm;
      Q[2][1] = 0.0;
    }
    else
    {
      norm = sqrt(1.0 / norm);
      for (j=0; j < 3; j++)
        Q[j][1] = Q[j][1] * norm;
    }
  }
  else
  {
      /* For degenerate eigenvalue, calculate second eigenvector according to
    //   v[1] = v[0] x (A - w[1]).e[i]
    //   
    // This would really get to complicated if we could not assume all of A to
    // contain meaningful values. */
    A[1][0]  = A[0][1];
    A[2][0]  = A[0][2];
    A[2][1]  = A[1][2];
    A[0][0] += w[0];
    A[1][1] += w[0];
    for (i=0; i < 3; i++)
    {
      A[i][i] -= w[1];
      n0       = SQR(A[0][i]) + SQR(A[1][i]) + SQR(A[2][i]);
      if (n0 > thresh)
      {
        Q[0][1]  = Q[1][0]*A[2][i] - Q[2][0]*A[1][i];
        Q[1][1]  = Q[2][0]*A[0][i] - Q[0][0]*A[2][i];
        Q[2][1]  = Q[0][0]*A[1][i] - Q[1][0]*A[0][i];
        norm     = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1]);
        if (norm > SQR(256.0 * DBL_EPSILON) * n0) /* Accept cross product only if the angle between */
        {                                         /* the two vectors was not too small */
          norm = sqrt(1.0 / norm);
          for (j=0; j < 3; j++)
            Q[j][1] = Q[j][1] * norm;
          break;
        }
      }
    }
    
    if (i == 3)    /* This means that any vector orthogonal to v[0] is an EV. */
    {
      for (j=0; j < 3; j++)
	  if (Q[j][0] != 0.0)                                   /* Find nonzero element of v[0] ... */
	  {                                                     /* ... and swap it with the next one */
          norm          = 1.0 / sqrt(SQR(Q[j][0]) + SQR(Q[(j+1)%3][0]));
          Q[j][1]       = Q[(j+1)%3][0] * norm;
          Q[(j+1)%3][1] = -Q[j][0] * norm;
          Q[(j+2)%3][1] = 0.0;
          break;
        }
    }
  }
      
  
  /* Calculate third eigenvector according to */
  /*   v[2] = v[0] x v[1] */
  Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1];
  Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1];
  Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1];
#endif

  return 0;
}

int dsyevd3(double A[3][3], double Q[3][3], double w[3])
/*< ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using Cuppen's Divide & Conquer algorithm.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   dsyev2(), slvsec3(), dsytrd3()
// ---------------------------------------------------------------------------->*/
{
  int i, j, k;
#define n 3
  double R[3][3];                /* Householder transformation matrix */
  double P[3][3];                /* Unitary transformation matrix which diagonalizes D + w w^T */
  double e[2];                   /* Off-diagonal elements after Householder transformation */
  double d[3];                   /* Eigenvalues of split matrix in the "divide" step) */
  double c, s;                   /* Eigenvector of 2x2 block in the "divide" step */
  double z[3];                   /* Numerators of secular equation / Updating vector */
  double t;                      /* Miscellaenous temporary stuff */

  /* Initialize Q */
#ifndef EVALS_ONLY
  memset(Q, 0.0, 9*sizeof(double));
#endif
  
  /* Transform A to real tridiagonal form by the Householder method */
  dsytrd3(A, R, w, e);
 
  
  /* "Divide"
  // --------
  
  // Detect matrices that factorize to avoid multiple eigenvalues in the Divide/Conquer algorithm */
  for (i=0; i < n-1; i++)
  {
    t = fabs(w[i]) + fabs(w[i+1]);
    if (fabs(e[i]) <= 8.0*DBL_EPSILON*t)
    {
      if (i == 0)
      {
        dsyev2(w[1], e[1], w[2], &d[1], &d[2], &c, &s);
        w[1] = d[1];
        w[2] = d[2];
#ifndef EVALS_ONLY
        Q[0][0] = 1.0;
        for (j=1; j < n; j++)
        {
          Q[j][1] = s*R[j][2] + c*R[j][1];
          Q[j][2] = c*R[j][2] - s*R[j][1];
        }
#endif
      }
      else
      {
        dsyev2(w[0], e[0], w[1], &d[0], &d[1], &c, &s);
        w[0] = d[0];
        w[1] = d[1];
#ifndef EVALS_ONLY
        Q[0][0]   = c;
        Q[0][1]   = -s;
        Q[1][0]   = R[1][1]*s;
        Q[1][1]   = R[1][1]*c;
        Q[1][2]   = R[1][2];
        Q[2][0]   = R[2][1]*s;
        Q[2][1]   = R[2][1]*c;
        Q[2][2]   = R[2][2];
#endif
      }

      return 0;
    }
  }
  
  /* Calculate eigenvalues and eigenvectors of 2x2 block */
  dsyev2(w[1]-e[0], e[1], w[2], &d[1], &d[2], &c, &s);
  d[0] = w[0] - e[0];

  
  /* "Conquer"
  // ---------

  // Determine coefficients of secular equation */
  z[0] = e[0];
  z[1] = e[0] * SQR(c);
  z[2] = e[0] * SQR(s);

  /* Call slvsec3 with d sorted in ascending order. We make
  // use of the fact that dsyev2 guarantees d[1] >= d[2]. */
  if (d[0] < d[2])
    slvsec3(d, z, w, P, 0, 2, 1);
  else if (d[0] < d[1])
    slvsec3(d, z, w, P, 2, 0, 1);
  else
    slvsec3(d, z, w, P, 2, 1, 0);

#ifndef EVALS_ONLY
  /* Calculate eigenvectors of matrix D + beta * z * z^t and store them in the
  // columns of P */
  z[0] = sqrt(fabs(e[0]));
  z[1] = c * z[0];
  z[2] = -s * z[0];

  /* Detect duplicate elements in d to avoid division by zero */
  t = 8.0*DBL_EPSILON*(fabs(d[0]) + fabs(d[1]) + fabs(d[2]));
  if (fabs(d[1] - d[0]) <= t)
  {
    for (j=0; j < n; j++)
    {
      if (P[0][j] * P[1][j] <= 0.0)
      {
        P[0][j] = z[1];
        P[1][j] = -z[0];
        P[2][j] = 0.0;
      }
      else
        for (i=0; i < n; i++)
          P[i][j] = z[i]/P[i][j];
    }
  }
  else if (fabs(d[2] - d[0]) <= t)
  {
    for (j=0; j < n; j++)
    {
      if (P[0][j] * P[2][j] <= 0.0)
      {
        P[0][j] = z[2];
        P[1][j] = 0.0;
        P[2][j] = -z[0];
      }
      else
        for (i=0; i < n; i++)
          P[i][j] = z[i]/P[i][j];
    }
  }
  else
  {
    for (j=0; j < n; j++)
      for (i=0; i < n; i++)
      {
        if (P[i][j] == 0.0)
        {
          P[i][j]       = 1.0;
          P[(i+1)%n][j] = 0.0;
          P[(i+2)%n][j] = 0.0;
          break;
        }
        else
          P[i][j] = z[i]/P[i][j];
      }
  }

  /* Normalize eigenvectors of D + beta * z * z^t */
  for (j=0; j < n; j++)
  {
    t = SQR(P[0][j]) + SQR(P[1][j]) + SQR(P[2][j]);
    t = 1.0 / sqrt(t);
    for (i=0; i < n; i++)
      P[i][j] *= t;
  }
  
  /* Undo diagonalization of 2x2 block */
  for (j=0; j < n; j++)
  {
    t       = P[1][j];
    P[1][j] = c*t - s*P[2][j];
    P[2][j] = s*t + c*P[2][j];
  }

  /* Undo Householder transformation */
  for (j=0; j < n; j++)
    for (k=0; k < n; k++)
    {
      t = P[k][j];
      for (i=0; i < n; i++)
        Q[i][j] += t * R[i][k];
    }
#endif

  return 0;
}

int dsyevq3(double A[3][3], double Q[3][3], double w[3])
/*< ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ---------------------------------------------------------------------------->*/
{
  int i, k, l;
#define n 3
  double e[3];                   /* The third element is used only as temporary workspace */
  double g, r, p, f, b, s, c, t; /* Intermediate storage */
  int nIter;
  int m;

  /* Transform A to real tridiagonal form by the Householder method */
  dsytrd3(A, Q, w, e);
  
  /* Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements */
  for (l=0; l < n-1; l++)
  {
    nIter = 0;
    while (1)
    {
	/* Check for convergence and exit iteration loop if off-diagonal
	// element e(l) is zero */
      for (m=l; m <= n-2; m++)
      {
        g = fabs(w[m])+fabs(w[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      /* Calculate g = d_m - k */
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        /* Form eigenvectors */
#ifndef EVALS_ONLY
        for (k=0; k < n; k++)
        {
          t = Q[k][i+1];
          Q[k][i+1] = s*Q[k][i] + c*t;
          Q[k][i]   = c*Q[k][i] - s*t;
        }
#endif 
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }

  return 0;
}


int dsyevh3(double A[3][3], double Q[3][3], double w[3])
/*< ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using Cardano's method for the eigenvalues and an analytical
// method based on vector cross products for the eigenvectors. However,
// if conditions are such that a large error in the results is to be
// expected, the routine falls back to using the slower, but more
// accurate QL algorithm. Only the diagonal and upper triangular parts of A need
// to contain meaningful values. Access to A is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   dsyevc3(), dsytrd3(), dsyevq3()
// ----------------------------------------------------------------------------
// Version history:
//   v1.1: Simplified fallback condition --> speed-up
//   v1.0: First released version
// ---------------------------------------------------------------------------->*/
{
#ifndef EVALS_ONLY
    double norm;          /* Squared norm or inverse norm of current eigenvector */
/*  double n0, n1;        // Norm of first and second columns of A */
    double error;         /* Estimated maximum roundoff error */
    double t, u;          /* Intermediate storage */
    int j;                /* Loop counter */
#endif

    /* Calculate eigenvalues */
  dsyevc3(A, w);

#ifndef EVALS_ONLY
/*  n0 = SQR(A[0][0]) + SQR(A[0][1]) + SQR(A[0][2]); */
/*  n1 = SQR(A[0][1]) + SQR(A[1][1]) + SQR(A[1][2]); */
  
  t = fabs(w[0]);
  if ((u=fabs(w[1])) > t)
    t = u;
  if ((u=fabs(w[2])) > t)
    t = u;
  if (t < 1.0)
    u = t;
  else
    u = SQR(t);
  error = 256.0 * DBL_EPSILON * SQR(u);

/*  error = 256.0 * DBL_EPSILON * (n0 + u) * (n1 + u); */

  Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
  Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
  Q[2][1] = SQR(A[0][1]);

  /* Calculate first eigenvector by the formula
  //   v[0] = (A - w[0]).e1 x (A - w[0]).e2 */
  Q[0][0] = Q[0][1] + A[0][2]*w[0];
  Q[1][0] = Q[1][1] + A[1][2]*w[0];
  Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1];
  norm    = SQR(Q[0][0]) + SQR(Q[1][0]) + SQR(Q[2][0]);

  /* If vectors are nearly linearly dependent, or if there might have
  // been large cancellations in the calculation of A[i][i] - w[0], fall
  // back to QL algorithm
  // Note that this simultaneously ensures that multiple eigenvalues do
  // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
  // i.e. all columns of A - w[0] * I are linearly dependent. */
  if (norm <= error)
    return dsyevq3(A, Q, w);
  else                      /* This is the standard branch */
  {
    norm = sqrt(1.0 / norm);
    for (j=0; j < 3; j++)
      Q[j][0] = Q[j][0] * norm;
  }
  
  /* Calculate second eigenvector by the formula
  //   v[1] = (A - w[1]).e1 x (A - w[1]).e2 */
  Q[0][1]  = Q[0][1] + A[0][2]*w[1];
  Q[1][1]  = Q[1][1] + A[1][2]*w[1];
  Q[2][1]  = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1];
  norm     = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1]);
  if (norm <= error)
    return dsyevq3(A, Q, w);
  else
  {
    norm = sqrt(1.0 / norm);
    for (j=0; j < 3; j++)
      Q[j][1] = Q[j][1] * norm;
  }
  
  /* Calculate third eigenvector according to
  //   v[2] = v[0] x v[1] */
  Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1];
  Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1];
  Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1];
#endif

  return 0;
}


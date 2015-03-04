/* generalized minimum residual method for general complex operators
 * Copyright (C) 2014 University of Texas at Austin
 * Inspired by:
 * Copyright (C) 1998-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.c,v 2.10 2006/10/10 19:53:27 ichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
/* Modified for inclusion with Madagascar. */
#include <rsf.h>

#include "cgmres.h"

static int m, n;
static float * c;
static sf_complex * v, * h, * g, * s; /* c can be complex as well, but not necessarily */
/* The triangular matrix R is stored in place of H */
/* The matrices are stored by row, except for vi (Q), which is store by column.*/

/* scales a complex vector
   m : dimension of vector x [] (m) */
static void
my_cscal (int m, float a, sf_complex *x)
{
  int i;

  for (i = 0; i < m; i ++)
    {
      //#ifdef SF_HAS_COMPLEX_H
      x[i] = x[i] * a;
      //#else
      //      x[i] = sf_crmul(x[i], a);
      //#end
    }

}

/* computes y = ax + y
   m : dimension of vector x [] (m) */
static void
my_caxpy (int m, float a, const sf_complex *x, sf_complex *y)
{
  int i;

  for (i = 0; i < m; i ++)
    {
      //#ifdef SF_HAS_COMPLEX_H
      y[i] += x[i] * a;
      //#else
      //      y[i] = sf_cadd(y[i], sf_crmul(x[i], a));
      //#end
    }

}

/* computes a*b (conjugating the first vector)
   m : dimension of vector x [] (m) */
static sf_complex
my_cdotc (int m,
	  const sf_complex *a, const sf_complex *b)
{
  int i;
  sf_complex res;
  
  res = sf_cmplx(0.0f, 0.0f);
  for (i = 0; i < m; i ++) {
    res += conjf(a[i]) * b[i];
  }
  return res;
}

/* computes ab (without conjugating the first vector)
   m : dimension of vector x [] (m) */
static sf_complex
my_cdotu (int m,
	  const sf_complex *a, const sf_complex *b)
{
  int i;
  sf_complex res;
  
  res = sf_cmplx(0.0f, 0.0f);
  for (i = 0; i < m; i ++) {
    res += a[i] * b[i];
  }
  return res;
}

/* double-precision L2 norm of a complex number */
static float
my_cnorm (int n, const sf_complex* x)
{
    double prod, xi, yi;
    int i;

    prod = 0.;
    for (i = 0; i < n; i++) {
	xi = (double) crealf(x[i]);
	yi = (double) cimagf(x[i]);
	prod += xi*xi + yi*yi;
    }
    prod = sqrt(prod);
    return (float) prod;
}

/* m  : number of iteration */
/* nn : dimension of matrix r [] (nnxnn) */
static void
back_sub (int m, int nn,
	  const sf_complex *r, const sf_complex *g,
	  sf_complex *y)
{
  int i, j, jj;

  /*for (j = m - 1;j >= 0;j --){*/
  /* above for-loop fail, because j is unsigned!! */
  for (j = m - 1, jj = 0; jj < m; j --, jj ++)
    {
      y [j] = g [j];
      for (i = j + 1; i < m; i ++)
	{
	  //#ifdef SF_HAS_COMPLEX_H
	  y [j] -= r [j * nn + i] * y [i];
	  //#else
	  //	  y [j] = sf_csub(y [j] , sf_cmul(r [j * nn + i] , y [i]));
	  //#endif
	}
      //#ifdef SF_HAS_COMPLEX_H
      y [j] /= r [j * nn + j]; /* what about division by zero? */
      //#else
      //      y [j] = sf_cdiv(y [j] , r [j * nn + j]);
      //#endif
    }
}

static void
showvec (int m, const sf_complex* v)
{
  int i;
  for (i=0; i < n; i ++) {
    printf ("(%5.2f, %5.2f)",crealf(v[i]),cimagf(v[i]));
  }
  printf("\n");
}

void cgmres_init(int nx      /* data size */, 
		int restart /* memory */) 
/*< initialize >*/
{
    m=restart;
    n=nx;

    sf_warning("Using complex GMRES(m) method...");

    v   = sf_complexalloc((m + 1) * n); /* orthogonal basis q */
    h   = sf_complexalloc(m * m);
    g   = sf_complexalloc(m + 1);
    c   = sf_floatalloc(m);
    s   = sf_complexalloc(m);
}

void cgmres_close(void)
/*< free allocated storage >*/
{   
  free (v);
  free (h);
  free (g);
  free (c);
  free (s);
}

void cgmres (const sf_complex *f                                              /* data */, 
	    sf_complex *x                                                    /* estimated model */,
	    void (*myatimes) (int, const sf_complex *, sf_complex *, void *) /* operator */,
	    void * user_data                                                 /* internal data */,
	    int itmax                                                        /* number of iterations */, 
	    float tol                                                        /* tolerance */, 
	    bool verb                                                        /* verbosity */)
/*< GMRES solver >*/
{
  int iter;
  double res = 0.0;

  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  /* m: # of iteration at once */
  int i, j, k;
  sf_complex rr, hv, r1, r2, g0, mu, tau;
  float tmpf, hh, rho;

  iter = 0;
  /* 1. start: */
  /* compute r0 */
  /* compute v1 */
  /* beta */
  myatimes (n, x, v + 0, user_data); /* use v [0] temporaliry */

  /* v = f - v */
  my_cscal (n, -1.0f, v); /* v = - v */
  my_caxpy (n, 1.0, f, v); /* v = f - v */

  /* g[0] = cblas_dnrm2 (n, v, 1); */
  tmpf = my_cnorm (n, v); /* sqrt(v'*v) */
  g[0] = sf_cmplx( tmpf, 0.0f ); /* beta = norm(r_0) */
  my_cscal (n, 1.0f / tmpf, v); /* v1 = r/beta */

  sf_warning("Original residual %e\n", tmpf*tmpf);

  /* main loop */
  while (iter <= itmax)
    {
      ++iter;
      /* 2. iterate: */
      for (j = 0; j < m; j ++) /* column */
	{
	  /* tmp = A.vj : use v [(j + 1) * n] directly */
	  myatimes (n, v + j * n, v + (j + 1) * n, user_data);

	  /* h_i,j (i=1,...,j) */
	  for (i = 0; i <= j; i ++)
	    {

	      h [i * m + j] =
		my_cdotc (n, v + i * n, v + (j + 1) * n);

	    }
	  /* vv_j+1 */
	  for (k = 0; k < n; k ++)
	    {
	      hv = sf_cmplx(0.0f, 0.0f);
	      for (i = 0; i <= j; i ++)
		{
		  hv += h [i * m + j] * v [i * n + k];
		}
	      v [(j + 1) * n + k] -= hv;
	    }

	  /* h_j+1,j */
	  /* v_j+1 */

	  /* hh = cblas_dnrm2 (n, v + (j + 1) * n, 1); */
	  hh = my_cnorm(n, v + (j + 1) * n); /* h_j+1,j */
	  my_cscal (n, 1.0 / hh, v + (j + 1) * n); /* v_j+1 */

	  /*--- Arnoldi Iteration ---*/

	  /* rotate */
	  /* apply previous givens rotation to new column h_:,j */
	  for (i = 0; i < j; i ++)
	    {
	      r1 = h [ i      * m + j];
	      r2 = h [(i + 1) * m + j];
	      h [ i      * m + j] =   c [i] * r1 + conjf(s [i]) * r2; /* Given's rotation */
	      h [(i + 1) * m + j] = - s [i] * r1 + c [i] * r2; /* Given's rotation */
	    }
	  rr = h [j * m + j]; /*h_j,j (complex)*/
	  tmpf = cabsf(rr);
	  rho = sqrtf (tmpf * tmpf + hh * hh); /* temporary variable */

	  if (1) { /*Fancier Givens rotation*/

	    if (tmpf < hh) {
	      mu = rr/hh;
	      tau = conjf(mu)/cabsf(mu);
	    } else {
	      mu = hh/rr;
	      tau = mu/cabsf(mu);
	    }
	  
	    if (tmpf < 1e-10) {
	      c[j] = 0.0f;
	      s[j] = sf_cmplx(1.0f,0.0f); // 90 degree rotation
	    } else {
	      /* c [j] = rr / rho, s [j] = hh / rho, h [j * m + j] = rho */
	      c [j] = tmpf / rho; /* new givens rotation matrix that zeros h_j+1,j */
	      s [j] = hh*tau / rho;
	    }

	  } else { /*Givens rotation*/

	    c [j] = tmpf/rho;
	    s [j] = c [j] * (hh/rr);

	  }

	  h [j * m + j] = c[j]*rr + conjf(s[j]) * hh; /* resultant (after rotated) element */

	  /* g is givens rotation iteratively applied to ||r_0||e_1 */
	  g0 = g [j];
	  g [j    ] =   c [j] * g0;  /* Given's rotation applied to b vector, and the other term is always zero */
	  g [j + 1] = - s [j] * g0;  /* Given's rotation applied to b vector, and the other term is always zero */
	}
      /* 3. form the approximate solution */
      /* solve for y_k, j=m, h is Hessenberg rotated to upper triangular matrix, and g is 
         givens rotations applied to the RHS (||r0||e_1), and s is the solution (y) */
      back_sub (j/*m*/, m, h, g, s); /* use s [] as y_k */
      /* x_m */
      /* x = Qy = Qs */
      for (i = 0; i < n; i ++)
	{
	  for (k = 0; k < j/*m*/; k ++)
	    {
	      x [i] += v [k * n + i] * s [k];
	    }
	}

      /* 4. restart */
      res = cabsf (g [j/*m*/]); /* residual */
      /*fprintf (stderr, "# iter %d res %e\n", iter, *res);*/
      /* if satisfied, */
      if (verb)
	  sf_warning("complex-gmres(%d) %d %d %e\n",
		     m, iter, j, res*res);
 
      if (res <= tol) break;
      /* else */
      /* compute r_m */
      myatimes (n, x, v + 0, user_data);
      /* r_m */

      /* compute v1 */
      /* v = f - v */
      my_cscal (n, -1.0f, v); /* v = - v */
      my_caxpy (n, 1.0, f, v); /* v = f - v */

      /* g[0] = cblas_dnrm2 (n, v, 1); */
      tmpf = my_cnorm (n, v); /* sqrt(v'*v) */
      g[0] = sf_cmplx( tmpf, 0.0f ); /* beta = norm(r_0) */
      my_cscal (n, 1.0f / tmpf, v); /* v1 = r/beta */

    }

    /* adjust iter */
  iter *= m;

  sf_warning("complex-gmres(%d) it= %d res^2= %e\n", m, iter, res*res);
}


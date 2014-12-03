/* generalized minimum residual method
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

#include "gmres.h"

#include "_bool.h"
/*^*/

#include "error.h"
#include "alloc.h"
#include "blas.h"

static int m, n;
static float * v, * h, * g, * c, * s;

/* m  : number of iteration */
/* nn : dimension of matrix r [] (nnxnn) */
static void
back_sub (int m, int nn,
	  const float *r, const float *g,
	  float *y)
{
  int i, j, jj;

  /*for (j = m - 1;j >= 0;j --){*/
  /* above for-loop fail, because j is unsigned!! */
  for (j = m - 1, jj = 0; jj < m; j --, jj ++)
    {
      y [j] = g [j];
      for (i = j + 1; i < m; i ++)
	{
	  y [j] -= r [j * nn + i] * y [i];
	}
      y [j] /= r [j * nn + j];
    }
}

void sf_gmres_init(int nx      /* data size */, 
		   int restart /* memory */) 
/*< initialize >*/
{
    m=restart;
    n=nx;
    
    v   = sf_floatalloc((m + 1) * n);
    h   = sf_floatalloc(m * m);
    g   = sf_floatalloc(m + 1);
    c   = sf_floatalloc(m);
    s   = sf_floatalloc(m);
}

void sf_gmres_close(void)
/*< free allocated storage >*/
{   
  free (v);
  free (h);
  free (g);
  free (c);
  free (s);
}

void sf_gmres (const float *f                                         /* data */, 
	       float *x                                               /* estimated model */,
	       void (*myatimes) (int, const float *, float *, void *) /* operator */,
	       void * user_data                                       /* internal data */,
	       int itmax                                              /* number of iterations */, 
	       float tol                                              /* tolerance */, 
	       bool verb                                              /* verbosity */)
/*< GMRES solver >*/
{
  int iter;
  double res = 0.0;

  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  /* m: # of iteration at once */
  int i, j, k;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;

  iter = 0;
  /* 1. start: */
  /* compute r0 */
  /* compute v1 */
  /* beta */
  myatimes (n, x, v + 0, user_data); /* use v [0] temporaliry */

  /* use ATLAS' CBLAS routines */

  /* v = f - v */
  cblas_sscal (n, -1.0, v, 1); /* v = - v */
  cblas_saxpy (n, 1.0, f, 1, v, 1); /* v = f - v */

  /* g[0] = cblas_dnrm2 (n, v, 1); */
  g[0] = sqrt (cblas_sdot (n, v, 1, v, 1));
  cblas_sscal (n, 1.0 / g[0], v, 1);

  /* main loop */
  while (iter <= itmax)
    {
      ++iter;
      /* 2. iterate: */
      for (j = 0; j < m; j ++)
	{
	  /* tmp = A.vj : use v [(j + 1) * n] directly */
	  myatimes (n, v + j * n, v + (j + 1) * n, user_data);
	  /* h_i,j (i=1,...,j) */
	  for (i = 0; i <= j; i ++)
	    {
	      /* use ATLAS' CBLAS routines */

	      h [i * m + j] =
		cblas_sdot (n, v + (j + 1) * n, 1,
			    v + i * n, 1);
	    }
	  /* vv_j+1 */
	  for (k = 0; k < n; k ++)
	    {
	      hv = 0.0;
	      for (i = 0; i <= j; i ++)
		{
		  hv += h [i * m + j] * v [i * n + k];
		}
	      v [(j + 1) * n + k] -= hv;
	    }
	  /* h_j+1,j */
	  /* v_j+1 */

	  /* use ATLAS' CBLAS routines */

	  /* hh = cblas_dnrm2 (n, v + (j + 1) * n, 1); */
	  hh = sqrt (cblas_sdot (n, v + (j + 1) * n, 1, v + (j + 1) * n, 1));
	  cblas_sscal (n, 1.0 / hh, v + (j + 1) * n, 1);

	  /* rotate */
	  for (i = 0; i < j; i ++)
	    {
	      r1 = h [ i      * m + j];
	      r2 = h [(i + 1) * m + j];
	      h [ i      * m + j] = c [i] * r1 - s [i] * r2;
	      h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
	    }
	  rr = h [j * m + j];
	  hv = sqrt (rr * rr + hh * hh); /* temporary variable */
	  c [j] =  rr / hv;
	  s [j] = -hh / hv;
	  h [j * m + j] = hv; /* resultant (after rotated) element */

	  g0 = g [j];
	  g [j    ] = c [j] * g0;
	  g [j + 1] = s [j] * g0;
	}
      /* 3. form the approximate solution */
      /* solve y_k */
      back_sub (j/*m*/, m, h, g, c); /* use c [] as y_k */
      /* x_m */
      for (i = 0; i < n; i ++)
	{
	  for (k = 0; k < j/*m*/; k ++)
	    {
	      x [i] += v [k * n + i] * c [k];
	    }
	}

      /* 4. restart */
      res = fabs (g [j/*m*/]); /* residual */
      /*fprintf (stderr, "# iter %d res %e\n", iter, *res);*/
      /* if satisfied, */
      if (verb)
	  sf_warning("libiter-gmres(%d) %d %d %e\n",
		     m, iter, j, res*res);
 
      if (res <= tol) break;
      /* else */
      /* compute r_m */
      myatimes (n, x, v + 0, user_data);
      /* r_m */
      /* compute v1 */

      /* use ATLAS' CBLAS routines */

      /* v = f - v */
      cblas_sscal (n, -1.0, v, 1); /* v = - v */
      cblas_saxpy (n, 1.0, f, 1, v, 1); /* v = f - v */

      /* g [0] = cblas_dnrm2 (n, v, 1); */
      g[0] = sqrt (cblas_sdot (n, v, 1, v, 1));
      cblas_sscal (n, 1.0 / g[0], v, 1);

    }

    /* adjust iter */
  iter *= m;

  sf_warning("libiter-gmres(%d) it= %d res^2= %e\n", m, iter, res*res);
}


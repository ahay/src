/* 2-D missing data interpolation by gradient or Laplacian regularization */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
/*^*/

#include "lapfill.h"

#include "laplac2.h"

static int n12;
static float* zero;
static bool grad, verb;

void lapfill_init (int m1, int m2 /* model size */, 
		   bool grad1     /* gradient (or Laplacian) */,
                   bool verb1     /* verbosity */)
/*< initialize >*/
{
    int i;
    int n;

    n12 = m1*m2;
    grad = grad1;
    verb = verb1;

    n = grad? 2*n12: n12;

    zero = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	zero[i] = 0.;
    }

    if (grad) {
	sf_igrad2_init (m1,m2);
    } else {
	laplac2_init (m1,m2);
    }
}

void lapfill_close (void)
/*< free allocated storage >*/
{
    free (zero);
}

void lapfill(int niter   /* number of iterations */, 
	     float* mm   /* model [m1*m2] */, 
	     bool *known /* mask for known data [m1*m2] */)
/*< interpolate >*/
{
    if (grad) {
	sf_solver (sf_igrad2_lop, sf_cgstep, n12, 2*n12, mm, zero, 
		   niter, "x0", mm, "known", known, 
		   "verb", verb, "end");
    } else {
	sf_solver (laplac2_lop, sf_cgstep, n12, n12, mm, zero, 
		   niter, "x0", mm, "known", known, 
		   "verb", verb, "end");
    }
    sf_cgstep_close ();
}

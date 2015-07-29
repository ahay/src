/* 2-D Missing data interpolation by gradient minimization. */
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

#include "grad2fill.h"

#include "_bool.h"
/*^*/

#include "alloc.h"
#include "igrad2.h"
#include "bigsolver.h"
#include "cgstep.h"

static int n12;
static float* zero;

void sf_grad2fill_init (int m1, int m2)
/*< Initialize with data dimensions >*/
{
    int i;

    n12 = m1*m2;

    zero = sf_floatalloc (2*n12);

    for (i=0; i < 2*n12; i++) {
	zero[i] = 0.;
    }

    sf_igrad2_init (m1,m2);
}

void sf_grad2fill_close (void)
/*< Free allocate storage >*/
{
    free (zero);
}

void sf_grad2fill(int niter   /* number of iterations */, 
	       float* mm   /* estimated model */, 
	       bool *known /* mask */)
/*< Run optimization >*/
{
    sf_solver (sf_igrad2_lop, sf_cgstep, n12, 2*n12, mm, zero, 
	       niter, "x0", mm, "known", known, "end");
    sf_cgstep_close ();
}

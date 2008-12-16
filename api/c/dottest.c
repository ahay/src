/* Dot product test for linear operators. */
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

#include <stdlib.h>

#include "dottest.h"
#include "alloc.h"
#include "randn.h"
#include "blas.h"

#include "_solver.h"
/*^*/

void sf_dot_test(sf_operator oper /* linear operator */, 
		 int nm           /* model size */, 
		 int nd           /* data size */, 
		 float* dot1      /* first output */, 
		 float* dot2      /* second output */) 
/*< The dot product test to see if the adjoint is coded correctly.
   In the output dot1[0] shpould be equal to dot1[1] 
   (within machine precision),
   and dot2[0] should be equal to dot2[1]. >*/
{
    float *mod1, *mod2, *dat1, *dat2;

    mod1 = sf_floatalloc (nm);
    mod2 = sf_floatalloc (nm);
    dat1 = sf_floatalloc (nd);
    dat2 = sf_floatalloc (nd);

    sf_random( nm, mod1);
    sf_random( nd, dat2);

    oper(false, false, nm, nd, mod1, dat1);
    dot1[0] = cblas_sdot( nd, dat1, 1, dat2, 1);

    oper(true, false, nm, nd, mod2, dat2);
    dot1[1] = cblas_sdot( nm, mod1, 1, mod2, 1);

    oper(false, true, nm, nd, mod1, dat1);
    dot2[0] = cblas_sdot( nd, dat1, 1, dat2, 1);

    oper(true, true, nm, nd, mod2, dat2);
    dot2[1] = cblas_sdot( nm, mod1, 1, mod2, 1);

    free (mod1);
    free (mod2);
    free (dat1);
    free (dat2);
}


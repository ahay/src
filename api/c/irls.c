/* Weighting for iteratively-reweighted least squares */
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

#include <math.h>

#include "irls.h"
#include "quantile.h"
#include "alloc.h"

static float *abs1;

void sf_irls_init(int n) 
/*< Initialize with data size >*/
{
    abs1 = sf_floatalloc(n);
}

void sf_irls_close(void) 
/*< free allocated storage >*/
{
    free (abs1);
}

void sf_l1 (int n, const float *res, float *weight)  
/*< weighting for L1 norm >*/
{
    float rbar;
    int i;

   /* take absolute value */
    for (i=0; i < n; i++) {
	abs1[i] = fabsf(res[i]);
    }

    /* find median (rbar) */
    rbar = sf_quantile(n/2,n,abs1);

    /* weight = 1/sqrt(sqrt(1+(res/rbar)^2)) */
    for (i=0; i < n; i++) {
	weight[i] = 1./sqrtf(hypotf(1.,res[i]/rbar));
    }
}

void sf_cauchy (int n, const float *res, float *weight)  
/*< weighting for Cauchy norm >*/
{
    float rbar;
    int i;

    /* take absolute value */
    for (i=0; i < n; i++) {
	abs1[i] = fabsf(res[i]);
    }

    /* find median (rbar) */
    rbar = sf_quantile(n/2,n,abs1);

    /* weight = 1/sqrt(1+(res/rbar)^2) */
    for (i=0; i < n; i++) {
	weight[i] = 1./hypotf(1.,res[i]/rbar);
    }
}

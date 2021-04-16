/* Blocky derivative by shaping regularization */
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

#include "blockder.h"
#include "doubint.h"

static float *tmp;

void blockder_init(int n1, int n2       /* data size */, 
                   float perc           /* sharpening percentage */,
                   const float *block   /* boundaries */,
                   float *weight)
/*< initialize >*/
{
    int n;

    n = n1*n2;
    
    sf_repeat_init(n1,n2,doubint_lop);

    doubint_init(n1);
    tmp = sf_floatalloc(n);
    sf_sharpen_init(n,perc,0.5);
    sf_sharpen(block);
    sf_weight2_init(1,n,weight);

    sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, false);    
}

void blockder_close(void)
/*< free allocated storage >*/
{
    free(tmp);

    doubint_close();
    sf_sharpen_close();
    sf_weight2_close();
    sf_conjgrad_close();
}

void blockder(int niter   /* number of iterations */, 
	      float* data /* input data */, 
	      float* der  /* output derivative */) 
/*< find the derivative >*/
{ 
    sf_conjgrad(sf_weight2_lop,sf_repeat_lop,sf_weight_lop,tmp,der,data,niter);
}

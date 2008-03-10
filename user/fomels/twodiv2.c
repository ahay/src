/* Smooth division with several components */
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

#include "twodiv2.h"
#include "gauss2.h"
#include "freqfilt2.h"
#include "triangle2.h"
#include "repeat.h"
#include "weight2.h"

static int n, niter;
static float *p;
static bool gauss;

void twodiv2_init(int nw             /* number of components */, 
		  int n1, int n2     /* data size */, 
		  float f1, float f2 /* smoothing radius */, 
		  int niter1         /* number of iteration */, 
		  bool gauss1        /* Gaussian or triangle smoothing */,
		  bool verb          /* verbosity flag */,
		  float* den         /* denominator */) 
/*< initialize >*/
{
    n = n1*n2;
    niter = niter1;
    gauss = gauss1;

    if (gauss) {
	gauss2_init(n1,n2,f1,f2);
	repeat_init(n,nw,freqfilt2_lop);
    } else {
	triangle2_init((int) f1, (int) f2, n1, n2, 1);
	repeat_init(n,nw,triangle2_lop);
    }
    sf_conjgrad_init(nw*n, nw*n, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (nw*n);
    weight2_init(nw,n,den);
}

void twodiv2_close (void)
/*< free allocated storage >*/
{
    if (gauss) {
	gauss2_close();
    } else { 
	triangle2_close();
    }
    sf_conjgrad_close();
    free (p);
    weight2_close();
}

void twodiv2 (float* num, float* rat)
/*< smoothly divide num/rat >*/
{
    sf_conjgrad(NULL,weight2_lop,repeat_lop,p,rat,num,niter);
}

/* 	$Id$	 */

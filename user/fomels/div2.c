/* 2-D smooth division */
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

#include "div2.h"
#include "gauss2.h"

static int n, niter;
static float *p;
static bool gauss;

void div2_init(int n1, int n2     /* data dimensions */, 
	       float f1, float f2 /* smoothing */, 
	       int niter1         /* number of iterations */, 
	       bool gauss1        /* if exact gaussian */,
	       bool verb          /* verbosity flag */)
/*< initialize >*/ 
{
    n = n1*n2;
    niter = niter1;
    gauss = gauss1;

    if (gauss) {
	gauss2_init(n1,n2,f1,f2);
    } else {
	sf_triangle2_init((int) f1, (int) f2, n1, n2, 1);
    }
    sf_conjgrad_init(n, n, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (n);
}

void div2_close (void)
/*< free allocated storage >*/
{
    if (gauss) {
	gauss2_close();
    } else { 
	sf_triangle2_close();
    }
    sf_conjgrad_close();
    free (p);
}

void div2 (float* num, float* den,  float* rat)
/*< smooth division: rat=num/den >*/
{
    sf_weight_init(den);
    if (gauss) {
	sf_conjgrad(NULL, sf_weight_lop,sf_freqfilt2_lop,p,rat,num,niter);
    } else {
	sf_conjgrad(NULL, sf_weight_lop,sf_triangle2_lop,p,rat,num,niter); 
    }
}

/* 	$Id$	 */

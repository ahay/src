/* 1-D smooth division */
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

#include "div1.h"
#include "gauss.h"
#include "weight.h"

static int niter;
static float *p;
static bool gauss;

void div1_init(int n1      /* data size */, 
	       float f1    /* smoothing */, 
	       int niter1  /* number of iterations */, 
	       bool gauss1 /* if use exact gaussian */) 
/*< initialize >*/
{
    niter = niter1;
    gauss = gauss1;

    if (gauss) {
	gauss_init(n1,f1);
    } else {
	sf_triangle1_init((int) f1, n1);
    }
    sf_conjgrad_init(n1, n1, n1, n1, 1., 1.e-6, true, false);
    p = sf_floatalloc (n1);
}

void div1_close (void)
/*< free allocated storage >*/
{
    if (gauss) {
	gauss_close();
    } else { 
	sf_triangle1_close();
    }
    sf_conjgrad_close();
    free (p);
}

void div1 (float* num, float* den,  float* rat)
/*< smooth division: rat=num/den >*/
{
    weight_init(den);
    if (gauss) {
	sf_conjgrad(NULL, weight_lop,sf_freqfilt_lop,p,rat,num,niter);
    } else {
	sf_conjgrad(NULL, weight_lop,sf_triangle1_lop,p,rat,num,niter); 
    }
}

/* 	$Id$	 */

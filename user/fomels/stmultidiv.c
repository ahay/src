/* Smooth division with several components and nonstationary smoothing */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static int n, nk;
static float lam, **bas, *prev;

void stmultidiv_init(int nw       /* number of components */, 
		     int nd       /* data size */,
		     float** den  /* denominator [nw*nd] */,
		     float lambda)
/*< initialize >*/
{
    n = nd;
    nk = nw;

    lam = lambda*lambda;
    
    bas = den;

    prev = sf_floatalloc(nk);
}

void stmultidiv_close(void)
{
    free(prev);
}

void stmultidiv (float* num  /* numerator */, 
		 float** rat /* ratio */)
/*< smoothly divide num/rat >*/
{
    int i, k;
    float res, eps;

    for (k=0; k < nk; k++) {
	prev[k] = 0.0;
    }
    for (i=0; i < n; i++) {
	res = num[i];
	eps = lam;
	for (k=0; k < nk; k++) {
	    res -= bas[k][i]*prev[k];
	    eps += bas[k][i]*bas[k][i];
	}
	for (k=0; k < nk; k++) {
	    prev[k] += res/eps*bas[k][i];
	    rat[k][i] = prev[k];
	}
    }
}

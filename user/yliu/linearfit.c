/* Linear fitting using least-squares regression. */
/*
  Copyright (C) 2011 Jilin University
  
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

#include "linearfit.h"

void linearfit (int nd        /* data length */,
		float *dat    /* input data */,
		float **func   /* fitting function */,
		float *sol    /* fitting coefficients */)
/*< solve a+bx ~= dat >*/
{
    int i1, ic, nc, id;
    float *rhs, **mat;

    nc = 2;
    sf_gaussel_init(nc);
    rhs = sf_floatalloc(nc);
    mat = sf_floatalloc2(nc,nc);

	
    /* compute A'A matrix */
    for (ic=0; ic < nc; ic++) {
	for (id=0; id <= ic; id++) {
	    mat[ic][id] = 0.;
	    for (i1=0; i1 < nd; i1++) {
		mat[ic][id] += func[ic][i1]*func[id][i1];
	    }
	    mat[id][ic] = mat[ic][id];
	}
    }
    
    /* compute A'd */
    for (ic=0; ic < nc; ic++) {
	rhs[ic] = 0.;
	for (i1=0; i1 < nd; i1++) {
	    rhs[ic] += func[ic][i1]*dat[i1];
	}
    }
    
    /* inversion */
    sf_gaussel_solve(mat,rhs,sol);
    
    /* compute Ac */
    for (i1=0; i1 < nd; i1++) {
	dat[i1] = 0.;
	for (ic=0; ic < nc; ic++) {
	    dat[i1] += func[ic][i1]*sol[ic];
	}
    }

}

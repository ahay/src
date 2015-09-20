/* Preconditioned conjugate gradients */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "pcg.h"

#ifndef _pcg_h

typedef void (*normal_operator)(int,const float*,float*);
/*^*/

#endif

void conjgrad(normal_operator oper /* operator to invert */,
	      normal_operator prec /* preconditioning */,
	      int nd /* data size */, 
	      const float *d /* data (right-hand side) */, 
	      const float *x0 /* starting model */, 
	      float *x /* estimated model */,
	      int niter /* maximum number of iterations */,
	      float tol /* tolerance */)
/*< Conjugate-gradient algorithm for solving B x = d >*/
{
    int ix, iter;
    float *g, *cg, *s, *cs, *z;
    double gn, gnp=0.0, beta;

    g = sf_floatalloc(nd);
    cg = sf_floatalloc(nd);
    s = sf_floatalloc(nd);
    cs = sf_floatalloc(nd);
    z = sf_floatalloc(nd);

    for (ix=0; ix < nd; ix++) {
	x[ix] = (NULL != x0)? x0[ix]: 0.0f;
    }
    oper(nd,x,g);
    for (ix=0; ix < nd; ix++) {
	g[ix] -= d[ix];
    }

    for (iter=0; iter < niter; iter++) {
	oper(nd,g,cg);
	if (NULL != prec) {
	    prec(nd,g,z);
	} else {
	    for (ix=0; ix < nd; ix++) {
		z[ix] = g[ix];
	    }
	}

	gn = cblas_dsdot(nd,g,1,z,1);
	
	if (gn < tol) {
	    sf_warning("Converged after %d iterations",iter);
	    return;
	}
	
	sf_warning("iteration %d: %g",iter+1,gn);
	
        if (0==iter) {
	    for (ix=0; ix < nd; ix++) {
		s[ix] = g[ix];
		cs[ix] = cg[ix];
	    }
        } else {
            beta = gn/gnp;

            /* cs = cg + cs*beta */
	    cblas_saxpy(nd,beta,cs,1,cg,1); 
	    cblas_sswap(nd,cs,1,cg,1);

	    cblas_saxpy(nd,beta,s,1,z,1); /* s = z + s*beta */
	    cblas_sswap(nd,s,1,z,1);

	}
        gnp = gn;

	beta = -gn/cblas_dsdot(nd,cs,1,s,1);

	cblas_saxpy(nd,beta,s,1,x,1);
	cblas_saxpy(nd,beta,cs,1,g,1);
    }

    free(g); free(cg);
    free(s); free(cs);
}


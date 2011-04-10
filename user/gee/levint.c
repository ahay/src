/* 1-D inverse interpolation with PEF estimation */
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

#include "levint.h"

#include "tcai2.h"
#include "peftc.h"
#include "lint1.h"

void levint1 (int niter       /* number of iterations */, 
	      int nm          /* model size */, 
	      int nr          /* model size + filter size */, 
	      int nd          /* data size */, 
	      float *coord    /* data coordinates [nd] */, 
	      const float *dd /* data [nd] */, 
	      float o1        /* model origin */, 
	      float d1        /* model sampling */, 
	      float *rr       /* initial model + filter [nr] */, 
	      float eps       /* regularization scaling */)
/*< interpolate >*/ 
{
    float *aa;
    bool *ma;
    int ia, ir, na;

    ma = sf_boolalloc(nr);
    na = nr-nm;
    aa = rr+nm;

    lint1_init (o1, d1, coord);
    peftc_init (na, nm, aa, rr);
    tcai2_init (na, nm, aa);

    /* starting guess */
    for (ir=0; ir < nr; ir++) {
	ma[ir] = (bool) (ir >= nm);
	rr[ir] = 0.;
    }
    aa[0] = 1.;
    aa[1] = -2.;
    aa[2] = 1.;
    for (ia=3; ia < na; ia++) {
	aa[ia] = 0.;
    }

    sf_solver_reg (lint1_lop, sf_cgstep, peftc_lop, nr, nr, nd, rr, dd, niter,
		   eps, "nlreg", tcai2_lop, "x0", rr, "known", ma, "end");

    /* free filter coefficients */
    for (ir=nm+1; ir < nr; ir++) {
	ma[ir] = false;
    }
    
    sf_cgstep_close();

    sf_solver_reg (lint1_lop, sf_cgstep, peftc_lop, nr, nr, nd, rr, dd, niter,
		   eps, "nlreg", tcai2_lop, "x0", rr, "known", ma, "end");

    sf_cgstep_close();

    free(ma);
}
